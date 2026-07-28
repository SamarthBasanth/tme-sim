"""Binary wire payloads for the tmesim Studio API (API_CONTRACT.md v1).

Two formats, both **little-endian**, both laid out so the front end can slice a
single ``ArrayBuffer`` into typed arrays with zero parsing and hand them straight
to WebGL2. JSON for 850k cells is not an option: the geometry payload for the
GSE324268 section (851,509 cells) is 9 bytes/cell = 7.7 MB of binary versus
~40 MB of text that then has to be parsed into JS numbers.

Geometry (``/api/runs/{id}/geometry``)::

    magic "TMSG" (4B) | version u32 | n_cells u32 | extent_x f32 | extent_y f32
    | x f32[n] | y f32[n] | type u8[n]

Field (``/api/runs/{id}/field/{name}``)::

    magic "TMSF" (4B) | version u32 | rows u32 | cols u32 | dx_um f32
    | min f32 | max f32 | values f32[rows*cols]

Mechanism (``/api/runs/{id}/mechanism``)::

    magic "TMSM" (4B) | version u32 | n_cells u32 | n_labels u32 | json_len u32
    | aggregate json[json_len] (UTF-8) | label u8[n_cells]

One byte per cell, in the SAME order as the geometry payload, so the front end
uploads it as a second vertex attribute against the buffer it already has. The
aggregate breakdown (counts, percentages, thresholds, the reconciliation against
``kill_fraction``) rides in the same response as a short JSON block, because a
legend without its counts is decoration.

Time-lapse (``/api/runs/{id}/timelapse``)::

    magic "TMSL" (4B) | version u32 | n_frames u32 | interval_s f32 | window_s f32
    | n_frames x [ frame_index u32 | time_s f32 | <a complete TMSF block> ]

The per-frame payload is **byte-identical to what** :func:`encode_field`
**produces**, prefixed with its frame index and time stamp — the existing TMSF
shape extended with a frame index, not a new grid format. The front end decodes
a frame with the TMSF reader it already ships, and each frame carries its own
min/max so a scrubber can show either a per-frame or a global colour range
without a second pass.

Field values are **row-major** (``values[r*cols + c]``), matching the numpy
``(ny, nx)`` grid the engine works on.

``min`` / ``max`` are computed over the FINITE entries only. A field that is
entirely non-finite reports ``nan`` for both rather than a fabricated range —
the front end must treat a non-finite range as "no data", not as ``0..1``.

Nothing here interprets the science: encoders take arrays, decoders return
arrays, and :func:`decode_geometry` / :func:`decode_field` exist so the byte
layout can be round-tripped in a test rather than trusted.
"""
from __future__ import annotations

import json
import struct
from typing import Dict, Optional, Sequence

import numpy as np

__all__ = [
    "MAGIC_GEOMETRY",
    "MAGIC_FIELD",
    "MAGIC_MECHANISM",
    "MAGIC_TIMELAPSE",
    "PAYLOAD_VERSION",
    "GEOMETRY_HEADER_STRUCT",
    "FIELD_HEADER_STRUCT",
    "MECHANISM_HEADER_STRUCT",
    "TIMELAPSE_HEADER_STRUCT",
    "TIMELAPSE_FRAME_STRUCT",
    "FIELD_NAMES",
    "MAX_TIMELAPSE_FRAMES",
    "encode_geometry",
    "decode_geometry",
    "encode_field",
    "decode_field",
    "encode_mechanism",
    "decode_mechanism",
    "encode_timelapse",
    "decode_timelapse",
]

MAGIC_GEOMETRY = b"TMSG"
MAGIC_FIELD = b"TMSF"
MAGIC_MECHANISM = b"TMSM"
MAGIC_TIMELAPSE = b"TMSL"

#: Bumped only if the byte layout changes. The contract is frozen at 1.
PAYLOAD_VERSION = 1

#: "<" = little-endian, no padding. 4s + u32 + u32 + f32 + f32 = 20 bytes.
GEOMETRY_HEADER_STRUCT = struct.Struct("<4sIIff")
#: 4s + u32 + u32 + u32 + f32 + f32 + f32 = 28 bytes.
FIELD_HEADER_STRUCT = struct.Struct("<4sIIIfff")
#: 4s + u32 + u32 + u32 + u32 = 20 bytes.
MECHANISM_HEADER_STRUCT = struct.Struct("<4sIIII")
#: 4s + u32 + u32 + f32 + f32 = 20 bytes.
TIMELAPSE_HEADER_STRUCT = struct.Struct("<4sIIff")
#: u32 + f32 = 8 bytes, prefixing each embedded TMSF block.
TIMELAPSE_FRAME_STRUCT = struct.Struct("<If")

#: The field names the contract exposes on ``/api/runs/{id}/field/{name}``.
FIELD_NAMES = ("vascular_density", "ecm_density", "oxygen", "drug_auc")

#: Hard ceiling on time-lapse frames. A frame is a full f32 grid, so the payload
#: grows linearly: at the flagship 677x425 grid one frame is 1.15 MB and 40 is
#: 46 MB. The cap is a wire-size decision, not a physical one.
MAX_TIMELAPSE_FRAMES = 40


def _finite_range(values: np.ndarray) -> tuple:
    """(min, max) over the finite entries; (nan, nan) if none are finite.

    Never substitutes a default range — a caller that gets ``nan`` knows the
    field carried no usable value, which is different from a field of zeros.
    """
    finite = np.isfinite(values)
    if not finite.any():
        return float("nan"), float("nan")
    usable = values[finite]
    return float(usable.min()), float(usable.max())


def encode_geometry(
    x: Sequence[float],
    y: Sequence[float],
    cell_type: Sequence[int],
    extent_um: Sequence[float],
    version: int = PAYLOAD_VERSION,
) -> bytes:
    """Pack per-cell positions + types into the TMSG payload.

    Parameters
    ----------
    x, y : array-like of float
        Cell centroids in micrometres, in the reconstruction's domain frame.
    cell_type : array-like of int
        ``tmesim.core.CellType`` codes, one per cell. Must fit in ``uint8``
        (the enum is small); a value outside 0..255 is a programming error and
        raises rather than wrapping silently.
    extent_um : (width_um, height_um)
        Domain extent, so the front end can set up its projection without a
        second request.
    """
    xs = np.ascontiguousarray(np.asarray(x, dtype="<f4"))
    ys = np.ascontiguousarray(np.asarray(y, dtype="<f4"))
    if xs.ndim != 1 or ys.ndim != 1:
        raise ValueError(f"x and y must be 1-D; got {xs.shape} and {ys.shape}")
    if xs.size != ys.size:
        raise ValueError(f"x ({xs.size}) and y ({ys.size}) length mismatch")

    types_i = np.asarray(cell_type)
    if types_i.size != xs.size:
        raise ValueError(
            f"cell_type ({types_i.size}) does not match n_cells ({xs.size})"
        )
    types_i = types_i.astype(np.int64, copy=False)
    if types_i.size and (types_i.min() < 0 or types_i.max() > 255):
        raise ValueError(
            f"cell_type values must fit in uint8; got range "
            f"[{int(types_i.min())}, {int(types_i.max())}]"
        )
    types = np.ascontiguousarray(types_i.astype(np.uint8))

    ext = np.asarray(extent_um, dtype=np.float64).ravel()
    if ext.size != 2:
        raise ValueError(f"extent_um must be (width, height); got {extent_um!r}")

    header = GEOMETRY_HEADER_STRUCT.pack(
        MAGIC_GEOMETRY, int(version), int(xs.size), float(ext[0]), float(ext[1])
    )
    return b"".join((header, xs.tobytes(), ys.tobytes(), types.tobytes()))


def decode_geometry(buf: bytes) -> Dict:
    """Inverse of :func:`encode_geometry`. Raises on a malformed payload."""
    if len(buf) < GEOMETRY_HEADER_STRUCT.size:
        raise ValueError(
            f"geometry payload truncated: {len(buf)} B < "
            f"{GEOMETRY_HEADER_STRUCT.size} B header"
        )
    magic, version, n, ex, ey = GEOMETRY_HEADER_STRUCT.unpack_from(buf, 0)
    if magic != MAGIC_GEOMETRY:
        raise ValueError(f"bad magic {magic!r}; expected {MAGIC_GEOMETRY!r}")

    off = GEOMETRY_HEADER_STRUCT.size
    expected = off + n * 4 + n * 4 + n
    if len(buf) != expected:
        raise ValueError(
            f"geometry payload is {len(buf)} B; expected {expected} B for "
            f"n_cells={n}"
        )
    xs = np.frombuffer(buf, dtype="<f4", count=n, offset=off)
    off += n * 4
    ys = np.frombuffer(buf, dtype="<f4", count=n, offset=off)
    off += n * 4
    types = np.frombuffer(buf, dtype=np.uint8, count=n, offset=off)
    return {
        "magic": magic,
        "version": int(version),
        "n_cells": int(n),
        "extent_um": (float(ex), float(ey)),
        "x": xs,
        "y": ys,
        "cell_type": types,
    }


def encode_field(
    values: np.ndarray,
    dx_um: float,
    version: int = PAYLOAD_VERSION,
) -> bytes:
    """Pack a 2-D ``(rows, cols)`` grid field into the TMSF payload.

    ``values`` is written row-major as ``float32``. ``min`` / ``max`` in the
    header are the finite range (see :func:`_finite_range`), so the front end
    can build its colour ramp without a second pass over 300k floats.
    """
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"field must be 2-D (rows, cols); got shape {arr.shape}")
    rows, cols = arr.shape
    vmin, vmax = _finite_range(arr)
    payload = np.ascontiguousarray(arr.astype("<f4"))
    header = FIELD_HEADER_STRUCT.pack(
        MAGIC_FIELD, int(version), int(rows), int(cols),
        float(dx_um), float(vmin), float(vmax),
    )
    return header + payload.tobytes()


def decode_field(buf: bytes) -> Dict:
    """Inverse of :func:`encode_field`. Raises on a malformed payload."""
    if len(buf) < FIELD_HEADER_STRUCT.size:
        raise ValueError(
            f"field payload truncated: {len(buf)} B < "
            f"{FIELD_HEADER_STRUCT.size} B header"
        )
    magic, version, rows, cols, dx, vmin, vmax = FIELD_HEADER_STRUCT.unpack_from(buf, 0)
    if magic != MAGIC_FIELD:
        raise ValueError(f"bad magic {magic!r}; expected {MAGIC_FIELD!r}")

    off = FIELD_HEADER_STRUCT.size
    n = rows * cols
    expected = off + n * 4
    if len(buf) != expected:
        raise ValueError(
            f"field payload is {len(buf)} B; expected {expected} B for "
            f"{rows}x{cols}"
        )
    vals = np.frombuffer(buf, dtype="<f4", count=n, offset=off).reshape(rows, cols)
    return {
        "magic": magic,
        "version": int(version),
        "rows": int(rows),
        "cols": int(cols),
        "dx_um": float(dx),
        "min": float(vmin),
        "max": float(vmax),
        "values": vals,
    }


# --------------------------------------------------------------------------- #
# Mechanism attribution
# --------------------------------------------------------------------------- #
def encode_mechanism(
    labels: Sequence[int],
    aggregate: Dict,
    version: int = PAYLOAD_VERSION,
) -> bytes:
    """Pack per-cell mechanism labels + the aggregate breakdown into TMSM.

    Parameters
    ----------
    labels : array-like of int
        One :mod:`tmesim.sim.mechanism` label code per cell, in the SAME order
        as the geometry payload. Must fit in ``uint8``; a value outside 0..255
        is a programming error and raises rather than wrapping silently.
    aggregate : dict
        The JSON-serialisable breakdown (counts, thresholds, reconciliation).
        It must NOT contain the label array — that is the binary half of this
        payload and duplicating it as JSON would defeat the point.
    """
    lab_i = np.asarray(labels)
    lab_i = lab_i.astype(np.int64, copy=False)
    if lab_i.ndim != 1:
        raise ValueError(f"labels must be 1-D; got shape {lab_i.shape}")
    if lab_i.size and (lab_i.min() < 0 or lab_i.max() > 255):
        raise ValueError(
            f"label codes must fit in uint8; got range "
            f"[{int(lab_i.min())}, {int(lab_i.max())}]"
        )
    lab = np.ascontiguousarray(lab_i.astype(np.uint8))

    body = json.dumps(aggregate, allow_nan=False,
                      separators=(",", ":")).encode("utf-8")
    n_labels = int(len(aggregate.get("label_names") or ()))
    header = MECHANISM_HEADER_STRUCT.pack(
        MAGIC_MECHANISM, int(version), int(lab.size), n_labels, len(body))
    return b"".join((header, body, lab.tobytes()))


def decode_mechanism(buf: bytes) -> Dict:
    """Inverse of :func:`encode_mechanism`. Raises on a malformed payload."""
    if len(buf) < MECHANISM_HEADER_STRUCT.size:
        raise ValueError(
            f"mechanism payload truncated: {len(buf)} B < "
            f"{MECHANISM_HEADER_STRUCT.size} B header"
        )
    magic, version, n, n_labels, json_len = MECHANISM_HEADER_STRUCT.unpack_from(buf, 0)
    if magic != MAGIC_MECHANISM:
        raise ValueError(f"bad magic {magic!r}; expected {MAGIC_MECHANISM!r}")
    off = MECHANISM_HEADER_STRUCT.size
    expected = off + json_len + n
    if len(buf) != expected:
        raise ValueError(
            f"mechanism payload is {len(buf)} B; expected {expected} B for "
            f"n_cells={n}, json_len={json_len}"
        )
    aggregate = json.loads(bytes(buf[off:off + json_len]).decode("utf-8"))
    labels = np.frombuffer(buf, dtype=np.uint8, count=n, offset=off + json_len)
    return {
        "magic": magic,
        "version": int(version),
        "n_cells": int(n),
        "n_labels": int(n_labels),
        "aggregate": aggregate,
        "labels": labels,
    }


# --------------------------------------------------------------------------- #
# Time-lapse
# --------------------------------------------------------------------------- #
def encode_timelapse(
    frames: np.ndarray,
    frame_times: Sequence[float],
    dx_um: float,
    interval_s: Optional[float] = None,
    version: int = PAYLOAD_VERSION,
) -> bytes:
    """Pack a ``(n_frames, rows, cols)`` concentration series into TMSL.

    Each frame is written as a complete, byte-identical :func:`encode_field`
    TMSF block prefixed by its frame index and time stamp, so the front end
    reuses its existing TMSF decoder per frame and every frame carries its own
    finite min/max.

    ``interval_s`` is the NOMINAL spacing the solver was asked for; the realised
    ``frame_times`` are authoritative, because the capture points are rounded
    onto the solver's own step grid.
    """
    arr = np.asarray(frames, dtype=np.float64)
    if arr.ndim != 3:
        raise ValueError(
            f"frames must be 3-D (n_frames, rows, cols); got shape {arr.shape}")
    times = np.asarray(frame_times, dtype=np.float64).ravel()
    if times.size != arr.shape[0]:
        raise ValueError(
            f"frame_times has {times.size} entries but there are "
            f"{arr.shape[0]} frames")
    if arr.shape[0] > MAX_TIMELAPSE_FRAMES:
        raise ValueError(
            f"{arr.shape[0]} frames exceeds MAX_TIMELAPSE_FRAMES="
            f"{MAX_TIMELAPSE_FRAMES}; the payload is n_frames full f32 grids")

    window = float(times[-1]) if times.size else 0.0
    nominal = float(interval_s) if interval_s is not None else (
        window / (times.size - 1) if times.size > 1 else 0.0)
    out = [TIMELAPSE_HEADER_STRUCT.pack(
        MAGIC_TIMELAPSE, int(version), int(arr.shape[0]), nominal, window)]
    for i in range(arr.shape[0]):
        out.append(TIMELAPSE_FRAME_STRUCT.pack(i, float(times[i])))
        out.append(encode_field(arr[i], dx_um, version=version))
    return b"".join(out)


def decode_timelapse(buf: bytes) -> Dict:
    """Inverse of :func:`encode_timelapse`. Raises on a malformed payload."""
    if len(buf) < TIMELAPSE_HEADER_STRUCT.size:
        raise ValueError(
            f"timelapse payload truncated: {len(buf)} B < "
            f"{TIMELAPSE_HEADER_STRUCT.size} B header"
        )
    magic, version, n_frames, interval, window = \
        TIMELAPSE_HEADER_STRUCT.unpack_from(buf, 0)
    if magic != MAGIC_TIMELAPSE:
        raise ValueError(f"bad magic {magic!r}; expected {MAGIC_TIMELAPSE!r}")

    off = TIMELAPSE_HEADER_STRUCT.size
    frames = []
    for expected_index in range(n_frames):
        if len(buf) < off + TIMELAPSE_FRAME_STRUCT.size + FIELD_HEADER_STRUCT.size:
            raise ValueError(
                f"timelapse truncated at frame {expected_index} of {n_frames}")
        index, t = TIMELAPSE_FRAME_STRUCT.unpack_from(buf, off)
        if int(index) != expected_index:
            raise ValueError(
                f"timelapse frame index out of order: got {index}, expected "
                f"{expected_index}")
        off += TIMELAPSE_FRAME_STRUCT.size
        fmagic, _v, rows, cols, _dx, _mn, _mx = \
            FIELD_HEADER_STRUCT.unpack_from(buf, off)
        if fmagic != MAGIC_FIELD:
            raise ValueError(
                f"timelapse frame {expected_index} is not a TMSF block "
                f"(magic {fmagic!r})")
        size = FIELD_HEADER_STRUCT.size + rows * cols * 4
        block = decode_field(bytes(buf[off:off + size]))
        block["frame"] = int(index)
        block["time_s"] = float(t)
        frames.append(block)
        off += size
    if off != len(buf):
        raise ValueError(
            f"timelapse payload has {len(buf) - off} trailing bytes after "
            f"{n_frames} frames")
    return {
        "magic": magic,
        "version": int(version),
        "n_frames": int(n_frames),
        "interval_s": float(interval),
        "window_s": float(window),
        "frames": frames,
    }
