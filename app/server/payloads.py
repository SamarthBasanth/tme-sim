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

import struct
from typing import Dict, Sequence

import numpy as np

__all__ = [
    "MAGIC_GEOMETRY",
    "MAGIC_FIELD",
    "PAYLOAD_VERSION",
    "GEOMETRY_HEADER_STRUCT",
    "FIELD_HEADER_STRUCT",
    "FIELD_NAMES",
    "encode_geometry",
    "decode_geometry",
    "encode_field",
    "decode_field",
]

MAGIC_GEOMETRY = b"TMSG"
MAGIC_FIELD = b"TMSF"

#: Bumped only if the byte layout changes. The contract is frozen at 1.
PAYLOAD_VERSION = 1

#: "<" = little-endian, no padding. 4s + u32 + u32 + f32 + f32 = 20 bytes.
GEOMETRY_HEADER_STRUCT = struct.Struct("<4sIIff")
#: 4s + u32 + u32 + u32 + f32 + f32 + f32 = 28 bytes.
FIELD_HEADER_STRUCT = struct.Struct("<4sIIIfff")

#: The field names the contract exposes on ``/api/runs/{id}/field/{name}``.
FIELD_NAMES = ("vascular_density", "ecm_density", "oxygen", "drug_auc")


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
