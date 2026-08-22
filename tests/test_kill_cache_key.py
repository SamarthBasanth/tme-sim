"""The kill-calibration cache key must cover every drug parameter.

A cache key that under-keys is a check that has stopped checking: the guard
still runs, still passes, and no longer constrains anything. That is exactly
what happened here. ``tmesim/validation/kill_calibration.py`` keyed its cached
unit-transport solves on

    f"{name}|{dx_um}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}"

— drug name and grid only, with **no drug parameter at all** — while the cached
object (``e_unit``, the per-cell unit-concentration exposure field) is the
output of a PDE solve that consumed ``diffusivity``, ``binding``, ``clearance``
and ``bbb_penetrance``. Any change to a descriptor, a partition constant or
``HBD_LOGP_DECREMENT`` was therefore silently ignored on a warm cache, and
``vivum validate`` published 3-D:2-D kill ratios that did not recompute (see
``results/hbd_convention_audit.md`` sections 6.2-6.3).

The fix derives the key from ``dataclasses.fields(DrugParams)`` so it cannot go
stale when a parameter is added years from now. These tests are what make that
structural rather than remembered: they perturb each field in turn and require
every perturbation to move the key. They would have FAILED on the old key for
every field, and on the breast twin's key
(``experiments/kill_calibration_breast.py``, which keys on ``ic50`` — the one
drug parameter that does *not* determine ``e_unit`` — and omits the four that
do) for every field except ``ic50``. Detecting the defect that was actually
present is the bar a regression test has to clear.
"""
from __future__ import annotations

import ast
import dataclasses
import pathlib
import re

import pytest

from tmesim.core import DRUG_LIBRARY, DrugParams
from tmesim.validation import kill_calibration as KC


def _perturb(value):
    """Return a value guaranteed different from ``value`` for key purposes."""
    if isinstance(value, bool):
        return not value
    if isinstance(value, (int, float)):
        return value * 1.5 + 1.0
    return "PERTURBED"


def _base() -> DrugParams:
    base = DRUG_LIBRARY.get("doxorubicin-like") or DRUG_LIBRARY.get("doxorubicin")
    assert base is not None, "no doxorubicin entry in DRUG_LIBRARY"
    return base


def test_kill_cache_key_covers_every_drug_parameter():
    """Every field of DrugParams must move the key. No exceptions, ever.

    Deliberately over-keys: ``ic50`` and ``hill`` do not affect ``e_unit``, so
    including them costs an occasional needless resolve. That is the safe
    direction -- a spurious recompute is cheap, a stale published number is not.
    """
    base = _base()
    k0 = KC.unit_transport_cache_key(base.name, base, 12.0)
    unmoved = []
    for f in dataclasses.fields(DrugParams):
        alt = _perturb(getattr(base, f.name))
        moved = KC.unit_transport_cache_key(base.name,
                         dataclasses.replace(base, **{f.name: alt}),
                         12.0)
        if moved == k0:
            unmoved.append(f.name)
    assert not unmoved, f"cache key ignores DrugParams fields: {unmoved}"


def test_kill_cache_key_covers_the_four_transport_determining_fields():
    """Named explicitly, because these four are what ``unit_transport`` reads.

    The generic test above subsumes this one today. It is kept separate so that
    if someone ever narrows the key "for performance", the failure message names
    the physics rather than a list of field names.
    """
    base = _base()
    k0 = KC.unit_transport_cache_key(base.name, base, 12.0)
    for field in ("diffusivity", "binding", "clearance", "bbb_penetrance"):
        alt = _perturb(getattr(base, field))
        assert KC.unit_transport_cache_key(base.name,
                        dataclasses.replace(base, **{field: alt}),
                        12.0) != k0, (
            f"cache key ignores {field!r}, which unit_transport consumes -- "
            "cached exposure fields will be served for the wrong drug")


def test_kill_cache_key_still_separates_name_and_grid():
    """The properties the OLD key did have must survive the fix."""
    base = _base()
    k0 = KC.unit_transport_cache_key(base.name, base, 12.0)
    assert KC.unit_transport_cache_key(base.name, base, 10.0) != k0
    assert KC.unit_transport_cache_key("other-drug", base, 12.0) != k0
    # deterministic across calls: a key that is not stable never hits
    assert KC.unit_transport_cache_key(base.name, base, 12.0) == k0


def test_kill_cache_key_is_not_the_old_unparameterised_format():
    """Guard against a revert to the hand-written f-string.

    The old key was fully determined by ``(name, dx_um)`` plus module
    constants, so two drugs differing only in their parameters collided.
    """
    base = _base()
    twin = dataclasses.replace(base, diffusivity=base.diffusivity * 2.0)
    assert KC.unit_transport_cache_key(base.name, base, 12.0) != KC.unit_transport_cache_key(base.name, twin, 12.0), (
        "two drugs with different diffusivity share a cache key -- this is the "
        "pre-fix behaviour that published stale kill ratios")


@pytest.mark.parametrize("name", sorted(KC.DRUGS))
def test_kill_cache_key_builds_for_every_published_drug(name):
    """The key must be constructible for each drug the scorecard publishes."""
    proto = _base()
    key = KC.unit_transport_cache_key(name, dataclasses.replace(proto, name=name), 12.0)
    assert key.startswith(f"{name}|12.0|")
    assert len(key.rsplit("|", 1)[-1]) == 16


# --------------------------------------------------------------------------- #
# One function, imported everywhere
# --------------------------------------------------------------------------- #
#: Every module that reads or writes ``data_cache/kill_calibration*_runs.json``.
#: Five sites used to hand-build this key independently and they did not agree
#: with each other, which is how the defect survived. The real bug was the
#: duplication, so this is the test that stops a sixth site reintroducing it.
_CACHE_KEY_CALL_SITES = (
    "tmesim/validation/kill_calibration.py",
    "experiments/kill_calibration_breast.py",
    "experiments/grid_convergence_gbm.py",
    "experiments/grid_convergence_breast.py",
    "experiments/dtissue_sensitivity.py",
)

#: The shape of the old hand-written key, as a regex over source text. Matches
#: an f-string interpolating a drug name straight into a ``|``-joined key.
_HANDBUILT = re.compile(
    r'f"\{\s*name\s*\}\|.*\|\{\s*(?:DT_S|dt_s)\s*\}')


#: The definition site, which is the ONE place allowed to build the string.
_DEFINITION_SITE = "tmesim/validation/kill_calibration.py"


def _repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def _code_lines(path: pathlib.Path):
    """Source lines with docstrings removed.

    Docstrings quote the old key on purpose — that is the record of what went
    wrong — and a guard that cannot tell documentation from code would force
    the explanation to be deleted. Deleting the explanation is the failure mode
    this whole exercise is about, so the guard is taught to parse instead.
    """
    src = path.read_text()
    doc_lines = set()
    for node in ast.walk(ast.parse(src)):
        if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef)):
            continue
        body = getattr(node, "body", None)
        if (body and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)):
            d = body[0]
            doc_lines.update(range(d.lineno, (d.end_lineno or d.lineno) + 1))
    return [(i, ln) for i, ln in enumerate(src.splitlines(), 1)
            if i not in doc_lines]


@pytest.mark.parametrize("relpath", _CACHE_KEY_CALL_SITES)
def test_no_site_hand_builds_the_unit_transport_cache_key(relpath):
    """No module may reconstruct the key as an f-string. Import it.

    A hand-built key is not merely duplication: it silently diverges. The GBM
    site omitted every drug parameter; the breast site included ``ic50``, the
    one that does not matter, and omitted the four that do; and two readers
    reconstructed the GBM shape to read fields they then re-scored. Only the
    definition site is allowed to build the string, and there it must fold in
    the parameter digest.
    """
    offenders = [(i, ln) for i, ln in _code_lines(_repo_root() / relpath)
                 if _HANDBUILT.search(ln)]
    if relpath == _DEFINITION_SITE:
        # allowed exactly once, and only in its digest-bearing form
        offenders = [(i, ln) for i, ln in offenders if "{digest}" not in ln]
    assert not offenders, (
        f"{relpath} hand-builds the unit-transport cache key:\n  "
        + "\n  ".join(f"line {i}: {ln.strip()}" for i, ln in offenders)
        + "\nImport kill_calibration.unit_transport_cache_key instead -- a key "
          "rebuilt by hand goes stale the next time DrugParams gains a field.")


def test_the_definition_site_builds_the_key_exactly_once():
    """One function, one construction. Two would already be a divergence."""
    built = [ln for _, ln in _code_lines(_repo_root() / _DEFINITION_SITE)
             if _HANDBUILT.search(ln)]
    assert len(built) == 1, (
        f"expected exactly one cache-key construction in {_DEFINITION_SITE}, "
        f"found {len(built)}: {built}")
    assert "{digest}" in built[0]


@pytest.mark.parametrize("relpath", _CACHE_KEY_CALL_SITES[1:])
def test_every_cache_site_imports_the_shared_key_function(relpath):
    """The counterpart: each site must actually import the shared helper."""
    src = (_repo_root() / relpath).read_text()
    assert ("unit_transport_cache_key" in src
            or "require_cached_unit_field" in src), (
        f"{relpath} touches the unit-transport cache but imports neither "
        "unit_transport_cache_key nor require_cached_unit_field")


# --------------------------------------------------------------------------- #
# A miss must be loud
# --------------------------------------------------------------------------- #
def test_cache_miss_raises_rather_than_returning_nothing():
    """A reader that cannot find its entry must fail, not publish a short table.

    This is the half of the defect that is easy to skip. A silent miss is
    indistinguishable from "nothing was cached": it does not look like a
    defect, it looks like an absence, and the table it produces still looks
    complete. Every reader now goes through ``require_cached_unit_field``.
    """
    base = _base()
    with pytest.raises(KeyError) as exc:
        KC.require_cached_unit_field({}, base.name, base, 12.0, source="(empty)")
    msg = str(exc.value)
    assert "no cached unit-transport field" in msg
    assert "wanted key" in msg
    assert "STALE" in msg


def test_cache_hit_returns_the_field_and_is_not_over_strict():
    """The loud path must not make a legitimate hit fail."""
    base = _base()
    key = KC.unit_transport_cache_key(base.name, base, 12.0)
    cache = {key: {"e_unit": [1.0, 2.0, 3.0]}}
    got = KC.require_cached_unit_field(cache, base.name, base, 12.0,
                                       source="(fixture)")
    assert got.tolist() == [1.0, 2.0, 3.0]


def test_cache_miss_message_names_the_stale_entries_it_found():
    """The message must be actionable: a key that exists under the OLD format
    is the single most likely cause, so say so and list what is present."""
    base = _base()
    stale = f"{base.name}|12.0|10.0|21600.0|10"
    with pytest.raises(KeyError) as exc:
        KC.require_cached_unit_field({stale: {"e_unit": []}}, base.name, base,
                                     12.0, source="(fixture)")
    assert stale in str(exc.value)
    assert "PRE-FIX" in str(exc.value)
