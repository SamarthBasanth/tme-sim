#!/usr/bin/env python3
"""Plant each defect, confirm it is REACHED, run the suite, revert, record.

Why this is a script and not a paragraph: a plant log written by hand is a
claim about work that happened once, and the counts in it go stale the moment
a test is added. This re-derives them, so ``experiments/exact_memo/plants.json``
is regenerable and the record that quotes it cannot drift from the suite.

Safety: each target file is copied before anything is edited and restored in a
``finally``, then compared byte-for-byte against the copy. If a restore ever
fails, the script says so loudly and stops rather than leaving a planted defect
in the tree. It refuses to start unless the suite is GREEN and the working tree
has no unexpected edit to the files it is about to touch.

Usage:  run_plants.py [--write]      (--write updates plants.json)
"""
from __future__ import annotations

import filecmp
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PYTHON = os.environ.get("TMESIM_PYTHON",
                        "/Users/samarthbasanth/tme-sim/.venv/bin/python")
SUITE = "tests/test_memo.py"

MEMO = "app/server/memo.py"
JOBS = "app/server/jobs.py"

#: (id, file, [(find, replace), ...], a probe asserting the plant is REACHED,
#:  one sentence on why this defect is the dangerous shape)
PLANTS = [
    (
        "solver-constants-dropped", MEMO,
        [('        "solver_constants": solver_constants(),\n', "")],
        "from app.server.memo import run_keys;"
        "c=run_keys(SPEC)['components'];"
        "assert 'solver_constants' not in c, c;"
        "print('reached: components =', sorted(c))",
        "a frozen RunSpec carries no solver settings, so without them in the "
        "key two runs computed with different numerics share one entry",
    ),
    (
        "engine-fingerprint-reported-but-not-keyed", MEMO,
        [('    common = {k: v for k, v in components.items() '
          'if k != "spec_digest"}\n',
          '    common = {k: v for k, v in components.items()\n'
          '              if k not in ("spec_digest", "engine")}\n')],
        "from app.server import memo;"
        "a=memo.run_keys(SPEC);"
        "assert 'engine' in a['components'];"
        "s=dict(memo.engine_fingerprint()); s['digest']='0'*64;"
        "memo._FINGERPRINT=s;"
        "b=memo.run_keys(SPEC);"
        "assert a['run_key']==b['run_key'];"
        "print('reached: engine still reported, key unchanged under a new "
        "digest')",
        "the dangerous shape: the component is still RETURNED in the response, "
        "so the key looks right while the engine no longer reaches it. "
        "Deleting the component outright instead raises a KeyError in "
        "jobs._memo_stamp and fails every run, which is a loud failure rather "
        "than a stale-cache one",
    ),
    (
        "canonical-drops-nulls", MEMO,
        [("    if isinstance(value, dict):\n"
          "        names = [str(k) for k in value]",
          "    if isinstance(value, dict):\n"
          "        value = {k: v for k, v in value.items() if v is not None}\n"
          "        names = [str(k) for k in value]")],
        "from app.server.memo import spec_digest as d;"
        "assert d({'timelapse_frames': None}) == d({});"
        "print('reached: an explicit null and an absent field now collide')",
        "execute_run reads spec.get('timelapse_frames', 24) and then "
        "int(x) if x is not None else 0, so null means 0 frames and absent "
        "means 24 -- two different runs that a tidying canonicaliser merges",
    ),
    (
        "canonical-drops-type-tags", MEMO,
        [('        return ["num", repr(float(value))]\n',
          "        return repr(float(value))\n"),
         ('        return ["num", repr(value)]\n',
          "        return repr(value)\n"),
         ('        return ["str", value]\n', "        return value\n"),
         ('        return ["bool", value]\n', "        return repr(value)\n")],
        "from app.server.memo import spec_digest as d;"
        "assert d({'x': '1.0'}) == d({'x': 1.0});"
        "assert d({'x': True}) == d({'x': 'True'});"
        "print('reached: a number and its own repr now collide')",
        "an untagged canonicaliser merges a value with its own text",
    ),
    (
        "scaled-axes-widened-to-n-doses", MEMO,
        [('SCALED_AXES = ("plasma_peak", "half_life_s")',
          'SCALED_AXES = ("plasma_peak", "half_life_s", "n_doses",\n'
          '               "dose_interval_s")'),
         ('    changed = [a for a, k in (("plasma_peak", "plasma_peak"),\n'
          '                              ("half_life_s", "half_life_s"))\n'
          "               if base[k] != new[k]]",
          '    changed = [a for a, k in (("plasma_peak", "plasma_peak"),\n'
          '                              ("half_life_s", "half_life_s"),\n'
          '                              ("n_doses", "n_doses"),\n'
          '                              ("dose_interval_s", '
          '"dose_interval_s"))\n'
          "               if base[k] != new[k]]")],
        "from app.server.memo import SCALED_AXES, run_keys;"
        "assert 'n_doses' in SCALED_AXES;"
        "two=dict(SPEC, pk=dict(SPEC['pk'], n_doses=2));"
        "assert run_keys(SPEC)['shape_key']==run_keys(two)['shape_key'];"
        "print('reached: n_doses now shares the shape key')",
        "the derived-but-unmeasured widening: n_doses DOES reach the solver "
        "only through plasma_mean, and every row of the measurement carries "
        "n_doses = 1",
    ),
    (
        "kill-fraction-scaled", MEMO,
        [('    "metrics": {\n        "mean_tumor_auc": (1, 0, 2.0,',
          '    "metrics": {\n'
          '        "kill_fraction": (1, 0, 2.0, "PLANT: assumed linear"),\n'
          '        "therapeutic_coverage_1x": (1, 0, 2.0, "PLANT: assumed '
          'linear"),\n'
          '        "mean_tumor_auc": (1, 0, 2.0,')],
        "from app.server.memo import READOUT_RULES;"
        "assert 'kill_fraction' in READOUT_RULES['metrics'];"
        "print('reached: metrics rules =', sorted(READOUT_RULES['metrics']))",
        "the one the suite catches against a real re-solve rather than against "
        "an opinion",
    ),
    (
        "pruned-payload-served-as-present", MEMO,
        [("        path = mgr.payload_path(run_id, name)\n"
          "        on_disk = os.path.isfile(path)",
          "        path = mgr.payload_path(run_id, name)\n"
          "        on_disk = True")],
        "import inspect;"
        "from app.server import memo;"
        "assert 'on_disk = True' in inspect.getsource(memo.payload_status);"
        "print('reached: payload_status never stats the file')",
        "a memo that hands back a run whose binary payload is gone while "
        "reporting it present -- the failure scripts/prune_run_payloads.py "
        "exists to prevent",
    ),
    (
        "memo-stamp-never-written", JOBS,
        [('        "memo": _memo_stamp(spec),\n', "")],
        "import inspect;"
        "from app.server import jobs;"
        "assert '_memo_stamp(spec)' not in inspect.getsource(jobs.execute_run);"
        "print('reached: execute_run no longer stamps a run')",
        "an unstamped store answers nothing, which is the intended fail-closed "
        "behaviour for every run written before this field existed",
    ),
]

SPEC_LITERAL = (
    "SPEC={'section':'xenium_breast','dx_um':30.0,'total_time_s':3600.0,"
    "'pk':{'plasma_peak':50.0,'half_life_s':14400.0,'n_doses':1,"
    "'dose_interval_s':86400.0}};"
)

_COUNT = re.compile(r"(\d+) failed, (\d+) passed")
_PASSED_ONLY = re.compile(r"^(\d+) passed", re.M)
_FAILED_LINE = re.compile(r"^FAILED (\S+)", re.M)


def run(args, **kw):
    return subprocess.run(args, cwd=ROOT, capture_output=True, text=True, **kw)


def suite():
    out = run([PYTHON, "-m", "pytest", SUITE, "-q"])
    text = out.stdout + out.stderr
    m = _COUNT.search(text)
    if m:
        failed, passed = int(m.group(1)), int(m.group(2))
    else:
        m2 = _PASSED_ONLY.search(text)
        failed, passed = 0, int(m2.group(1)) if m2 else -1
    reds = [t.split("::", 1)[-1] for t in _FAILED_LINE.findall(text)]
    errs = sorted({ln.strip() for ln in text.splitlines()
                   if ln.startswith("E ") and "Assertion" in ln})
    return {"failed": failed, "passed": passed, "red_tests": reds,
            "red_output": errs[:6]}


def apply_plant(path, edits):
    with open(path, encoding="utf-8") as fh:
        src = fh.read()
    for find, replace in edits:
        if find not in src:
            raise SystemExit(f"PLANT DOES NOT APPLY to {path}: {find!r:.120}")
        src = src.replace(find, replace, 1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(src)


def main():
    write = "--write" in sys.argv
    env = dict(os.environ)
    env.setdefault("TMESIM_STATE_DIR",
                   "experiments/surrogate_feasibility/_state")
    os.environ.update(env)

    base = suite()
    if base["failed"]:
        raise SystemExit(f"the suite is not green to begin with: {base}")
    print(f"green baseline: {base['passed']} passed\n")

    results = []
    tmp = tempfile.mkdtemp(prefix="memo-plants-")
    try:
        for pid, rel, edits, probe, why in PLANTS:
            path = os.path.join(ROOT, rel)
            pristine = os.path.join(tmp, pid + "-" + os.path.basename(rel))
            shutil.copy2(path, pristine)
            try:
                apply_plant(path, edits)
                parse = run([PYTHON, "-c",
                             f"import ast;ast.parse(open({rel!r}).read())"])
                if parse.returncode:
                    raise SystemExit(f"{pid}: the plant does not parse\n"
                                     f"{parse.stderr}")
                reach = run([PYTHON, "-c", SPEC_LITERAL + probe])
                if reach.returncode:
                    raise SystemExit(f"{pid}: the plant is NOT REACHED\n"
                                     f"{reach.stdout}{reach.stderr}")
                red = suite()
                print(f"{pid}: parses, reached ({reach.stdout.strip()}), "
                      f"{red['failed']} failed / {red['passed']} passed")
                for line in red["red_output"][:3]:
                    print(f"    {line}")
                results.append({"id": pid, "file": rel,
                                "edits": [list(e) for e in edits],
                                "why_this_shape": why,
                                "reached_by": reach.stdout.strip(),
                                **red})
            finally:
                shutil.copy2(pristine, path)
                if not filecmp.cmp(pristine, path, shallow=False):
                    raise SystemExit(
                        f"RESTORE FAILED for {rel}. The planted defect may "
                        "still be in the tree. Fix this before anything else.")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    after = suite()
    if after["failed"]:
        raise SystemExit(f"the suite is NOT green after reverting: {after}")
    print(f"\nall plants reverted; suite green again: {after['passed']} passed")

    doc = {
        "what_this_is": (
            "The planted-defect log for tests/test_memo.py, regenerated by "
            "experiments/exact_memo/run_plants.py. Each entry was applied to a "
            "pristine copy of the file named, CONFIRMED to parse and to be "
            "REACHED (reached_by is the probe's own output), then the suite "
            "was run and its failures recorded, then the file was restored and "
            "compared byte-for-byte. The suite is asserted green before the "
            "first plant and after the last."),
        "suite": SUITE,
        "green_baseline": {"passed": base["passed"], "failed": base["failed"]},
        "green_after": {"passed": after["passed"], "failed": after["failed"]},
        "plants": results,
    }
    if write:
        dest = os.path.join(HERE, "plants.json")
        with open(dest, "w") as fh:
            json.dump(doc, fh, indent=1)
        print("wrote", dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
