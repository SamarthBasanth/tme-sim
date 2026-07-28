"""Generate the REAL runs that back results/studio_export.html.

Runs the production engine through the same entry points the Studio backend
uses (``app.server.jobs.JobManager`` -> ``execute_run``), so these are ordinary
tmesim runs that happen to have been submitted from a script rather than from
the browser. Nothing here is a mock.

Three runs, identical in every RunSpec field except ``drug``, on the flagship
Xenium 319-panel section at the production dx = 12 um:

    temozolomide | cisplatin | doxorubicin

Holding the PK fixed is what makes the comparison a comparison: the only thing
that differs between the three is the molecule, so the delta is a delta in
transport + potency, not in dosing.
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.server.api import default_spec, validate_spec  # noqa: E402
from app.server.jobs import JobManager  # noqa: E402

DRUGS = ["temozolomide", "cisplatin", "doxorubicin"]


def main() -> int:
    base = default_spec()
    base["dx_um"] = 12.0
    base["sensitivity_mode"] = "uniform"

    specs = []
    for drug in DRUGS:
        spec = json.loads(json.dumps(base))
        spec["drug"] = drug
        verdict = validate_spec(spec)
        if not verdict["ok"]:
            print(f"REFUSED {drug}: {verdict['errors']}", flush=True)
            return 2
        specs.append(spec)

    mgr = JobManager()
    ids = []
    for spec in specs:
        job = mgr.submit(spec)
        ids.append(job.run_id)
        print(f"submitted {job.run_id}  drug={spec['drug']}", flush=True)

    last = {}
    while True:
        jobs = [mgr.get(i) for i in ids]
        for j in jobs:
            key = (j.run_id, j.status, j.stage)
            if last.get(j.run_id) != key:
                last[j.run_id] = key
                print(f"  {time.strftime('%H:%M:%S')} {j.run_id} "
                      f"{j.status}/{j.stage} {j.progress:.2f}", flush=True)
        if all(j.status in ("done", "failed", "cancelled") for j in jobs):
            break
        time.sleep(2.0)

    mgr.shutdown()

    ok = True
    for j in (mgr.get(i) for i in ids):
        print(f"\n=== {j.run_id} {j.spec['drug']} -> {j.status} "
              f"({j.elapsed_s:.1f}s)", flush=True)
        if j.status != "done":
            ok = False
            print(j.error, flush=True)
            print("\n".join(j.log_tail(30)), flush=True)
            continue
        print(json.dumps(j.body.get("metrics"), indent=1), flush=True)

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "_export_run_ids.json"), "w") as fh:
        json.dump(ids, fh)
    print("\nrun ids:", ids, flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
