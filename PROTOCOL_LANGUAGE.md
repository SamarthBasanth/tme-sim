# TPL — Tumour Protocol Language (v1 spec)

A small domain-specific language for defining a drug, an administration protocol, and a
parameter sweep, then running it against a reconstructed specimen. Frozen interface: the
parser/executor (`tmesim/protocol/`) and the editor (`app/web/js/ui/protocol.js`) both
implement exactly this.

**Design stance.** MATLAB-adjacent in feel — assignment, ranges, loops over parameters, an
interactive run — but domain-specific, declarative where possible, and impossible to express
something the engine cannot honour. It is a *protocol* language, not a general programming
language: no arbitrary control flow, no user functions, no I/O. Every construct maps onto a
real quantity in the model, and anything the engine cannot back is a parse error, not a
silent default.

---

## Units are mandatory

Every physical quantity carries a unit. `dose 150` is an error; `dose 150 mg/m2` is not.
Recognised: `mg/m2`, `mg/kg`, `mg`, `uM`, `nM`, `mM`, `um`, `mm`, `s`, `min`, `h`, `d`,
`um2/s`, `%`. Conversions are explicit in the executor and reported in provenance.

---

## 1. Defining a drug

Three sources, in decreasing order of evidence. The executor records which was used.

```tpl
# (a) from a structure — compiled through the molecular layer
drug TMZ from smiles "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N"
  ic50 from gdsc                  # measured; errors if the drug is absent from the screen
end

# (b) from the built-in library
drug CIS from library "cisplatin"
end

# (c) from an uploaded file — .smi/.mol/.sdf structure, or .json/.yaml parameters
drug MYCPD from file "uploads/compound_47.sdf"
  ic50 12 uM                      # explicit override; recorded as user-supplied
  bbb_penetrance 0.31             # any DrugParams field may be overridden
end
```

Overridable fields: `ic50`, `hill`, `diffusivity`, `binding`, `clearance`,
`bbb_penetrance`, `molecular_weight`. Anything not given is derived by the molecular layer
and labelled `computed`; anything given is labelled `user-supplied`. `ic50 from gdsc` is
labelled `measured` with its n.

---

## 2. Administration protocol

```tpl
protocol stupp_adjuvant
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 5
  then rest 23 d
  cycles 6
end

protocol infusion_24h
  route iv
  dose 100 mg/m2 over 24 h        # `over` = infusion duration, not a bolus
  repeat 1
end
```

`dose ... over T` produces a sustained plasma profile; without `over` it is a bolus. The
executor converts dose → `PKConfig(plasma_peak, half_life_s, n_doses, dose_interval_s)`
and **must state the conversion and its assumptions in provenance** (body-surface area,
bioavailability, volume of distribution). Requires `PDEConfig.time_resolved = true`;
the executor sets it and says so.

---

## 3. Sweeps

```tpl
sweep dose from 50 to 400 mg/m2 step 50
sweep every from 6 to 48 h in 4 steps        # log or linear, `in N steps` = linear
sweep dx 12 um                                # a fixed override, not a sweep
```

Sweeps are a cross-product. The executor reports the run count and **refuses to launch
more than `max_runs` (default 64) without an explicit `allow runs N`** — a two-axis sweep
is trivially a hundred PDE solves.

---

## 4. Running

```tpl
run TMZ with stupp_adjuvant on xenium_319
  dx 12 um
  state_mode markers            # errors if the specimen's panel cannot support it
  sensitivity uniform
  report coverage, kill, mean_auc, penetration_length
end
```

`on <section>` accepts a built-in section id or an ingested sample id. If the specimen's
capability report says a requested mode is unsupported, that is an **error at validation
time**, before anything runs.

---

## 5. Comparison and assertions

```tpl
sweep dose from 100 to 400 mg/m2 step 100

compare TMZ with stupp_adjuvant on xenium_319   # 4 RUNS of ONE compound
  dx 12 um
  report coverage, kill
end

expect kill > 0.05                # recorded as pass/fail; never alters the result
```

`compare` takes **two or more RUNS of ONE compound** on one section. The runs come from a
`sweep` — a dose ladder, a schedule interval, a grid pitch. Three refusals, all at
**validation time, before anything runs**:

| program | verdict |
|---|---|
| `compare TMZ, CIS with ...` (two distinct compounds) | **error** — `withdrawn-ranking` |
| `compare TMZ, TMZ with ...` | **error** — names one compound twice, so every run is the same run |
| `compare TMZ with ...` and no `sweep` | **error** — plans one run, and one run is not a comparison |
| ANY program naming two compounds across SEPARATE statements — two `run` blocks, `compare` plus `run`, two `compare` blocks, or a `sweep` over either | **error** — `withdrawn-ranking`. Added 2026-09-03: the three rows above are scoped to ONE STATEMENT, and a plan is one job with one `grid`, one `rows` list and one colour scale spanning every cell, so the compounds land on a shared axis however they got there. Measured: all four rewrites planned cleanly, one of them an eight-point two-compound dose ladder. |

`rank by <metric>` is **gone from the language**. A program containing it is a
validation error naming the withdrawal, not a silent no-op.

Assertions are *reported*, not enforced — a failed `expect` marks the run, it does not
error. The point is to state a hypothesis before running it. `expect` and `report` still
take any metric in the metric table; both are single-run constructs and neither orders
anything.

---

### CHANGE — 2026-08-21: cross-compound ordering withdrawn

**This section previously read `compare TMZ, CIS with stupp_adjuvant on xenium_319` /
`rank by coverage`, and the shipped editor template called that ranking "the defensible
output". Both are withdrawn.** The old form is recorded here rather than edited away,
because programs written against it exist and their authors are owed the reason.

**What was withdrawn:** ordering two or more compounds by predicted effect —
`compare A, B ...` and the `rank by` option, in the engine, the editor, the wire format
and the result payload (the `ranking` key is gone from `GET /api/protocol/runs/{id}`).

**What was kept, and why it is a different question:** comparing two or more RUNS of the
same compound on the same tissue. *Ordering compounds by predicted effect goes; comparison
of runs stays.* One compound at two doses on one section is a question this map answers
honestly.

**Why.** The ordering was measured this week and does not hold.

* `results/pk_ranking_robustness.md` — per-compound plasma exposure spans **113.1x**
  across the compounds the panel ranked, while the published GBM order first breaks at a
  spread of only **1.53x**. The input spread is ~74x wider than the margin the output
  survives.
* `results/capacity_ranking_robustness.md` — the breast **top-1 pick changes** under this
  repo's own structural capacity rule (`two_route.binding_capacity_from_logp`), and the
  full GBM order breaks at a uniform capacity of R = 36.3, inside the measured
  brain-slice range of 3.2–4184.

**Why the fix is arity in drug-space, and NOT a narrower metric.** Someone will propose
keeping `rank by` and restricting it to a metric that survives. It has been looked at, and
there is no such metric. `penetration_length` is a *shape* quantity, so it structurally
survives both failures above — the plasma factor is an exact linear rescale of the whole
field (1.41e-14 over a 5,200x span, `results/pk_ranking_robustness.md` §1.2) and the
capacity term leaves the converged gradient shape exactly unmoved (1e-9,
`results/capacity_ranking_robustness.md` §2b). It then falls straight into the descriptor
gap: shape is set by the transport law, whose descriptor inputs have never been scored
(`results/descriptor_falsification.md`). The repo already carries the verdict on the
predecessor of that metric — `docs/architecture.md:361` says "do not rank drugs by it"
against `penetration_index`, and `docs/VERIFICATION_METHODOLOGY.md:413-417` records that
that ranking was **measured to be noise** (grid drift 13x the entire across-drug spread).
So: no safe metric to narrow to, and the restriction is on how many compounds a statement
may name instead.

**What a within-compound comparison is and is not invariant to.** A ratio taken inside one
map is invariant to the dose/plasma **scale** exactly, and to the omitted tissue capacity
only where the field has **converged** — capacity enters as a time rescale, so over a
finite unconverged window that second invariance is not exact. Neither is an invariance to
error in the **transport law**: a wrong effective diffusivity moves the *shape* of the
gradient, not just its level. The perturbations actually swept are the three cell-state
routes (`results/cell_state_impact.md`), tissue capacity uniform and per-compound over
3.2–4184 (`results/capacity_ranking_robustness.md`), and the per-compound plasma factor
(`results/pk_ranking_robustness.md`). That list is a list; an axis missing from it was not
swept, and silence is not coverage.

**Migration.** Replace `compare A, B with P on S` with one `run` statement per compound,
**each in its own program**, read each map on its own, and delete the `rank by` line. There
is no flag that restores the old behaviour.

> **"Each in its own program" was added on 2026-09-03, and it is the whole
> instruction.** Without it this sentence was the documented route around the
> prohibition it exists to serve: two `run` statements in ONE program parsed,
> planned and executed, and put both compounds in one grid on one colour scale.
> The executor now refuses that, so a program following the old wording gets
> `withdrawn-ranking` rather than a ranking. Reading each map on its own means
> running each on its own.

---

## Provenance and refusal (non-negotiable)

Every executed run emits, alongside its metrics:
- the resolved `DrugParams` with **per-field provenance**: `measured` / `computed` /
  `user-supplied` / `flagged-default`;
- the dose→PK conversion with its assumptions written out;
- the specimen capability report;
- the standing caveats — `sensitivity uniform` yields a **no-resistance upper bound**;
  absolute kill is not defensible on panels that cannot support cell state; nothing here
  is validated against a patient outcome.

The executor must **refuse rather than default**. Unknown drug, unsupported mode for the
specimen, missing unit, a sweep exceeding `max_runs`, an `ic50 from gdsc` for a drug absent
from the screen — all are errors carrying the reason and the line number.

---

## Grammar (informal)

```
program     := (drug_def | protocol_def | sweep | run_stmt | compare_stmt | expect_stmt)*
drug_def    := "drug" IDENT "from" source NEWLINE field* "end"
source      := "smiles" STRING | "library" STRING | "file" STRING
protocol_def:= "protocol" IDENT NEWLINE proto_stmt* "end"
proto_stmt  := "route" ("oral"|"iv"|"ip") | "dose" NUM UNIT ["over" NUM UNIT]
             | "every" NUM UNIT | "repeat" INT | "then" "rest" NUM UNIT | "cycles" INT
sweep       := "sweep" IDENT ("from" NUM "to" NUM UNIT ("step" NUM | "in" INT "steps") | NUM UNIT)
run_stmt    := "run" IDENT "with" IDENT "on" IDENT NEWLINE run_opt* "end"
compare_stmt:= "compare" IDENT "with" IDENT "on" IDENT NEWLINE run_opt* "end"
run_opt     := "dx" NUM UNIT | "state_mode" IDENT | "sensitivity" IDENT
             | "report" IDENT ("," IDENT)*
expect_stmt := "expect" IDENT COMPARATOR NUM [UNIT]
```

`compare_stmt` names exactly ONE compound; its two-or-more runs come from a `sweep`. The
parser still *accepts* the old comma list `compare A, B` so that the planner can answer
with the withdrawal diagnostic instead of a bare syntax error, but no program with two
distinct compounds validates. `run_opt` has no `rank by` — withdrawn 2026-08-21, §5.

Comments are `#` to end of line. Blocks close with `end`. Indentation is not significant.
