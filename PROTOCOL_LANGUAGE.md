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
compare TMZ, CIS with stupp_adjuvant on xenium_319
  rank by coverage
end

expect kill > 0.05                # recorded as pass/fail; never alters the result
```

Assertions are *reported*, not enforced — a failed `expect` marks the run, it does not
error. The point is to state a hypothesis before running it.

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
compare_stmt:= "compare" IDENT ("," IDENT)* "with" IDENT "on" IDENT NEWLINE run_opt* "end"
expect_stmt := "expect" IDENT COMPARATOR NUM [UNIT]
```

Comments are `#` to end of line. Blocks close with `end`. Indentation is not significant.
