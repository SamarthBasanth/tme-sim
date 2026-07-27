#!/usr/bin/env bash
#
# reproduce_real.sh — one-command regeneration of every REAL result in tmesim.
#
# This runs the real-data pipeline end to end:
#
#   data  ->  reconstruct  ->  flagship simulation  ->  drug panel  ->  cross-tissue
#         ->  internal validation  ->  external validation   [ ->  deep experiments ]
#
# Nothing synthetic is touched: the synthetic cohort generator and the
# fabricated-geometry loader are deleted from the package (preflight asserts it),
# and the synthetic-era scripts live retired under archive/synthetic_era/.
#
# Everything is deterministic (seed 0) — re-running reproduces the same numbers.
#
# Usage:
#   ./scripts/reproduce_real.sh                 # core pipeline (~2h15m)
#   ./scripts/reproduce_real.sh --with-deep     # + the deep experiments (~3h45m)
#   ./scripts/reproduce_real.sh --list          # print the stages and exit
#   ./scripts/reproduce_real.sh --dry-run       # print the commands, run nothing
#   ./scripts/reproduce_real.sh --only panel    # run a single stage
#   ./scripts/reproduce_real.sh --from cross    # resume from a stage onwards
#
# Environment:
#   PY   python interpreter (default: <repo>/.venv/bin/python)
#
# See docs/REPRODUCE.md for setup, data acquisition and expected outputs, and
# docs/DATA_INVENTORY.md for the re-download record of every cached dataset.
#
# EVERY command this script RUNS has been executed in this repository. Modules
# that look runnable but are NOT (their docstrings advertise commands that are
# currently no-ops) are emitted as explicit [skip] lines with the reason — see
# stage_validate() and stage_external(). Do not "fix" them by adding them back.

set -euo pipefail

# --------------------------------------------------------------------------
# Paths / interpreter
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

PY="${PY:-${REPO_ROOT}/.venv/bin/python}"
CACHE="${REPO_ROOT}/data_cache"
HUNT="${CACHE}/hunt"
RESULTS="${REPO_ROOT}/results"

DRY_RUN=0
WITH_DEEP=0
ONLY=""
FROM=""

ALL_STAGES=(preflight data reconstruct flagship panel cross validate external)
DEEP_STAGES=(deep)

# --------------------------------------------------------------------------
# Pretty printing
# --------------------------------------------------------------------------
_t0_total=$(date +%s)

hr()   { printf '%s\n' "------------------------------------------------------------------------"; }
say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '\033[33m    WARN: %s\033[0m\n' "$*" >&2; }
die()  { printf '\n\033[31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# run <human label> <approx runtime> -- <command...>
run() {
  local label="$1"; shift
  local approx="$1"; shift
  [ "${1:-}" = "--" ] && shift
  printf '\n  \033[36m[run]\033[0m %s   (approx %s)\n' "${label}" "${approx}"
  printf '        $ %s\n' "$*"
  if [ "${DRY_RUN}" -eq 1 ]; then return 0; fi
  local t0 t1
  t0=$(date +%s)
  "$@"
  t1=$(date +%s)
  printf '        done in %ds\n' "$((t1 - t0))"
}

# skip <label> <reason> — record a step that is deliberately NOT run.
skip() {
  printf '\n  \033[33m[skip]\033[0m %s\n' "$1"
  printf '        reason: %s\n' "$2"
}

usage() {
  sed -n '3,35p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  exit 0
}

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --with-deep) WITH_DEEP=1; shift ;;
    --dry-run)   DRY_RUN=1; shift ;;
    --only)      ONLY="${2:-}"; shift 2 ;;
    --from)      FROM="${2:-}"; shift 2 ;;
    --list)
      printf 'stages: %s\n' "${ALL_STAGES[*]}"
      printf 'extra : %s (enable with --with-deep)\n' "${DEEP_STAGES[*]}"
      exit 0 ;;
    -h|--help)   usage ;;
    *)           die "unknown argument '$1' (try --help)" ;;
  esac
done

STAGES=("${ALL_STAGES[@]}")
[ "${WITH_DEEP}" -eq 1 ] && STAGES+=("${DEEP_STAGES[@]}")

stage_enabled() {
  local want="$1" s seen=0
  if [ -n "${ONLY}" ]; then
    # preflight always runs; it is cheap and catches a broken venv early.
    [ "${want}" = "preflight" ] && return 0
    [ "${want}" = "${ONLY}" ] && return 0
    return 1
  fi
  if [ -n "${FROM}" ]; then
    [ "${want}" = "preflight" ] && return 0
    for s in "${STAGES[@]}"; do
      [ "${s}" = "${FROM}" ] && seen=1
      [ "${s}" = "${want}" ] && { [ "${seen}" -eq 1 ] && return 0 || return 1; }
    done
    return 1
  fi
  for s in "${STAGES[@]}"; do [ "${s}" = "${want}" ] && return 0; done
  return 1
}

# ==========================================================================
# STAGE 0 — preflight
# ==========================================================================
stage_preflight() {
  say "STAGE 0/7  PREFLIGHT — interpreter, package, metric self-test   (approx 10 s)"

  [ -x "${PY}" ] || die "python interpreter not found at '${PY}'.
  Create it with:   python3 -m venv .venv && ./.venv/bin/pip install -e .
  Or point PY at your own interpreter:   PY=\$(which python) $0"

  info "interpreter : ${PY}"
  info "repo root   : ${REPO_ROOT}"
  info "results dir : ${RESULTS}"

  if [ "${DRY_RUN}" -eq 0 ]; then
    "${PY}" - <<'PYCHECK' || die "tmesim is not importable — run: ./.venv/bin/pip install -e ."
import sys
import tmesim
import numpy, scipy, pandas, matplotlib, sklearn, skimage, requests   # noqa: F401
print(f"    python      : {sys.version.split()[0]}")
print(f"    tmesim      : {tmesim.__version__}")

# Optional-but-needed extras. NOT in pyproject.toml — see docs/REPRODUCE.md section 1.
missing = []
for mod, why in (("rdkit", "molecular descriptors (falls back to a table without it)"),
                 ("openpyxl", "validation.gbm_tmz_exvivo reads .xlsx"),
                 ("statsmodels", "validation.moncayo_explant mixed-effects fit")):
    try:
        __import__(mod)
    except Exception:
        missing.append(f"{mod} ({why})")
if missing:
    print("    optional extras MISSING:")
    for m in missing:
        print(f"      - {m}")
else:
    print("    optional extras: rdkit, openpyxl, statsmodels all present")
PYCHECK
    mkdir -p "${RESULTS}"
  fi

  # Fail loudly if a synthetic-era module ever reappears on the import path.
  if [ -f "${REPO_ROOT}/tmesim/data/synthetic.py" ] || [ -f "${REPO_ROOT}/tmesim/data/loaders.py" ]; then
    die "a synthetic-era module is back in tmesim/data/ — this harness is REAL-DATA ONLY"
  fi
  info "synthetic modules absent (tmesim/data/synthetic.py, tmesim/data/loaders.py) — OK"

  # Cheap known-answer checks on the metric definitions. No simulation, no data.
  # This is what proves penetration_index is scale-invariant (and therefore
  # retired) while therapeutic_coverage is not.
  run "delivery-metric known-answer self-test" "1 s" -- \
      "${PY}" -m tmesim.validation.delivery_metrics --selftest
}

# ==========================================================================
# STAGE 1 — data
# ==========================================================================
# Three classes of input (full record: docs/DATA_INVENTORY.md):
#   AUTO     — fetched by `python -m tmesim.cli acquire` (TCGA, Ivy GAP zip, GDSC2)
#   MANUAL   — large open files fetched once with curl (both Xenium sections,
#              Visium x4, GSE84465, HGCC)
#   HUNT     — data_cache/hunt/*, the external-validation corpus (5.3 GB)
stage_data() {
  say "STAGE 1/7  DATA — real open datasets in data_cache/   (approx 0 s cached, ~20 s if acquiring)"

  local auto_missing=0
  [ -f "${CACHE}/gdsc_glioma_dose_response.csv" ] || auto_missing=1
  [ -f "${CACHE}/ivygap/fpkm_table.csv" ] || [ -f "${CACHE}/ivy_gap/gene_expression_matrix.zip" ] || auto_missing=1

  if [ "${auto_missing}" -eq 1 ]; then
    info "auto-acquirable data missing — downloading (GDSC2 glioma, Ivy GAP, TCGA-GBM)"
    run "acquire open datasets" "20 s" -- "${PY}" -m tmesim.cli acquire --out "${CACHE}"
  else
    info "auto-acquirable data present (GDSC2, Ivy GAP) — skipping download"
  fi

  # Ivy GAP: acquire caches the zip; the loader needs it extracted.
  if [ ! -f "${CACHE}/ivygap/fpkm_table.csv" ] && [ -f "${CACHE}/ivy_gap/gene_expression_matrix.zip" ]; then
    info "extracting Ivy GAP gene_expression_matrix.zip -> data_cache/ivygap/"
    if [ "${DRY_RUN}" -eq 0 ]; then
      mkdir -p "${CACHE}/ivygap"
      unzip -o -q "${CACHE}/ivy_gap/gene_expression_matrix.zip" -d "${CACHE}/ivygap"
    fi
  fi

  # Manual (large) open datasets — fail with an exact pointer, never fabricate.
  local missing=()
  [ -f "${CACHE}/xenium_gbm/Xenium_GBM_cells.csv.gz" ]                 || missing+=("data_cache/xenium_gbm/Xenium_GBM_cells.csv.gz")
  [ -f "${CACHE}/xenium_gbm/Xenium_GBM_cell_feature_matrix.tar.gz" ]   || missing+=("data_cache/xenium_gbm/Xenium_GBM_cell_feature_matrix.tar.gz")
  [ -f "${CACHE}/xenium_gbm/Xenium_GBM_gene_panel.json" ]              || missing+=("data_cache/xenium_gbm/Xenium_GBM_gene_panel.json")
  [ -d "${CACHE}/visium_gbm/spatial" ]                                 || missing+=("data_cache/visium_gbm/spatial/")
  [ -d "${CACHE}/visium_gbm/filtered_feature_bc_matrix" ]              || missing+=("data_cache/visium_gbm/filtered_feature_bc_matrix/")
  [ -f "${CACHE}/scref_gbm/GSE84465_GBM_All_data.csv.gz" ]             || missing+=("data_cache/scref_gbm/GSE84465_GBM_All_data.csv.gz")
  [ -f "${CACHE}/drugresponse/hgcc_drugresponse.json" ]                || missing+=("data_cache/drugresponse/hgcc_drugresponse.json")

  if [ "${#missing[@]}" -gt 0 ]; then
    printf '\n' >&2
    for m in "${missing[@]}"; do printf '    MISSING: %s\n' "${m}" >&2; done
    die "the manually-fetched open datasets above are not in the cache.
  They are open (no login, no DUA) but too large for the acquire step.
  Fetch them with the curl commands in docs/DATA_INVENTORY.md (Tier A),
  then re-run this script."
  fi
  info "Xenium 319-panel / Visium / GSE84465 / HGCC caches present — OK"

  # The 480-gene Immuno-Oncology section: optional, warn only.
  if [ -f "${CACHE}/xenium_gbm_io/Xenium_GBM_cells.csv.gz" ]; then
    info "Xenium 480-gene Immuno-Oncology section present (816,769 cells) — OK"
  else
    warn "data_cache/xenium_gbm_io/ absent — the IO-panel reconstruct step will be skipped (see docs/DATA_INVENTORY.md A2)"
  fi

  # Cross-tissue needs the three GEO Visium sections; warn (not fail) if absent.
  local n_visium=0
  for s in GBM1 GBM2 GBM3; do [ -d "${CACHE}/visium_gbm/${s}" ] && n_visium=$((n_visium + 1)); done
  if [ "${n_visium}" -lt 3 ]; then
    warn "only ${n_visium}/3 GSE194329 Visium sections (GBM1/2/3) are cached — the cross-tissue stage needs all 3"
  else
    info "4 Visium sections present (10x parent + GBM1/2/3) — OK"
  fi

  # Hunt caches drive the EXTERNAL validation stage; warn (not fail) if absent.
  local n_hunt=0
  for d in autoradiography_microdialysis maldi_msi_drug verified_0; do
    [ -d "${HUNT}/${d}" ] && n_hunt=$((n_hunt + 1))
  done
  if [ "${n_hunt}" -lt 3 ]; then
    warn "${n_hunt}/3 required data_cache/hunt/ caches present — the 'external' stage will fail. See docs/DATA_INVENTORY.md Tier B."
  else
    info "external-validation hunt caches present (regional PK, GIST MSI) — OK"
  fi
}

# ==========================================================================
# STAGE 2 — reconstruct
# ==========================================================================
stage_reconstruct() {
  say "STAGE 2/7  RECONSTRUCT — real tissue -> real whole-niche microenvironment   (approx 5 min)"
  info "PRODUCTION RESOLUTION IS dx = 12 um. Coarser grids do not resolve the"
  info "~45 um perivascular boundary layer; the CLI default of 30 um is NOT quotable."
  info "The tissue mask is 'geometric' (union of 30 um discs around measured"
  info "centroids), so the simulated domain does not move with dx."

  run "reconstruct Xenium 319-panel @ dx=12 um (40,887 cells)" "10 s" -- \
      "${PY}" -m tmesim.cli reconstruct --source xenium --dx 12 --out "${RESULTS}"

  # The 480-gene Immuno-Oncology section. Optional — 20x the cells, and it carries
  # VEGFA/HIF1A (which the 319 panel lacks) but loses PECAM1.
  if [ -f "${CACHE}/xenium_gbm_io/Xenium_GBM_cells.csv.gz" ] || [ "${DRY_RUN}" -eq 1 ]; then
    run "reconstruct Xenium 480-gene IO panel @ dx=12 um (816,769 cells)" "45 s" -- \
        "${PY}" -m tmesim.cli reconstruct --source xenium \
                --sample "${CACHE}/xenium_gbm_io" --dx 12 --out "${RESULTS}"
  else
    skip "reconstruct Xenium 480-gene IO panel" "data_cache/xenium_gbm_io/ not present"
  fi

  # NOTE: the CLI Visium path snaps dx to the ~100 um spot pitch (dx < 50 is
  # rejected in reconstruct.real_pipeline._reconstruct_visium). The cross-tissue
  # experiment re-runs the SAME Visium tissues on the matched dx = 12 um grid via
  # its own reconstruction path, and reports the paired resolution control.
  run "reconstruct Visium (10x parent, native ~100 um spot pitch)" "1 min" -- \
      "${PY}" -m tmesim.cli reconstruct --source visium --out "${RESULTS}"

  run "reconstruction fidelity (cross-platform + permutation tests)" "3 min" -- \
      "${PY}" -m tmesim.validation.reconstruction_fidelity
}

# ==========================================================================
# STAGE 3 — flagship simulation
# ==========================================================================
stage_flagship() {
  say "STAGE 3/7  FLAGSHIP — validated physics on the real microenvironment   (approx 8 min)"
  info "Real geometry + real drug physics. The kill number is a PREDICTION with no"
  info "outcome label to check it against — this section is an untreated resection."
  info "NOTE: the CLI still PRINTS penetration_index. It is RETIRED (see stage 4);"
  info "ignore that line."

  run "flagship real simulation (temozolomide, dx=12 um)" "4 min" -- \
      "${PY}" experiments/real_simulation.py

  run "CLI single-drug simulate (smoke check of the documented command)" "4 min" -- \
      "${PY}" -m tmesim.cli simulate --source xenium --dx 12 --drug temozolomide --out "${RESULTS}"
}

# ==========================================================================
# STAGE 4 — drug panel
# ==========================================================================
stage_panel() {
  say "STAGE 4/7  DRUG PANEL — 5 real drugs, one real tumor   (approx 21 min)"
  info "HEADLINE METRIC: therapeutic_coverage @ 1x IC50 (fraction of tumour cells"
  info "whose time-averaged exposure clears the drug's measured IC50). Quote it at"
  info "dx = 12 um with the 0.5x/2x band beside it and a +-0.03 error bar."
  info "Do NOT rank drugs by penetration_index: it is RETIRED. It is a"
  info "scale-invariant SHAPE statistic of the vascular field and spans just 0.0014"
  info "across a panel whose delivered exposure spans 120x."

  run "real drug panel @ dx=12 um" "10 min" -- \
      "${PY}" experiments/real_drug_panel.py

  run "defensible delivery metrics (+ grid ladder)" "11 min" -- \
      "${PY}" -m tmesim.validation.delivery_metrics
}

# ==========================================================================
# STAGE 5 — cross tissue
# ==========================================================================
stage_cross() {
  say "STAGE 5/7  CROSS-TISSUE — same drug, 5 real tumors, matched dx = 12 um   (approx 20 min)"
  info "Includes the paired dx=100 um resolution control that retracted the earlier"
  info "'kill spans 2%->24%' headline (that spread was a grid artifact)."

  run "cross-tissue at matched resolution" "20 min" -- \
      "${PY}" experiments/cross_tissue.py
}

# ==========================================================================
# STAGE 6 — internal validation
# ==========================================================================
stage_validate() {
  say "STAGE 6/7  INTERNAL VALIDATION — numerics, calibration, and the model's own geometry   (approx 85 min)"
  info "alpha is fit ONCE on doxorubicin's real perivascular anchor using the"
  info "E-FOLD convention L = d_half/ln2 (alpha ~= 148.6). The earlier 309.6 was a"
  info "convention error (fitting the raw half-distance inflates alpha by 1/ln^2(2))."

  run "calibration (single alpha, out-of-sample re-validation)" "2 min" -- \
      "${PY}" -m tmesim.validation.calibration

  run "multi-drug tissue-PK anchors (PMID-verified anchors, 6 drugs)" "3 min" -- \
      "${PY}" -m tmesim.validation.pk_multidrug

  run "pharmacology vs real patient-derived screens (GDSC2 vs HGCC)" "1 min" -- \
      "${PY}" -m tmesim.validation.pharmacology_real

  # This is the check that found the vasculature is barely distinguishable from a
  # random vessel field at the same density. It is the evidence for correction 5.
  run "vascular morphometry vs published human GBM (is the vasculature real?)" "2 min" -- \
      "${PY}" -m tmesim.validation.vascular_validation

  run "vessel profile (AUC vs distance-from-vessel + SOURCE NULL)" "6 min" -- \
      "${PY}" -m tmesim.validation.vessel_profile

  run "numerical rigor (mass conservation, dx ladder, dt ladder, elasticities)" "45 min" -- \
      "${PY}" -m tmesim.validation.numerical_rigor

  run "uncertainty quantification (80 Latin-hypercube draws)" "25 min" -- \
      "${PY}" -m tmesim.validation.uncertainty

  # Print-only modules. They regenerate NOTHING under results/ — run them for the
  # console table, not for a file.
  run "tissue-PK anchor registry (prints only, writes nothing)" "5 s" -- \
      "${PY}" -m tmesim.validation.pk_anchors

  skip "validation.pk_reality" \
       "its __main__ only PRINTS. results/pk_reality.md has no producer in the current module and is stale — do not expect it to refresh."
}

# ==========================================================================
# STAGE 7 — external validation (against other people's measurements)
# ==========================================================================
stage_external() {
  say "STAGE 7/7  EXTERNAL VALIDATION — the model vs measured drug concentrations   (approx 12 min)"
  info "These are the only checks scored against drug concentrations measured in"
  info "real human tissue by someone else. Both are allowed to fail, and one does."

  if [ ! -d "${HUNT}/autoradiography_microdialysis" ] && [ "${DRY_RUN}" -eq 0 ]; then
    die "data_cache/hunt/autoradiography_microdialysis/ is missing.
  The external stage needs the hunt caches. See docs/DATA_INVENTORY.md Tier B."
  fi

  # Validation run #1: partition magnitude PASSES (1.01x bias excluding one
  # outlier), cross-drug ordering PASSES (3.65x predicted vs 4.66x measured),
  # REGIONAL ordering FAILS (model predicts exactly 1.00x — no barrier term).
  run "human glioma regional PK (125 samples / 42 patients + 27 PMID anchors)" "5 min" -- \
      "${PY}" -m tmesim.validation.glioma_regional_pk

  # Validation run #2: transport-only is FALSIFIED as a sufficient explanation
  # (7.9x geometric-mean fold error vs 13.0x for the no-partitioning null), and
  # the failure is diagnostic — 4 of 5 liver pairs lie below the model's
  # structural floor, so no penetration length can reproduce them.
  run "GIST imatinib MSI (6 matched human tumour/normal pairs)" "5 min" -- \
      "${PY}" -m tmesim.validation.gist_imatinib

  # Right disease, right drug, measured response: GSE232173 patient-derived GSC
  # cultures under temozolomide (Ntafoulis 2023, PMID 37620410). Potency axis only
  # — no drug concentration was measured, and dissociated cultures have no
  # microenvironment.
  if [ -d "${HUNT}/geo_organoid_screens" ] || [ "${DRY_RUN}" -eq 1 ]; then
    run "GBM + temozolomide measured ex-vivo response (GSE232173)" "3 min" -- \
        "${PY}" -m tmesim.validation.gbm_tmz_exvivo
  else
    skip "validation.gbm_tmz_exvivo" "data_cache/hunt/geo_organoid_screens/ not present"
  fi

  # Moncayo HGSOC explants: the only open asset with a per-pixel measured drug
  # field AND a molecular state on the same tissue. Ovarian, bath-dosed ex vivo —
  # it validates diffusion + uptake + sequestration, never vascular delivery.
  if [ -d "${HUNT}/geomx_dsp" ] || [ "${DRY_RUN}" -eq 1 ]; then
    run "HGSOC explant MALDI-MSI (376,736 pixels + 114 GeoMx ROIs)" "15 min" -- \
        "${PY}" -m tmesim.validation.moncayo_explant
  else
    skip "validation.moncayo_explant" "data_cache/hunt/geomx_dsp/ not present"
  fi

  # ---- NOT RUN, and why. Do not "fix" these by adding them: the modules have no
  # ---- entry point and exit without writing anything. Their docstrings advertise
  # ---- commands that are currently no-ops. Check for `if __name__ == "__main__":`
  # ---- before trusting a docstring here.
  skip "validation.cellstate_validation (oxygen-banded states vs measured MKI67/TOP2A)" \
       "INCOMPLETE — run()/make_figure() exist; no write_report, no __main__."
  skip "validation.ranking_robustness  (is the top-1 drug call decision-grade?)" \
       "INCOMPLETE — only run_ranking_mc(); no writer, no __main__. Until it lands, report cisplatin/doxorubicin as a TIE."
  skip "validation.vessel_network_impact (KDE source vs discrete vessel network)" \
       "INCOMPLETE — declares vessel_network_impact.{md,png} but has no __main__ yet."
  skip "hypoxia vs measured VEGFA/HIF1A on the 480-gene IO panel" \
       "NOT WRITTEN — marker_panels points at validation/hypoxia_validation.py, which does not exist."

  # The discrete vessel network that would fix the KDE confound. Built and
  # self-tested, but DORMANT: nothing imports it, and delivery_source still uses
  # the 40 um KDE. Running the self-test documents that it works.
  run "vessel-network self-test (DORMANT module — not wired into delivery)" "2 min" -- \
      "${PY}" -m tmesim.reconstruct.vessel_network
}

# ==========================================================================
# OPTIONAL — deep experiments
# ==========================================================================
stage_deep() {
  say "OPTIONAL  DEEP EXPERIMENTS — niche, combination, dose/schedule, patient report   (approx 90 min)"
  warn "niche_response and dose_schedule both rest on vessel-relative geometry, which"
  warn "is CONFOUNDED by the KDE vascular field (README correction 5). Read their"
  warn "reports for the confound, not for the headline."

  run "response by anatomic niche region (Visium + Ivy GAP signatures)" "15 min" -- \
      "${PY}" experiments/niche_response.py

  run "combination therapy (spatial complementarity)" "15 min" -- \
      "${PY}" experiments/combination.py

  run "dose & schedule sweep" "30 min" -- \
      "${PY}" experiments/dose_schedule.py

  run "per-tumor virtual drug panel report (Xenium)" "12 min" -- \
      "${PY}" experiments/patient_report.py --source xenium

  run "per-tumor virtual drug panel report (Visium)" "12 min" -- \
      "${PY}" experiments/patient_report.py --source visium
}

# ==========================================================================
# Drive
# ==========================================================================
hr
printf 'tmesim — REAL-DATA reproducibility harness\n'
printf 'stages: %s' "${STAGES[*]}"
[ -n "${ONLY}" ] && printf '   (--only %s)' "${ONLY}"
[ -n "${FROM}" ] && printf '   (--from %s)' "${FROM}"
[ "${DRY_RUN}" -eq 1 ] && printf '   [DRY RUN]'
printf '\n'
hr

for stage in preflight data reconstruct flagship panel cross validate external deep; do
  case "${stage}" in
    deep) [ "${WITH_DEEP}" -eq 1 ] || { [ "${ONLY}" = "deep" ] || continue; } ;;
  esac
  if stage_enabled "${stage}"; then
    "stage_${stage}"
  fi
done

_t1_total=$(date +%s)
hr
printf '\nDONE in %d min. Regenerated artifacts are under results/.\n' "$(( (_t1_total - _t0_total) / 60 ))"
printf 'Headline write-ups: results/SESSION_AUDIT.md (authoritative current state),\n'
printf '                    results/real_drug_panel.md, results/delivery_metrics.md,\n'
printf '                    results/glioma_regional_pk.md, results/gist_imatinib.md,\n'
printf '                    results/numerical_rigor.md, results/vascular_validation.md.\n'
printf '\nReminders before you quote anything:\n'
printf '  * penetration_index is RETIRED. Use therapeutic_coverage @ 1x IC50.\n'
printf '  * The vascular field is a 40 um KDE, not vessels — every vessel-relative\n'
printf '    spatial claim is confounded (results/vascular_validation.md).\n'
printf '  * The regional CE-vs-non-enhancing axis FAILS; the additive-leak fix\n'
printf '    (L ~ 0.24) exists in reconstruct/barrier.py but is OFF by default.\n'
printf '  * Predicted drug RESPONSE is unvalidated and cannot be validated with\n'
printf '    open data (results/DATA_HUNT.md).\n'
hr
