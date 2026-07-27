#!/usr/bin/env bash
#
# reproduce_real.sh — one-command regeneration of every REAL result in tmesim.
#
# This runs the real-data pipeline end to end:
#
#   data  ->  reconstruct  ->  flagship simulation  ->  drug panel
#         ->  cross-tissue  ->  validations   [ ->  deep experiments ]
#
# Nothing synthetic is touched: the synthetic cohort generator and the
# fabricated-geometry loader are deleted from the package, and the synthetic-era
# scripts live (retired) under archive/synthetic_era/.
#
# Everything is deterministic (seed 0) — re-running reproduces the same numbers.
#
# Usage:
#   ./scripts/reproduce_real.sh                 # core pipeline (~1h45m)
#   ./scripts/reproduce_real.sh --with-deep     # + the deep experiments (~3h)
#   ./scripts/reproduce_real.sh --list          # print the stages and exit
#   ./scripts/reproduce_real.sh --dry-run       # print the commands, run nothing
#   ./scripts/reproduce_real.sh --only panel    # run a single stage
#   ./scripts/reproduce_real.sh --from cross    # resume from a stage onwards
#
# Environment:
#   PY   python interpreter (default: <repo>/.venv/bin/python)
#
# See docs/REPRODUCE.md for setup, data acquisition and expected outputs.

set -euo pipefail

# --------------------------------------------------------------------------
# Paths / interpreter
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

PY="${PY:-${REPO_ROOT}/.venv/bin/python}"
CACHE="${REPO_ROOT}/data_cache"
RESULTS="${REPO_ROOT}/results"

DRY_RUN=0
WITH_DEEP=0
ONLY=""
FROM=""

ALL_STAGES=(preflight data reconstruct flagship panel cross validate)
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

usage() {
  sed -n '3,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
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
  say "STAGE 0/6  PREFLIGHT — interpreter, package, output dir   (approx 5 s)"

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
PYCHECK
    mkdir -p "${RESULTS}"
  fi

  # Fail loudly if a synthetic-era module ever reappears on the import path.
  if [ -f "${REPO_ROOT}/tmesim/data/synthetic.py" ] || [ -f "${REPO_ROOT}/tmesim/data/loaders.py" ]; then
    die "a synthetic-era module is back in tmesim/data/ — this harness is REAL-DATA ONLY"
  fi
  info "synthetic modules absent (tmesim/data/synthetic.py, tmesim/data/loaders.py) — OK"
}

# ==========================================================================
# STAGE 1 — data
# ==========================================================================
# Two classes of input:
#   AUTO   — fetched by `python -m tmesim.cli acquire` (TCGA, Ivy GAP zip, GDSC)
#   MANUAL — large open files fetched once with curl (Xenium, Visium, GSE84465,
#            HGCC). URLs and byte sizes: docs/REPRODUCE.md + docs/real_data/.
stage_data() {
  say "STAGE 1/6  DATA — real open datasets in data_cache/   (approx 0 s cached, ~2 min if acquiring)"

  local auto_missing=0
  [ -f "${CACHE}/gdsc_glioma_dose_response.csv" ] || auto_missing=1
  [ -f "${CACHE}/ivygap/fpkm_table.csv" ] || [ -f "${CACHE}/ivy_gap/gene_expression_matrix.zip" ] || auto_missing=1

  if [ "${auto_missing}" -eq 1 ]; then
    info "auto-acquirable data missing — downloading (GDSC glioma, Ivy GAP, TCGA-GBM)"
    run "acquire open datasets" "2 min" -- "${PY}" -m tmesim.cli acquire --out "${CACHE}"
  else
    info "auto-acquirable data present (GDSC, Ivy GAP) — skipping download"
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
  Fetch them with the curl commands in docs/REPRODUCE.md section 'Manual downloads',
  then re-run this script."
  fi
  info "Xenium / Visium / GSE84465 / HGCC caches present — OK"

  # Cross-tissue needs the three GEO Visium sections; warn (not fail) if absent.
  local n_visium=0
  for s in GBM1 GBM2 GBM3; do [ -d "${CACHE}/visium_gbm/${s}" ] && n_visium=$((n_visium + 1)); done
  if [ "${n_visium}" -lt 3 ]; then
    warn "only ${n_visium}/3 GSE194329 Visium sections (GBM1/2/3) are cached — the cross-tissue stage needs all 3"
  else
    info "4 Visium sections present (10x parent + GBM1/2/3) — OK"
  fi
}

# ==========================================================================
# STAGE 2 — reconstruct
# ==========================================================================
stage_reconstruct() {
  say "STAGE 2/6  RECONSTRUCT — real tissue -> real whole-niche microenvironment   (approx 4 min)"
  info "40,887-cell 10x Xenium GBM section at the PRODUCTION resolution dx = 12 um."
  info "Coarser grids do NOT resolve the ~45 um perivascular boundary layer."

  run "reconstruct Xenium @ dx=12 um" "20 s" -- \
      "${PY}" -m tmesim.cli reconstruct --source xenium --dx 12 --out "${RESULTS}"

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
  say "STAGE 3/6  FLAGSHIP — validated physics on the real microenvironment   (approx 5 min)"
  info "Real geometry + real drug physics. The kill number is a PREDICTION with no"
  info "outcome label to check it against — this section is an untreated resection."

  run "flagship real simulation (temozolomide, dx=12 um)" "4 min" -- \
      "${PY}" experiments/real_simulation.py

  run "CLI single-drug simulate (smoke check of the documented command)" "2 min" -- \
      "${PY}" -m tmesim.cli simulate --source xenium --dx 12 --drug temozolomide --out "${RESULTS}"
}

# ==========================================================================
# STAGE 4 — drug panel
# ==========================================================================
stage_panel() {
  say "STAGE 4/6  DRUG PANEL — 5 real drugs, one real tumor   (approx 12 min)"
  info "Compare drugs by mean tissue AUC / therapeutic (IC50) coverage."
  info "Do NOT rank drugs by penetration_index: it is a scale-invariant SHAPE"
  info "statistic of the vascular field and collapses to ~0.195 across the panel."

  run "real drug panel @ dx=12 um" "10 min" -- \
      "${PY}" experiments/real_drug_panel.py

  run "defensible delivery metrics (+ grid ladder)" "11 min" -- \
      "${PY}" -m tmesim.validation.delivery_metrics
}

# ==========================================================================
# STAGE 5 — cross tissue
# ==========================================================================
stage_cross() {
  say "STAGE 5/6  CROSS-TISSUE — same drug, 5 real tumors, matched dx = 12 um   (approx 20 min)"
  info "Includes the paired dx=100 um resolution control that retracted the earlier"
  info "'kill spans 2%->24%' headline (that spread was a grid artifact)."

  run "cross-tissue at matched resolution" "20 min" -- \
      "${PY}" experiments/cross_tissue.py
}

# ==========================================================================
# STAGE 6 — validations
# ==========================================================================
stage_validate() {
  say "STAGE 6/6  VALIDATION — physics vs real measured PK, and honest error bars   (approx 70 min)"
  info "alpha is fit ONCE on doxorubicin's real perivascular anchor using the"
  info "E-FOLD convention L = d_half/ln2 (alpha ~= 148.6). The earlier 309.6 was a"
  info "convention error (fitting to the raw half-distance inflates alpha by 1/ln^2(2))."

  run "calibration (single alpha, out-of-sample re-validation)" "2 min" -- \
      "${PY}" -m tmesim.validation.calibration

  run "multi-drug tissue-PK anchors (7 PMID-verified anchors, 6 drugs)" "3 min" -- \
      "${PY}" -m tmesim.validation.pk_multidrug

  run "PK reality check (predicted vs published tissue penetration)" "3 min" -- \
      "${PY}" -m tmesim.validation.pk_reality

  run "pharmacology vs real patient-derived screens (GDSC vs HGCC)" "1 min" -- \
      "${PY}" -m tmesim.validation.pharmacology_real

  run "vessel profile (AUC vs distance-from-vessel + source null)" "6 min" -- \
      "${PY}" -m tmesim.validation.vessel_profile

  run "numerical rigor (mass conservation, dx ladder, dt ladder, elasticities)" "45 min" -- \
      "${PY}" -m tmesim.validation.numerical_rigor

  run "uncertainty quantification (80 Latin-hypercube draws @ dx=15 um)" "25 min" -- \
      "${PY}" -m tmesim.validation.uncertainty
}

# ==========================================================================
# OPTIONAL — deep experiments
# ==========================================================================
stage_deep() {
  say "OPTIONAL  DEEP EXPERIMENTS — niche, combination, dose/schedule, patient report   (approx 90 min)"

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

for stage in preflight data reconstruct flagship panel cross validate deep; do
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
printf 'Headline write-ups: results/REAL_VALIDATION.md, results/real_drug_panel.md,\n'
printf '                    results/cross_tissue.md, results/numerical_rigor.md.\n'
printf '\nReminder — what is NOT validated: predicted drug RESPONSE against real\n'
printf 'patient outcomes. No open dataset pairs spatial GBM tissue with measured\n'
printf 'response. See results/REAL_VALIDATION.md section 6.\n'
hr
