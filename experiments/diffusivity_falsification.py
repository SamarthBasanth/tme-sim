"""FALSIFICATION TEST -- score the engine's compiled tissue diffusivity against
the 71-compound Pruijn 2008 multicellular-layer dataset.

Nothing here is registered in ``tmesim/validation/scorecard.py`` and nothing in
``tmesim/`` is modified or fitted. This is a one-shot falsification test of a
structural claim, run through the unmodified molecular pipeline.

THE CLAIM UNDER TEST
--------------------
``MolecularProfile.effective_diffusivity`` (the ``D_eff`` that feeds the
transport PDE) is computed by this chain:

    mw -> descriptors.hydrodynamic_radius_nm  = 0.066 * mw**(1/3)     [MW ONLY]
       -> stokes_einstein.free_diffusivity(radius, T)                 [radius ONLY]
       -> brownian.simulate_effective_diffusivity(free_D, ecm_frac)   [radius + ECM]
       -> ecm_mesh.ogston_renkin_hindrance(radius)                    [radius ONLY]

Every input is either molecular weight or the tissue ECM volume fraction.
``logP``, hydrogen-bond donors/acceptors, TPSA and charge are computed by
``descriptors`` but NEVER reach ``effective_diffusivity``. Verified empirically
in ``_verify_mw_only_path()`` below: compiling tirapazamine from its SMILES and
compiling a bare ``Molecule`` carrying only ``mw=178.151`` give the same
``effective_diffusivity`` to within 0.01%.

Therefore, at fixed ECM fraction, the engine MUST predict a near-constant
``D_eff`` across any set of compounds of similar molecular weight. A measured
dataset that spans a wide range at fixed MW falsifies the model outright.

THE MEASUREMENT
---------------
Pruijn FB, Patel K, Hay MP, Wilson WR, Hicks KO. "Prediction of Tumour Tissue
Diffusion Coefficients of Hypoxia-Activated Prodrugs from Physicochemical
Parameters." Aust J Chem 2008;61(9):687-693. doi:10.1071/CH08240

  * 67 benzotriazine di-N-oxides (tirapazamine analogues) + 4 miscellaneous
    compounds = 71, all measured in HT29 human colon carcinoma multicellular
    layers, ONE method, ONE cell line.
  * ``Dmcl`` obtained by fitting flux between two stirred compartments across
    the MCL to Fick's second law -- i.e. a tissue-scale effective diffusion
    coefficient in exactly the engine's functional role.

PROVENANCE OF THE NUMBERS BELOW -- READ THIS
--------------------------------------------
The publisher (CSIRO Publishing / connectsci.au) returns HTTP 403 to automated
retrieval. That 403 is a Cloudflare bot challenge ("Enable JavaScript and
cookies to continue"), NOT a paywall: Unpaywall and OpenAlex both record this
DOI as bronze open access with the publisher PDF as the OA location.

Both the article PDF and its Accessory Publication (Table S1) were retrieved
from the Internet Archive Wayback Machine, which holds HTTP-200 snapshots of
the publisher's own OA PDF URLs:

  article : http://web.archive.org/web/20240530035922id_/https://www.publish.csiro.au/CH/pdf/CH08240
  Table S1: http://web.archive.org/web/2015/https://www.publish.csiro.au/?act=view_file&file_id=CH08240_AC.pdf

These are the publisher's own files, archived. The per-compound values below
were extracted from the Table S1 PDF by pdfplumber word-coordinate column
assignment (not by text-flow parsing, which mis-assigns the wrapped cells) and
spot-checked against an independent ``pdftotext -layout`` rendering.

Every ``d_meas`` value below is READ FROM Table S1. None is derived, inferred,
imputed, or predicted. The engine values are computed at run time.

Table S1 header, verbatim:
  "Table S1. Structures, physicochemical parameters, calculated Dmcl (from Eq. 1
   using the coefficients for HT29 in Study A as listed in Table 1) and measured
   Dmcl for HT29 MCLs. Errors are SEM and units of Dmcl are cm2s-1."

Run:  ./.venv/bin/python experiments/diffusivity_falsification.py
"""
from __future__ import annotations

import math
from typing import List, Tuple

from tmesim.molecular import compile_drug, descriptors
from tmesim.molecular.types import Molecule

# Table S1 reports Dmcl as (value x 10^6) in cm^2/s. 1 cm^2/s = 1e8 um^2/s,
# so a tabulated value v means v * 1e-6 cm^2/s = v * 100 um^2/s.
TABLE_UNITS_TO_UM2_S = 100.0

#: Pruijn 2008 Table S1, all 71 compounds, ordered by molecular weight.
#: (name, MW_daltons, HD, HA, logP7.4, logP_source, Dmcl_x1e6_cm2_s,
#:  Dmcl_SEM_x1e6_cm2_s, Dmcl_calculated_by_the_paper_x1e6_cm2_s)
#:
#: logP_source: "meas" = measured logP7.4 column, "calc" = calculated logP7.4
#: column (used only where the measured column is empty / "np" = not published).
#: The final field is the paper's OWN Eqn-1 prediction; it is carried only so
#: the authors' published QSAR can be shown alongside the engine's, and is
#: never mixed with the measured column.
PRUIJN_2008_TABLE_S1: List[Tuple] = [
    ("Urea", 60.1, 4, 3, -2.01, "meas", 0.450, 0.023, 0.498),
    ("SN26434", 162.0, 2, 5, 1.06, "meas", 3.04, 0.1, 1.8),
    ("SN28287", 177.0, 0, 5, -0.22, "meas", 1.92, 0.28, 2.69),
    ("SN29047", 177.0, 0, 5, 1.09, "calc", 3.8, 0.24, 3.39),
    ("SN24280 (TPZ)", 178.0, 2, 6, -0.38, "meas", 0.404, 0.013, 0.417),
    ("Mannitol", 182.0, 6, 6, -2.62, "meas", 0.211, 0.0002, 0.245),
    ("SN28266", 191.0, 0, 5, 0.22, "meas", 2.61, 0.011, 2.86),
    ("SN27758", 192.0, 2, 6, 0.127, "meas", 0.662, 0.022, 0.579),
    ("SN27759", 192.0, 2, 6, 0.068, "meas", 0.543, 0.025, 0.553),
    ("SN27769", 192.0, 2, 6, 0.241, "meas", 0.75, 0.018, 0.632),
    ("acridan-9(10H)-one", 195.0, 1, 2, 2.96, "meas", 3.17, 0.007, 3.29),
    ("SN28293", 207.0, 1, 6, -0.91, "calc", 0.441, 0.054, 0.584),
    ("SN27791", 213.0, 2, 6, 0.441, "meas", 0.87, 0.022, 0.699),
    ("SN27800", 213.0, 2, 6, 0.166, "meas", 0.669, 0.024, 0.56),
    ("SN27891", 213.0, 2, 6, -0.1, "meas", 0.354, 0.006, 0.46),
    ("SN29053", 219.0, 2, 6, -0.83, "calc", 0.519, 0.015, 0.307),
    ("SN28246", 221.0, 1, 6, -0.47, "meas", 0.494, 0.018, 0.807),
    ("SN28400", 221.0, 0, 6, -0.187, "meas", 2.09, 0.036, 2.09),
    ("SN28205", 231.0, 1, 7, -0.53, "calc", 0.503, 0.049, 0.409),
    ("SN29051", 233.0, 1, 6, 0.39, "calc", 1.83, 0.31, 1.52),
    ("SN28614", 235.0, 0, 6, -0.02, "calc", 2.13, 0.174, 2.143),
    ("SN29039", 235.0, 0, 6, 0.2, "calc", 1.74, 0.114, 2.3),
    ("SN28875", 236.0, 2, 7, None, "calc", 0.248, 0.004, 0.363),
    ("SN27849", 246.0, 2, 6, 0.815, "meas", 0.793, 0.054, 0.876),
    ("SN27916", 246.0, 2, 6, 0.479, "meas", 0.515, 0.012, 0.659),
    ("SN28019", 249.0, 2, 7, 0.83, "meas", 0.453, 0.018, 0.654),
    ("SN28436", 249.0, 3, 7, -1.48, "meas", 0.325, 0.021, 0.207),
    ("SN28870", 250.0, 1, 7, 0.356, "meas", 1.32, 0.108, 1.11),
    ("SN28098", 254.0, 1, 6, 1.43, "meas", 2.7, 0.058, 2.265),
    ("SN27820", 256.0, 2, 8, -1.15, "meas", 0.213, 0.031, 0.219),
    ("SN27866", 256.0, 2, 8, -1.38, "meas", 0.319, 0.034, 0.213),
    ("SN28402", 263.0, 1, 7, -2.01, "meas", 0.187, 0.012, 0.24),
    ("SN28427", 263.0, 1, 7, -0.74, "meas", 0.509, 0.028, 0.438),
    ("SN28433", 263.0, 2, 7, -2.62, "meas", 0.139, 0.001, 0.1996),
    ("SN28682", 263.0, 1, 7, -0.71, "meas", 0.399, 0.013, 0.448),
    ("SN28823", 279.0, 1, 8, -0.854, "meas", 0.232, 0.008, 0.315),
    ("SN28648", 288.0, 0, 6, -1.27, "meas", 0.65, 0.01, 0.887),
    ("SN28620", 289.0, 1, 7, -0.471, "meas", 0.723, 0.045, 0.509),
    ("SN29143", 289.0, 1, 7, -0.071, "meas", 0.59, 0.033, 0.712),
    ("SN29147", 289.0, 1, 7, 0.179, "meas", 0.835, 0.046, 0.878),
    ("SN28627", 290.0, 0, 7, -1.86, "calc", 0.747, 0.064, 0.41),
    ("SN28606", 291.0, 1, 8, -0.35, "meas", 0.363, 0.004, 0.424),
    ("SN29650", 291.0, 1, 8, -1.05, "meas", 0.427, 0.007, 0.279),
    ("SN28864", 303.0, 1, 7, 0.001, "meas", 0.749, 0.071, 0.735),
    ("SN28643", 305.0, 1, 7, 0.488, "meas", 1.71, 0.226, 1.082),
    ("SN29102", 305.0, 1, 7, 0.655, "meas", 0.863, 0.035, 1.22),
    ("SN29267", 305.0, 1, 8, 0.072, "meas", 0.452, 0.006, 0.577),
    ("SN29662", 305.0, 1, 8, -0.54, "meas", 0.454, 0.061, 0.36),
    ("SN28361", 312.0, 1, 7, 1.51, "meas", 2.4, 0.001, 1.84),
    ("SN28852", 317.0, 1, 7, -0.21, "calc", 0.681, 0.051, 0.597),
    ("SN28425", 319.0, 1, 8, 0.223, "meas", 0.77, 0.087, 0.638),
    ("SN29270", 319.0, 1, 7, 0.95, "calc", 1.48, 0.07, 1.42),
    ("SN29390", 324.0, 1, 7, 0.44, "calc", 1.1, 0.08, 1.01),
    ("SN29146", 329.0, 1, 7, 0.619, "meas", 1.19, 0.05, 1.13),
    ("SN28890", 331.0, 1, 7, 0.252, "meas", 0.945, 0.059, 0.856),
    ("NLCQ-1", 332.0, 1, 7, 0.22, "meas", 0.606, 0.003, 0.834),
    ("SN29347", 345.0, 1, 7, 1.23, "calc", 1.5, 0.07, 1.55),
    ("SN29375", 345.0, 1, 7, 1.61, "calc", 1.93, 0.14, 1.79),
    ("SN29366", 348.0, 1, 7, 1.82, "calc", 2.17, 0.17, 1.47),
    ("SN29273", 361.0, 1, 7, 2.63, "calc", 1.75, 0.11, 2.09),
    ("SN26955", 455.0, 2, 8, 1.69, "meas", 0.617, 0.125, 0.691),
    ("SN28565", 462.0, 2, 10, 0.148, "meas", 0.128, 0.041, 0.1693),
    ("SN27140", 470.0, 2, 10, 1.84, "meas", 0.541, 0.008, 0.426),
    ("SN28660", 470.0, 3, 10, 1.36, "meas", 0.192, 0.016, 0.172),
    ("SN28860", 484.0, 2, 10, 2.15, "meas", 0.399, 0.034, 0.545),
    ("SN28713", 497.0, 1, 9, 1.88, "calc", 1.47, 0.006, 1.21),
    ("SN27149", 498.0, 3, 10, 0.88, "meas", 0.152, 0.0074, 0.15),
    ("SN28874", 498.0, 2, 10, None, "calc", 0.734, 0.183, 0.776),
    ("SN28141", 512.0, 2, 10, 1.53, "meas", 0.382, 0.004, 0.315),
    ("SN28588", 539.0, 2, 11, 0.665, "meas", 0.976, 0.046, 0.161),
    ("SN28871", 553.0, 2, 11, None, "calc", 0.225, 0.035, 0.188),
]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def _spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return _pearson(rank(xs), rank(ys))


def _ols(xs, ys):
    """Least-squares slope and intercept of ys on xs."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    return slope, my - slope * mx


def _verify_mw_only_path() -> None:
    """Prove D_eff is a pure function of MW before trusting the MW-only compile."""
    tpz_smiles = "C1=CC=C2C(=C1)[N+](=NC(=[N+]2[O-])N)[O-]"
    from_smiles = compile_drug.build_profile(
        descriptors.from_smiles("tirapazamine", tpz_smiles)
    )
    mw_only = compile_drug.build_profile(
        Molecule(name="tpz_mw_only", smiles="", descriptors={"mw": 178.151})
    )
    a, b = from_smiles.effective_diffusivity, mw_only.effective_diffusivity
    rel = abs(a - b) / a
    print("PRE-FLIGHT: is D_eff a pure function of MW?")
    print(f"  tirapazamine via descriptors.from_smiles(SMILES) : D_eff = {a:.4f} um^2/s")
    print(f"  tirapazamine via Molecule(mw=178.151) only       : D_eff = {b:.4f} um^2/s")
    print(f"  relative difference                              : {rel:.2%}")
    print(f"  -> logP / HBD / HBA / TPSA / charge are computed but DO NOT reach D_eff."
          f" {'CONFIRMED' if rel < 1e-3 else 'NOT CONFIRMED'}\n")


def compile_engine_deff(mw: float) -> float:
    """Compile one molecular weight through the UNMODIFIED pipeline -> D_eff."""
    return compile_drug.build_profile(
        Molecule(name=f"mw{mw}", smiles="", descriptors={"mw": float(mw)})
    ).effective_diffusivity


def fold(vals) -> float:
    return max(vals) / min(vals)


# --------------------------------------------------------------------------- #
def main() -> None:
    print("=" * 78)
    print("DIFFUSIVITY FALSIFICATION TEST")
    print("Engine D_eff vs Pruijn 2008 (doi:10.1071/CH08240) Table S1, n=71,")
    print("HT29 multicellular layers, one method, one cell line.")
    print("=" * 78 + "\n")

    _verify_mw_only_path()

    rows = []
    print(f"Compiling {len(PRUIJN_2008_TABLE_S1)} compounds through "
          f"descriptors -> compile_drug.build_profile (ECM reference 0.4)...")
    for (name, mw, hd, ha, logp, src, d_meas, d_sem, d_paper) in PRUIJN_2008_TABLE_S1:
        rows.append({
            "name": name, "mw": mw, "hd": hd, "ha": ha, "logp": logp, "src": src,
            "meas": d_meas * TABLE_UNITS_TO_UM2_S,
            "sem": d_sem * TABLE_UNITS_TO_UM2_S,
            "paper": d_paper * TABLE_UNITS_TO_UM2_S,
            "eng": compile_engine_deff(mw),
        })
    print("done.\n")

    meas = [r["meas"] for r in rows]
    eng = [r["eng"] for r in rows]
    paper = [r["paper"] for r in rows]

    # ---------------- per-compound table ----------------
    print("-" * 78)
    print("PER-COMPOUND: measured vs engine (um^2/s)")
    print("-" * 78)
    print(f"{'compound':<20}{'MW':>7}{'HD':>4}{'HA':>4}{'logP':>7}"
          f"{'measured':>11}{'engine':>10}{'eng/meas':>10}")
    for r in sorted(rows, key=lambda r: r["meas"]):
        lp = "  n/a " if r["logp"] is None else f"{r['logp']:>6.2f}"
        print(f"{r['name']:<20}{r['mw']:>7.1f}{r['hd']:>4}{r['ha']:>4}{lp}"
              f"{r['meas']:>11.1f}{r['eng']:>10.1f}{r['eng']/r['meas']:>10.2f}")

    # ---------------- headline: spread ----------------
    print("\n" + "=" * 78)
    print("RESULT 1 -- THE SPREAD (the decisive statistic)")
    print("=" * 78)
    print(f"  measured Dmcl : {min(meas):8.1f} -> {max(meas):8.1f} um^2/s"
          f"   = {fold(meas):6.1f}-fold")
    print(f"  engine  D_eff : {min(eng):8.1f} -> {max(eng):8.1f} um^2/s"
          f"   = {fold(eng):6.1f}-fold")
    print(f"  paper's Eqn 1 : {min(paper):8.1f} -> {max(paper):8.1f} um^2/s"
          f"   = {fold(paper):6.1f}-fold")
    print(f"\n  The engine reproduces {fold(eng)/fold(meas):.1%} of the measured spread.")
    print(f"  MW range spanned: {min(r['mw'] for r in rows):.1f} - "
          f"{max(r['mw'] for r in rows):.1f} Da")

    # ---------------- correlation ----------------
    lm = [math.log10(v) for v in meas]
    le = [math.log10(v) for v in eng]
    lp_ = [math.log10(v) for v in paper]
    print("\n" + "=" * 78)
    print("RESULT 2 -- CORRELATION WITH MEASUREMENT")
    print("=" * 78)
    print(f"  engine   vs measured : Pearson r(log10) = {_pearson(le, lm):+.3f}"
          f"   R^2 = {_pearson(le, lm)**2:.3f}   Spearman = {_spearman(eng, meas):+.3f}")
    print(f"  paper Eqn1 vs measured: Pearson r(log10) = {_pearson(lp_, lm):+.3f}"
          f"   R^2 = {_pearson(lp_, lm)**2:.3f}   Spearman = {_spearman(paper, meas):+.3f}")
    print("  (the paper reports R^2 = 0.907, Q^2 = 0.879 for its own fit)")

    # ---------------- MW exponent ----------------
    lmw = [math.log10(r["mw"]) for r in rows]
    b_meas, _ = _ols(lmw, lm)
    b_eng, _ = _ols(lmw, le)
    print("\n" + "=" * 78)
    print("RESULT 3 -- THE MOLECULAR-WEIGHT EXPONENT")
    print("=" * 78)
    print(f"  d log10(Dmcl) / d log10(MW), MEASURED (this dataset) : {b_meas:+.3f}")
    print(f"  d log10(D_eff)/ d log10(MW), ENGINE                  : {b_eng:+.3f}")
    print( "  Pruijn 2008 Table 1, coefficient b (fitted, n=69)    :  -0.62 +/- 0.14")
    print(f"  variance in measured log10(Dmcl) explained by MW alone: "
          f"{_pearson(lmw, lm)**2:.3f}")

    # ---------------- the kill shot: narrow MW bands ----------------
    print("\n" + "=" * 78)
    print("RESULT 4 -- SPREAD WITHIN NARROW MOLECULAR-WEIGHT BANDS")
    print("=" * 78)
    print("A MW-only model cannot produce ANY spread among compounds of equal MW.")
    print("Bands are chosen to hold MW nearly constant; +-5% or exact ties.\n")
    print(f"{'MW band':<14}{'n':>3}{'MW spread':>11}{'measured range':>22}"
          f"{'meas fold':>11}{'engine fold':>13}")
    bands = [(160, 200), (200, 225), (240, 265), (285, 310), (310, 350), (450, 500)]
    for lo, hi in bands:
        sub = [r for r in rows if lo <= r["mw"] < hi]
        if len(sub) < 3:
            continue
        m_ = [r["meas"] for r in sub]
        e_ = [r["eng"] for r in sub]
        mws = [r["mw"] for r in sub]
        print(f"{f'{lo}-{hi}':<14}{len(sub):>3}{f'{min(mws):.0f}-{max(mws):.0f}':>11}"
              f"{f'{min(m_):.1f} - {max(m_):.1f}':>22}{fold(m_):>10.1f}x"
              f"{fold(e_):>12.2f}x")

    print("\n  EXACT-MW TIES (identical molecular weight -> engine is exactly flat):")
    by_mw = {}
    for r in rows:
        by_mw.setdefault(r["mw"], []).append(r)
    for mw, sub in sorted(by_mw.items()):
        if len(sub) < 3:
            continue
        m_ = [r["meas"] for r in sub]
        e_ = [r["eng"] for r in sub]
        names = ", ".join(r["name"] for r in sub)
        print(f"    MW {mw:.0f} Da, n={len(sub)}: measured {min(m_):.1f} - {max(m_):.1f} "
              f"um^2/s ({fold(m_):.1f}x), engine {fold(e_):.2f}x")
        print(f"        {names}")

    # ---------------- what the spread tracks instead ----------------
    have = [r for r in rows if r["logp"] is not None]
    llp = [r["logp"] for r in have]
    lmh = [math.log10(r["meas"]) for r in have]
    print("\n" + "=" * 78)
    print("RESULT 5 -- WHAT THE MEASURED SPREAD ACTUALLY TRACKS")
    print("=" * 78)
    print(f"  log10(Dmcl) vs logP7.4      (n={len(have)}): "
          f"Pearson r = {_pearson(llp, lmh):+.3f}  R^2 = {_pearson(llp, lmh)**2:.3f}")
    print(f"  log10(Dmcl) vs log10(MW)    (n={len(rows)}): "
          f"Pearson r = {_pearson(lmw, lm):+.3f}  R^2 = {_pearson(lmw, lm)**2:.3f}")
    hd_ = [float(r["hd"]) for r in rows]
    print(f"  log10(Dmcl) vs H-bond donors(n={len(rows)}): "
          f"Pearson r = {_pearson(hd_, lm):+.3f}  R^2 = {_pearson(hd_, lm)**2:.3f}")
    print("\n  logP7.4 and HBD are computed by descriptors.compute_descriptors but")
    print("  neither reaches effective_diffusivity. That is the structural gap.")

    # ---------------- accuracy ----------------
    ratios = sorted(r["eng"] / r["meas"] for r in rows)
    n = len(ratios)
    med = ratios[n // 2] if n % 2 else (ratios[n // 2 - 1] + ratios[n // 2]) / 2
    within3 = sum(1 for x in ratios if 1 / 3 <= x <= 3)
    print("\n" + "=" * 78)
    print("RESULT 6 -- ABSOLUTE ACCURACY (secondary to the spread result)")
    print("=" * 78)
    print(f"  engine/measured ratio: min {min(ratios):.2f}x, median {med:.2f}x, "
          f"max {max(ratios):.2f}x")
    print(f"  within 3x of measurement: {within3}/{n} = {within3/n:.0%}")
    print(f"  geometric mean bias    : "
          f"{10**(sum(math.log10(x) for x in ratios)/n):.2f}x")

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  Measured spread {fold(meas):.1f}x across 71 compounds, one method, one")
    print(f"  cell line. Engine spread {fold(eng):.1f}x, driven entirely by MW.")
    print(f"  Within compounds of IDENTICAL molecular weight the measurement still")
    print(f"  spans several fold and the engine is flat by construction.")
    print("  The MW-only diffusivity model is FALSIFIED as a predictor of")
    print("  compound-to-compound variation in tissue diffusivity.")


if __name__ == "__main__":
    main()
