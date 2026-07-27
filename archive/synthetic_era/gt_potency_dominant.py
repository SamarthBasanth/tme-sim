"""
Ground-truth motif #8: POTENCY_DOMINANT (honesty / null case).

Biological assumption
---------------------
For this patient cohort and drug panel, intrinsic pharmacodynamic POTENCY is the
overwhelming determinant of tumour-cell kill, and drug delivery to the tumour is
only a WEAK modulator. Biologically this is the regime where the drugs are so
different in on-target potency (IC50 spanning ~3 orders of magnitude, from
doxorubicin at 0.5 to temozolomide at 606) that even large differences in
delivery cannot reorder them: a highly potent agent that reaches the tumour at a
somewhat reduced concentration still out-kills a weak agent delivered perfectly.
This corresponds to a well-perfused, thin-stroma tumour where *all* drugs arrive
at broadly therapeutic fractions, so the delivery term saturates near 1 and acts
only as a small (~0.8-1.0) multiplier.

    kill[p,d] = sens[p] * Hill(dose, ic50[d]) * delivery(p,d) * noise

1. Hill(dose, ic50[d]) = dose / (dose + ic50[d])   (Hill coeff h = 1)
   The dominant term. dose is fixed (30). Because IC50 spans 3 orders of
   magnitude, the Hill values are widely separated (0.98, 0.71, 0.55, 0.38,
   0.05), so the drug ranking within every patient is set almost entirely here.

2. delivery(p,d) in [0.8, 1.0]  --  MINOR multiplier.
   Built from real transport biology (perfusion from vasc, BBB permeability
   drug_bbb, ECM sieving of large MW through dense ECM) but COMPRESSED to a
   narrow band: delivery = 0.8 + 0.2 * raw_norm, where raw_norm is the min-max
   normalised transport score. The full swing (0.8 -> 1.0) is only a 25%
   relative effect -- deliberately too small to overturn the Hill ordering. The
   most potent drug (doxorubicin) is also the worst delivered (bbb=0.05, largest
   MW), yet even at delivery=0.8 its kill (0.98*0.8=0.79) exceeds the next drug's
   best case (methotrexate 0.71*1.0=0.71). So the globally-most-potent drug is
   almost always the per-patient best.

3. sens[p]  --  patient-level, drug-agnostic chemosensitivity. Scales overall
   magnitude but does NOT change the per-patient drug ranking.

Expected outcome (the point of this null motif): the pure POTENCY shortcut should
be near-oracle, and the mechanistic transport sim should NOT beat it -- physics
buys nothing when potency dominates. We report that honestly.

Everything is deterministic given seed=0.
"""
import numpy as np
import selection_lib as sl


def build_gt_kill(cache):
    vasc = cache["vasc"]            # (P,)
    ecm = cache["ecm"]             # (P,)
    sens = cache["sens"]            # (P,)
    mw = cache["drug_mw"]        # (D,)
    ic50 = cache["drug_ic50"]      # (D,)
    bbb = cache["drug_bbb"]       # (D,)
    dose = float(cache["dose"])

    # ---- dominant potency term (Hill, h=1) ----
    hill = dose / (dose + ic50)                          # (D,)

    # ---- raw transport score from real biology ----
    perfusion = np.clip(vasc / 0.45, 0.0, 1.0)          # (P,)
    mw_norm = mw / mw.max()                              # (D,)
    beta = 2.2
    ecm_pen = np.exp(-beta * np.outer(ecm, mw_norm))     # (P,D)
    raw = perfusion[:, None] * bbb[None, :] * ecm_pen    # (P,D)

    # ---- compress delivery into a NARROW minor band [0.8, 1.0] ----
    raw_norm = (raw - raw.min()) / (raw.max() - raw.min())
    delivery = 0.8 + 0.2 * raw_norm                      # (P,D) in [0.8, 1.0]

    # ---- combine: potency dominates, delivery barely modulates ----
    kill = sens[:, None] * hill[None, :] * delivery      # (P,D)

    # mild multiplicative log-normal noise, fixed seed
    rng = np.random.default_rng(0)
    noise = np.exp(rng.normal(0.0, 0.05, size=kill.shape))
    kill = kill * noise

    return np.clip(kill, 0.0, 1.0)


def main():
    cache = sl.load_cache()
    gt_kill = build_gt_kill(cache)
    result = sl.score_policies(gt_kill, seed=0)
    sl.print_json(result)


if __name__ == "__main__":
    main()
