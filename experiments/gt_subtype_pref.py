"""GROUND-TRUTH MODEL -- Motif #6: subtype_pref.

Biological assumption
---------------------
Each GBM transcriptional subtype has an INTRINSIC pharmacologic preference:
its driver pathways make it depend on (and respond to) some drug classes more
than others, independent of how much drug reaches the tumor. We encode this as a
subtype x drug affinity matrix A[s,d] (target/pathway engagement, in [0,1]):

  proneural  (PDGFRA/IDH-ish): alkylator-responsive -> temozolomide, lomustine
  classical  (EGFR-amp)      : crosslinker/antimetabolite -> cisplatin, methotrexate
  mesenchymal(NF1, aggressive): anthracycline/nitrosourea -> doxorubicin, lomustine

Intrinsic affinity is then MODULATED by DELIVERY (does the drug physically reach
viable tumor) and by POTENCY (once there, does the dose kill). So a drug a subtype
"prefers" can still lose if it cannot be delivered -- e.g. doxorubicin is the most
potent agent (IC50 0.5) and mesenchymal's top intrinsic partner, yet it barely
crosses the BBB (bbb 0.05) and is a large molecule (MW 543) that penetrates dense
ECM poorly, so it is frequently NOT the realized best drug.

  kill[p,d] = A[sub(p),d] * delivery[p,d] * potency[p,d] * (1 + noise)

delivery[p,d] = bbb[d] * vasc[p]^0.5 * exp(-beta * ecm[p] * mw_norm[d])
  (BBB permeability gates entry; perfusion feeds convective delivery; dense ECM
   blocks diffusion, worse for large molecules)
potency[p,d]  = sens[p] * dose/(dose + ic50[d])
  (Hill/Emax occupancy scaled by patient chemo-sensitivity)

Key params: beta=4.0 (ECM/size penetration penalty), vasc exponent 0.5,
affinity matrix below, multiplicative-lognormal noise sigma=0.05, seed=0.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selection_lib

c = selection_lib.load_cache()
panel = c["panel"]
vasc, ecm, sens = c["vasc"], c["ecm"], c["sens"]
mw, ic50, bbb = c["drug_mw"], c["drug_ic50"], c["drug_bbb"]
dose = float(c["dose"])
feat, feat_names = c["feat"], list(c["feat_names"])
P, D = len(vasc), len(panel)

# subtype index per patient from one-hots (proneural,classical,mesenchymal)
si = {n: feat_names.index(n) for n in
      ["subtype_proneural", "subtype_classical", "subtype_mesenchymal"]}
oneh = np.stack([feat[:, si["subtype_proneural"]],
                 feat[:, si["subtype_classical"]],
                 feat[:, si["subtype_mesenchymal"]]], axis=1)
subtype = oneh.argmax(1)  # 0=proneural 1=classical 2=mesenchymal

# --- intrinsic subtype x drug affinity (motif core) ---
# drug order: temozolomide, lomustine, cisplatin, methotrexate, doxorubicin
# drug order: temozolomide, lomustine, cisplatin, methotrexate, doxorubicin
# proneural  : lomustine-dominant (clear intrinsic winner)
# classical  : cisplatin-dominant (MTX close in affinity but poorly delivered)
# mesenchymal: doxorubicin is the intrinsic favorite but is delivery-dead
#              (bbb 0.05, MW 543); lomustine & cisplatin are a near-tie, so the
#              REALIZED best flips patient-by-patient with delivery -> subtype
#              alone cannot predict it, only a delivery-aware model can.
A = np.array([
    [0.60, 0.95, 0.40, 0.35, 0.25],   # proneural  -> lomustine
    [0.45, 0.55, 0.95, 0.85, 0.30],   # classical  -> cisplatin
    [0.35, 0.82, 0.80, 0.30, 0.90],   # mesenchymal-> dox(intrinsic) / lom~cis tie
])
aff = A[subtype]                       # (P,D)

# --- delivery term ---
mw_norm = mw / mw.max()                # in (0,1]
beta = 4.0
vasc_f = (vasc / vasc.max()) ** 0.5    # perfusion, gently compressed
ecm_pen = np.exp(-beta * ecm[:, None] * mw_norm[None, :])
delivery = bbb[None, :] * vasc_f[:, None] * ecm_pen   # (P,D)

# --- potency term ---
potency = sens[:, None] * (dose / (dose + ic50[None, :]))  # (P,D)

# --- combine + mild seeded noise ---
rng = np.random.default_rng(0)
raw = aff * delivery * potency
noise = rng.lognormal(mean=0.0, sigma=0.05, size=raw.shape)
raw = raw * noise
gt_kill = raw / raw.max()              # scale into [0,1]
gt_kill = np.clip(gt_kill, 0.0, 1.0)

result = selection_lib.score_policies(gt_kill, seed=0)
selection_lib.print_json(result)
