"""ADVERSARIAL GROUND TRUTH -- make the mechanistic SIM LOSE to feature-ML.

Blind spot exploited
--------------------
The PDE sim is a pure DELIVERY/PK model: it predicts kill from how much drug
physically reaches viable tumor (vasc perfusion, ECM diffusion, BBB penetrance,
molecular size). It has NO notion of the tumor's molecular target dependence.
On this cohort the sim's realized recommendation is degenerate -- it always
picks cisplatin (best-delivered mid-size agent) -- and potency always picks
doxorubicin (lowest IC50).

Real GBM response, however, is dominated by TRANSCRIPTIONAL-SUBTYPE pathway
addiction: which driver pathway the tumor depends on decides which drug class
actually kills it, largely independent of delivery once a therapeutic amount
gets in (all agents here are dosed to CNS-relevant exposure). Subtype is an
OBSERVABLE biomarker (RNA panel), so a cheap classifier that sees subtype +
perfusion can learn the response map -- but the delivery sim, blind to target
dependence, cannot.

Ground truth
------------
  kill[p,d] = A[sub(p),d] * potency_mod[p,d] * delivery_mod[p,d] * noise

  A[s,d]  intrinsic subtype x drug pathway-engagement (the DOMINANT term):
    proneural   -> lomustine   (alkylator/nitrosourea-responsive)
    classical   -> cisplatin   (EGFR-amp, crosslinker-responsive)
    mesenchymal -> doxorubicin (anthracycline-responsive, aggressive)
  potency_mod = 0.85 + 0.30 * sens[p] * hill(dose, ic50[d])   (mild)
  delivery_mod= 0.85 + 0.30 * bbb[d]*sqrt(vasc)*exp(-b*ecm*mw) (mild)

Delivery/potency only MODULATE (band ~[0.85,1.15]); the affinity gap
(0.95 vs <=0.5) is wide enough that they cannot flip the per-subtype winner.
So the true best drug tracks subtype -> ML (which observes subtype) recovers it,
while sim (cisplatin-always) is right only for classical patients.

Params: affinity matrix below, potency/delivery modulation band 0.85..1.15,
beta=4.0 ECM penalty, lognormal noise sigma=0.04, seed=0.
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

si = {n: feat_names.index(n) for n in
      ["subtype_proneural", "subtype_classical", "subtype_mesenchymal"]}
oneh = np.stack([feat[:, si["subtype_proneural"]],
                 feat[:, si["subtype_classical"]],
                 feat[:, si["subtype_mesenchymal"]]], axis=1)
subtype = oneh.argmax(1)  # 0=proneural 1=classical 2=mesenchymal

# drug order: temozolomide, lomustine, cisplatin, methotrexate, doxorubicin
# DOMINANT intrinsic subtype->drug engagement (target dependence the sim ignores)
A = np.array([
    [0.50, 0.95, 0.40, 0.35, 0.35],   # proneural   -> lomustine
    [0.40, 0.50, 0.95, 0.45, 0.35],   # classical   -> cisplatin
    [0.35, 0.50, 0.45, 0.40, 0.95],   # mesenchymal -> doxorubicin
])
aff = A[subtype]  # (P,D)

def hill(e, k, h=1.3):
    e = np.maximum(e, 0.0)
    return e**h / (e**h + k**h)

# --- mild potency modulation (in-band, cannot flip affinity winner) ---
pot = sens[:, None] * hill(dose, ic50[None, :])
pot_mod = 0.85 + 0.30 * (pot / pot.max())

# --- mild delivery modulation ---
mw_norm = mw / mw.max()
beta = 4.0
vasc_f = np.sqrt(vasc / vasc.max())
delivery = bbb[None, :] * vasc_f[:, None] * np.exp(-beta * ecm[:, None] * mw_norm[None, :])
del_mod = 0.85 + 0.30 * (delivery / delivery.max())

rng = np.random.default_rng(0)
noise = rng.lognormal(mean=0.0, sigma=0.04, size=(P, D))
raw = aff * pot_mod * del_mod * noise
gt_kill = np.clip(raw / raw.max(), 0.0, 1.0)

result = selection_lib.score_policies(gt_kill, seed=0)
# extra diagnostics for the report
tb = gt_kill.argmax(1)
result["_extra"] = {
    "true_best_dist": [int((tb == d).sum()) for d in range(D)],
    "best_by_subtype": {int(s): int(np.bincount(tb[subtype == s], minlength=D).argmax())
                        for s in range(3)},
}
selection_lib.print_json(result)
