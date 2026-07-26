"""Summary figure for the drug-selection rigor study (15 independent ground truths)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# (name, acc_sim, acc_pot, acc_ml, frac_delivlim, rank_corr)  -- from the synthesis scoreboard
R = [
    ("multiplicative",1.00,0.00,1.00,1.00,0.40),
    ("liebig_min",0.57,0.00,0.77,1.00,0.48),
    ("threshold",0.77,0.00,0.80,1.00,0.42),
    ("hypoxia_resist",0.33,0.03,0.63,0.97,0.42),
    ("ecm_sieving",0.30,0.00,0.73,1.00,0.56),
    ("bbb_gated",0.87,0.03,0.87,0.97,0.79),
    ("subtype_pref",0.43,0.00,0.90,1.00,0.38),
    ("saturating",0.33,0.00,0.80,1.00,0.77),
    ("potency_dominant",0.00,1.00,1.00,0.00,0.51),
    ("delivery_dominant",0.00,0.00,0.87,1.00,-0.11),
    ("noisy_realistic",0.37,0.10,0.50,0.90,0.50),
    ("anti_delivery*",0.00,0.33,0.97,0.67,0.02),
    ("deliv×sens-gate*",0.73,0.10,0.63,0.90,0.88),
    ("subtype_blindspot*",0.33,0.33,1.00,0.67,0.35),
    ("potency_dominant2*",0.00,1.00,1.00,0.00,0.63),
]
names=[r[0] for r in R]; sim=np.array([r[1] for r in R]); pot=np.array([r[2] for r in R])
ml=np.array([r[3] for r in R]); frac=np.array([r[4] for r in R]); rc=np.array([r[5] for r in R])

fig=plt.figure(figsize=(14,5.5))
# Panel A: per-motif accuracy, sim vs potency vs ML
ax=fig.add_subplot(1,2,1); x=np.arange(len(R)); w=0.27
ax.bar(x-w,sim,w,label="mechanistic sim",color="#2ca02c")
ax.bar(x,pot,w,label="potency shortcut",color="#ff7f0e")
ax.bar(x+w,ml,w,label="feature-ML (CV)",color="#1f77b4")
ax.axhline(0.2,ls="--",c="gray",lw=1,label="random (0.20)")
ax.set_xticks(x); ax.set_xticklabels(names,rotation=55,ha="right",fontsize=7.5)
ax.set_ylabel("top-1 drug-selection accuracy"); ax.set_ylim(0,1.05)
ax.set_title("Per-motif accuracy — ML (blue) dominates the sim (green) in 14/15")
ax.legend(fontsize=8,loc="upper center",ncol=2)

# Panel B: circularity test — sim's edge over potency vs delivery-limitedness (genuine)
ax2=fig.add_subplot(1,2,2)
gap_pot=sim-pot
sc=ax2.scatter(frac,gap_pot,c=rc,cmap="viridis",s=90,edgecolor="k",linewidth=0.5,vmin=-0.1,vmax=0.9)
ax2.axhline(0,c="gray",lw=1)
ax2.set_xlabel("delivery-limitedness of ground truth\n(frac patients where best ≠ most-potent)")
ax2.set_ylabel("sim − potency accuracy gap")
ax2.set_title("Sim beats potency when delivery-limited (r=0.89)\ncolor = rank-corr to sim's own physics (circularity)")
cb=plt.colorbar(sc,ax=ax2); cb.set_label("sim_vs_gt_rank_corr")
ax2.annotate("genuine wins:\nhigh delivery-limit,\nMODERATE rank-corr",(0.55,0.55),fontsize=8,
             bbox=dict(boxstyle="round",fc="white",ec="gray",alpha=0.8))
fig.suptitle("Does spatial-transport physics help drug selection? Sim clears the potency floor but not the cheap-ML bar (15 independent ground truths)",fontsize=11)
plt.tight_layout(); fig.savefig("results/drug_selection_rigor.png",dpi=130)
print("wrote results/drug_selection_rigor.png")
print(f"sim beats potency: {int((gap_pot>0).sum())}/15 ; sim beats ML: {int(((sim-ml)>0).sum())}/15")
