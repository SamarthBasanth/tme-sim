"""Shared data objects for the molecular-scale layer.

This layer is the multiscale BRIDGE: it takes a real drug *molecule* (structure)
and computes the physical parameters the tissue-scale transport engine needs,
via physics-based estimators and a Brownian-dynamics particle simulation. Its
final product is a `core.DrugParams` (so the existing engine consumes it with no
changes) plus a richer `MolecularProfile` that records provenance.

Honesty note carried in-code: this is NOT atomistic molecular dynamics of a whole
tissue (that is not computable — ~10^15 atoms). It is (a) cheminformatics/
physicochemical estimation from structure, and (b) genuine overdamped Langevin /
Brownian dynamics at the micro-scale to derive an effective diffusivity in a
crowded medium. Target binding affinity/potency is the weakest link from
structure alone and is treated as an input or a flagged estimate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class Molecule:
    """A drug molecule and its computed physicochemical descriptors.

    descriptors keys (filled by molecular.descriptors), all optional but expected:
      mw (g/mol), logp, tpsa (A^2), hbd, hba, rotatable_bonds, aromatic_rings,
      formal_charge, radius_hydrodynamic_nm, is_antibody (bool-ish).
    """

    name: str
    smiles: str
    descriptors: Dict[str, float] = field(default_factory=dict)
    source: str = "rdkit"  # "rdkit" or "fallback_table"

    def get(self, key: str, default=None):
        return self.descriptors.get(key, default)


@dataclass
class MolecularProfile:
    """The molecular layer's full physics output for one drug in one tissue context.

    All transport-relevant quantities the tissue engine cares about, plus enough
    provenance to defend every number.

    free_diffusivity : Stokes-Einstein diffusivity in dilute interstitial fluid, um^2/s.
    effective_diffusivity : diffusivity after Brownian-dynamics hindrance through the
        crowded ECM/cell medium at `ecm_reference`, um^2/s. THIS feeds DrugParams.diffusivity.
    hindrance : effective_diffusivity / free_diffusivity in [0,1].
    hindrance_curve : {ecm_fraction: D_eff} sampled from the Brownian sim (for plots).
    membrane_permeability : passive permeability, um/s (cell/BBB-relevant).
    bbb_penetrance : [0,1] fraction of vascular conc crossing into brain tissue.
    free_fraction : [0,1] unbound fraction in plasma (1 - protein binding).
    tissue_binding_rate : first-order retention/binding in tissue, 1/s.
    clearance : first-order clearance/metabolism in tissue, 1/s.
    ic50 : potency (uM). Provided or flagged-estimate (see provenance).
    ecm_reference : the ECM volume fraction the effective_diffusivity was computed at.
    provenance : {quantity: short string on method + confidence}.
    """

    molecule: Molecule
    free_diffusivity: float
    effective_diffusivity: float
    hindrance: float
    membrane_permeability: float
    bbb_penetrance: float
    free_fraction: float
    tissue_binding_rate: float
    clearance: float
    ic50: float
    ecm_reference: float = 0.4
    hindrance_curve: Dict[float, float] = field(default_factory=dict)
    provenance: Dict[str, str] = field(default_factory=dict)

    def summary(self) -> Dict:
        return {
            "name": self.molecule.name,
            "mw": self.molecule.get("mw"),
            "logp": self.molecule.get("logp"),
            "free_D_um2s": round(self.free_diffusivity, 2),
            "eff_D_um2s": round(self.effective_diffusivity, 2),
            "hindrance": round(self.hindrance, 3),
            "bbb_penetrance": round(self.bbb_penetrance, 3),
            "free_fraction": round(self.free_fraction, 3),
            "ic50_uM": self.ic50,
        }
