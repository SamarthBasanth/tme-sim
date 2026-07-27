"""Load the REAL Ivy Glioblastoma Atlas Project (Ivy GAP) anatomic-structure RNA-seq.

Source: Puchalski et al. 2018, *Science* 360(6389):660-663,
"An anatomic transcriptional atlas of human glioblastoma."  Allen Institute for
Brain Science, Ivy Glioblastoma Atlas Project (Ivy GAP).  Fully OPEN access
(http://glioblastoma.alleninstitute.org; data released under the Allen Institute
terms of use, freely downloadable, no registration/controlled tier).

Ivy GAP dissected GBM tissue by laser microdissection into named **anatomic
structures** (the real regional architecture of the tumour) and profiled each by
RNA-seq (RSEM -> gene-level FPKM, cross-sample normalised).  Two study designs
share one FPKM matrix (270 RNA samples, 37 tumours):

  * ``reference_histology`` -- 122 samples of the 5 canonical anatomic structures
    identified by H&E in 10 tumours:
        LE    Leading Edge                              (infiltrative margin / brain interface)
        IT    Infiltrating Tumor                        (infiltrative margin)
        CT    Cellular Tumor                            (tumour core)
        CTmvp Microvascular Proliferation               (angiogenic / perivascular)
        CTpan Pseudopalisading cells around necrosis    (perinecrotic border)
  * ``csc_cluster`` -- 148 samples of putative cancer-stem-cell / niche clusters
    identified by ISH probe expression (an 18-gene reference set), sampled by
    HIGH ("-<GENE>") or LOW ("-control-...") expression, under seven base
    structures.  These add two further vascular / necrotic niches:
        CThbv Hyperplastic blood vessels in cellular tumor (perivascular)
        CTpnz Perinecrotic zone                            (perinecrotic / hypoxic)

So the FULL base-structure vocabulary present in the RNA-seq (``base_structure``)
is: LE, IT, CT, CTmvp, CThbv, CTpan, CTpnz.

Separately, Ivy GAP publishes an **11-feature "tumour feature" histology ontology**
as machine-classified fractional areas per tumour sub-block (``sub_block_details``):
le, lehbv, it, ithbv, ct, ctpnz, ctpnn, ctpan, cthbv, ctmvp, ctne (necrosis).
These are the whole-slide region composition used to describe each block.

Cached files (see ``docs/real_data/ivygap.md``), extracted from the open
``gene_expression_matrix.zip`` plus the four open summary CSVs:
  * ``fpkm_table.csv``              - genes x samples normalised FPKM matrix
                                      (index = gene_id, columns = rna_well_id).
  * ``rows-genes.csv``             - gene_id -> gene_symbol / entrez / chromosome.
  * ``columns-samples.csv``        - rna_well_id -> structure label + tumour/block.
  * ``rna_seq_samples_details.csv``- per-sample structure ontology (acronym, name,
                                      parent) + donor clinical fields.
  * ``sub_block_details.csv``      - tumour-feature fractional areas per sub-block.
  * ``tumor_details.csv``          - per-donor molecular subtype, MGMT, survival.

Only pandas/numpy/stdlib are used.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import pandas as pd

# Default cache location (repo-root/data_cache/ivygap).
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_DEFAULT_CACHE = os.path.join(_REPO_ROOT, "data_cache", "ivygap")
# The acquisition step (tmesim/data/acquire.py) writes the zip under ivy_gap/;
# accept that spelling too, transparently.
_ALT_CACHE = os.path.join(_REPO_ROOT, "data_cache", "ivy_gap")

_FPKM_FILE = "fpkm_table.csv"
_GENES_FILE = "rows-genes.csv"
_SAMPLES_FILE = "columns-samples.csv"
_RNASEQ_DETAILS_FILE = "rna_seq_samples_details.csv"
_SUBBLOCK_FILE = "sub_block_details.csv"
_TUMOR_FILE = "tumor_details.csv"

# Map each base anatomic structure to a coarse microenvironment REGION so the
# reconstruction can label margin / core / perivascular / necrotic zones from
# real per-structure signatures.  (Derived from Puchalski 2018 histology.)
REGION_MAP: Dict[str, str] = {
    "LE": "margin",         # leading edge: tumour/normal-brain interface
    "IT": "margin",         # infiltrating tumour: infiltrative margin
    "CT": "core",           # cellular tumour: dense proliferative core
    "CTmvp": "perivascular",  # microvascular proliferation (angiogenic)
    "CThbv": "perivascular",  # hyperplastic blood vessels
    "CTpan": "necrotic",    # pseudopalisading cells around necrosis
    "CTpnz": "necrotic",    # perinecrotic zone (hypoxic)
}

# The 11 tumour-feature (H&E histology) fractional-area columns in
# sub_block_details.csv, with human-readable names.
TUMOR_FEATURES: Dict[str, str] = {
    "le": "Leading Edge",
    "lehbv": "Leading Edge - hyperplastic blood vessels",
    "it": "Infiltrating Tumor",
    "ithbv": "Infiltrating Tumor - hyperplastic blood vessels",
    "ct": "Cellular Tumor",
    "ctpnz": "Cellular Tumor - perinecrotic zone",
    "ctpnn": "Cellular Tumor - perinecrotic (palisading) - narrow",
    "ctpan": "Cellular Tumor - pseudopalisading cells around necrosis",
    "cthbv": "Cellular Tumor - hyperplastic blood vessels",
    "ctmvp": "Cellular Tumor - microvascular proliferation",
    "ctne": "Cellular Tumor - necrosis",
}


def _base_structure(acronym: str) -> str:
    """Base anatomic structure = the token before the first '-' (e.g. CTpan)."""
    return str(acronym).split("-", 1)[0]


def _sampling_family(acronym: str) -> str:
    """How the sample was dissected: reference histology vs. ISH gene cluster."""
    a = str(acronym)
    if "reference-histology" in a:
        return "reference_histology"
    if "-control-" in a:
        return "csc_cluster_control"
    return "csc_cluster_gene"


def _resolve_cache(cache_dir: Optional[str]) -> str:
    if cache_dir is not None:
        return cache_dir
    if os.path.exists(os.path.join(_DEFAULT_CACHE, _FPKM_FILE)):
        return _DEFAULT_CACHE
    if os.path.exists(os.path.join(_ALT_CACHE, _FPKM_FILE)):
        return _ALT_CACHE
    return _DEFAULT_CACHE


@dataclass
class IvyGapStructures:
    """Real Ivy GAP anatomic-structure RNA-seq, joined and ready to use.

    Attributes
    ----------
    samples : pandas.DataFrame
        One row per RNA sample (index = ``rna_well_id`` as str). Columns include
        ``structure_abbreviation``, ``structure_name``, ``base_structure``,
        ``region`` (margin/core/perivascular/necrotic), ``sampling_family``,
        ``tumor_id``, ``tumor_name`` and, where joined, donor clinical fields
        (``molecular_subtype``, ``mgmt_methylation``, ``survival_days``).
    expression : pandas.DataFrame
        genes x samples FPKM matrix (index = gene symbol, columns aligned to
        ``samples.index``). Normalised FPKM as released by Ivy GAP.
    genes : pandas.DataFrame
        Gene annotation (index = gene symbol; columns gene_id, entrez, chromosome).
    """

    samples: pd.DataFrame
    expression: pd.DataFrame
    genes: pd.DataFrame

    def region_for(self, sample_id: str) -> str:
        return self.samples.loc[str(sample_id), "region"]


def load_structures(
    cache_dir: Optional[str] = None,
    reference_histology_only: bool = False,
) -> IvyGapStructures:
    """Load Ivy GAP samples x anatomic-structure labels joined to FPKM expression.

    Parameters
    ----------
    cache_dir : str, optional
        Directory with the extracted Ivy GAP files. Defaults to
        ``data_cache/ivygap`` (falls back to ``data_cache/ivy_gap``).
    reference_histology_only : bool
        If True, keep only the 122 canonical H&E anatomic-structure samples
        (LE/IT/CT/CTmvp/CTpan). If False (default), return all 270 samples
        including the CSC-cluster niches (adds CThbv, CTpnz).

    Returns
    -------
    IvyGapStructures
        ``.samples`` (sample x structure), ``.expression`` (genes x samples FPKM),
        ``.genes``. See the dataclass docstring.

    Raises
    ------
    FileNotFoundError
        If the cache files are missing (run acquisition first).
    """
    cache_dir = _resolve_cache(cache_dir)
    fpkm_path = os.path.join(cache_dir, _FPKM_FILE)
    genes_path = os.path.join(cache_dir, _GENES_FILE)
    samples_path = os.path.join(cache_dir, _SAMPLES_FILE)
    for p in (fpkm_path, genes_path, samples_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Missing cached Ivy GAP file: {p}. Extract "
                "gene_expression_matrix.zip into data_cache/ivygap/ first "
                "(see docs/real_data/ivygap.md)."
            )

    # --- gene annotation: gene_id -> gene_symbol -------------------------------
    genes = pd.read_csv(genes_path)
    genes["gene_id"] = genes["gene_id"].astype(int)
    id2symbol = dict(zip(genes["gene_id"], genes["gene_symbol"].astype(str)))

    # --- sample metadata -------------------------------------------------------
    samples = pd.read_csv(samples_path)
    samples["rna_well_id"] = samples["rna_well_id"].astype(int).astype(str)
    samples = samples.set_index("rna_well_id")
    samples["base_structure"] = samples["structure_abbreviation"].map(_base_structure)
    samples["sampling_family"] = samples["structure_abbreviation"].map(_sampling_family)
    samples["region"] = samples["base_structure"].map(REGION_MAP)

    # Join donor clinical / molecular fields from the RNA-seq details table.
    details_path = os.path.join(cache_dir, _RNASEQ_DETAILS_FILE)
    if os.path.exists(details_path):
        det = pd.read_csv(details_path)
        det["sample_id"] = det["sample_id"].astype(int).astype(str)
        keep = [
            "sample_id", "donor_id", "molecular_subtype", "mgmt_methylation",
            "survival_days", "egfr_amplification", "age_in_years",
            "structure_parent_acronym",
        ]
        det = det[[c for c in keep if c in det.columns]].drop_duplicates("sample_id")
        det = det.set_index("sample_id")
        samples = samples.join(det, how="left")

    if reference_histology_only:
        samples = samples[samples["sampling_family"] == "reference_histology"]

    # --- FPKM matrix (genes x samples) -----------------------------------------
    fpkm = pd.read_csv(fpkm_path, index_col=0)
    fpkm.columns = [str(int(c)) for c in fpkm.columns]  # rna_well_id headers
    fpkm.index = fpkm.index.astype(int).map(id2symbol)  # gene_id -> gene symbol
    fpkm.index.name = "gene_symbol"
    # Collapse any duplicate symbols by mean (a handful of gene_ids share a symbol).
    if fpkm.index.duplicated().any():
        fpkm = fpkm.groupby(level=0).mean()

    # Align expression columns to the (possibly filtered) sample set.
    common = [s for s in samples.index if s in fpkm.columns]
    samples = samples.loc[common]
    expression = fpkm[common]

    genes_idx = genes.drop_duplicates("gene_symbol").set_index("gene_symbol")
    return IvyGapStructures(samples=samples, expression=expression, genes=genes_idx)


def structure_signatures(
    data: Optional[IvyGapStructures] = None,
    cache_dir: Optional[str] = None,
    by: str = "base_structure",
    reference_histology_only: bool = True,
) -> pd.DataFrame:
    """Mean FPKM expression per anatomic region -- the real per-region signatures.

    Averages the FPKM matrix over the samples of each anatomic structure so the
    reconstruction can score/label a region (margin, core, perivascular,
    necrotic) against a genuine measured transcriptional signature.

    Parameters
    ----------
    data : IvyGapStructures, optional
        Pre-loaded data; if None, loads via :func:`load_structures`.
    cache_dir : str, optional
        Passed to :func:`load_structures` when ``data`` is None.
    by : {"base_structure", "region", "structure_abbreviation"}
        Grouping key. ``base_structure`` -> the 5-7 anatomic structures;
        ``region`` -> the 4 coarse microenvironment zones (margin/core/
        perivascular/necrotic); ``structure_abbreviation`` -> every fine label.
    reference_histology_only : bool
        If True (default) signatures use only the clean H&E anatomic-structure
        samples. Set False to include CSC-cluster niches (needed to populate the
        CThbv/CTpnz perivascular/necrotic groups when grouping by base/region).

    Returns
    -------
    pandas.DataFrame
        genes x groups mean FPKM (index = gene symbol, columns = the grouping
        values). Column order follows a fixed anatomic ordering where applicable.
    """
    if data is None:
        # When grouping by region/base and asking for the vascular-necrotic
        # niches (CThbv, CTpnz), those only exist among CSC-cluster samples, so
        # default to reference_histology_only for the clean 5-structure case.
        data = load_structures(cache_dir=cache_dir, reference_histology_only=reference_histology_only)

    if by not in data.samples.columns:
        raise ValueError(f"'by' must be a column of samples; got {by!r}")

    groups = data.samples[by]
    # mean over samples within each group -> genes x groups
    sig = data.expression.T.groupby(groups).mean().T
    sig.columns.name = by

    # Stable, anatomically-ordered columns where we know the order.
    if by == "base_structure":
        order = ["LE", "IT", "CT", "CTmvp", "CThbv", "CTpan", "CTpnz"]
    elif by == "region":
        order = ["margin", "core", "perivascular", "necrotic"]
    else:
        order = None
    if order:
        cols = [c for c in order if c in sig.columns] + \
               [c for c in sig.columns if c not in order]
        sig = sig[cols]
    return sig


if __name__ == "__main__":
    d = load_structures()
    print(f"samples: {d.samples.shape[0]}  expression: {d.expression.shape[0]} genes "
          f"x {d.expression.shape[1]} samples")
    print("\nbase-structure sample counts:")
    print(d.samples["base_structure"].value_counts().to_string())
    print("\nmicroenvironment-region sample counts:")
    print(d.samples["region"].value_counts(dropna=False).to_string())
    print("\nsampling-family counts:")
    print(d.samples["sampling_family"].value_counts().to_string())

    sig = structure_signatures(by="base_structure", reference_histology_only=False)
    print(f"\nsignature matrix: {sig.shape[0]} genes x {sig.shape[1]} structures")
    for g in ("VEGFA", "PECAM1", "CA9", "MKI67"):
        if g in sig.index:
            print(f"  {g:7s} mean FPKM by structure: "
                  + ", ".join(f"{c}={sig.loc[g, c]:.1f}" for c in sig.columns))
