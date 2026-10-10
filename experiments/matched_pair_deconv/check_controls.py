"""Are all 180 panel entries missing from the scRNA index control probes?"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
import io_loaders as L  # noqa: E402
from run_test import dedup_first  # noqa: E402

panel = L.xenium_panel_genes(
    "/tmp/vmp/dl_xenium/GSM9809297_WU1221_Edge_cell_feature_matrix.zarr.zip")
_sc, sc_genes, _ = L.read_10x("/tmp/vmp/ex_scrna")
sc_dedup, _ = dedup_first(sc_genes)
known = set(sc_dedup)
absent = [g for g in panel if g not in known]
prefixes = ("NegControlProbe", "NegControlCodeword", "UnassignedCodeword",
            "DeprecatedCodeword", "antisense", "BLANK", "Intergenic")
ctrl = [g for g in absent if g.startswith(prefixes)]
other = [g for g in absent if not g.startswith(prefixes)]
print(f"panel {len(panel)}  absent {len(absent)}  control-prefixed {len(ctrl)}  "
      f"real-gene-looking {len(other)}")
print("non-control absentees:", other)
real = [g for g in panel if not g.startswith(prefixes)]
print(f"real gene targets in panel: {len(real)};  of those absent from scRNA: "
      f"{len([g for g in real if g not in known])}")
