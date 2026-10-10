"""What feature_types does the GSE335275 panel declare, and for which keys?"""
import json
import sys
import zipfile
from collections import Counter

z = zipfile.ZipFile(sys.argv[1] if len(sys.argv) > 1 else
                    "/tmp/vmp/dl_xenium/"
                    "GSM9809297_WU1221_Edge_cell_feature_matrix.zarr.zip")
att = json.loads(z.read("cell_features/.zattrs").decode())
keys = [str(k) for k in att["feature_keys"]]
types = [str(t) for t in att["feature_types"]]
print(Counter(types))
for t in sorted(set(types)):
    ex = [k for k, tt in zip(keys, types) if tt == t]
    print(f"  {t!r}: n={len(ex)}  e.g. {ex[:4]} ... {ex[-2:]}")
