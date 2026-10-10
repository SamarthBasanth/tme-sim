"""Dump the keys and shapes of a Xenium cell_feature_matrix.zarr.zip .zattrs."""
import json
import sys
import zipfile

z = zipfile.ZipFile(sys.argv[1])
att = json.loads(z.read("cell_features/.zattrs").decode())
print("keys:", list(att.keys()))
for k, v in att.items():
    if isinstance(v, list):
        print(f"{k}: list len {len(v)} head {v[:8]}")
    else:
        print(f"{k}: {v!r}")
for n in z.namelist():
    if n.endswith(".zarray"):
        a = json.loads(z.read(n).decode())
        print(n, "shape", a.get("shape"), "dtype", a.get("dtype"))
