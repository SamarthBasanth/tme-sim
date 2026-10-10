"""Probe a Xenium cell_feature_matrix.zarr.zip: structure, panel genes, shape.

Read-only. Takes the archive path as argv[1]; never runs anything from the
download directory.
"""
import json
import sys
import zipfile

z = zipfile.ZipFile(sys.argv[1])
names = z.namelist()
print("n entries", len(names))
for n in names[:50]:
    print("  ", n)
print("--- metadata ---")
for n in names:
    base = n.rsplit("/", 1)[-1]
    if base in (".zattrs", ".zgroup", ".zarray", "zarr.json") and n.count("/") < 5:
        try:
            print("##", n)
            print(z.read(n)[:1500].decode())
        except Exception as exc:  # noqa: BLE001
            print("##", n, "ERR", exc)
