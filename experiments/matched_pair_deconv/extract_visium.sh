#!/bin/bash
# Extract every Visium tarball in $1 into its own directory under $2.
set -eu
SRC="$1"; DST="$2"
mkdir -p "$DST"
for f in "$SRC"/*.tar.gz; do
  b="$(basename "$f" .tar.gz)"
  mkdir -p "$DST/$b"
  /usr/bin/tar -xzf "$f" -C "$DST/$b"
  echo "extracted $b"
done
