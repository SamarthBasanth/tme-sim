#!/bin/bash
# Verify every downloaded file's byte count against the series' own
# suppl/filelist.txt (which carries an authoritative Size column).
set -u
TMP="$(/usr/bin/mktemp -d)"
for g in GSE318564 GSE318562 GSE335275; do
  /usr/bin/curl -sS --fail --max-time 120 \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${g:0:6}nnn/$g/suppl/filelist.txt" \
    >> "$TMP/manifest.tsv" || echo "MANIFEST FETCH FAILED $g"
done
BAD=0
for d in /tmp/vmp/dl_scrna /tmp/vmp/dl_visium /tmp/vmp/dl_xenium; do
  for f in "$d"/*; do
    [ -f "$f" ] || continue
    n="$(basename "$f")"
    have="$(/usr/bin/wc -c < "$f" | /usr/bin/tr -d ' ')"
    want="$(/usr/bin/awk -F'\t' -v N="$n" '$2==N {print $4}' "$TMP/manifest.tsv" | /usr/bin/head -1)"
    if [ -z "$want" ]; then
      echo "NO MANIFEST ENTRY  $n  ($have bytes on disk)"; BAD=$((BAD+1))
    elif [ "$have" = "$want" ]; then
      echo "OK    $n  $have"
    else
      echo "SIZE MISMATCH  $n  disk=$have manifest=$want"; BAD=$((BAD+1))
    fi
  done
done
echo "---"
echo "problems: $BAD"
/bin/rm -rf "$TMP"
