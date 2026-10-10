#!/bin/bash
# Fetch GEO supplementary files for the matched-pair deconvolution test.
# Usage: fetch.sh <outdir> <gsm> <filename> [<gsm> <filename> ...]
set -u
OUT="$1"; shift
mkdir -p "$OUT"
while [ "$#" -ge 2 ]; do
  GSM="$1"; FN="$2"; shift 2
  PRE="${GSM:0:7}nnn"
  URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/${PRE}/${GSM}/suppl/${FN}"
  if [ -s "$OUT/$FN" ]; then echo "have   $FN"; continue; fi
  if curl -sS --fail --max-time 1200 --retry 3 --retry-delay 5 -o "$OUT/$FN" "$URL"; then
    echo "ok     $FN  $(wc -c < "$OUT/$FN") bytes"
  else
    echo "FAILED $FN  ($URL)"
  fi
done
