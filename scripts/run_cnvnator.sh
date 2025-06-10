#!/bin/bash
# Usage: run_cnvnator.sh <input_cram> <output_prefix> <reference.fa> <chrom_fa_dir>
set -euo pipefail

CRAM=$1
OUT=$2
REF=$3
CHROM_DIR=$4

ROOT=${OUT}.root

cnvnator -root "$ROOT" -tree "$CRAM"
cnvnator -root "$ROOT" -his 100 -d "$CHROM_DIR"
cnvnator -root "$ROOT" -stat 100
cnvnator -root "$ROOT" -partition 100
cnvnator -root "$ROOT" -call 100 > "$OUT"
