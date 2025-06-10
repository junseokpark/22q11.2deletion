#!/bin/bash
# Usage: run_meerkat.sh <input_cram> <output_txt> <reference.fa> <sample>
set -euo pipefail

CRAM=$1
OUT=$2
REF=$3
SAMPLE=$4
OUTDIR="results/meerkat/$SAMPLE"

mkdir -p "$OUTDIR"

# Convert to BAM, sort, and index
BAM="$OUTDIR/$SAMPLE.bam"
SORTED="$OUTDIR/$SAMPLE.sorted.bam"
samtools view -b "$CRAM" > "$BAM"
samtools sort "$BAM" -o "$SORTED"
samtools index "$SORTED"

cd "$OUTDIR"
meerkat.pl config "$SAMPLE" "$SORTED" "$REF"
meerkat.pl all "$SAMPLE"

cp "$OUTDIR/$SAMPLE.final_fusion.txt" "$OUT"

