#!/bin/bash
# Usage: run_canvas.sh <input_cram> <output_vcf> <reference.fa> <sample>

set -euo pipefail

CRAM=$1
OUT=$2
REF=$3
SAMPLE=$4
OUTDIR="results/canvas/$SAMPLE"

canvas germline --bam "$CRAM" --reference "$REF" --output-type vcf --output-dir "$OUTDIR"
mv "$OUTDIR"/*.vcf "$OUT"
