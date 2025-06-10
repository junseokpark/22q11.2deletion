#!/bin/bash
# Usage: run_delly.sh <input_cram> <output_vcf> <reference.fa>

set -euo pipefail

CRAM=$1
OUT=$2
REF=$3
BAM=${OUT%.vcf}.bam
BCF=${OUT%.vcf}.bcf

samtools view -b "$CRAM" | samtools sort -o "$BAM"
samtools index "$BAM"
delly call -g "$REF" -o "$BCF" "$BAM"
bcftools view "$BCF" > "$OUT"
