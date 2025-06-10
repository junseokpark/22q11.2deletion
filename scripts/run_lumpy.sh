#!/bin/bash
# Usage: run_lumpy.sh <input_cram> <output_vcf>
set -euo pipefail

CRAM=$1
OUT=$2
BAM=${OUT%.vcf}.bam
SORTED_BAM=${OUT%.vcf}.sorted.bam

samtools view -b "$CRAM" > "$BAM"
samtools sort "$BAM" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"
lumpyexpress -B "$SORTED_BAM" -o "$OUT"
