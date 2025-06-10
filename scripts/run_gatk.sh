#!/bin/bash
# Usage: run_gatk.sh <reference.fa> <ploidy_dir> <output_dir>

set -euo pipefail

REF=$1
PLOIDY_DIR=$2
OUTDIR=$3
mkdir -p "$OUTDIR"

INTERVALS="$OUTDIR/intervals.interval_list"
gatk PreprocessIntervals -R "$REF" --bin-length 1000 --padding 0 -O "$INTERVALS"

for CRAM in source/*.cram; do
  SAMPLE=$(basename "$CRAM" .cram | sed 's/MSSM_//' | sed 's/_NeuN_pl//')
  OUT_HDF5="$OUTDIR/$SAMPLE.hdf5"
  gatk CollectReadCounts -R "$REF" -L "$INTERVALS" -I "$CRAM" -O "$OUT_HDF5"
done

SCATTER_DIR="$OUTDIR/scatter"
gatk SplitIntervals -R "$REF" -L "$INTERVALS" -scatter-count 50 -O "$SCATTER_DIR"

gatk GermlineCNVCaller \
  --run-mode COHORT \
  --contig-ploidy-calls "$PLOIDY_DIR" \
  --input $OUTDIR/*.hdf5 \
  -L "$SCATTER_DIR" \
  -O "$OUTDIR/gcnv_output" --output-prefix cohort

gatk PostprocessGermlineCNVCalls \
  -O "$OUTDIR/gcnv_cohort_calls.vcf" \
  -run-directory "$OUTDIR/gcnv_output" \
  -model-shard-index 0 --model-shard-count 1
