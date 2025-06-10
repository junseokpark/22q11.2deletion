#!/bin/bash

set -euo pipefail

CRAM=$1
OUT=$2
REF=$3
SAMPLE=$4
THREADS=$5

RUN_DIR="results/manta/$SAMPLE"

configManta.py --bam "$CRAM" --referenceFasta "$REF" --runDir "$RUN_DIR"
"$RUN_DIR/runWorkflow.py" -j "$THREADS"
