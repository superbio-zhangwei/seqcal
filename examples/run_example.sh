#!/usr/bin/env bash
set -euo pipefail

seqcal \
  -f examples/example.fna \
  -o examples/out \
  --forward-primer GATGAAGARCGYAGYRAA \
  --reverse-primer TCCTCCGCTTATTGATATGC \
  -t 2
