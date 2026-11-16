#!/usr/bin/env bash
# Run the A. thaliana example on the CPU with the precomputed PlantDNA16K-FM attention.
#   bash examples/run_example.sh [workers]
# 1. genome-wide mode: weighted Bonferroni eGenes and SuSiE, lambda cross-fitted over chromosomes 1-3 / 4-5
# 2. single-gene mode: FLC (AT5G10140), Bonferroni over its own 242 variants, lambda = 1.65
# To recompute the attention on a GPU first, see examples/README.md.
set -euo pipefail
workers=${1:-4}
cd "$(dirname "$0")/.."
d=examples/arabidopsis
plantdna-finemap finemap --tests $d/tests_B2.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 \
  --attention $d/attention_PlantDNA16K-FM.parquet --out results/arabidopsis --workers $workers
plantdna-finemap finemap --tests $d/tests_B2.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 \
  --attention $d/attention_PlantDNA16K-FM.parquet --single-gene AT5G10140 --lambda 1.65 --out results/arabidopsis_FLC --workers 1
