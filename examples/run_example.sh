#!/usr/bin/env bash
# Run the A. thaliana example (chromosomes 2 and 4, and FLC) on the CPU with the precomputed PlantDNA16K-FM attention.
#   bash examples/run_example.sh [workers]
# 1. published association tests: weighted Bonferroni eGenes and SuSiE, lambda cross-fitted between chromosomes 2 and 4
# 2. single-gene mode: FLC (AT5G10140), Bonferroni over its own 242 variants, lambda = 1.65
# 3. the whole workflow from genotypes and expression: genes -> assoc -> finemap
# To recompute the attention on a GPU, see examples/README.md.
set -euo pipefail
workers=${1:-4}
cd "$(dirname "$0")/.."
d=examples/arabidopsis
plantdna-finemap finemap --tests $d/tests_B2.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 \
  --attention $d/attention_PlantDNA16K-FM.parquet --folds 2:4,5 --out results/arabidopsis --workers $workers
plantdna-finemap finemap --tests $d/tests_B2.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 \
  --attention $d/attention_PlantDNA16K-FM.parquet --single-gene AT5G10140 --lambda 1.65 --out results/arabidopsis_FLC --workers 1
plantdna-finemap genes --gff $d/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap assoc --genotypes $d/genotypes_B2 --expression $d/expression_B2.tsv.gz --genes results/genes.tsv \
  --out results/assoc.tsv.gz --workers $workers
plantdna-finemap finemap --tests results/assoc.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 \
  --attention $d/attention_PlantDNA16K-FM.parquet --folds 2:4,5 --out results/arabidopsis_from_genotypes --workers $workers
