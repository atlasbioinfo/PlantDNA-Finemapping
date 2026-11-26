#!/usr/bin/env bash
# Reproduce the analysis of the paper for A. thaliana chromosomes 2 and 4 (and FLC) on the CPU, with the precomputed
# PlantDNA-FM attention. No GPU and no PyTorch are needed.
#   bash examples/run_example.sh [workers]
# 1. check the input files
# 2. fine-mapping with the association tests, M and lambda of the paper -> results/arabidopsis/finemapping.csv
# 3. compare the result with the table of the paper
# 4. single-gene mode: FLC (AT5G10140), Bonferroni over its own 242 variants, lambda = 1.65
# To recompute the attention on a GPU, see examples/README.md.
set -euo pipefail
workers=${1:-4}
cd "$(dirname "$0")/.."
mkdir -p results
d=examples/arabidopsis
in="--tests $d/tests_B2.tsv.gz --expression $d/expression_B2.tsv.gz --genotypes $d/genotypes_B2 --attention $d/attention_PlantDNA-FM.parquet"
paper="--folds 2:4,5 --n-tests 4297315 --lambda 1.35,1.65"     # genome-wide M and cross-fitted lambda of the paper

plantdna-finemap genes --gff $d/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap check $in --genes results/genes.tsv
plantdna-finemap finemap $in --genes results/genes.tsv $paper --out results/arabidopsis --workers $workers --no-check
plantdna-finemap compare --result results/arabidopsis/finemapping.csv \
  --expected examples/expected/arabidopsis_chr2_chr4_finemapping_paper.csv.gz
plantdna-finemap finemap $in --genes results/genes.tsv --single-gene AT5G10140 --lambda 1.65 --out results/arabidopsis_FLC --workers 1 --no-check
echo "Example finished: results/arabidopsis/finemapping.csv reproduces the table of the paper for chromosomes 2 and 4."
