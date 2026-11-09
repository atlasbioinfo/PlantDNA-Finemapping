#!/usr/bin/env bash
# Run the example for one species on the CPU with the precomputed attention.
#   bash examples/run_example.sh arabidopsis [workers]
#   bash examples/run_example.sh maize [workers]
# To recompute the attention on a GPU first, see examples/README.md.
set -euo pipefail
sp=${1:?usage: run_example.sh arabidopsis|maize [workers]}; workers=${2:-4}
cd "$(dirname "$0")/.."
case $sp in
  arabidopsis) d=examples/arabidopsis; c=5;  gff=$d/TAIR10_chr5.gff3.gz;     att=attention/arabidopsis_TAIR10_PlantDNA-FM.parquet ;;
  maize)       d=examples/maize;       c=10; gff=$d/AGPv3.31_chr10.gff3.gz; att=attention/maize_AGPv3_PlantDNA-FM.parquet ;;
  *) echo "unknown species: $sp" >&2; exit 1 ;;
esac
out=results/$sp; mkdir -p $out
plantdna-finemap windows --gff $gff --out $out/windows.tsv
plantdna-finemap map --vcf $d/genotypes_chr$c.vcf.gz --expression $d/expression_chr$c.tsv.gz --windows $out/windows.tsv \
  --kinship $d/kinship.tsv.gz --attention $att --out $out --workers $workers
