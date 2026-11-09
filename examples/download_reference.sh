#!/usr/bin/env bash
# Reference chromosome for the example (needed only to recompute the attention):
#   bash examples/download_reference.sh arabidopsis   # TAIR10 chromosome 5 (Ensembl Plants release 31), 8 MB
#   bash examples/download_reference.sh maize         # B73 AGPv3 chromosome 10 (Ensembl Plants release 31), 45 MB
set -euo pipefail
cd "$(dirname "$0")"
base=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-31/fasta
case ${1:?usage: download_reference.sh arabidopsis|maize} in
  arabidopsis) url=$base/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.31.dna.chromosome.5.fa.gz; fa=arabidopsis/TAIR10_chr5.fa ;;
  maize)       url=$base/zea_mays/dna/Zea_mays.AGPv3.31.dna.chromosome.10.fa.gz;               fa=maize/AGPv3_chr10.fa ;;
  *) echo "unknown species: $1" >&2; exit 1 ;;
esac
curl -fL "$url" | gunzip -c > "$fa"
"$(command -v python3 || command -v python)" -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$fa"
echo "wrote examples/$fa"
