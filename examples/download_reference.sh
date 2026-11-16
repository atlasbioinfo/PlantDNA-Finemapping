#!/usr/bin/env bash
# Reference genome for the example (needed only to recompute the attention):
#   bash examples/download_reference.sh     # TAIR10 (Ensembl Plants release 31), 36 MB download, 120 MB unpacked
set -euo pipefail
cd "$(dirname "$0")"
url=https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-31/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa.gz
fa=arabidopsis/TAIR10.fa
curl -fL "$url" | gunzip -c > "$fa"
"$(command -v python3 || command -v python)" -c "import pysam, sys; pysam.faidx(sys.argv[1])" "$fa"
echo "wrote examples/$fa"
