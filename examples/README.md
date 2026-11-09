# Examples

One chromosome per species, with the expression data, genotypes and genome-wide kinship of the paper.

| | *A. thaliana* | maize |
|---|---|---|
| chromosome | 5 | 10 |
| accessions | 471 | 263 |
| genes with expression | 3,772 | 1,155 |
| SNVs in the VCF | 320,096 | 220,668 |
| kinship | 755,957 genome-wide SNVs | 19,878,839 genome-wide SNVs |

## Run (CPU)

```bash
bash examples/run_example.sh arabidopsis 8
bash examples/run_example.sh maize 8
```

The script builds the windows from the GFF3 and runs `map` with the precomputed attention in `attention/`. Expected `summary.tsv` with 8 workers (times on a 32-core server):

*A. thaliana* (2 min 41 s):

| prior | genes tested | tests | eGenes | eGenes with eQTL | eQTLs | single-variant eQTLs | variants PIP > 0.9 | mean eQTL size |
|---|---|---|---|---|---|---|---|---|
| uniform | 3,771 | 367,954 | 1,114 | 1,082 | 1,654 | 559 | 681 | 4.70 |
| PlantDNA-FM | 3,771 | 367,954 | 1,119 | 1,088 | 1,684 | 586 | 707 | 4.48 |

Maize (1 min 28 s):

| prior | genes tested | tests | eGenes | eGenes with eQTL | eQTLs | single-variant eQTLs | variants PIP > 0.9 | mean eQTL size |
|---|---|---|---|---|---|---|---|---|
| uniform | 1,155 | 229,004 | 625 | 606 | 999 | 247 | 308 | 10.66 |
| PlantDNA-FM | 1,155 | 229,004 | 627 | 609 | 1,006 | 254 | 332 | 9.64 |

The eGene threshold uses the number of tests on this chromosome only, so these counts are not a subset of the genome-wide results in the paper.

## Recompute the attention (GPU)

```bash
bash examples/download_reference.sh arabidopsis        # TAIR10 chromosome 5 from Ensembl Plants
plantdna-finemap windows --gff examples/arabidopsis/TAIR10_chr5.gff3.gz --out results/arabidopsis/windows.tsv
plantdna-finemap attention --model <PlantDNA-FM directory or Hugging Face ID> --windows results/arabidopsis/windows.tsv \
    --fasta examples/arabidopsis/TAIR10_chr5.fa --vcf examples/arabidopsis/genotypes_chr5.vcf.gz --out results/arabidopsis/attention
plantdna-finemap map --vcf examples/arabidopsis/genotypes_chr5.vcf.gz --expression examples/arabidopsis/expression_chr5.tsv.gz \
    --windows results/arabidopsis/windows.tsv --kinship examples/arabidopsis/kinship.tsv.gz \
    --attention results/arabidopsis/attention/attention.parquet --out results/arabidopsis --workers 8
```

For maize, use `download_reference.sh maize`, `AGPv3.31_chr10.gff3.gz`, `AGPv3_chr10.fa` and the chromosome 10 files. On an RTX 5090, one gene takes about 1.5 s.

## Files

| file | content |
|---|---|
| `genotypes_chr*.vcf.gz` (+ `.tbi`) | SNVs in the 8-kb gene windows with MAF ≥ 0.03 in these accessions (`map` then applies MAF ≥ 0.05 per gene); GT only |
| `expression_chr*.tsv.gz` | expression of the genes on this chromosome, as used in the paper |
| `kinship.tsv.gz` | genome-wide relatedness matrix (output of `kinship` on all SNVs; MAF ≥ 0.05, call rate ≥ 0.9) |
| `*.gff3.gz` | genome annotation of this chromosome |

## Sources

Please cite the original data when you use these files.

- *A. thaliana* genotypes: 1001 Genomes Consortium. 1,135 genomes reveal the global pattern of polymorphism in *Arabidopsis thaliana*. *Cell* 166, 481–491 (2016); `1001genomes_snp-short-indel_only_ACGTN.vcf.gz`. Annotation: TAIR10 / Araport11 genes from Ensembl Plants.
- *A. thaliana* expression: leaf RNA-seq of 1001 Genomes accessions (Kawakatsu T. *et al.* Epigenomic diversity in a global collection of *Arabidopsis thaliana* accessions. *Cell* 166, 492–505 (2016)), reprocessed by pseudogenome quantification (Voichek Y. *et al.*, 2024); batch B2, 471 accessions, as used in the paper.
- Maize genotypes and expression: Kremling K. A. G. *et al.* Dysregulation of expression correlates with rare-allele burden and fitness loss in maize. *Nature* 555, 520–523 (2018). Genotypes: HapMap 3.2.1 (KNN-imputed, MAF > 0.05), REF set to the B73 AGPv3 base. Expression: third leaf base (L3Base). Annotation: Ensembl Plants release 31 (AGPv3.31).
