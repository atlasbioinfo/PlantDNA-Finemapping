# Example

*A. thaliana* batch B2 (471 accessions): the 4,872 tested genes of chromosomes 2 and 4 plus *FLC* (AT5G10140, chromosome 5).

| file | content |
|---|---|
| `tests_B2.tsv.gz` | 1,307,160 GEMMA tests of Voichek *et al.* (gene, chrom, pos, ref, alt, maf, logp = −log₁₀*P*); variants in gene body ± 10 kb. These are the tests used in the paper. |
| `expression_B2.tsv.gz` | expression of the 4,873 genes in the 471 accessions, as used for the tests |
| `genotypes_B2.bed.gz/.bim.gz/.fam` | PLINK 1 fileset (`.bed` and `.bim` gzipped; read directly) of the 385,382 tested positions that are biallelic SNVs in the 1001 Genomes VCF (A1 = ALT), 471 accessions |
| `attention_PlantDNA-FM.parquet` | PlantDNA-FM attention score of every test, as used in the paper (same as the `attention` output) |
| `TAIR10_example.gff3.gz` | TAIR10 / Araport11 annotation (Ensembl Plants) of chromosomes 2 and 4 and of *FLC*; input of `genes` |

## Run (CPU)

```bash
bash examples/run_example.sh 8
```

No GPU is needed: the script uses the CPU environment and the attention scores shipped here. It checks the input files, runs two analyses (about 3 min with 8 workers on a 32-core server) and compares the result with the table of the paper (`examples/expected/arabidopsis_chr2_chr4_finemapping_paper.csv.gz`). The run is correct if the log contains

```
  check: all checks passed
  compare: PASS, the result reproduces the expected table
```


**1. The analysis of the paper for chromosomes 2 and 4.** The published association tests, the genome-wide number of tests of the paper (*M* = 4,297,315) and the λ of the paper (1.35 for chromosomes 1–3, 1.65 for chromosomes 4–5):

```bash
plantdna-finemap genes --gff examples/arabidopsis/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap finemap --tests examples/arabidopsis/tests_B2.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA-FM.parquet \
    --genes results/genes.tsv --folds 2:4,5 --n-tests 4297315 --lambda 1.35,1.65 --out results/arabidopsis --workers 8
```

The result is `results/arabidopsis/finemapping.csv`: 52,683 variants that are significant with the uniform prior or with the PlantDNA-FM prior, one row each, in the format of the table of the paper. It is the table of the paper for chromosomes 2 and 4: the same variants, with identical prior weights, thresholds, significance and eGene calls, and PIPs within 4 × 10⁻⁵.

```
gene,gene_name,strand,gene_start,gene_end,rs,ref,alt,maf,raw_pval (-log10),n_tested_in_gene,FM_prior_weight,threshold_FM (-log10),sig_raw,eGene_raw,sig_FM,eGene_FM,FM_new_sig,PIP_raw,PIP_FM
AT2G01090,,+,79779,81710,2_80921,T,C,0.136,12.6677,275,5.5189,7.1924,True,True,True,True,False,0.55547,0.79387
```

`results/arabidopsis/summary.tsv` gives the counts, identical to the results of the paper for these two chromosomes:

| prior | λ | genes tested | *M* | eGenes | eGenes with credible set | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|
| uniform | – | 4,873 | 4,297,315 | 1,314 | 1,307 | 797 | 904 | 4.65 |
| PlantDNA-FM | 1.35 / 1.65 | 4,873 | 4,297,315 | 1,365 | 1,362 | 817 | 924 | 4.26 |

**2. Single-gene mode for *FLC*** (`--single-gene AT5G10140 --lambda 1.65`, Bonferroni over the 242 variants of *FLC*). Expected: `results/arabidopsis_FLC/finemapping.csv` lists the significant variants of *FLC*, and `signals.tsv` has one credible set with each prior, 7 variants (uniform; lead 3,178,749, PIP 0.70) and 6 variants (PlantDNA-FM; lead 3,178,749, PIP 0.81).

On your own data, leave out `--n-tests`, `--lambda` and `--folds`: *M* is then the number of tests in your table and λ is learned by cross-fitting.

## Recompute the attention (GPU)

```bash
bash examples/download_reference.sh        # TAIR10 from Ensembl Plants -> examples/arabidopsis/TAIR10.fa
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA-FM --genes results/genes.tsv \
    --fasta examples/arabidopsis/TAIR10.fa --tests examples/arabidopsis/tests_B2.tsv.gz --out results/attention
```

`results/genes.tsv` is written by `run_example.sh` (`genes --all-biotypes`, which keeps the few tested genes that are not protein-coding). Then pass `--attention results/attention/attention.parquet` to `finemap`. `--model` also takes a local model directory. On an RTX 5090 the 4,873 genes (10,544 windows) take 3 h 58 min; add `--only genes.txt` (one gene ID per line) to score a few genes first. The result matches the shipped scores (the paper's): the largest relative difference over the 1,307,160 scores is 2.6 × 10⁻⁷.

## Sources

Please cite the original data when you use these files.

- Expression and association tests: leaf RNA-seq of 1001 Genomes accessions (Kawakatsu T. *et al.* Epigenomic diversity in a global collection of *Arabidopsis thaliana* accessions. *Cell* 166, 492–505 (2016)), reprocessed by pseudogenome quantification, and GEMMA LMM tests per gene (`-lmm 2 -maf 0.05 -miss 0.5`) by Voichek Y. *et al.* (2024); batch B2, 471 accessions.
- Genotypes: 1001 Genomes Consortium. 1,135 genomes reveal the global pattern of polymorphism in *Arabidopsis thaliana*. *Cell* 166, 481–491 (2016); `1001genomes_snp-short-indel_only_ACGTN.vcf.gz`.
- Annotation and reference: TAIR10 / Araport11 from Ensembl Plants.
