# PlantDNA-Finemapping

expGWAS fine-mapping in plant populations with PlantDNA-FM attention priors: weighted Bonferroni eGene detection and SuSiE fine-mapping.

Association panels in plants are small and highly structured, so an expression GWAS (expGWAS) signal usually points to a block of variants in near-perfect linkage disequilibrium rather than to one variant. This code uses attention from [PlantDNA-FM](https://huggingface.co/ATLASBIOINFO/PlantDNA-FM), a plant DNA foundation model with a 16-kb context, as a prior that ranks the variants within each gene, and applies the same prior to eGene detection and to fine-mapping.

Definitions used throughout:

- **expGWAS**: association of the expression of each gene with the variants in and around it (gene body ± 10 kb).
- **eGene**: a gene with at least one significant variant.
- **expGWAS signal**: a 95% credible set from SuSiE fine-mapping of an eGene. A gene can have several signals.

## Contents

- [Installation](#installation)
- [Quick start: reproduce the example](#quick-start-reproduce-the-example-cpu-about-3-min)
- [User guide: your own data](#user-guide-your-own-data)
  - [What you need](#what-you-need)
  - [Step 1: `genes`](#step-1-genes--gene-table-from-the-annotation)
  - [Step 2: association tests (GEMMA)](#step-2-association-tests-gemma)
  - [Step 3: `attention`](#step-3-attention--plantdna-fm-attention-scores-gpu)
  - [Step 4: `finemap`](#step-4-finemap--egenes-and-expgwas-signals)
  - [Check your inputs](#check-your-inputs)
  - [Candidate gene](#candidate-gene-single-gene-mode)
- [Troubleshooting](#troubleshooting)
- [Method](#method)
- [Results in the paper](#results-in-the-paper)
- [Precomputed attention](#precomputed-attention)
- [Reproducibility](#reproducibility)

## Installation

There are two environments. Most users need only the first one.

**CPU environment** (`genes`, `check`, `finemap`, `compare`). No GPU and no PyTorch. It is enough to reproduce the example and to analyse your own data when you already have an attention file.

```bash
git clone https://github.com/atlasbioinfo/PlantDNA-Finemapping.git
cd PlantDNA-Finemapping
conda env create -f environment.yml        # Python 3.11, numpy, scipy, pandas, pyarrow, pysam
conda activate plantdna-finemap
```

Without conda, `pip install .` in any Python ≥ 3.9 does the same.

**GPU environment** (adds `attention`). Needed only to compute attention scores for a new panel or species. It can be on another machine: compute `attention.parquet` on a GPU server, copy the file, and run `finemap` in the CPU environment.

```bash
conda env create -f environment-gpu.yml    # the CPU environment plus PyTorch and transformers
conda activate plantdna-finemap-gpu
plantdna-finemap check --tests examples/arabidopsis/tests_B2.tsv.gz --gpu     # reports PyTorch, transformers and the GPU
```

The attention step was tested with PyTorch 2.10 and transformers 5.9 on an RTX 5090; a 16-kb window needs about 12 GB of GPU memory and 1.4 s. The model is downloaded from Hugging Face ([`ATLASBIOINFO/PlantDNA-FM`](https://huggingface.co/ATLASBIOINFO/PlantDNA-FM), 1 GB) the first time `attention` runs. For the *A. thaliana* and maize panels of the paper the attention files are already in this repository ([Precomputed attention](#precomputed-attention)).

## Quick start: reproduce the example (CPU, about 3 min)

In the CPU environment:

```bash
bash examples/run_example.sh 8      # 8 worker processes
```

The example is *A. thaliana* batch B2 (471 accessions) for the 4,872 tested genes of chromosomes 2 and 4 plus *FLC*, with the association tests of the paper (the GEMMA results of Voichek *et al.*), the attention scores of the paper, and the settings of the paper. The script

1. checks the input files;
2. runs `finemap`, which writes `results/arabidopsis/finemapping.csv` (52,683 significant variants; 1,314 eGenes with the uniform prior and 1,365 with the PlantDNA-FM prior);
3. compares this table with the table of the paper;
4. runs `finemap` in single-gene mode for *FLC*.

The run is correct if step 3 ends with

```
  compare: PASS, the result reproduces the expected table
```

See [examples/README.md](examples/README.md) for the files, the full expected output and for recomputing the attention on a GPU.

## User guide: your own data

### What you need

| you provide | file | used by |
|---|---|---|
| reference genome | FASTA with a `.fai` index (plain or bgzipped) | `attention` |
| genome annotation | GFF3 (`.gz` allowed) | `genes` |
| genotypes | VCF (bgzipped, with a tabix index) or a PLINK 1 fileset | GEMMA, `finemap` |
| expression | table of genes × samples (TSV, `.gz` allowed) | GEMMA, `finemap` |
| association tests | one row per gene and variant, from GEMMA (TSV, `.gz` allowed) | `attention`, `finemap` |

The steps run in this order, each one reading the output of the previous ones:

```bash
plantdna-finemap genes     --gff genes.gff3.gz --out genes.tsv
# associations.tsv.gz: one test per gene and variant, from GEMMA as in the paper (step 2)
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA-FM --genes genes.tsv --fasta genome.fa \
                           --tests associations.tsv.gz --out attention/
plantdna-finemap finemap   --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
                           --attention attention/attention.parquet --genes genes.tsv --out results/ --workers 8
```

The result is one table, `results/finemapping.csv`, with one row per significant variant ([step 4](#step-4-finemap--egenes-and-expgwas-signals)).

`genes` and `finemap` run in the CPU environment; `attention` needs the GPU environment. Without `--attention`, `finemap` reports the uniform prior only. `attention` and `finemap` first check their input files and stop with a message if something does not fit ([Check your inputs](#check-your-inputs)).

Before you start:

- **Sample IDs** must be the same in the genotypes and in the header of the expression table. Samples present in only one of them are ignored.
- **Chromosome names** must match between the FASTA, the GFF3 and the genotypes, apart from a `chr`/`Chr` prefix.
- **Coordinates**: the genotypes and the annotation must be on the same assembly as the FASTA.
- **Expression values** are used as given. Normalise them as for any GWAS (for example log-transformed values, or a rank-based inverse normal transformation per gene). Missing values are allowed.

### Step 1: `genes` — gene table from the annotation

```bash
plantdna-finemap genes --gff genes.gff3.gz --out genes.tsv
```

**Input.** A GFF3 with `gene`, `mRNA` (or `transcript`) and `exon` features in the Ensembl style (`ID=gene:AT1G01010`, `Parent=transcript:AT1G01010.1`):

```
1  araport11  gene  3631  5899  .  +  .  ID=gene:AT1G01010;biotype=protein_coding
1  araport11  mRNA  3631  5899  .  +  .  ID=transcript:AT1G01010.1;Parent=gene:AT1G01010
1  araport11  exon  3631  3913  .  +  .  Parent=transcript:AT1G01010.1
```

**Options.** `--all-biotypes` keeps all genes instead of only `biotype=protein_coding`; `--chroms 1,2,3` keeps some chromosomes.

**Output.** `genes.tsv`, one row per gene:

| column | content |
|---|---|
| `gene`, `chrom`, `strand` | gene ID, chromosome, strand |
| `gene_start`, `gene_end` | gene body (1-based, inclusive) |
| `gene_name` | `Name` attribute of the gene, if any |
| `exons` | exons of the representative transcript (`<gene>.1` or `<gene>_T01`), `start-end,start-end,…`; empty if the gene has no such transcript |

### Step 2: association tests (GEMMA)

One test per gene and variant: the expression of the gene against each variant in gene body ± 10 kb. As in the paper, which uses the GEMMA results of Voichek *et al.* for *A. thaliana*, the tests come from GEMMA: a centred relatedness matrix, and for each gene a mixed-model test of the variants in gene body ± 10 kb with MAF ≥ 0.05 and call rate ≥ 0.5.

```bash
gemma -bfile genotypes -gk 1 -o kinship                                # once: output/kinship.cXX.txt
# per gene: put its expression in column 6 of the .fam file and list its variants (gene body +/- 10 kb) in snps.txt
gemma -bfile genotypes_gene -k output/kinship.cXX.txt -lmm 2 -snps snps.txt -maf 0.05 -miss 0.5 -o AT1G01010
```

Add a `gene` column to each GEMMA output (`output/AT1G01010.assoc.txt`) and concatenate the files into one table. This association table (`--tests` of steps 3 and 4) is a TSV (`.gz` allowed) or parquet with one row per variant–gene test:

| column | content |
|---|---|
| `gene` | gene ID |
| `chrom` (or `chr`) | chromosome |
| `pos` (or `ps`) | position of the variant |
| `logp`, or `p` / `p_wald` / `p_lrt` / `p_score` | −log₁₀*P*, or the *P* value |
| `maf` (or `af`), `ref`, `alt` (optional) | minor allele frequency and alleles; copied to the result table |

The GEMMA columns `chr`, `ps`, `af` and `p_lrt` are read as they are. Every row counts towards *M*, the number of tests of the Bonferroni threshold. The tests of any other GWAS tool can be used if they are given in this format.

Example (`examples/arabidopsis/tests_B2.tsv.gz`):

```
gene        chrom  pos     ref   alt   maf     logp
AT2G01060   2      68105   C     G     0.07    0.10141
AT2G01060   2      68151   C     T     0.314   2.38087
```

### Step 3: `attention` — PlantDNA-FM attention scores (GPU)

```bash
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA-FM --genes genes.tsv --fasta genome.fa \
    --tests associations.tsv.gz --out attention/
```

**Inputs.** `--model`: the Hugging Face ID, or a local directory with the model; `--genes`: output of step 1; `--fasta`: reference genome with a `.fai` index (`samtools faidx genome.fa`); `--tests`: the association table, which lists the variants to score for each gene.

**Options.** `--only genes.txt` (one gene ID per line) scores a subset; `--window` (16000), `--step` (8000), `--flank` (10000); `--device` (cuda). Rerunning the same command resumes an interrupted run.

**Time.** About 1.4 s per 16-kb window on an RTX 5090; a gene needs about two windows, so 10,000 genes take about 8 h.

**Outputs** in `attention/`:

| file | content |
|---|---|
| `attention.parquet` | one row per variant–gene test: `gene`, `pos`, `attention` (the score), `attention_source` (`exon windows` or `windows without exon`) |
| `chunks/` | the score of every variant in every window (kept for resuming) |

For the *A. thaliana* and maize panels of the paper, the scores are already in this repository ([Precomputed attention](#precomputed-attention)).

### Step 4: `finemap` — eGenes and expGWAS signals

```bash
plantdna-finemap finemap --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
    --attention attention/attention.parquet --genes genes.tsv --out results/ --workers 8
```

Detects eGenes and fine-maps them, once with the uniform prior and once with the PlantDNA-FM prior, so the two can be compared.

**Inputs.**

| option | file |
|---|---|
| `--tests` | the association table of step 2 |
| `--expression` | TSV (`.gz` allowed): first column gene ID, then one column per sample; the values that went into the association tests |
| `--genotypes` | one of: a bgzipped VCF with a tabix index (`bgzip file.vcf; tabix -p vcf file.vcf.gz`; biallelic SNVs with a `GT` field are used); a PLINK 1 fileset, given as the prefix of `.bed`/`.bim`/`.fam` (samples = IID; `.bed` and `.bim` may be gzipped); or a numpy `.npz` with arrays `samples`, `chrom`, `pos` and `dosage` (samples × variants; NaN or negative = missing) |
| `--attention` | `attention.parquet` from step 3 (or a TSV with `gene`, `pos`, `attention`); tests without a score get a neutral weight |
| `--genes` (optional) | `genes.tsv` from step 1; adds the gene name, strand and coordinates to the result table |

Expression table:

```
gene        108       265       403       430
AT2G01060   1.74637   1.67785   1.92481   2.0227
AT2G01070   1.08573   2.20241   1.02561   1.73898
```

Tested variants that are missing from the genotypes are left out of fine-mapping.

**Options.**

| option | default | meaning |
|---|---|---|
| `--lambda` | cross-fitted | fixed prior strength λ: one value, or one per fold (`1.35,1.65`) |
| `--folds` | chromosome halves | cross-fitting folds, `1,2,3:4,5` or a TSV with columns `chrom`, `fold` (0/1) |
| `--alpha` | 0.05 | family-wise error rate |
| `--n-tests` | rows in `--tests` | *M* of the Bonferroni threshold, e.g. the genome-wide number of tests when `--tests` holds a subset |
| `--single-gene` | off | candidate-gene mode for one gene (see below) |
| `--susie-all` | off | fine-map all tested variants of each eGene instead of the significant ones |
| `--L` | 10 | maximum number of SuSiE effects per gene |
| `--workers` | 4 | processes |

By default λ is learned from your data by cross-fitting between two groups of chromosomes, which needs eGenes on at least two chromosomes.

**Output: `results/finemapping.csv`.** One row per variant that is significant with the uniform prior or with the PlantDNA-FM prior, in the format of the table of the paper:

| column | content |
|---|---|
| `gene`, `gene_name`, `strand`, `gene_start`, `gene_end` | the gene (name, strand and coordinates only with `--genes`) |
| `rs` | variant, `chromosome_position` |
| `ref`, `alt` | alleles (from the association table if it has them, otherwise from the genotypes) |
| `maf` | minor allele frequency (from the association table) |
| `raw_pval (-log10)` | −log₁₀*P* of the association test |
| `n_tested_in_gene` | number of variants tested for the gene |
| `FM_prior_weight` | prior weight *w* of the variant (mean 1 within the gene) |
| `threshold_FM (-log10)` | weighted Bonferroni threshold of the variant, −log₁₀(α*w*/*M*) |
| `sig_raw`, `sig_FM` | variant significant with the uniform prior / with the PlantDNA-FM prior |
| `eGene_raw`, `eGene_FM` | gene is an eGene with the uniform prior / with the PlantDNA-FM prior |
| `FM_new_sig` | variant significant only with the PlantDNA-FM prior |
| `PIP_raw`, `PIP_FM` | posterior inclusion probability from SuSiE with each prior (empty if the variant is not significant with that prior) |

```
gene,gene_name,strand,gene_start,gene_end,rs,ref,alt,maf,raw_pval (-log10),n_tested_in_gene,FM_prior_weight,threshold_FM (-log10),sig_raw,eGene_raw,sig_FM,eGene_FM,FM_new_sig,PIP_raw,PIP_FM
AT2G01090,,+,79779,81710,2_80921,T,C,0.136,12.6677,275,5.5189,7.1924,True,True,True,True,False,0.55547,0.79387
```

Variants with a high `PIP_FM` are the candidate causal variants. Without `--attention`, the table has the uniform-prior columns only.

**Other files** in `results/`, for more detail:

| file | content |
|---|---|
| `summary.tsv` | one row per prior: λ, genes tested, tests, *M*, eGenes, eGenes with a signal, expGWAS signals, single-variant signals, variants with PIP > 0.9, mean signal size (variants per credible set) |
| `signals.tsv` | one row per expGWAS signal (95% credible set): `gene`, `prior`, `n_variants`, `purity`, lead variant (`lead_pos`, `lead_pip`, `lead_logp`) and all variant positions |
| `genes.tsv` | one row per tested gene: number of tests, lead variant, eGene status and number of signals with each prior |
| `pip.tsv.gz` | every fine-mapped variant and prior: `pip`, `cs` (credible set of the gene, −1 = none), `purity`, `n_samples`, `converged` |
| `tests.tsv.gz` | every test, significant or not: `logp`, `attention`, `fold`, prior `weight`, threshold and significance with each prior |
| `lambda.tsv`, `lambda_curve.tsv` | λ of each fold; number of eGenes in each fold for every λ of the grid (cross-fitting only) |

### Check your inputs

A wrong input file usually does not crash the analysis; it gives an empty or a wrong result. `attention` and `finemap` therefore check the inputs before they start, and `check` runs the same checks alone, in a few seconds:

```bash
plantdna-finemap check --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
    --attention attention/attention.parquet --genes genes.tsv --fasta genome.fa
```

Give the files you have; `--gpu` also checks the environment of the attention step. Each check prints one line:

```
  OK    association tests: 1,307,160 tests, 4,873 genes (median 243 variants per gene), chromosomes 2, 4, 5
  OK    genes of the association tests with expression: 100.0%
  OK    samples with expression and genotypes (471 of 471): 100.0%
  OK    tested variants found in the genotypes (20 genes): 94.8%
  OK    tests with an attention score: 100.0%
  OK    genes of the association tests in the gene table: 100.0%
  OK    tested variants within 100 kb of their gene: 100.0%
  OK    REF alleles of the tests that match the FASTA (2,000 variants): 100.0%
  check: all checks passed
```

`FAIL` stops the run: the files do not fit together. `WARN` lets it continue, but look at the line. What is checked:

| check | what a failure means |
|---|---|
| files and indexes exist (VCF `.tbi`, FASTA `.fai`, complete PLINK fileset) | a path is wrong or an index is missing; the message gives the command to make it |
| association table has `gene`, `chrom`, `pos` and `logp` or `p`; −log₁₀*P* is finite and not negative | wrong column names, or *P* values in a column named `logp` |
| genes of the tests have expression | gene IDs differ between the files, or the expression table is transposed |
| samples are shared by expression and genotypes | sample IDs differ between the files |
| chromosomes of the tests are in the genotypes; tested variants are found there | chromosome names or the assembly differ |
| tests have an attention score | the attention file was computed for other tests, gene IDs or another assembly |
| genes of the tests are in the gene table; variants lie near their gene | another annotation version or another assembly |
| chromosomes of the tests are in the FASTA; REF alleles match the FASTA (if the tests have a `ref` column) | the FASTA is not the assembly on which the variants were called |
| at least two chromosomes for cross-fitting | give a fixed `--lambda` for a single chromosome |

To compare two result tables, for example your rerun of the example with the table of the paper:

```bash
plantdna-finemap compare --result results/arabidopsis/finemapping.csv \
    --expected examples/expected/arabidopsis_chr2_chr4_finemapping_paper.csv.gz
```

The variants, the alleles and all calls must be identical; numbers may differ by 0.001 (`--tolerance`).

### Candidate gene (single-gene mode)

```bash
plantdna-finemap finemap --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
    --attention attention/attention.parquet --single-gene AT5G10140 --lambda 1.65 --out results_FLC/
```

Only the tests of one gene are used and *M* is its number of variants (Bonferroni over the gene's own variants). λ must be given, for example the cross-fitted λ of the gene's chromosome from a genome-wide run. `--susie-all` fine-maps all tested variants of the gene instead of only the significant ones. The outputs are the same files as above.

## Troubleshooting

| message or symptom | cause and fix |
|---|---|
| `plantdna-finemap: command not found` | the environment is not active: `conda activate plantdna-finemap` |
| `FAIL  samples with expression and genotypes: 0 of …` | sample IDs differ between the two files (the message shows three of each); rename the columns of the expression table or the samples of the VCF (`bcftools reheader -s`) |
| `FAIL  genes of the association tests with expression: 0.0%` | gene IDs differ, or the expression table has samples in rows: it must have genes in rows |
| `FAIL  genotypes: the index … .tbi is missing` | `bgzip file.vcf && tabix -p vcf file.vcf.gz`; a VCF compressed with plain `gzip` cannot be indexed |
| `FAIL  tests with an attention score: 0.0%` | the attention file belongs to other data; rerun `attention` with this association table, or use the precomputed file of the same panel |
| `FAIL  REF alleles of the tests that match the FASTA` | the FASTA is another assembly (or version) than the one of the variants |
| `finemapping.csv` has few rows or `PIP_FM` is empty | few variants pass the Bonferroni threshold; look at `summary.tsv`. `PIP_FM` is empty for variants that are significant only with the uniform prior |
| `attention` runs out of GPU memory | a 16-kb window needs about 12 GB; use a larger GPU, or `--device cpu` for a few genes |
| cross-fitted λ is 0 | the attention prior did not increase eGenes on the other fold of your data; the result then equals the uniform prior |

## Method

**Association test.** Not part of this package. In the paper, the *A. thaliana* tests are the GEMMA results of Voichek *et al.* (mixed model with a relatedness matrix, likelihood-ratio test; variants in gene body ± 10 kb with MAF ≥ 0.05 and call rate ≥ 0.5).

**Attention score.** For each gene, the region gene body ± 10 kb, extended to cover every tested variant of the gene, is tiled with 16-kb windows at a step of 8 kb. Each window is read in the direction of transcription and run through the model once. The query positions are the exon positions of the gene (all positions if the window contains no exon of the gene, or if the gene has no `<gene>.1` / `<gene>_T01` transcript); attention from these positions is averaged over all layers and heads (36 layers × 12 heads = 432 attention matrices), giving one score per position. The score of a variant is the mean over the windows that contain an exon of the gene, otherwise over the other windows.

**Prior weights.** *q* is the attention score rank-normalised within each gene (ties get the average rank; variants without a score get *q* = 0, a neutral weight). The weight is *w* = exp(λ*q*), scaled to a mean of 1 within each gene.

**λ by two-fold cross-fitting.** The chromosomes are split into two folds (default: in chromosome order, first half and second half, e.g. *A. thaliana* 1–3 / 4–5 and maize 1–5 / 6–10; or give your own with `--folds`). For each fold, λ is chosen from 0, 0.05, …, 4 to maximise the number of eGenes on the *other* fold and is then applied to the tests of this fold, so no test's weight depends on its own data. If several λ give the same maximum, the median of the tied values is used (the upper middle value for an even number of ties). A fixed λ can be given with `--lambda`.

**eGene detection.** Uniform prior: −log₁₀*P* > −log₁₀(α/*M*). Attention prior: −log₁₀*P* > −log₁₀(α*w*/*M*), the weighted Bonferroni procedure, with *M* the number of tests (all variant–gene pairs) and α = 0.05. Because the weights average to 1 within each gene, the family-wise error rate stays at α.

**Fine-mapping.** For each eGene, SuSiE (*L* = 10) runs on its significant variants, with mean-imputed genotypes and the expression values as given. The prior inclusion probabilities are π = 0.8*w*/Σ*w* + 0.2/*m* over these *m* variants (uniform prior: 1/*m*). Signals are 95% credible sets with purity ≥ 0.5. SuSiE is a numpy port of `susieR::susie` with its defaults.

## Results in the paper

| | uniform prior | PlantDNA-FM prior | cross-fitted λ (fold 1 / fold 2) |
|---|---|---|---|
| *A. thaliana*, batch B2, 471 accessions: eGenes | 4,283 | 4,479 | 1.35 / 1.65 (chromosomes 1–3 / 4–5) |
| *A. thaliana*: expGWAS signals (95% credible sets) | 6,810 | 7,009 | |
| maize, 263 lines: eGenes | 8,557 | 8,758 | 1.25 / 1.85 (chromosomes 1–5 / 6–10) |
| maize: expGWAS signals (95% credible sets) | 14,424 | 14,824 | |

## Precomputed attention

`attention/` holds the PlantDNA-FM scores used in the paper, so `finemap` can run without a GPU on these panels:

| file | genome | genes | variant–gene pairs | tests in the paper |
|---|---|---|---|---|
| `arabidopsis_TAIR10_PlantDNA-FM.parquet` | TAIR10 | 16,053 | 4,297,315 | GEMMA tests of Voichek *et al.*, batch B2 (gene body ± 10 kb) |
| `maize_AGPv3_PlantDNA-FM.parquet` | B73 AGPv3 | 17,974 | 8,881,040 | HapMap 3.2.1 SNVs in gene body ± 10 kb |

Tests that are not in the file get a neutral weight, and `finemap` reports the coverage. For another panel, compute the scores with `attention`.

## Reproducibility

Checked against the analysis in the paper:

| step | check | result |
|---|---|---|
| `finemap` | the example (chromosomes 2 and 4) with the paper's settings (`--n-tests 4297315 --lambda 1.35,1.65`) | `finemapping.csv` has the same 52,683 variants as the table of the paper for these chromosomes, with identical weights, thresholds, significance and eGene calls (1,314 / 1,365 eGenes, 2,114 / 2,158 expGWAS signals with the uniform / PlantDNA-FM prior); PIPs within 4 × 10⁻⁵ |
| `finemap --single-gene` | *FLC* | the same credible sets as the candidate-gene analysis for the paper (7 / 6 variants) |
| `attention` | all 4,873 example genes (1,307,160 variant–gene pairs, 10,544 windows) recomputed with this code and the weights downloaded from Hugging Face, on an RTX 5090 (3 h 58 min) | largest relative difference to the paper's scores 2.6 × 10⁻⁷ |
| installation | new conda environment from `environment.yml`, `examples/run_example.sh` | all input checks pass and `compare` reports PASS against the table of the paper |

## Citation

Yu H., Korte A., Voichek Y., Menon G., Sun W., Dean C., Nordborg M., Zhang H. & Ding Y. DNA foundation model priors improve fine-mapping in plant populations (submitted).

SuSiE: Wang G., Sarkar A., Carbonetto P. & Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J. R. Stat. Soc. B* 82, 1273–1300 (2020).

Cross-weighting: Ignatiadis N. & Huber W. Covariate powered cross-weighted multiple testing. *J. R. Stat. Soc. B* 83, 720–751 (2021).

## License

MIT (see `LICENSE`). The example data come from published sources and keep their original terms; see [examples/README.md](examples/README.md).
