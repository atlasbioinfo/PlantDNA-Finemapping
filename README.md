# PlantDNA-Finemapping

expGWAS fine-mapping in plant populations with PlantDNA-FM attention priors: association tests, weighted Bonferroni eGene detection and SuSiE fine-mapping, starting from a VCF and an expression table.

Association panels in plants are small and highly structured, so an expression GWAS (expGWAS) signal usually points to a block of variants in near-perfect linkage disequilibrium rather than to one variant. This code uses attention from [PlantDNA-FM](https://huggingface.co/ATLASBIOINFO/PlantDNA-FM), a plant DNA foundation model with a 16-kb context, as a prior that ranks the variants within each gene, and applies the same prior to eGene detection and to fine-mapping.

Definitions used throughout:

- **expGWAS**: association of the expression of each gene with the variants in and around it (gene body ± 10 kb).
- **eGene**: a gene with at least one significant variant.
- **expGWAS signal**: a 95% credible set from SuSiE fine-mapping of an eGene. A gene can have several signals.

## Contents

- [Installation](#installation)
- [Quick start: the bundled example](#quick-start-the-bundled-example-cpu-about-5-min)
- [User guide: your own data](#user-guide-your-own-data)
  - [What you need](#what-you-need)
  - [Step 1: `genes`](#step-1-genes--gene-table-from-the-annotation)
  - [Step 2: `assoc`](#step-2-assoc--association-tests)
  - [Step 3: `attention`](#step-3-attention--plantdna-fm-attention-scores-gpu)
  - [Step 4: `finemap`](#step-4-finemap--egenes-and-expgwas-signals)
  - [Candidate gene](#candidate-gene-single-gene-mode)
- [Method](#method)
- [Results in the paper](#results-in-the-paper)
- [Precomputed attention](#precomputed-attention)
- [Reproducibility](#reproducibility)

## Installation

```bash
git clone https://github.com/atlasbioinfo/PlantDNA-Finemapping.git
cd PlantDNA-Finemapping
pip install .                  # CPU steps: genes, assoc, finemap
pip install ".[attention]"     # also the GPU step (PyTorch, transformers)
```

Requires Python ≥ 3.9 and numpy, scipy, pandas, pyarrow and pysam. The attention step was tested with PyTorch 2.10 and transformers 5.9 on an RTX 5090 (about 1.4 s per 16-kb window). The model is downloaded from Hugging Face ([`ATLASBIOINFO/PlantDNA-FM`](https://huggingface.co/ATLASBIOINFO/PlantDNA-FM)) the first time `attention` runs.

## Quick start: the bundled example (CPU, about 5 min)

```bash
bash examples/run_example.sh 8      # 8 worker processes
```

The example is *A. thaliana* batch B2 (471 accessions) for the 4,872 tested genes of chromosomes 2 and 4 plus *FLC*. The script runs

1. `finemap` with the association tests and the settings of the paper; the output is the paper's result for these two chromosomes (1,314 eGenes with the uniform prior and 1,365 with the PlantDNA-FM prior);
2. `finemap` in single-gene mode for *FLC*;
3. the whole workflow from the genotypes and the expression table (`genes`, `assoc`, `finemap`), as you would run it on your own data.

Results go to `results/`. See [examples/README.md](examples/README.md) for the files, the expected results and for recomputing the attention on a GPU.

## User guide: your own data

### What you need

| you provide | file | used by |
|---|---|---|
| reference genome | FASTA with a `.fai` index (plain or bgzipped) | `attention` |
| genome annotation | GFF3 (`.gz` allowed) | `genes` |
| genotypes | VCF (bgzipped, with a tabix index) or a PLINK 1 fileset | `assoc`, `finemap` |
| expression | table of genes × samples (TSV, `.gz` allowed) | `assoc`, `finemap` |

The four commands run in this order, each one reading the output of the previous ones:

```bash
plantdna-finemap genes     --gff genes.gff3.gz --out genes.tsv
plantdna-finemap assoc     --genotypes genotypes.vcf.gz --expression expression.tsv.gz --genes genes.tsv \
                           --out associations.tsv.gz --workers 8
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA-FM --genes genes.tsv --fasta genome.fa \
                           --tests associations.tsv.gz --out attention/
plantdna-finemap finemap   --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
                           --attention attention/attention.parquet --out results/ --workers 8
```

`genes`, `assoc` and `finemap` run on CPUs; `attention` needs a GPU. Without `--attention`, `finemap` reports the uniform prior only, so the three CPU steps already give a standard expGWAS with fine-mapping.

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
| `exons` | exons of the representative transcript (`<gene>.1` or `<gene>_T01`), `start-end,start-end,…`; empty if the gene has no such transcript |

### Step 2: `assoc` — association tests

```bash
plantdna-finemap assoc --genotypes genotypes.vcf.gz --expression expression.tsv.gz --genes genes.tsv \
    --out associations.tsv.gz --workers 8
```

For each gene, the variants in gene body ± 10 kb are tested against the expression of the gene in a linear mixed model that accounts for population structure ([Method](#method)).

**Inputs.**

| option | file |
|---|---|
| `--genotypes` | one of: a bgzipped VCF with a tabix index (`bgzip file.vcf; tabix -p vcf file.vcf.gz`; biallelic SNVs with a `GT` field are used); a PLINK 1 fileset, given as the prefix of `.bed`/`.bim`/`.fam` (samples = IID; `.bed` and `.bim` may be gzipped); or a numpy `.npz` with arrays `samples`, `chrom`, `pos` and `dosage` (samples × variants; NaN or negative = missing) |
| `--expression` | TSV: first column gene ID, then one column per sample |
| `--genes` | output of step 1 |
| `--kinship` (optional) | relatedness matrix as a TSV: square, sample IDs in the first column and in the header. By default the matrix is computed from all variants of `--genotypes`. |

Expression table:

```
gene        108       265       403       430
AT2G01060   1.74637   1.67785   1.92481   2.0227
AT2G01070   1.08573   2.20241   1.02561   1.73898
```

**Options.** `--flank` (10000, nt on either side of the gene body), `--maf` (0.05) and `--min-call` (0.5): minimum minor allele frequency and call rate of a variant among the samples with expression of the gene; `--min-samples` (30): genes with expression in fewer samples are skipped; `--kinship-every` (1): use every *n*-th variant for the relatedness matrix, for very large VCFs; `--workers` (4).

The relatedness matrix should reflect the whole genome: give `--genotypes` a genome-wide file, or pass your own matrix with `--kinship` (for example one from GEMMA or PLINK, with sample IDs added).

**Output.** `associations.tsv.gz`, one row per variant–gene test:

| column | content |
|---|---|
| `gene`, `chrom`, `pos` | gene, chromosome and position of the variant |
| `maf` | minor allele frequency among the samples used for this gene |
| `logp` | −log₁₀*P* |
| `h2` | heritability explained by relatedness (REML estimate for the gene) |
| `n` | number of samples |

**Using another GWAS tool instead.** Steps 3 and 4 take any table (TSV or parquet) with one row per variant–gene test and the columns `gene`, `chrom` (or `chr`), `pos` (or `ps`), and `logp` (−log₁₀*P*) or `p` (also `p_wald`, `p_lrt`, `p_score`). In the paper, the *A. thaliana* tests are the GEMMA results of Voichek *et al.*

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
    --attention attention/attention.parquet --out results/ --workers 8
```

Detects eGenes and fine-maps them, once with the uniform prior and once with the PlantDNA-FM prior, so the two can be compared.

**Inputs.** `--tests`, `--expression` and `--genotypes` as in step 2; `--attention`: `attention.parquet` from step 3 (or a TSV with `gene`, `pos`, `attention`). Tested variants that are missing from the genotypes are left out of fine-mapping; tests without an attention score get a neutral weight.

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

**Outputs** in `results/`:

| file | content |
|---|---|
| `summary.tsv` | one row per prior: λ, genes tested, tests, *M*, eGenes, eGenes with a signal, expGWAS signals, single-variant signals, variants with PIP > 0.9, mean signal size (variants per credible set) |
| `genes.tsv` | one row per gene: number of tests, lead variant (`lead_pos`, `lead_logp`), and eGene status and number of signals with each prior |
| `signals.tsv` | one row per expGWAS signal: `gene`, `prior`, `signal` (ID), `n_variants`, `purity`, lead variant (`lead_pos`, `lead_pip`, `lead_logp`) and all variant positions (`variants`) |
| `pip.tsv.gz` | one row per fine-mapped variant and prior: `pip` (posterior inclusion probability), `cs` (credible set of the gene, −1 = none), `purity`, `n_samples`, `converged` |
| `tests.tsv.gz` | every test: `logp`, `maf`, `attention`, `fold`, prior `weight`, and the threshold and significance with each prior |
| `lambda.tsv`, `lambda_curve.tsv` | λ of each fold; number of eGenes in each fold for every λ of the grid (cross-fitting only) |

Start with `summary.tsv` for the overall numbers and with `signals.tsv` for the candidate causal variants: a signal with `n_variants` = 1 is resolved to a single variant, and `pip.tsv.gz` ranks the variants of the others.

### Candidate gene (single-gene mode)

```bash
plantdna-finemap finemap --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes.vcf.gz \
    --attention attention/attention.parquet --single-gene AT5G10140 --lambda 1.65 --out results_FLC/
```

Only the tests of one gene are used and *M* is its number of variants (Bonferroni over the gene's own variants). λ must be given, for example the cross-fitted λ of the gene's chromosome from a genome-wide run. `--susie-all` fine-maps all tested variants of the gene instead of only the significant ones. The outputs are the same files as above.

## Method

**Association test.** For each gene, `assoc` tests the variants in gene body ± 10 kb (MAF ≥ 0.05 and call rate ≥ 0.5 among the samples with expression of that gene) in a linear mixed model, *y* = *x*β + *g* + *e* with *g* ~ N(0, σ²*h*²**K**): *h*² is estimated per gene by REML, genotypes and expression are whitened, and each variant gets a *t* test with *n* − 2 degrees of freedom. **K** is the centred relatedness matrix.

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
| `finemap` | the example (chromosomes 2 and 4) with the paper's settings (`--n-tests 4297315 --lambda 1.35,1.65`) | the paper's results for these chromosomes: 1,314 / 1,365 eGenes and 2,114 / 2,158 expGWAS signals with the uniform / PlantDNA-FM prior, the same credible sets; PIPs within 4 × 10⁻⁵ |
| `finemap --single-gene` | *FLC* | the same credible sets as the candidate-gene analysis for the paper (7 / 6 variants) |
| `attention` | 892 genes (248,022 variant–gene pairs, 1,915 windows) recomputed with this code on an RTX 5090; 8 example genes recomputed with the weights downloaded from Hugging Face | largest relative difference to the paper's scores 2.4 × 10⁻⁷ in both checks |
| `assoc` | the example, against the GEMMA tests of Voichek *et al.* | −log₁₀*P* correlated at *r* = 0.997 with the genome-wide relatedness matrix of the example; VCF and PLINK input give identical tests |
| installation | clean environment (Python 3.11), `pip install .`, `examples/run_example.sh` | the numbers in [examples/README.md](examples/README.md) |

## Citation

Yu H., Korte A., Voichek Y., Menon G., Sun W., Dean C., Nordborg M., Zhang H. & Ding Y. DNA foundation model priors improve fine-mapping in plant populations (in preparation).

SuSiE: Wang G., Sarkar A., Carbonetto P. & Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J. R. Stat. Soc. B* 82, 1273–1300 (2020).

Cross-weighting: Ignatiadis N. & Huber W. Covariate powered cross-weighted multiple testing. *J. R. Stat. Soc. B* 83, 720–751 (2021).

## License

MIT (see `LICENSE`). The example data come from published sources and keep their original terms; see [examples/README.md](examples/README.md).
