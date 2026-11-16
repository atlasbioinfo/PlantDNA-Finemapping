# PlantDNA-Finemapping

expGWAS fine-mapping in plant populations with PlantDNA16K-FM attention priors: weighted Bonferroni eGene detection and SuSiE fine-mapping.

Association panels in plants are small and highly structured, so an expression GWAS (expGWAS) signal usually points to a block of variants in near-perfect linkage disequilibrium rather than to one variant. This code uses attention from PlantDNA16K-FM, a plant DNA foundation model with a 16-kb context, as a prior that ranks the variants within each gene, and applies the same prior to eGene detection and to fine-mapping.

Definitions used throughout:

- **expGWAS**: association of the expression of each gene with the variants in and around it (gene body ± 10 kb).
- **eGene**: a gene with at least one significant variant.
- **expGWAS signal**: a 95% credible set from SuSiE fine-mapping of an eGene. A gene can have several signals.

## Method

| step | command | hardware |
|---|---|---|
| 1. Gene table (gene body, strand, exons of a representative transcript) from a GFF3 | `genes` | CPU |
| 2. PlantDNA16K-FM attention score per tested variant | `attention` | GPU (or use the [precomputed scores](#precomputed-attention)) |
| 3. Prior weights, λ by cross-fitting, eGene detection and SuSiE fine-mapping, under the uniform and the attention prior | `finemap` | CPU |

The association test itself is not part of this package: `finemap` takes one row per variant–gene test from any GWAS tool, for example GEMMA run per gene on the variants in gene body ± 10 kb, as in Voichek *et al.*

**Attention score.** For each gene, the region gene body ± 10 kb, extended to cover every tested variant of the gene, is tiled with 16-kb windows at a step of 8 kb. Each window is read in the direction of transcription and run through the model once. The query positions are the exon positions of the gene (all positions if the window contains no exon of the gene, or if the gene has no `<gene>.1` / `<gene>_T01` transcript); attention from these positions is averaged over all layers and heads (36 layers × 12 heads = 432 attention matrices), giving one score per position. The score of a variant is the mean over the windows that contain an exon of the gene, otherwise over the other windows.

**Prior weights.** *q* is the attention score rank-normalised within each gene (ties get the average rank; variants without a score get *q* = 0, a neutral weight). The weight is *w* = exp(λ*q*), scaled to a mean of 1 within each gene.

**λ by two-fold cross-fitting.** The chromosomes are split into two folds (default: in chromosome order, first half and second half, e.g. *A. thaliana* 1–3 / 4–5 and maize 1–5 / 6–10; or give your own with `--folds`). For each fold, λ is chosen from 0, 0.05, …, 4 to maximise the number of eGenes on the *other* fold and is then applied to the tests of this fold, so no test's weight depends on its own data. If several λ give the same maximum, the median of the tied values is used (the upper middle value for an even number of ties). A fixed λ can be given with `--lambda`.

**eGene detection.** Uniform prior: −log₁₀*P* > −log₁₀(α/*M*). Attention prior: −log₁₀*P* > −log₁₀(α*w*/*M*), the weighted Bonferroni procedure, with *M* the number of tests (all variant–gene pairs) and α = 0.05. Because the weights average to 1 within each gene, the family-wise error rate stays at α.

**Fine-mapping.** For each eGene, SuSiE (*L* = 10) runs on its significant variants, with mean-imputed genotypes and the expression values as given. The prior inclusion probabilities are π = 0.8*w*/Σ*w* + 0.2/*m* over these *m* variants (uniform prior: 1/*m*). Signals are 95% credible sets with purity ≥ 0.5. SuSiE is a numpy port of `susieR::susie` with its defaults.

**Single-gene mode** (candidate-gene analysis, `--single-gene`): only the tests of one gene are used and *M* is its number of variants (Bonferroni over the gene's own variants), with a fixed λ, for example the cross-fitted λ of the gene's chromosome from a genome-wide run. `--susie-all` fine-maps all tested variants of the gene instead of only the significant ones.

## Results in the paper

| | uniform prior | PlantDNA16K-FM prior | cross-fitted λ (fold 1 / fold 2) |
|---|---|---|---|
| *A. thaliana*, batch B2, 471 accessions: eGenes | 4,283 | 4,479 | 1.35 / 1.65 (chromosomes 1–3 / 4–5) |
| *A. thaliana*: expGWAS signals (95% credible sets) | 6,810 | 7,009 | |
| maize, 263 lines: eGenes | 8,557 | 8,758 | 1.25 / 1.85 (chromosomes 1–5 / 6–10) |
| maize: expGWAS signals (95% credible sets) | 14,424 | 14,824 | |

## Installation

```bash
git clone https://github.com/atlasbioinfo/PlantDNA-Finemapping.git
cd PlantDNA-Finemapping
pip install .                  # CPU steps: genes, finemap
pip install ".[attention]"     # also the GPU step (PyTorch, transformers)
```

Requires Python ≥ 3.9 and numpy, scipy, pandas, pyarrow and pysam. The attention step was tested with PyTorch 2.10 and transformers 5.9 on an RTX 5090 (about 1.4 s per 16-kb window).

**Model.** `attention --model` takes a local directory or a Hugging Face ID. The commands below use `ATLASBIOINFO/PlantDNA16K-FM`; this ID is a placeholder until the 16K weights are released on Hugging Face, so please check the final ID there. The 8-kb model of the same family is [`ATLASBIOINFO/PlantDNA-FM`](https://huggingface.co/ATLASBIOINFO/PlantDNA-FM); with it, use `--window 8000 --step 4000`.

## Quick start: the bundled example (CPU, about 1 min)

```bash
bash examples/run_example.sh 8      # 8 worker processes
```

The example is *A. thaliana* batch B2 (471 accessions) for the 892 genes in the first megabase of each chromosome plus *FLC*, with the GEMMA tests of Voichek *et al.* and the precomputed attention. The script runs `finemap` genome-wide (λ cross-fitted over chromosomes 1–3 / 4–5) and in single-gene mode for *FLC*; results go to `results/arabidopsis/` and `results/arabidopsis_FLC/`. See [examples/README.md](examples/README.md) for the expected results and for recomputing the attention on a GPU.

## Your own data

```bash
# 1. gene table from the genome annotation
plantdna-finemap genes --gff genes.gff3.gz --out genes.tsv

# 2. attention scores for every tested variant (GPU); rerunning the same command resumes
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA16K-FM --genes genes.tsv --fasta genome.fa \
    --tests associations.tsv.gz --out attention/

# 3. eGenes and expGWAS signals under the uniform and the attention prior
plantdna-finemap finemap --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes \
    --attention attention/attention.parquet --out results/ --workers 8

# candidate gene: Bonferroni over the gene's own variants, fixed lambda
plantdna-finemap finemap --tests associations.tsv.gz --expression expression.tsv.gz --genotypes genotypes \
    --attention attention/attention.parquet --single-gene AT5G10140 --lambda 1.65 --out results_FLC/
```

Without `--attention`, `finemap` reports the uniform prior only.

### Inputs

| input | format |
|---|---|
| association tests (`--tests`) | TSV (optionally gzipped) or parquet, one row per variant–gene test: `gene`, `chrom` (or `chr`), `pos` (or `ps`), and `logp` (−log₁₀*P*) or `p` (also `p_wald`, `p_lrt`, `p_score`). Every row counts towards *M* (or set *M* with `--n-tests`). |
| expression (`--expression`) | TSV: first column gene ID, then one column per accession; missing values allowed. Use the values that went into the association test. |
| genotypes (`--genotypes`) | PLINK 1 fileset (prefix of `.bed`/`.bim`/`.fam`; samples = IID), a bgzipped VCF with a tabix index (GT, biallelic SNVs only) or a numpy `.npz` with arrays `samples`, `chrom`, `pos` and `dosage` (samples × variants; NaN or negative = missing). Variants are matched to the tests by chromosome and position; tested variants missing from the genotypes are left out of fine-mapping. |
| annotation (`genes --gff`) | GFF3 with `gene`, `mRNA`/`transcript` and `exon` features (Ensembl style: `ID=gene:…`, `Parent=transcript:…`). Genes with `biotype=protein_coding` are kept (`--all-biotypes` keeps all). The query exons are those of the transcript `<gene>.1` or `<gene>_T01`; genes without such a transcript use all positions as queries in every window. |
| reference (`attention --fasta`) | FASTA with a `.fai` index (plain or bgzipped) |
| attention (`finemap --attention`) | parquet or TSV with columns `gene`, `pos`, `attention` (output of `attention`) |

Chromosome names may differ by a `chr`/`Chr` prefix between files.

### Main options of `finemap`

| option | default | meaning |
|---|---|---|
| `--lambda` | cross-fitted | fixed λ: one value, or one per fold (`1.35,1.65`) |
| `--folds` | chromosome halves | cross-fitting folds, `1,2,3:4,5` or a TSV with columns `chrom`, `fold` (0/1) |
| `--alpha` | 0.05 | family-wise error rate |
| `--n-tests` | rows in `--tests` | *M*, e.g. the genome-wide number of tests when `--tests` holds a subset |
| `--single-gene` | off | candidate-gene mode for one gene (*M* = its number of tests; needs `--lambda`) |
| `--susie-all` | off | fine-map all tested variants of each eGene instead of the significant ones |
| `--L` | 10 | maximum number of SuSiE effects per gene |
| `--workers` | 4 | processes |

Options of `attention`: `--window` (16000), `--step` (8000), `--flank` (10000), `--only` (file of gene IDs), `--device` (cuda).

### Outputs of `finemap`

| file | content |
|---|---|
| `summary.tsv` | per prior: λ, genes tested, tests, *M*, eGenes, eGenes with a signal, expGWAS signals, single-variant signals, variants with PIP > 0.9, mean signal size (variants per credible set) |
| `genes.tsv` | per gene: tests, lead variant, eGene status and number of signals under each prior |
| `signals.tsv` | one row per expGWAS signal: gene, prior, number of variants, purity, lead variant (position, PIP, −log₁₀*P*), all variant positions |
| `pip.tsv.gz` | PIP and credible set (`cs`, −1 = none) of every fine-mapped variant, per prior |
| `tests.tsv.gz` | every test: −log₁₀*P*, attention, fold, prior weight, threshold and significance under each prior |
| `lambda.tsv`, `lambda_curve.tsv` | λ per fold; eGenes per fold for every λ of the grid (cross-fitting only) |

### Outputs of `attention`

`attention.parquet` (`gene`, `pos`, `attention`, `attention_source` = `exon windows` or `windows without exon`) and `chunks/` with the score of every variant in every window.

## Precomputed attention

`attention/` holds the PlantDNA16K-FM scores used in the paper, so `finemap` can run without a GPU on these panels:

| file | genome | genes | variant–gene pairs | tests in the paper |
|---|---|---|---|---|
| `arabidopsis_TAIR10_PlantDNA16K-FM.parquet` | TAIR10 | 16,053 | 4,297,315 | GEMMA tests of Voichek *et al.*, batch B2 (gene body ± 10 kb) |
| `maize_AGPv3_PlantDNA16K-FM.parquet` | B73 AGPv3 | 17,974 | 8,881,040 | HapMap 3.2.1 SNVs in gene body ± 10 kb |

Tests that are not in the file get a neutral weight, and `finemap` reports the coverage. For another panel, compute the scores with `attention`.

## Reproducibility

Checked against the analysis in the paper:

| step | check | result |
|---|---|---|
| `attention` | the 892 example genes (248,022 variant–gene pairs, 1,915 windows), recomputed with this code on an RTX 5090 (44 min) | largest relative difference to the paper's scores 2.4 × 10⁻⁷ (AT4G00520, the one example gene without a `.1` transcript, all positions as queries: 2.2 × 10⁻⁷) |
| `finemap` | the 892 example genes with the paper's settings (`--n-tests 4297320 --lambda 1.35,1.65`) | the same eGenes (259 / 270 under the uniform / PlantDNA16K-FM prior) and the same 417 / 421 credible sets; PIPs within 4 × 10⁻⁵ |
| `finemap --single-gene` | *FLC* | the same credible sets as the candidate-gene analysis for the paper at λ = 1.6 (7 / 6 variants; with `--susie-all` 4, 3, 2 / 4, 2, 2); λ = 1.65 gives the same sets |

## Citation

Yu H., Korte A., Voichek Y., Menon G., Sun W., Dean C., Nordborg M., Zhang H. & Ding Y. DNA foundation model priors improve fine-mapping in plant populations (in preparation).

SuSiE: Wang G., Sarkar A., Carbonetto P. & Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J. R. Stat. Soc. B* 82, 1273–1300 (2020).

Cross-weighting: Ignatiadis N. & Huber W. Covariate powered cross-weighted multiple testing. *J. R. Stat. Soc. B* 83, 720–751 (2021).

## License

MIT (see `LICENSE`). The example data come from published sources and keep their original terms; see [examples/README.md](examples/README.md).
