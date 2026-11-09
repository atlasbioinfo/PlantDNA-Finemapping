# PlantDNA-Finemapping

Fine-mapping cis-eQTLs in plant populations with PlantDNA-FM attention priors: weighted Bonferroni eGene detection and SuSiE fine-mapping.

Association panels in plants are small and highly structured, so a cis-eQTL usually points to a block of variants in near-perfect linkage disequilibrium rather than to one variant. This pipeline uses attention from PlantDNA-FM, a DNA foundation model pretrained on 260 plant species, as a prior that ranks variants within each gene, and applies the same prior to eGene detection and to fine-mapping.

## What it does

| step | command | hardware |
|---|---|---|
| 1. 8-kb window per gene (TSS −3 kb … +5 kb, direction of transcription) and the exons of one representative transcript | `windows` | CPU |
| 2. Relatedness matrix from genome-wide SNVs | `kinship` | CPU |
| 3. PlantDNA-FM attention score per variant | `attention` | GPU (or use the [precomputed scores](#precomputed-attention)) |
| 4. Association, eGene detection and SuSiE fine-mapping, under the uniform prior and the attention prior | `map` | CPU |

Definitions used throughout:

- **eGene**: a gene with at least one significant cis association.
- **eQTL**: a 95% credible set from SuSiE fine-mapping of an eGene. A gene can have several eQTLs.
- **single-variant eQTL**: an eQTL that contains one variant.

**Attention score.** Each window is read in the direction of transcription and run through the model once. In every layer, attention from the exon positions is averaged over heads and query positions. The score of a variant is the attention to its position, averaged over all layers (36 layers × 12 heads = 432 attention matrices).

**Association.** For each gene, the pipeline uses the accessions with expression and the SNVs in the window with MAF ≥ 0.05 and call rate ≥ 0.5 in those accessions. Missing genotypes are mean-imputed. Expression and genotypes are whitened with the kinship matrix (h² by REML for each gene), and each variant is tested by marginal regression. Every tested variant–gene pair counts towards *M*.

**Prior weights.** The attention score is rank-normalised within each gene (*q*). The prior weight is *w* = exp(λ*q*), scaled to a mean of 1 within each gene; the default is λ = 1. Variants without a score get *q* = 0, which is a neutral weight.

**eGene detection.** Uniform prior: −log₁₀*P* > −log₁₀(α/*M*). Attention prior: −log₁₀*P* > −log₁₀(α*w*/*M*), the weighted Bonferroni procedure. Because the weights average to 1 within each gene, every gene keeps its share of α and the family-wise error rate stays at α (default 0.05). An uninformative prior can lower power but cannot create false positives.

**Fine-mapping.** SuSiE with *L* = 10 runs on the same whitened data. The prior inclusion probabilities are 0.8*w*/Σ*w* + 0.2/*m* for a gene with *m* variants, so no variant falls below 0.2/*m*; under the uniform prior every variant gets 1/*m*. Credible sets are 95% sets with purity ≥ 0.5.

## Installation

```bash
git clone https://github.com/atlasbioinfo/PlantDNA-Finemapping.git
cd PlantDNA-Finemapping
pip install .                  # CPU steps: windows, kinship, map
pip install ".[attention]"     # also the GPU step (PyTorch, transformers)
```

Requires Python ≥ 3.9 and numpy, scipy, pandas, pyarrow and pysam. The attention step was tested with PyTorch 2.10 and transformers 5.9 on an RTX 5090.

## Quick start: run an example (CPU, no GPU needed)

```bash
bash examples/run_example.sh arabidopsis 8    # A. thaliana chromosome 5, 471 accessions
bash examples/run_example.sh maize 8          # maize chromosome 10, 263 accessions
```

Results are written to `results/<species>/`. See [examples/README.md](examples/README.md) for the data, the expected results and how to recompute the attention on a GPU.

## Your own data

```bash
# 1. windows from the genome annotation
plantdna-finemap windows --gff genes.gff3.gz --out windows.tsv

# 2. kinship from genome-wide SNVs of the accessions with expression
plantdna-finemap kinship --vcf genome.vcf.gz --expression expression.tsv --out kinship.tsv --workers 8

# 3. attention scores (GPU); rerunning the same command resumes
plantdna-finemap attention --model <PlantDNA-FM directory or Hugging Face ID> \
    --windows windows.tsv --fasta genome.fa --vcf genome.vcf.gz --out attention/

# 4. eGenes and eQTLs under the uniform and the attention prior
plantdna-finemap map --vcf genome.vcf.gz --expression expression.tsv --windows windows.tsv \
    --kinship kinship.tsv --attention attention/attention.parquet --out results/ --workers 8
```

Without `--attention`, `map` reports the uniform prior only.

### Inputs

| input | format |
|---|---|
| genotypes | bgzipped VCF with a tabix index (`.tbi`). GT is read; only biallelic SNVs (REF and ALT in A/C/G/T) are used; the dosage is the number of ALT alleles; haploid calls count as homozygous. |
| expression | TSV: first column gene ID (as in the GFF3), then one column per accession (VCF sample names); missing values allowed. Normalise expression before use, for example with a rank-based inverse normal transformation per gene. |
| annotation | GFF3 with `gene`, `mRNA`/`transcript` and `exon` features (Ensembl style: `ID=gene:…`, `Parent=transcript:…`). Genes with `biotype=protein_coding` are kept (`--all-biotypes` keeps all). The representative transcript is `<gene>.1` or `<gene>_T01`, otherwise the one with the longest total exon length. |
| reference | FASTA with a `.fai` index (plain or bgzipped); needed only for `attention` |
| kinship | TSV with accessions as row and column names (output of `kinship`) |
| attention | parquet or TSV with columns `gene`, `pos`, `attention` |

Chromosome names may differ by a `chr`/`Chr` prefix between files.

### Main options of `map`

| option | default | meaning |
|---|---|---|
| `--lambda` | 1 | prior strength λ |
| `--alpha` | 0.05 | family-wise error rate |
| `--maf` | 0.05 | minimum MAF in the accessions with expression, per gene |
| `--min-call` | 0.5 | minimum call rate in the accessions with expression, per gene |
| `--min-samples` | 30 | genes with fewer accessions with expression are not tested |
| `--scored-only` | off | test only variants with an attention score (otherwise unscored variants are tested with a neutral weight) |
| `--L` | 10 | maximum number of SuSiE effects per gene |
| `--genes` | all | file with one gene ID per line |
| `--workers` | 4 | processes |

### Outputs of `map`

| file | content |
|---|---|
| `summary.tsv` | per prior: genes tested, tests (*M*), eGenes, eGenes with an eQTL, eQTLs, single-variant eQTLs, variants with PIP > 0.9, mean eQTL size (variants per credible set) |
| `genes.tsv` | per gene: accessions, h², variants tested, lead variant, eGene status and number of eQTLs under each prior |
| `eqtls.tsv` | one row per eQTL: gene, prior, number of variants, purity, lead variant (position, PIP, −log₁₀*P*), all variant positions |
| `pip.tsv.gz` | PIP and credible set of every variant in every eGene, per prior |
| `associations.tsv.gz` | every variant–gene test: MAF, −log₁₀*P*, attention, prior weight and the threshold under each prior |

## Precomputed attention

`attention/` holds the attention scores used in the paper, so `map` can run without a GPU:

| file | genome | genes | variants |
|---|---|---|---|
| `arabidopsis_TAIR10_PlantDNA-FM.parquet` | TAIR10 | 15,910 | 1,399,706 |
| `maize_AGPv3_PlantDNA-FM.parquet` | B73 AGPv3 | 17,967 | 3,556,844 |

They cover the variants tested in the paper: 1001 Genomes SNVs in *A. thaliana* and HapMap 3.2.1 SNVs in maize. Variants of your panel that are not in the file get a neutral weight, and `map` reports the coverage. For another panel, compute the scores with `attention`.

## Reproducibility

The code was checked against the analysis in the paper:

| step | check | result |
|---|---|---|
| `windows` | TSS and window of 22,623 *A. thaliana* and 39,148 maize genes | identical; exons identical except 28 *A. thaliana* genes without a `.1` transcript |
| `kinship` | *A. thaliana*, 471 accessions, 1001 Genomes VCF | identical (755,957 SNVs; largest difference 0) |
| `attention` | two runs of 20 *A. thaliana* genes (3,737 variants), one from the example files | largest relative difference 3.4 × 10⁻⁷ |
| `map`, association | 1,397,910 *A. thaliana* tests; 231,287 maize tests (chromosome 10) | −log₁₀*P* within 5 × 10⁻⁴ |
| `map`, eGenes | *A. thaliana*, uniform / PlantDNA-FM prior | 4,293 / 4,334 (paper 4,298 / 4,340; see below) |
| `map`, SuSiE | maize chromosome 10, 534 / 540 eGenes | PIPs within 6 × 10⁻⁷; the same 893 / 903 credible sets |

Settings that reproduce the paper's sets of tested variants:

- *A. thaliana*: `windows --all-biotypes` (the paper also tested 16 genes annotated as `nontranslating_CDS`) and `map --scored-only`. The remaining 407 tests (0.03%) differ because the paper computed the call rate over all 471 accessions rather than over the accessions with expression of each gene.
- maize: `map --scored-only --min-call 0` (the paper applied no call-rate filter to the imputed HapMap genotypes).

## Citation

Yu H., Korte A., Voichek Y., Menon G., Sun W., Dean C., Nordborg M., Zhang H. & Ding Y. DNA foundation model priors improve fine-mapping in plant populations (in preparation).

SuSiE: Wang G., Sarkar A., Carbonetto P. & Stephens M. A simple new approach to variable selection in regression, with application to genetic fine mapping. *J. R. Stat. Soc. B* 82, 1273–1300 (2020).

## License

MIT (see `LICENSE`). The example data come from published sources and keep their original terms; see [examples/README.md](examples/README.md).
