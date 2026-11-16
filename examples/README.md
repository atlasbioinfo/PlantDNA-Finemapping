# Example

*A. thaliana* batch B2 (471 accessions): the 892 genes whose tested variants all lie in the first megabase of a chromosome, plus *FLC* (AT5G10140).

| file | content |
|---|---|
| `tests_B2.tsv.gz` | 248,022 GEMMA tests of Voichek *et al.* (gene, chrom, pos, maf, logp = −log₁₀*P*); variants in gene body ± 10 kb |
| `expression_B2.tsv.gz` | expression of the 892 genes in the 471 accessions, as used for the tests |
| `genotypes_B2.bed/.bim/.fam` | PLINK 1 fileset of the 60,582 tested positions that are biallelic SNVs in the 1001 Genomes VCF (A1 = ALT), 471 accessions |
| `attention_PlantDNA16K-FM.parquet` | PlantDNA16K-FM attention score of every test, as used in the paper (same as the `attention` output) |
| `TAIR10_example.gff3.gz` | TAIR10 / Araport11 annotation (Ensembl Plants) of the example regions; input of `genes` |

## Run (CPU)

```bash
bash examples/run_example.sh 8
```

This runs

```bash
plantdna-finemap finemap --tests examples/arabidopsis/tests_B2.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA16K-FM.parquet \
    --out results/arabidopsis --workers 8
plantdna-finemap finemap ... --single-gene AT5G10140 --lambda 1.65 --out results/arabidopsis_FLC --workers 1
```

Expected `results/arabidopsis/summary.tsv` (55 s with 8 workers on a 32-core server; λ cross-fitted over chromosomes 1–3 / 4–5):

| prior | λ | genes tested | tests (*M*) | eGenes | eGenes with signal | expGWAS signals | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|---|
| uniform | – | 892 | 248,022 | 300 | 296 | 483 | 163 | 188 | 5.09 |
| PlantDNA16K-FM | 2.15 / 2.4 | 892 | 248,022 | 325 | 324 | 501 | 196 | 228 | 3.97 |

*M* here is the number of tests in the example, so the threshold is less strict than genome-wide and λ is learned on 892 genes only; these counts are not a subset of the genome-wide results. With the paper's genome-wide *M* and λ (`--n-tests 4297320 --lambda 1.35,1.65`), the example genes give the paper's results exactly: 259 / 270 eGenes and 417 / 421 expGWAS signals.

Expected `results/arabidopsis_FLC/signals.tsv` (single-gene mode, Bonferroni over the 242 variants of *FLC*): one signal under each prior, 7 variants (uniform; lead 3,178,749, PIP 0.70) and 6 variants (PlantDNA16K-FM, λ = 1.65; lead 3,178,749, PIP 0.81).

## Recompute the attention (GPU)

```bash
bash examples/download_reference.sh        # TAIR10 from Ensembl Plants -> examples/arabidopsis/TAIR10.fa
plantdna-finemap genes --gff examples/arabidopsis/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA16K-FM --genes results/genes.tsv \
    --fasta examples/arabidopsis/TAIR10.fa --tests examples/arabidopsis/tests_B2.tsv.gz --out results/attention
```

Then pass `--attention results/attention/attention.parquet` to `finemap`. `--all-biotypes` keeps the few tested genes that are not protein-coding. The model ID is a placeholder until the 16K weights are on Hugging Face (see the main README); a local model directory works too. On an RTX 5090 the 892 genes (1,915 windows) take 44 min. The result matches the shipped scores (the paper's) to a relative difference of 2.4 × 10⁻⁷, including AT4G00520, which has no `.1` transcript and therefore uses all positions as queries.

## Sources

Please cite the original data when you use these files.

- Expression and association tests: leaf RNA-seq of 1001 Genomes accessions (Kawakatsu T. *et al.* Epigenomic diversity in a global collection of *Arabidopsis thaliana* accessions. *Cell* 166, 492–505 (2016)), reprocessed by pseudogenome quantification, and GEMMA LMM tests per gene (`-lmm 2 -maf 0.05 -miss 0.5`) by Voichek Y. *et al.* (2024); batch B2, 471 accessions.
- Genotypes: 1001 Genomes Consortium. 1,135 genomes reveal the global pattern of polymorphism in *Arabidopsis thaliana*. *Cell* 166, 481–491 (2016); `1001genomes_snp-short-indel_only_ACGTN.vcf.gz`.
- Annotation and reference: TAIR10 / Araport11 from Ensembl Plants.
