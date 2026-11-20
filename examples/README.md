# Example

*A. thaliana* batch B2 (471 accessions): the 4,872 tested genes of chromosomes 2 and 4 (one chromosome for each cross-fitting fold) plus *FLC* (AT5G10140, chromosome 5).

| file | content |
|---|---|
| `tests_B2.tsv.gz` | 1,307,160 GEMMA tests of Voichek *et al.* (gene, chrom, pos, maf, logp = −log₁₀*P*); variants in gene body ± 10 kb |
| `expression_B2.tsv.gz` | expression of the 4,873 genes in the 471 accessions, as used for the tests |
| `genotypes_B2.bed/.bim/.fam` | PLINK 1 fileset of the 385,382 tested positions that are biallelic SNVs in the 1001 Genomes VCF (A1 = ALT), 471 accessions |
| `attention_PlantDNA16K-FM.parquet` | PlantDNA16K-FM attention score of every test, as used in the paper (same as the `attention` output) |
| `TAIR10_example.gff3.gz` | TAIR10 / Araport11 annotation (Ensembl Plants) of chromosomes 2 and 4 and of *FLC*; input of `genes` |

## Run (CPU)

```bash
bash examples/run_example.sh 8
```

This runs three analyses (about 5 min with 8 workers on a 32-core server).

**1. Fine-mapping with the published association tests**, λ cross-fitted between the two chromosomes:

```bash
plantdna-finemap finemap --tests examples/arabidopsis/tests_B2.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA16K-FM.parquet \
    --folds 2:4,5 --out results/arabidopsis --workers 8
```

Expected `results/arabidopsis/summary.tsv`:

| prior | λ | genes tested | tests (*M*) | eGenes | eGenes with signal | expGWAS signals | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|---|
| uniform | – | 4,873 | 1,307,160 | 1,401 | 1,395 | 2,247 | 822 | 943 | 4.78 |
| PlantDNA16K-FM | 2.25 / 2.2 | 4,873 | 1,307,160 | 1,487 | 1,482 | 2,305 | 917 | 1,044 | 3.93 |

*M* here is the number of tests in the example, so the threshold is less strict than genome-wide and λ is learned on two chromosomes only. With the paper's genome-wide *M* and λ (`--n-tests 4297315 --lambda 1.35,1.65`), the example gives the paper's results for chromosomes 2 and 4 exactly: 1,314 / 1,365 eGenes and 2,114 / 2,158 expGWAS signals.

**2. Single-gene mode for *FLC*** (`--single-gene AT5G10140 --lambda 1.65`, Bonferroni over the 242 variants of *FLC*). Expected `results/arabidopsis_FLC/signals.tsv`: one signal under each prior, 7 variants (uniform; lead 3,178,749, PIP 0.70) and 6 variants (PlantDNA16K-FM; lead 3,178,749, PIP 0.81).

**3. The whole workflow from genotypes and expression**, as for your own data:

```bash
plantdna-finemap genes --gff examples/arabidopsis/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap assoc --genotypes examples/arabidopsis/genotypes_B2 --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genes results/genes.tsv --out results/assoc.tsv.gz --workers 8
plantdna-finemap finemap --tests results/assoc.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA16K-FM.parquet \
    --folds 2:4,5 --out results/arabidopsis_from_genotypes --workers 8
```

`assoc` takes 15 s and gives 1,260,663 tests in 4,873 genes. Expected `results/arabidopsis_from_genotypes/summary.tsv`:

| prior | λ | genes tested | tests (*M*) | eGenes | eGenes with signal | expGWAS signals | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|---|
| uniform | – | 4,873 | 1,260,663 | 1,376 | 1,370 | 2,223 | 813 | 936 | 4.74 |
| PlantDNA16K-FM | 1.8 / 1.7 | 4,873 | 1,260,663 | 1,470 | 1,470 | 2,317 | 893 | 1,022 | 4.07 |

These numbers differ a little from analysis 1 because the association tests differ: `assoc` uses its own mixed model with a relatedness matrix from the two example chromosomes, and it tests only the positions in the PLINK fileset. Its −log₁₀*P* values correlate at *r* = 0.987 with the GEMMA tests (0.998 with a genome-wide relatedness matrix). A VCF (`--genotypes file.vcf.gz`, bgzipped with a tabix index) works in the same way.

## Recompute the attention (GPU)

```bash
bash examples/download_reference.sh        # TAIR10 from Ensembl Plants -> examples/arabidopsis/TAIR10.fa
plantdna-finemap genes --gff examples/arabidopsis/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap attention --model ATLASBIOINFO/PlantDNA-FM --genes results/genes.tsv \
    --fasta examples/arabidopsis/TAIR10.fa --tests examples/arabidopsis/tests_B2.tsv.gz --out results/attention
```

Then pass `--attention results/attention/attention.parquet` to `finemap`. `--all-biotypes` keeps the few tested genes that are not protein-coding. `--model` also takes a local model directory. On an RTX 5090 a 16-kb window takes about 1.4 s, so the 4,873 genes need about 4 h; add `--only genes.txt` (one gene ID per line) to score a few genes first. The result matches the shipped scores (the paper's) to a relative difference of 2.4 × 10⁻⁷.

## Sources

Please cite the original data when you use these files.

- Expression and association tests: leaf RNA-seq of 1001 Genomes accessions (Kawakatsu T. *et al.* Epigenomic diversity in a global collection of *Arabidopsis thaliana* accessions. *Cell* 166, 492–505 (2016)), reprocessed by pseudogenome quantification, and GEMMA LMM tests per gene (`-lmm 2 -maf 0.05 -miss 0.5`) by Voichek Y. *et al.* (2024); batch B2, 471 accessions.
- Genotypes: 1001 Genomes Consortium. 1,135 genomes reveal the global pattern of polymorphism in *Arabidopsis thaliana*. *Cell* 166, 481–491 (2016); `1001genomes_snp-short-indel_only_ACGTN.vcf.gz`.
- Annotation and reference: TAIR10 / Araport11 from Ensembl Plants.
