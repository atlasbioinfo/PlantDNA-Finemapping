# Example

*A. thaliana* batch B2 (471 accessions): the 4,872 tested genes of chromosomes 2 and 4 plus *FLC* (AT5G10140, chromosome 5).

| file | content |
|---|---|
| `tests_B2.tsv.gz` | 1,307,160 GEMMA tests of Voichek *et al.* (gene, chrom, pos, maf, logp = −log₁₀*P*); variants in gene body ± 10 kb. These are the tests used in the paper. |
| `expression_B2.tsv.gz` | expression of the 4,873 genes in the 471 accessions, as used for the tests |
| `genotypes_B2.bed.gz/.bim.gz/.fam` | PLINK 1 fileset (`.bed` and `.bim` gzipped; read directly) of the 385,382 tested positions that are biallelic SNVs in the 1001 Genomes VCF (A1 = ALT), 471 accessions |
| `attention_PlantDNA-FM.parquet` | PlantDNA-FM attention score of every test, as used in the paper (same as the `attention` output) |
| `kinship_B2.tsv.gz` | genome-wide relatedness matrix of the 471 accessions (755,957 SNPs of all five chromosomes); input of `assoc --kinship` |
| `TAIR10_example.gff3.gz` | TAIR10 / Araport11 annotation (Ensembl Plants) of chromosomes 2 and 4 and of *FLC*; input of `genes` |

## Run (CPU)

```bash
bash examples/run_example.sh 8
```

This runs three analyses (about 5 min with 8 workers on a 32-core server).

**1. The analysis of the paper for chromosomes 2 and 4.** The published association tests, the genome-wide number of tests of the paper (*M* = 4,297,315) and the λ of the paper (1.35 for chromosomes 1–3, 1.65 for chromosomes 4–5):

```bash
plantdna-finemap finemap --tests examples/arabidopsis/tests_B2.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA-FM.parquet \
    --folds 2:4,5 --n-tests 4297315 --lambda 1.35,1.65 --out results/arabidopsis --workers 8
```

Expected `results/arabidopsis/summary.tsv`, identical to the results of the paper for these two chromosomes:

| prior | λ | genes tested | *M* | eGenes | eGenes with signal | expGWAS signals | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|---|
| uniform | – | 4,873 | 4,297,315 | 1,314 | 1,307 | 2,114 | 797 | 904 | 4.65 |
| PlantDNA-FM | 1.35 / 1.65 | 4,873 | 4,297,315 | 1,365 | 1,362 | 2,158 | 817 | 924 | 4.26 |

**2. Single-gene mode for *FLC*** (`--single-gene AT5G10140 --lambda 1.65`, Bonferroni over the 242 variants of *FLC*). Expected `results/arabidopsis_FLC/signals.tsv`: one signal with each prior, 7 variants (uniform; lead 3,178,749, PIP 0.70) and 6 variants (PlantDNA-FM; lead 3,178,749, PIP 0.81).

**3. The whole workflow from genotypes and expression**, as for your own data:

```bash
plantdna-finemap genes --gff examples/arabidopsis/TAIR10_example.gff3.gz --all-biotypes --out results/genes.tsv
plantdna-finemap assoc --genotypes examples/arabidopsis/genotypes_B2 --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genes results/genes.tsv --kinship examples/arabidopsis/kinship_B2.tsv.gz --out results/assoc.tsv.gz --workers 8
plantdna-finemap finemap --tests results/assoc.tsv.gz --expression examples/arabidopsis/expression_B2.tsv.gz \
    --genotypes examples/arabidopsis/genotypes_B2 --attention examples/arabidopsis/attention_PlantDNA-FM.parquet \
    --folds 2:4,5 --n-tests 4297315 --lambda 1.35,1.65 --out results/arabidopsis_from_genotypes --workers 8
```

Here the association tests are computed by `assoc` (15 s; 1,260,663 tests in 4,873 genes) instead of taken from the published GEMMA run. The two sets of tests agree closely (−log₁₀*P* correlated at *r* = 0.997), and so do the results. Expected `results/arabidopsis_from_genotypes/summary.tsv`:

| prior | λ | genes tested | *M* | eGenes | eGenes with signal | expGWAS signals | single-variant signals | variants PIP > 0.9 | mean signal size |
|---|---|---|---|---|---|---|---|---|---|
| uniform | – | 4,873 | 4,297,315 | 1,306 | 1,303 | 2,107 | 795 | 915 | 4.63 |
| PlantDNA-FM | 1.35 / 1.65 | 4,873 | 4,297,315 | 1,355 | 1,354 | 2,143 | 798 | 913 | 4.30 |

1,281 of the 1,314 eGenes of the paper (uniform prior) and 1,333 of the 1,365 (PlantDNA-FM prior) are found again. A VCF (`--genotypes file.vcf.gz`, bgzipped with a tabix index) works in the same way as the PLINK fileset.

On your own data, leave out `--n-tests`, `--lambda` and `--folds`: *M* is then the number of tests in your table and λ is learned by cross-fitting.

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
