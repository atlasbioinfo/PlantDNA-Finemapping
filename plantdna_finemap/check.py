"""Checks of the input files, run by `plantdna-finemap check` and at the start of `attention` and `finemap`.

Each check prints one line: OK, WARN (the run continues, but look at it) or FAIL (the inputs do not fit together and
the results would be wrong or empty). The checks look for the usual mistakes: a missing index, a table without the
required columns, sample IDs or chromosome names that differ between files, gene IDs that differ between files, a
genome assembly that does not match the variants, or an attention file computed for other tests.
"""
import os
import numpy as np
import pandas as pd
from .genotypes import open_genotypes, resolve_contig


class Report:
    def __init__(self, log=print):
        self.log, self.fails, self.warns = log, [], []

    def ok(self, msg):
        self.log(f"  OK    {msg}")

    def warn(self, msg):
        self.warns.append(msg)
        self.log(f"  WARN  {msg}")

    def fail(self, msg):
        self.fails.append(msg)
        self.log(f"  FAIL  {msg}")

    def rate(self, frac, what, fail_below, warn_below, hint):
        """One line for a fraction: FAIL below fail_below, WARN below warn_below."""
        msg = f"{what}: {frac:.1%}"
        if frac < fail_below:
            self.fail(f"{msg}. {hint}")
        elif frac < warn_below:
            self.warn(f"{msg}. {hint}")
        else:
            self.ok(msg)

    def finish(self, what="inputs"):
        if self.fails:
            raise SystemExit(f"\n{what}: {len(self.fails)} check(s) failed. Fix the FAIL lines above; "
                             "`plantdna-finemap check` repeats the checks without running the analysis.")
        self.log(f"  {what}: all checks passed" + (f" ({len(self.warns)} warning(s))" if self.warns else ""))


def _exists(R, path, what, extra=()):
    if not path or not os.path.exists(path):
        R.fail(f"{what}: file not found: {path}")
        return False
    for e, hint in extra:
        if not os.path.exists(e):
            R.fail(f"{what}: {e} is missing. {hint}")
            return False
    return True


def _genotype_files(R, path):
    if path.endswith(".npz"):
        return _exists(R, path, "genotypes")
    pre = path
    for ext in (".bed.gz", ".bim.gz", ".bed", ".bim", ".fam"):
        if pre.endswith(ext):
            pre = pre[:-len(ext)]
            break
    has = lambda e: os.path.exists(pre + e) or os.path.exists(pre + e + ".gz")
    if has(".bed") or has(".bim") or os.path.exists(pre + ".fam"):
        miss = [e for e in (".bed", ".bim") if not has(e)] + ([] if os.path.exists(pre + ".fam") else [".fam"])
        if miss:
            R.fail(f"genotypes: PLINK fileset {pre} is incomplete, missing {', '.join(miss)}")
            return False
        return True
    if not _exists(R, path, "genotypes"):
        return False
    if not path.endswith((".vcf.gz", ".vcf.bgz")):
        R.fail(f"genotypes: {path} is not a .vcf.gz, a PLINK prefix or a .npz. Compress a VCF with `bgzip file.vcf`.")
        return False
    if not (os.path.exists(path + ".tbi") or os.path.exists(path + ".csi")):
        R.fail(f"genotypes: the index {path}.tbi is missing. Make it with `tabix -p vcf {path}` (the file must be bgzipped, not gzipped).")
        return False
    return True


def check_tests(R, tests):
    """Read and check the association table. Returns the table or None."""
    from .finemap import read_tests
    if not _exists(R, tests, "association tests"):
        return None
    try:
        T = read_tests(tests)
    except Exception as e:
        R.fail(f"association tests: {e}")
        return None
    if not len(T):
        R.fail("association tests: the table is empty")
        return None
    bad = int((~np.isfinite(T.logp) | (T.logp < 0)).sum())
    if bad:
        R.fail(f"association tests: {bad:,} values of -log10 P are negative or not finite. `logp` must be -log10 P; give P values in a column `p`.")
    if T.logp.max() <= 1:
        R.warn("association tests: all values of `logp` are <= 1. If the column holds P values, name it `p`.")
    per = T.groupby("gene").size()
    R.ok(f"association tests: {len(T):,} tests, {len(per):,} genes (median {int(per.median())} variants per gene), "
         f"chromosomes {', '.join(sorted(T.chrom.unique(), key=lambda c: (len(c), c))[:12])}")
    if per.median() < 5:
        R.warn("association tests: fewer than 5 variants per gene in most genes; the prior ranks variants within a gene and needs several")
    return T


def check_finemap(R, T, expression, genotypes, attention=None, genes=None, single_gene=None, lam=None, folds=None):
    """Checks for `finemap`: expression, genotypes, attention and gene table against the association tests."""
    from .finemap import read_expression, _table
    tg = set(T.gene)
    # expression
    E = None
    if _exists(R, expression, "expression"):
        try:
            E = read_expression(expression)
        except Exception as e:
            R.fail(f"expression: cannot read the table ({e}). It must be tab-separated: first column gene ID, then one column per sample.")
    if E is not None:
        if E.shape[1] < 2:
            R.fail("expression: fewer than 2 sample columns. The table must be tab-separated, genes in rows and samples in columns.")
        if E.notna().sum().sum() == 0:
            R.fail("expression: no numeric values")
        share = len(tg & set(E.index)) / len(tg)
        hint = (f"Gene IDs differ between the files: tests e.g. {sorted(tg)[0]}, expression e.g. {E.index[0]}. "
                "If the expression table has samples in rows, transpose it.")
        R.rate(share, "genes of the association tests with expression", 0.01, 0.9, hint)
        if E.index.duplicated().any():
            R.warn(f"expression: {int(E.index.duplicated().sum())} duplicated gene IDs; the first row of each is used")
    # genotypes
    geno = None
    if _genotype_files(R, genotypes):
        try:
            geno = open_genotypes(genotypes)
        except Exception as e:
            R.fail(f"genotypes: cannot open {genotypes} ({e})")
    acc = []
    if geno is not None and E is not None:
        acc = [s for s in E.columns if s in set(geno.samples)]
        hint = (f"Sample IDs must be the same in both files: expression e.g. {list(E.columns[:3])}, genotypes e.g. {geno.samples[:3]}.")
        if len(acc) < 10:
            R.fail(f"samples with expression and genotypes: {len(acc)} of {E.shape[1]}. {hint}")
        else:
            R.rate(len(acc) / E.shape[1], f"samples with expression and genotypes ({len(acc)} of {E.shape[1]})", 0.0, 0.8, hint)
    if geno is not None:
        chroms = sorted(T.chrom.unique())
        miss = [c for c in chroms if resolve_contig(c, geno.contigs) is None]
        if len(miss) == len(chroms):
            R.fail(f"genotypes: none of the chromosomes of the tests ({', '.join(chroms[:6])}) is in the genotypes ({', '.join(map(str, geno.contigs[:6]))})")
        elif miss:
            R.warn(f"genotypes: chromosomes of the tests that are not in the genotypes: {', '.join(miss[:10])}")
        if len(miss) < len(chroms) and acc:
            rng = np.random.default_rng(0)
            gs = T[~T.chrom.isin(miss)].gene.unique()
            n = found = 0
            for g in rng.choice(gs, min(20, len(gs)), replace=False):
                d = T[T.gene == g]
                P, _ = geno.fetch(d.chrom.iloc[0], d.pos.values, acc[:5])
                n, found = n + len(d), found + len(P)
            R.rate(found / max(n, 1), "tested variants found in the genotypes (20 genes)", 0.5, 0.9,
                   "Positions of the tests and of the genotypes must be on the same assembly; only biallelic SNVs are read from a VCF.")
    # attention
    if attention is not None and _exists(R, attention, "attention"):
        try:
            S = _table(attention)
            miss = [c for c in ("gene", "pos", "attention") if c not in S.columns]
            if miss:
                R.fail(f"attention: column(s) {', '.join(miss)} missing; the file must be the attention.parquet written by `attention`")
            else:
                S = S.assign(gene=S.gene.astype(str), pos=S.pos.astype(np.int64))
                cov = T[["gene", "pos"]].merge(S[["gene", "pos"]].drop_duplicates(), on=["gene", "pos"], how="left", indicator=True)._merge.eq("both").mean()
                R.rate(cov, "tests with an attention score", 0.01, 0.95,
                       "The attention file must be computed for these tests (same gene IDs, same assembly); tests without a score get a neutral weight.")
                if (S.attention <= 0).any() or S.attention.isna().any():
                    R.warn("attention: some scores are missing, zero or negative")
        except Exception as e:
            R.fail(f"attention: cannot read {attention} ({e})")
    # gene table
    if genes is not None and _exists(R, genes, "gene table"):
        check_genes(R, T, genes)
    # folds
    if single_gene:
        if single_gene not in tg:
            R.fail(f"single-gene mode: gene {single_gene} is not in the association tests")
    elif attention is not None and lam is None and folds is None and T.chrom.nunique() < 2:
        R.fail("cross-fitting needs tests on at least two chromosomes; give a fixed --lambda for a single chromosome")
    return geno


def check_genes(R, T, genes, flank=10000):
    from .genes import read_genes
    try:
        G = read_genes(genes)
    except Exception as e:
        R.fail(f"gene table: cannot read {genes} ({e})")
        return None
    miss = [c for c in ("gene", "chrom", "strand", "gene_start", "gene_end") if c not in G.columns]
    if miss:
        R.fail(f"gene table: column(s) {', '.join(miss)} missing; use the output of `plantdna-finemap genes`")
        return None
    G = G.drop_duplicates("gene").set_index("gene")
    tg = pd.Index(T.gene.unique())
    share = tg.isin(G.index).mean()
    R.rate(share, "genes of the association tests in the gene table", 0.01, 0.95,
           f"Gene IDs differ: tests e.g. {tg[0]}, gene table e.g. {G.index[0]}. `genes --all-biotypes` keeps genes that are not protein-coding.")
    if share > 0:
        d = T[T.gene.isin(G.index)]
        s, e = d.gene.map(G.gene_start).astype(int), d.gene.map(G.gene_end).astype(int)
        near = ((d.pos >= s - 10 * flank) & (d.pos <= e + 10 * flank)).mean()
        R.rate(near, "tested variants within 100 kb of their gene", 0.5, 0.95,
               "The annotation and the variant positions seem to be on different assemblies, or the tests are not cis tests.")
        same = (d.chrom.map(lambda c: c[3:] if c.lower().startswith("chr") else c)
                == d.gene.map(G.chrom).map(lambda c: c[3:] if str(c).lower().startswith("chr") else str(c))).mean()
        if same < 0.99:
            R.warn(f"gene table: {1 - same:.1%} of the tests are on another chromosome than their gene")
    return G


def check_attention_inputs(R, T, genes, fasta):
    """Checks for `attention`: gene table and reference genome against the association tests."""
    import pysam
    G = check_genes(R, T, genes) if _exists(R, genes, "gene table") else None
    if not _exists(R, fasta, "reference genome", [(fasta + ".fai", f"Make it with `samtools faidx {fasta}`.")]):
        return
    try:
        fa = pysam.FastaFile(fasta)
    except Exception as e:
        R.fail(f"reference genome: cannot open {fasta} ({e}); it must be plain or bgzipped FASTA with a .fai index")
        return
    L = {c: fa.get_reference_length(c) for c in fa.references}
    chroms = sorted(T.chrom.unique())
    miss = [c for c in chroms if resolve_contig(c, L) is None]
    if len(miss) == len(chroms):
        R.fail(f"reference genome: none of the chromosomes of the tests ({', '.join(chroms[:6])}) is in the FASTA ({', '.join(list(L)[:6])})")
        return
    if miss:
        R.warn(f"reference genome: chromosomes of the tests that are not in the FASTA: {', '.join(miss[:10])}")
    d = T[~T.chrom.isin(miss)]
    over = (d.pos > d.chrom.map(lambda c: L[resolve_contig(c, L)])).mean()
    if over > 0:
        R.fail(f"reference genome: {over:.1%} of the tested positions lie beyond the end of their chromosome; the FASTA is another assembly")
    if "ref" in T.columns and over == 0:
        s = d[d.ref.astype(str).str.len() == 1].sample(min(2000, len(d)), random_state=0)
        hit = np.mean([fa.fetch(resolve_contig(c, L), p - 1, p).upper() == r.upper() for c, p, r in zip(s.chrom, s.pos, s.ref)]) if len(s) else 1.0
        R.rate(hit, "REF alleles of the tests that match the FASTA (2,000 variants)", 0.9, 0.99,
               "The FASTA must be the assembly on which the variants were called.")
    else:
        R.ok(f"reference genome: {len(L)} sequences, all chromosomes of the tests present")


def check_gpu(R):
    """Environment of the attention step."""
    try:
        import torch
        import transformers
    except ImportError as e:
        R.fail(f"attention environment: {e}. Install with `pip install \".[attention]\"`.")
        return
    R.ok(f"torch {torch.__version__}, transformers {transformers.__version__}")
    if not torch.cuda.is_available():
        R.warn("no CUDA GPU is visible; `attention --device cpu` works but takes minutes per window")
        return
    gb = torch.cuda.get_device_properties(0).total_memory / 2 ** 30
    (R.ok if gb >= 14 else R.warn)(f"GPU {torch.cuda.get_device_name(0)}, {gb:.0f} GB" + ("" if gb >= 14 else "; a 16-kb window needs about 12 GB"))


def run(tests, expression=None, genotypes=None, attention=None, genes=None, fasta=None, gpu=False, log=print):
    """`plantdna-finemap check`: all checks for the files that are given."""
    import sys, numpy, scipy, pysam
    log(f"environment: python {sys.version.split()[0]}, numpy {numpy.__version__}, pandas {pd.__version__}, scipy {scipy.__version__}, pysam {pysam.__version__}")
    R = Report(log)
    if gpu:
        check_gpu(R)
    T = check_tests(R, tests)
    if T is not None:
        if expression and genotypes:
            check_finemap(R, T, expression, genotypes, attention, None)
        elif expression or genotypes:
            R.warn("give both --expression and --genotypes to check them")
        if genes and fasta:
            check_attention_inputs(R, T, genes, fasta)
        elif genes:
            check_genes(R, T, genes)
    R.finish("check")


def compare(result, expected, tol=1e-3, log=print):
    """`plantdna-finemap compare`: is a finemapping.csv the same as an expected table?
    Rows are matched by gene and variant. Calls (True / False columns) and alleles must be identical; numbers may
    differ by `tol` (rounding and floating-point differences between machines). Gene annotation columns are not compared."""
    A, B = pd.read_csv(result), pd.read_csv(expected)
    R = Report(log)
    key = ["gene", "rs"]
    for name, D in (("result", A), ("expected table", B)):
        if any(c not in D.columns for c in key):
            raise SystemExit(f"compare: {name} has no columns gene and rs; give a finemapping.csv")
    M = A.merge(B, on=key, how="outer", suffixes=("_a", "_b"), indicator=True)
    only_a, only_b = int((M._merge == "left_only").sum()), int((M._merge == "right_only").sum())
    if only_a or only_b:
        R.fail(f"variants: {only_a:,} only in the result, {only_b:,} only in the expected table ({len(B):,} expected)")
    else:
        R.ok(f"variants: the same {len(B):,} significant variants in {B.gene.nunique():,} genes")
    both = M[M._merge == "both"]
    skip = set(key) | {"gene_name", "strand", "gene_start", "gene_end"}
    for c in [c for c in B.columns if c not in skip]:
        if c not in A.columns:
            R.fail(f"{c}: column missing in the result")
            continue
        a, b = both[c + "_a"], both[c + "_b"]
        if a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
            nan = int((a.isna() != b.isna()).sum())
            d = float((a - b).abs().max()) if len(a) else 0.0
            d = 0.0 if np.isnan(d) else d
            if nan or d > tol:
                R.fail(f"{c}: largest difference {d:.2g}; {nan:,} values present in only one table")
            else:
                R.ok(f"{c}: identical" if d == 0 else f"{c}: largest difference {d:.0e}")
        else:
            n = int((a.fillna("").astype(str) != b.fillna("").astype(str)).sum())
            (R.fail if n else R.ok)(f"{c}: {n:,} rows differ" if n else f"{c}: identical")
    if R.fails:
        raise SystemExit(f"\ncompare: the result differs from the expected table ({len(R.fails)} difference(s) above).")
    log("  compare: PASS, the result reproduces the expected table")
