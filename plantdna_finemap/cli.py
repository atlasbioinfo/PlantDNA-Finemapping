"""Command line: plantdna-finemap {genes, assoc, attention, finemap}."""
import argparse
import os
import sys
import time


def _log(s):
    print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)


def _lines(path):
    return [l.strip() for l in open(path) if l.strip()] if path else None


def _folds(s):
    """'1,2,3:4,5' -> {chrom: fold}; or a TSV file with columns chrom and fold (0 / 1)."""
    if s is None:
        return None
    if os.path.exists(s):
        import pandas as pd
        F = pd.read_csv(s, sep="\t", dtype=str)
        return {c: int(f) for c, f in zip(F.chrom, F.fold)}
    parts = s.split(":")
    if len(parts) != 2:
        raise SystemExit("--folds: give two comma-separated chromosome lists separated by ':' (e.g. 1,2,3:4,5) or a TSV file")
    return {c.strip(): f for f, p in enumerate(parts) for c in p.split(",") if c.strip()}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="plantdna-finemap", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("genes", help="gene table and representative-transcript exons from a GFF3")
    p.add_argument("--gff", required=True, help="GFF3 (.gz allowed)")
    p.add_argument("--out", required=True, help="output gene table (TSV)")
    p.add_argument("--chroms", help="comma-separated chromosomes to keep")
    p.add_argument("--all-biotypes", action="store_true", help="keep all genes, not only biotype=protein_coding")

    p = sub.add_parser("assoc", help="expGWAS association tests from genotypes and expression, linear mixed model (CPU)")
    p.add_argument("--genotypes", required=True, help="VCF (.vcf.gz + .tbi), PLINK 1 fileset (prefix or .bed) or .npz")
    p.add_argument("--expression", required=True, help="TSV: first column gene ID, one column per accession")
    p.add_argument("--genes", required=True, help="gene table (output of `genes`)")
    p.add_argument("--out", required=True, help="output association table (TSV, .gz allowed)")
    p.add_argument("--flank", type=int, default=10000, help="test the variants in gene body +/- flank (default 10000)")
    p.add_argument("--maf", type=float, default=0.05, help="minimum minor allele frequency (default 0.05)")
    p.add_argument("--min-call", type=float, default=0.5, help="minimum call rate of a variant (default 0.5)")
    p.add_argument("--min-samples", type=int, default=30, help="skip genes with expression in fewer samples (default 30)")
    p.add_argument("--kinship", help="relatedness matrix (TSV, sample IDs in the first column and the header); "
                                     "default: computed from the genotype file")
    p.add_argument("--kinship-every", type=int, default=1, help="use every n-th variant for the relatedness matrix (default 1)")
    p.add_argument("--workers", type=int, default=4)

    p = sub.add_parser("attention", help="per-variant PlantDNA16K-FM attention scores, sliding windows (GPU)")
    p.add_argument("--model", required=True, help="model directory or Hugging Face model ID")
    p.add_argument("--genes", required=True, help="gene table (output of `genes`)")
    p.add_argument("--fasta", required=True, help="reference genome (indexed FASTA, plain or bgzipped)")
    p.add_argument("--tests", required=True, help="association table; the variants of each gene to score (columns gene, chrom, pos)")
    p.add_argument("--out", required=True, help="output directory; rerunning resumes")
    p.add_argument("--window", type=int, default=16000, help="window length in nt (default 16000)")
    p.add_argument("--step", type=int, default=8000, help="window step in nt (default 8000)")
    p.add_argument("--flank", type=int, default=10000, help="region = gene body +/- flank, extended to all tested variants (default 10000)")
    p.add_argument("--only", help="file with one gene ID per line")
    p.add_argument("--device", default="cuda")

    p = sub.add_parser("finemap", help="weighted Bonferroni eGene detection and SuSiE fine-mapping (CPU)")
    p.add_argument("--tests", required=True, help="association table (TSV or parquet): gene, chrom, pos, and logp (-log10 P) or p")
    p.add_argument("--expression", required=True, help="TSV: first column gene ID, one column per accession")
    p.add_argument("--genotypes", required=True, help="VCF (.vcf.gz + .tbi), PLINK 1 fileset (prefix or .bed) or .npz")
    p.add_argument("--attention", help="attention scores (parquet or TSV with gene, pos, attention); omit for the uniform prior only")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--lambda", dest="lam", help="fixed prior strength: one value, or two values 'fold0,fold1' (default: cross-fitted)")
    p.add_argument("--folds", help="cross-fitting folds: '1,2,3:4,5' or a TSV with columns chrom, fold "
                                   "(default: chromosomes in order, split into two halves)")
    p.add_argument("--alpha", type=float, default=0.05, help="family-wise error rate (default 0.05)")
    p.add_argument("--n-tests", type=int, help="M for the Bonferroni threshold (default: number of rows in --tests)")
    p.add_argument("--single-gene", help="candidate-gene mode: only this gene, M = its number of tests (needs --lambda)")
    p.add_argument("--susie-all", action="store_true", help="fine-map all tested variants of each eGene, not only the significant ones")
    p.add_argument("--L", type=int, default=10, help="maximum number of SuSiE effects per gene")
    p.add_argument("--workers", type=int, default=4)

    a = ap.parse_args(argv)
    # BLAS threads per worker process, set before numpy is imported
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1" if a.cmd in ("finemap", "assoc") else str(os.cpu_count() or 1))
    if a.cmd == "genes":
        from .genes import make_genes
        G = make_genes(a.gff, None if a.all_biotypes else "protein_coding", a.chroms.split(",") if a.chroms else None)
        G.to_csv(a.out, sep="\t", index=False)
        _log(f"genes: {len(G)} genes, {(G.exons != '').sum()} with exons -> {a.out}")
    elif a.cmd == "assoc":
        from .assoc import run
        run(a.genotypes, a.expression, a.genes, a.out, a.flank, a.maf, a.min_call, a.min_samples, a.kinship, a.kinship_every,
            a.workers, log=_log)
    elif a.cmd == "attention":
        from .attention import run
        run(a.model, a.genes, a.fasta, a.tests, a.out, a.window, a.step, a.flank, _lines(a.only), a.device, log=_log)
    elif a.cmd == "finemap":
        from .finemap import run
        lam = [float(x) for x in a.lam.split(",")] if a.lam else None
        if lam is not None and len(lam) not in (1, 2):
            raise SystemExit("--lambda: one value or two values 'fold0,fold1'")
        run(a.tests, a.expression, a.genotypes, a.out, a.attention, lam, _folds(a.folds), a.alpha, a.n_tests, a.L,
            a.single_gene, a.susie_all, a.workers, log=_log)


if __name__ == "__main__":
    sys.exit(main())
