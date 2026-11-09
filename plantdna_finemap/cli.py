"""Command line: plantdna-finemap {windows, kinship, attention, map}."""
import argparse
import os
import sys
import time


def _log(s):
    print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)


def _lines(path):
    return [l.strip() for l in open(path) if l.strip()] if path else None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="plantdna-finemap", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("windows", help="8-kb gene windows and representative-transcript exons from a GFF3")
    p.add_argument("--gff", required=True, help="GFF3 (.gz allowed)")
    p.add_argument("--out", required=True, help="output windows TSV")
    p.add_argument("--chroms", help="comma-separated chromosomes to keep")
    p.add_argument("--all-biotypes", action="store_true", help="keep all genes, not only biotype=protein_coding")

    p = sub.add_parser("kinship", help="relatedness matrix from genome-wide SNVs")
    p.add_argument("--vcf", required=True, help="bgzipped, tabix-indexed VCF")
    p.add_argument("--out", required=True, help="output kinship TSV")
    p.add_argument("--expression", help="use the accessions (columns) of this expression table")
    p.add_argument("--samples", help="file with one accession per line (default: all VCF samples)")
    p.add_argument("--chroms", help="comma-separated chromosomes (default: all VCF contigs)")
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--min-call", type=float, default=0.9)
    p.add_argument("--workers", type=int, default=4)

    p = sub.add_parser("attention", help="per-variant PlantDNA-FM attention scores (GPU recommended)")
    p.add_argument("--model", required=True, help="model directory or Hugging Face model ID")
    p.add_argument("--windows", required=True)
    p.add_argument("--fasta", required=True, help="reference genome (indexed FASTA, plain or bgzipped)")
    p.add_argument("--vcf", required=True, help="variants to score (all SNVs in each window)")
    p.add_argument("--out", required=True, help="output directory; rerunning resumes")
    p.add_argument("--genes", help="file with one gene ID per line")
    p.add_argument("--device", default="cuda")

    p = sub.add_parser("map", help="association, eGene detection and SuSiE fine-mapping (CPU)")
    p.add_argument("--vcf", required=True)
    p.add_argument("--expression", required=True, help="TSV: first column gene ID, one column per accession")
    p.add_argument("--windows", required=True)
    p.add_argument("--kinship", required=True)
    p.add_argument("--attention", help="attention scores (parquet or TSV with gene, pos, attention); omit for the uniform prior only")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--lambda", dest="lam", type=float, default=1.0, help="prior strength (default 1)")
    p.add_argument("--alpha", type=float, default=0.05, help="family-wise error rate (default 0.05)")
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--min-call", type=float, default=0.5, help="minimum call rate per variant in the accessions with expression")
    p.add_argument("--min-samples", type=int, default=30)
    p.add_argument("--scored-only", action="store_true", help="test only variants with an attention score")
    p.add_argument("--L", type=int, default=10, help="maximum number of SuSiE effects per gene")
    p.add_argument("--genes", help="file with one gene ID per line")
    p.add_argument("--workers", type=int, default=4)

    a = ap.parse_args(argv)
    # BLAS threads per worker process, set before numpy is imported
    threads = 1 if a.cmd == "map" else max(1, (os.cpu_count() or 1) // getattr(a, "workers", 1))
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, str(threads))
    if a.cmd == "windows":
        from .windows import make_windows
        W = make_windows(a.gff, None if a.all_biotypes else "protein_coding", a.chroms.split(",") if a.chroms else None)
        W.to_csv(a.out, sep="\t", index=False)
        _log(f"windows: {len(W)} genes, {(W.exons != '').sum()} with exons -> {a.out}")
    elif a.cmd == "kinship":
        from .kinship import kinship
        samples = _lines(a.samples)
        if a.expression:
            from .mapping import read_expression
            samples = list(read_expression(a.expression).columns)
        if samples is not None:
            from .genotypes import VCF
            vs = set(VCF(a.vcf).samples)
            n0, samples = len(samples), [s for s in samples if s in vs]
            if len(samples) < n0:
                _log(f"kinship: {n0 - len(samples)} accessions are not in the VCF and are skipped")
        K, _ = kinship(a.vcf, samples, a.chroms.split(",") if a.chroms else None, a.maf, a.min_call, a.workers, log=_log)
        K.to_csv(a.out, sep="\t")
    elif a.cmd == "attention":
        from .attention import run
        run(a.model, a.windows, a.fasta, a.vcf, a.out, _lines(a.genes), a.device, log=_log)
    elif a.cmd == "map":
        from .mapping import run
        run(a.vcf, a.expression, a.windows, a.kinship, a.out, a.attention, a.lam, a.alpha, a.maf, a.min_call, a.min_samples, a.L,
            a.workers, _lines(a.genes), a.scored_only, log=_log)


if __name__ == "__main__":
    sys.exit(main())
