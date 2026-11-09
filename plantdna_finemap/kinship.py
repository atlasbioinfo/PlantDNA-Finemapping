"""Centred relatedness matrix K = Gc Gc^T / p from genome-wide SNVs.

SNVs with MAF >= 0.05 and call rate >= 0.9 in the given samples are used; missing calls are
mean-imputed and each SNV is centred. Chromosomes are processed in parallel.
"""
from multiprocessing import Pool
import numpy as np
import pandas as pd
from .genotypes import VCF


def _block(X, maf, min_call):
    called = ~np.isnan(X)
    cnt = called.sum(0)
    mu = np.where(called, X, 0).sum(0) / np.maximum(cnt, 1)
    af = mu / 2
    keep = (np.minimum(af, 1 - af) >= maf) & (cnt >= min_call * X.shape[0])
    Xk = np.where(called, X, mu[None, :])[:, keep]
    Xk -= Xk.mean(0, keepdims=True)
    return Xk.astype(np.float64) @ Xk.T.astype(np.float64), int(keep.sum())


def _chrom(args):
    path, chrom, samples, maf, min_call, chunk = args
    v = VCF(path)
    cols = v.sample_index(samples)
    n = len(samples)
    K, m, buf = np.zeros((n, n)), 0, []
    for _, _, _, d in v.records(chrom, cols=cols):
        buf.append(d)
        if len(buf) == chunk:
            k, c = _block(np.stack(buf, 1), maf, min_call)
            K += k
            m += c
            buf = []
    if buf:
        k, c = _block(np.stack(buf, 1), maf, min_call)
        K += k
        m += c
    return chrom, K, m


def kinship(vcf, samples=None, chroms=None, maf=0.05, min_call=0.9, workers=4, chunk=20000, log=print):
    v = VCF(vcf)
    samples = list(samples) if samples is not None else v.samples
    missing = set(samples) - set(v.samples)
    if missing:
        raise ValueError(f"{len(missing)} samples are not in the VCF, e.g. {sorted(missing)[:3]}")
    chroms = chroms or v.contigs
    K, M = np.zeros((len(samples), len(samples))), 0
    with Pool(min(workers, len(chroms))) as pool:
        for c, k, m in pool.imap_unordered(_chrom, [(vcf, c, samples, maf, min_call, chunk) for c in chroms]):
            K += k
            M += m
            log(f"kinship: chromosome {c}, {m:,} SNVs")
    if M == 0:
        raise ValueError("no SNV passed the filters")
    log(f"kinship: {len(samples)} samples, {M:,} SNVs")
    return pd.DataFrame((K / M).astype(np.float32), index=samples, columns=samples), M
