"""expGWAS association tests from genotypes and an expression table (linear mixed model).

For each gene, the variants in gene body +/- flank are tested against the expression of the gene:
    y = x_j b_j + g + e,   g ~ N(0, s2 h2 K),   e ~ N(0, s2 (1 - h2) I).
K is the centred relatedness matrix G_c G_c' / p (G_c: mean-imputed, column-centred dosages of p variants), computed
from the genotype file or read from --kinship. Per gene: h2 by REML; genotypes and expression are whitened with
(h2 K + (1 - h2) I)^-1/2 and the intercept is projected out; each variant is tested by its correlation r with the
expression, t = r sqrt((n - 2) / (1 - r2)), two-sided P from a t distribution with n - 2 d.f.
Variants are filtered per gene, in the samples with expression of that gene: MAF >= maf, call rate >= min_call.
The output is the association table that `attention` and `finemap` take (one row per variant-gene test).
Any other GWAS tool can be used instead (e.g. GEMMA per gene); `finemap` only needs gene, chrom, pos and P.
"""
import time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from .finemap import read_expression
from .genes import read_genes
from .genotypes import open_genotypes

_G = {}


def clean(X, maf=0.05, min_call=0.5):
    """Filter variants (columns) by call rate and MAF and mean-impute. Returns kept mask, imputed X, MAF."""
    ok = ~np.isnan(X)
    call = ok.mean(0)
    af = np.where(call > 0, np.nansum(X, 0) / np.maximum(ok.sum(0), 1) / 2, 0)
    m = np.minimum(af, 1 - af)
    keep = (call >= min_call) & (m >= maf)
    X = X[:, keep]
    X = np.where(np.isnan(X), 2 * af[keep][None, :], X)
    return keep, X.astype(np.float64), m[keep]


def kinship(geno, cols, maf=0.05, min_call=0.9, every=1, log=print):
    """Centred relatedness matrix G_c G_c' / p from the variants of the genotype file (every `every`-th one)."""
    n = len(cols)
    K, p, t0 = np.zeros((n, n)), 0, time.time()
    for chrom, pos, X in geno.chunks(cols):
        if every > 1:
            X = X[:, ::every]
        _, X, _ = clean(X, maf, min_call)
        X = X - X.mean(0)
        K += X @ X.T
        p += X.shape[1]
    if p == 0:
        raise ValueError("kinship: no variant passed the filters")
    log(f"assoc: kinship from {p:,} variants ({time.time() - t0:.0f} s)")
    return K / p


def read_kinship(path, samples):
    """Square matrix with sample IDs in the first column and in the header (TSV); reordered to `samples`."""
    K = pd.read_csv(path, sep="\t", index_col=0)
    K.index, K.columns = K.index.astype(str), K.columns.astype(str)
    miss = [s for s in samples if s not in K.index]
    if miss:
        raise ValueError(f"{path}: {len(miss)} samples are missing, e.g. {miss[:3]}")
    return K.loc[samples, samples].values.astype(np.float64)


def lmm_logp(X, y, lam, U):
    """-log10 P of each column of X in the mixed model, given the eigendecomposition K = U diag(lam) U'.
    Returns -log10 P and the REML estimate of h2."""
    n = len(y)
    lam = np.maximum(lam, 0)
    yr, one = U.T @ y, U.T @ np.ones(n)

    def negll(h2):
        w = 1.0 / (h2 * lam + (1 - h2))
        b = (one * w) @ yr / ((one * w) @ one)
        r = yr - one * b
        s2 = (r * r * w).sum() / (n - 1)
        return 0.5 * (-np.log(w).sum() + (n - 1) * np.log(s2) + np.log((one * w) @ one))

    grid = np.linspace(0.01, 0.99, 50)
    h2 = grid[int(np.argmin([negll(h) for h in grid]))]
    h2 = float(minimize_scalar(negll, bounds=(max(0.001, h2 - 0.03), min(0.999, h2 + 0.03)), method="bounded").x)
    sw = np.sqrt(1.0 / (h2 * lam + (1 - h2)))
    ys, Xs, o = sw * yr, sw[:, None] * (U.T @ X), sw * one
    Xs = Xs - np.outer(o, (o @ Xs) / (o @ o))             # project out the (whitened) intercept
    ys = ys - o * ((o @ ys) / (o @ o))
    r = (Xs * ys[:, None]).sum(0) / np.sqrt((Xs ** 2).sum(0) * (ys ** 2).sum() + 1e-300)
    r = np.clip(r, -0.999999999, 0.999999999)
    t = r * np.sqrt((n - 2) / (1 - r ** 2))
    return -(stats.t.logsf(np.abs(t), n - 2) + np.log(2)) / np.log(10), h2


def _init(genotypes, Y, K, samples, par):
    _G.update(geno=open_genotypes(genotypes), Y=Y, K=K, samples=samples, par=par, eig={})
    _G["cols"] = np.array(_G["geno"].sample_index(samples), int)


def _gene(job):
    gene, chrom, start, end = job
    par = _G["par"]
    y = _G["Y"].loc[gene].values.astype(float)
    ok = ~np.isnan(y)
    if ok.sum() < par["min_samples"]:
        return None
    pos, X = _G["geno"].region(chrom, max(1, start - par["flank"]), end + par["flank"], list(_G["cols"][ok]))
    if not len(pos):
        return None
    keep, X, maf = clean(X, par["maf"], par["min_call"])
    if not keep.any():
        return None
    key = ok.tobytes()
    if key not in _G["eig"]:                              # most genes share the same set of samples
        if len(_G["eig"]) > 8:
            _G["eig"].clear()
        _G["eig"][key] = np.linalg.eigh(_G["K"][np.ix_(ok, ok)])
    lam, U = _G["eig"][key]
    logp, h2 = lmm_logp(X, y[ok], lam, U)
    return pd.DataFrame({"gene": gene, "chrom": chrom, "pos": pos[keep], "maf": np.round(maf, 4), "logp": np.round(logp, 4),
                         "h2": round(h2, 3), "n": int(ok.sum())})


def run(genotypes, expression, genes, out, flank=10000, maf=0.05, min_call=0.5, min_samples=30, kinship_file=None,
        kinship_every=1, workers=4, log=print):
    t0 = time.time()
    geno = open_genotypes(genotypes)
    Y = read_expression(expression)
    samples = [s for s in Y.columns if s in set(geno.samples)]
    if len(samples) < min_samples:
        raise ValueError(f"only {len(samples)} samples are in both the expression table and the genotypes")
    Y = Y[samples]
    G = read_genes(genes)
    G = G[G.gene.isin(Y.index)]
    log(f"assoc: {len(samples)} samples, {len(G)} genes with expression (of {len(Y)} in the expression table)")
    if kinship_file:
        K = read_kinship(kinship_file, samples)
    else:
        K = kinship(geno, geno.sample_index(samples), maf, 0.9, kinship_every, log)
    jobs = [(g, str(c), int(s), int(e)) for g, c, s, e in zip(G.gene, G.chrom, G.gene_start, G.gene_end)]
    par = dict(flank=flank, maf=maf, min_call=min_call, min_samples=min_samples)
    res, first = 0, True
    import gzip
    opener = gzip.open if out.endswith(".gz") else open
    with opener(out, "wt") as fh, Pool(workers, initializer=_init, initargs=(genotypes, Y, K, samples, par)) as pool:
        tested = 0
        for i, d in enumerate(pool.imap(_gene, jobs, chunksize=16), 1):
            if d is not None:
                d.to_csv(fh, sep="\t", index=False, header=first)
                first, res, tested = False, res + len(d), tested + 1
            if i % 2000 == 0:
                log(f"assoc: {i}/{len(jobs)} genes, {res:,} tests")
    log(f"assoc: {res:,} tests in {tested} genes -> {out} ({time.time() - t0:.0f} s)")
