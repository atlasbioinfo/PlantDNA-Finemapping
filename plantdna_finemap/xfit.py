"""Two-fold cross-fitting of the prior strength lambda (w = exp(lambda * q), mean 1 within each gene).

Folds are sets of chromosomes (default: the chromosomes in order, split into two halves; A. thaliana 1-3 / 4-5,
maize 1-5 / 6-10). For each fold f, the lambda in GRID that maximises the number of eGenes on the other fold is
applied to the tests of fold f (independent hypothesis weighting with cross-weighting, Ignatiadis et al.), so
no test's weight depends on its own data. Ties: the median of the tied lambda values (upper middle value for an
even number of ties). M = all tests (threshold -log10(alpha * w / M)).
"""
import numpy as np
import pandas as pd

GRID = np.round(np.arange(0, 4.0001, 0.05), 2)


def gene_codes(gene):
    return pd.factorize(np.asarray(gene))[0]


def weights(q, g, lam):
    w = np.exp(lam * q)
    return w / (np.bincount(g, w) / np.bincount(g))[g]


def egenes_by_fold(logp, q, g, fold, lam, M=None, alpha=0.05):
    M = len(logp) if M is None else M
    hit = logp > -np.log10(alpha * weights(q, g, lam) / M)
    gf = np.zeros(g.max() + 1, int)
    gf[g] = fold
    eg = np.zeros(g.max() + 1, bool)
    eg[g[hit]] = True
    return np.array([int((eg & (gf == f)).sum()) for f in (0, 1)])


def xfit(logp, q, g, fold, M=None, alpha=0.05, grid=GRID):
    """Returns (w, lam_by_fold, curve). w: per-test weights with each fold's lambda learned on the other fold.
    curve: DataFrame lambda x (eGenes fold0, eGenes fold1). Genes must not span folds (fold is constant within a gene)."""
    C = np.stack([egenes_by_fold(logp, q, g, fold, l, M, alpha) for l in grid])
    lam = {}
    for f in (0, 1):
        tr = C[:, 1 - f]
        best = grid[tr == tr.max()]
        lam[f] = float(best[len(best) // 2])
    w = np.empty(len(logp))
    for f in (0, 1):
        wf = weights(q, g, lam[f])
        w[fold == f] = wf[fold == f]
    return w, lam, pd.DataFrame({"lambda": grid, "eGenes_fold0": C[:, 0], "eGenes_fold1": C[:, 1]})


def default_folds(chroms):
    """Chromosomes in natural order, first half (rounded up) -> fold 0, the rest -> fold 1."""
    def key(c):
        s = str(c)
        b = s[3:] if s.lower().startswith("chr") else s
        return (0, int(b), s) if b.isdigit() else (1, 0, s)
    u = sorted({str(c) for c in chroms}, key=key)
    h = (len(u) + 1) // 2
    return {c: int(i >= h) for i, c in enumerate(u)}
