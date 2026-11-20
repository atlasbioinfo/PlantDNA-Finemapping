"""Attention scores -> prior weights.

q    : within-gene rank-based normal quantile of the score (ties get the average rank; missing -> 0, neutral)
w    : exp(lambda * q), scaled to mean 1 within each gene (lambda: see xfit.py). The weighted Bonferroni threshold
       -log10(alpha * w / M) keeps every gene's share of alpha, so the family-wise error rate stays at alpha.
pi   : SuSiE prior inclusion probabilities 0.8 * w / sum(w) + 0.2 / m over the m variants fitted for a gene,
       so no variant falls below 0.2 / m.
"""
import numpy as np
import pandas as pd
from scipy.stats import norm


def rank_q(v, group):
    v = np.asarray(v, float)
    q = np.zeros(len(v))
    ok = ~np.isnan(v)
    if ok.sum() < 2:
        return q
    s = pd.Series(v[ok])
    g = np.asarray(group)[ok]
    r = s.groupby(g).rank(method="average").values
    n = s.groupby(g).transform("size").values
    q[ok] = norm.ppf(r / (n + 1))
    return q


def gene_weights(score, gene, lam=1.0):
    """Prior weights with mean 1 within each gene."""
    x = lam * rank_q(score, gene)
    e = np.exp(x - x.max())
    return e / pd.Series(e).groupby(np.asarray(gene)).transform("mean").values


def thresholds(w, M, alpha=0.05):
    """Per-test -log10 P thresholds of the weighted Bonferroni procedure."""
    return -np.log10(alpha * np.asarray(w, float) / M)


def susie_pi(w):
    w = np.asarray(w, float)
    pi = 0.8 * w / w.sum() + 0.2 / len(w)
    return pi / pi.sum()
