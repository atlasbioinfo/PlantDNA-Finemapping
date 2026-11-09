"""Linear mixed model by whitening (FaST-LMM style) and the marginal association test."""
import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar


def lmm_transform(X, y, K):
    """Cov(y) = h2 K + (1 - h2) I with h2 by REML; whiten X and y with Cov^-1/2 and project out the intercept.
    Returns X*, y*, h2."""
    n = len(y)
    lam, U = np.linalg.eigh(K)
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
    ys, Xs, os_ = sw * yr, sw[:, None] * (U.T @ X), sw * one
    Xs = Xs - np.outer(os_, (os_ @ Xs) / (os_ @ os_))
    ys = ys - os_ * ((os_ @ ys) / (os_ @ os_))
    return Xs, ys, h2


def marginal_logp(Xt, yt):
    """-log10 P of the per-variant regression on whitened data (t test, n - 2 d.f.)."""
    n = len(yt)
    Xc, yc = Xt - Xt.mean(0), yt - yt.mean()
    r = (Xc * yc[:, None]).sum(0) / np.sqrt((Xc ** 2).sum(0) * (yc ** 2).sum() + 1e-12)
    r = np.clip(r, -0.999999, 0.999999)
    t = r * np.sqrt((n - 2) / (1 - r ** 2))
    return -np.log10(np.maximum(2 * stats.t.sf(np.abs(t), n - 2), 1e-300))
