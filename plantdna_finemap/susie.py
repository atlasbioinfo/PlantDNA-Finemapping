"""Sum of Single Effects regression (SuSiE; Wang et al. 2020, J. R. Stat. Soc. B) in numpy.

Defaults follow susieR::susie: standardised genotypes, scaled prior variance 0.2 re-estimated per
effect (effects with zero prior variance are dropped), estimated residual variance, 95% credible
sets with purity (minimum absolute correlation) >= 0.5, tolerance 1e-3 on the ELBO, at most 100
iterations. prior_weights are the prior inclusion probabilities. PIPs use effects with prior
variance > 1e-9, as in susieR::susie_get_pip.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


def _lbf(betahat, shat2, V):
    if V <= 0:
        return np.zeros_like(betahat)
    lbf = 0.5 * np.log(shat2 / (shat2 + V)) + 0.5 * betahat ** 2 / shat2 * V / (V + shat2)
    lbf[~np.isfinite(shat2)] = 0.0
    return lbf


def _ser_loglik(logV, betahat, shat2, logpi):
    return logsumexp(_lbf(betahat, shat2, np.exp(logV)) + logpi)


def single_effect(Xr, d, sigma2, V, logpi, optimize_V=True):
    with np.errstate(divide="ignore", invalid="ignore"):
        betahat, shat2 = Xr / d, sigma2 / d
    betahat = np.where(np.isfinite(betahat), betahat, 0.0)
    if optimize_V:
        v0 = max(float(np.sum(betahat ** 2 * np.exp(logpi))), 1e-6)
        res = minimize_scalar(lambda lv: -_ser_loglik(lv, betahat, shat2, logpi), bounds=(-30, 15), method="bounded", options={"xatol": 1e-4})
        cand = [np.log(max(V, 1e-12)), np.log(v0), res.x]
        ll = [_ser_loglik(c, betahat, shat2, logpi) for c in cand]
        V = float(np.exp(cand[int(np.argmax(ll))]))
        if _ser_loglik(-700.0, betahat, shat2, logpi) >= max(ll):
            V = 0.0
    lbf = _lbf(betahat, shat2, V)
    logw = lbf + logpi
    lbf_model = logsumexp(logw)
    alpha = np.exp(logw - lbf_model)
    if V > 0:
        post_var = np.where(np.isfinite(1.0 / (1.0 / shat2 + 1.0 / V)), 1.0 / (1.0 / shat2 + 1.0 / V), 0.0)
        mu = np.where(np.isfinite(post_var / shat2 * betahat), post_var / shat2 * betahat, 0.0)
    else:
        post_var, mu = np.zeros_like(betahat), np.zeros_like(betahat)
    return alpha, mu, post_var + mu ** 2, V, lbf_model, lbf


def susie(X, y, L=10, prior_weights=None, scaled_prior_variance=0.2, estimate_prior_variance=True,
          estimate_residual_variance=True, coverage=0.95, min_abs_corr=0.5, tol=1e-3, max_iter=100,
          intercept=True, standardize=True):
    X, y = np.asarray(X, np.float64), np.asarray(y, np.float64)
    n, p = X.shape
    cm = X.mean(0) if intercept else np.zeros(p)
    csd = X.std(0, ddof=1) if standardize else np.ones(p)
    csd[csd == 0] = 1.0
    Xs = (X - cm) / csd
    yc = y - (y.mean() if intercept else 0.0)
    d = (Xs ** 2).sum(0)
    pi = np.full(p, 1.0 / p) if prior_weights is None else np.asarray(prior_weights, np.float64) / np.sum(prior_weights)
    logpi = np.log(np.maximum(pi, 1e-300))
    sigma2 = yc.var(ddof=1) if n > 1 else 1.0
    V = np.full(L, scaled_prior_variance * sigma2)
    alpha, mu, mu2 = np.full((L, p), 1.0 / p), np.zeros((L, p)), np.zeros((L, p))
    lbf_model, lbf_var = np.zeros(L), np.zeros((L, p))
    elbo = [-np.inf]
    for it in range(max_iter):
        b = alpha * mu
        r = yc - Xs @ b.sum(0)
        for l in range(L):
            r_l = r + Xs @ b[l]
            alpha[l], mu[l], mu2[l], V[l], lbf_model[l], lbf_var[l] = single_effect(Xs.T @ r_l, d, sigma2, V[l], logpi, estimate_prior_variance)
            b[l] = alpha[l] * mu[l]
            r = r_l - Xs @ b[l]
        XB = Xs @ b.T
        postb2 = alpha * mu2
        er2 = np.sum((yc - XB.sum(1)) ** 2) - np.sum(XB ** 2) + np.sum(d * postb2)
        kl = 0.0
        for l in range(L):
            r_l = yc - XB.sum(1) + XB[:, l]
            ser = -0.5 * n * np.log(2 * np.pi * sigma2) - 0.5 / sigma2 * (r_l @ r_l - 2 * (r_l @ (Xs @ b[l])) + np.sum(d * postb2[l]))
            kl += -lbf_model[l] + ser
        elbo.append(-n / 2 * np.log(2 * np.pi * sigma2) - er2 / (2 * sigma2) - kl)
        if estimate_residual_variance:
            sigma2 = max(er2 / n, 1e-8)
        if abs(elbo[-1] - elbo[-2]) < tol:
            break
    keep = V > 1e-9
    pip = 1 - np.prod(1 - alpha[keep], axis=0) if keep.any() else np.zeros(p)
    sets, seen = [], set()
    if keep.any():
        Xc = Xs / np.sqrt(np.maximum(d, 1e-12))
        for l in np.where(keep)[0]:
            o = np.argsort(-alpha[l])
            c = np.cumsum(alpha[l][o])
            k = int(np.searchsorted(c, coverage)) + 1
            idx = np.sort(o[:k])
            sub = idx if len(idx) <= 200 else np.random.default_rng(0).choice(idx, 200, replace=False)
            purity = float(np.abs(Xc[:, sub].T @ Xc[:, sub]).min()) if len(sub) > 1 else 1.0
            if purity >= min_abs_corr and tuple(idx) not in seen:
                seen.add(tuple(idx))
                sets.append({"effect": int(l), "idx": idx, "purity": purity, "coverage": float(c[k - 1]), "lbf": float(lbf_model[l])})
    return {"alpha": alpha, "mu": mu / csd, "V": V, "sigma2": sigma2, "pip": pip, "lbf_variable": lbf_var, "lbf": lbf_model,
            "elbo": np.array(elbo[1:]), "niter": it + 1, "sets": sets, "converged": abs(elbo[-1] - elbo[-2]) < tol}
