"""cis-eQTL mapping: association, weighted Bonferroni eGene detection and SuSiE fine-mapping.

1. Association. For each gene, accessions with expression (at least --min-samples), SNVs in the
   gene window with MAF >= --maf and call rate >= --min-call in these accessions (with
   --scored-only, only variants with an attention score), missing genotypes mean-imputed. Expression
   and genotypes are whitened with the kinship matrix (h2 by REML per gene), then each variant is
   tested by marginal regression. Every tested variant-gene pair counts towards M.
2. eGenes. Uniform prior: a gene is an eGene if any variant has -log10 P > -log10(alpha / M).
   Attention prior: thresholds -log10(alpha * w / M) with w = exp(lambda * q), mean 1 within each
   gene (q: within-gene rank-normal attention; variants without a score get q = 0).
3. eQTLs. Each eGene is fine-mapped with SuSiE (L effects) on the same whitened data, with prior
   inclusion probabilities 0.8 * w / sum(w) + 0.2 / m (uniform: 1 / m). An eQTL is a 95% credible
   set; a gene can have several eQTLs, and an eQTL with one variant is a single-variant eQTL.
"""
import os
import time
import warnings
from multiprocessing import Pool
import numpy as np
import pandas as pd
from .genotypes import VCF
from .lmm import lmm_transform, marginal_logp
from .prior import gene_weights, susie_pi, thresholds
from .susie import susie
from .windows import read_windows

_G = {}


def read_expression(path):
    E = pd.read_csv(path, sep="\t", index_col=0)
    E.index = E.index.astype(str)
    E.columns = E.columns.astype(str)
    return E.apply(pd.to_numeric, errors="coerce")


def read_kinship(path):
    K = pd.read_csv(path, sep="\t", index_col=0)
    K.index = K.index.astype(str)
    K.columns = K.columns.astype(str)
    return K


def _init(vcf, Y, KA, cols, maf, min_call, min_samples, L):
    _G.update(vcf=VCF(vcf), Y=Y, KA=KA, cols=cols, maf=maf, min_call=min_call, min_samples=min_samples, L=L)


def _prepare(gene, chrom, ws, we, only=None):
    """Whitened genotypes and expression of one gene, or None. only: restrict to these positions."""
    yv = _G["Y"].loc[gene].values
    sel = np.where(~np.isnan(yv))[0]
    if len(sel) < _G["min_samples"]:
        return None
    pos, X = _G["vcf"].region(chrom, ws, we, [_G["cols"][i] for i in sel])
    if not len(pos):
        return None
    m = np.isnan(X)
    called = 1 - m.mean(0)
    if m.any():
        with np.errstate(invalid="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            X = np.where(m, np.nanmean(X, 0)[None, :], X)
    af = np.nan_to_num(X.mean(0) / 2)              # all-missing variants get MAF 0 and are dropped
    keep = (np.minimum(af, 1 - af) >= _G["maf"]) & (called >= _G["min_call"])
    if only is not None:
        keep &= np.isin(pos, only)
    if not keep.any():
        return None
    Xt, yt, h2 = lmm_transform(X[:, keep], yv[sel], _G["KA"][np.ix_(sel, sel)])
    return pos[keep], np.minimum(af, 1 - af)[keep], Xt, yt, h2, len(sel)


def _assoc(job):
    gene, chrom, ws, we = job
    r = _prepare(gene, chrom, ws, we)
    if r is None:
        return None
    pos, maf, Xt, yt, h2, n = r
    return pd.DataFrame({"gene": gene, "chrom": chrom, "pos": pos, "maf": maf.astype(np.float32),
                         "logp": marginal_logp(Xt, yt).astype(np.float32), "n_samples": n, "h2": np.float32(h2)})


def _finemap(job):
    gene, chrom, ws, we, tested, priors = job      # priors: {name: weights aligned to the tested positions}
    pos, _, Xt, yt, _, _ = _prepare(gene, chrom, ws, we, tested)
    assert np.array_equal(pos, tested), gene
    out = []
    for name, w in priors.items():
        res = susie(Xt, yt, L=_G["L"], intercept=False, prior_weights=None if w is None else susie_pi(w))
        cs, purity = np.full(len(pos), -1), np.full(len(pos), np.nan)
        for c in res["sets"]:
            cs[c["idx"]] = c["effect"]
            purity[c["idx"]] = c["purity"]
        out.append(pd.DataFrame({"gene": gene, "prior": name, "pos": pos, "pip": res["pip"].astype(np.float32), "cs": cs,
                                 "purity": purity.astype(np.float32)}))
    return pd.concat(out, ignore_index=True)


def _eqtls(P, A):
    """One row per 95% credible set."""
    rows = []
    lp = A.set_index(["gene", "pos"]).logp
    for (g, prior, c), d in P[P.cs >= 0].groupby(["gene", "prior", "cs"], sort=False):
        d = d.sort_values("pip", ascending=False)
        lead = d.iloc[0]
        rows.append(dict(gene=g, prior=prior, eqtl=f"{g}:{int(c)}", n_variants=len(d), purity=round(float(lead.purity), 4),
                         lead_pos=int(lead.pos), lead_pip=round(float(lead.pip), 4), lead_logp=round(float(lp[(g, int(lead.pos))]), 3),
                         variants=";".join(str(p) for p in sorted(d.pos.astype(int)))))
    cols = ["gene", "prior", "eqtl", "n_variants", "purity", "lead_pos", "lead_pip", "lead_logp", "variants"]
    return pd.DataFrame(rows, columns=cols)


def run(vcf, expression, windows, kinship, out, attention=None, lam=1.0, alpha=0.05, maf=0.05, min_call=0.5, min_samples=30,
        L=10, workers=4, genes=None, scored_only=False, log=print):
    t0 = time.time()
    os.makedirs(out, exist_ok=True)
    E, K, W = read_expression(expression), read_kinship(kinship), read_windows(windows)
    v = VCF(vcf)
    acc = sorted(set(E.columns) & set(K.index) & set(v.samples))
    if len(acc) < min_samples:
        raise ValueError(f"only {len(acc)} accessions are in the expression table, kinship matrix and VCF")
    W = W[W.gene.isin(E.index)]
    if genes:
        W = W[W.gene.isin(set(genes))]
    log(f"map: {len(acc)} accessions, {len(W)} genes with a window and expression")
    Y = E.loc[W.gene.unique(), acc]
    init = (vcf, Y, K.loc[acc, acc].values.astype(np.float64), v.sample_index(acc), maf, min_call, min_samples, L)
    jobs = [(r.gene, r.chrom, int(r.win_start), int(r.win_end)) for r in W.itertuples(index=False)]
    with Pool(workers, _init, init) as pool:
        A = pd.concat([d for d in pool.imap(_assoc, jobs, chunksize=8) if d is not None], ignore_index=True)
    if attention is not None:
        S = pd.read_parquet(attention) if attention.endswith(".parquet") else pd.read_csv(attention, sep="\t")
        S = S[["gene", "pos", "attention"]].drop_duplicates(["gene", "pos"])
        S["gene"] = S.gene.astype(str)
        A = A.merge(S, on=["gene", "pos"], how="left")
        log(f"map: attention for {A.attention.notna().mean():.1%} of variant-gene pairs"
            + (" (testing only these)" if scored_only else " (missing = neutral weight)"))
        if scored_only:
            A = A[A.attention.notna()].reset_index(drop=True)
    M = len(A)
    log(f"map: {A.gene.nunique()} genes tested, M = {M:,} variant-gene tests ({time.time() - t0:.0f} s)")

    priors = {"uniform": np.ones(M)}
    if attention is not None:
        priors["PlantDNA-FM"] = gene_weights(A.attention.values, A.gene.values, lam)
        A["weight"] = priors["PlantDNA-FM"].astype(np.float32)
    egenes = {}
    for name, w in priors.items():
        thr = thresholds(w, M, alpha)
        A[f"threshold_{name}"] = thr.astype(np.float32)
        egenes[name] = set(A.gene[A.logp.values > thr])
        log(f"map: {name}: {len(egenes[name])} eGenes")

    starts = A.groupby("gene", sort=False).indices
    win = W.drop_duplicates("gene").set_index("gene")
    fm_jobs = []
    for g in sorted(set().union(*egenes.values())):
        ix = starts[g]
        fm_jobs.append((g, win.at[g, "chrom"], int(win.at[g, "win_start"]), int(win.at[g, "win_end"]), A.pos.values[ix],
                        {name: (None if name == "uniform" else priors[name][ix]) for name in priors if g in egenes[name]}))
    with Pool(workers, _init, init) as pool:
        P = pd.concat(pool.imap(_finemap, fm_jobs, chunksize=4), ignore_index=True)
    Q = _eqtls(P, A)

    G = A.groupby("gene", sort=False).agg(chrom=("chrom", "first"), n_samples=("n_samples", "first"), h2=("h2", "first"),
                                          n_variants=("pos", "size"))
    lead = A.loc[A.groupby("gene", sort=False).logp.idxmax()].set_index("gene")
    G["lead_pos"], G["lead_logp"] = lead.pos, lead.logp.round(3)
    for name in priors:
        G[f"egene_{name}"] = G.index.isin(egenes[name])
        n_eq = Q[Q.prior == name].groupby("gene").size()
        G[f"eqtls_{name}"] = n_eq.reindex(G.index).fillna(0).astype(int).where(G[f"egene_{name}"], 0)
    rows = []
    for name in priors:
        q = Q[Q.prior == name]
        p = P[P.prior == name]
        rows.append(dict(prior=name, genes_tested=A.gene.nunique(), tests=M, eGenes=len(egenes[name]),
                         eGenes_with_eQTL=q.gene.nunique(), eQTLs=len(q), single_variant_eQTLs=int((q.n_variants == 1).sum()),
                         variants_PIP_gt_0_9=int((p.pip > 0.9).sum()), mean_eQTL_size=round(float(q.n_variants.mean()), 3) if len(q) else np.nan))
    Sm = pd.DataFrame(rows)

    A.to_csv(f"{out}/associations.tsv.gz", sep="\t", index=False, float_format="%.6g")
    G.reset_index().to_csv(f"{out}/genes.tsv", sep="\t", index=False)
    P.to_csv(f"{out}/pip.tsv.gz", sep="\t", index=False, float_format="%.6g")
    Q.to_csv(f"{out}/eqtls.tsv", sep="\t", index=False)
    Sm.to_csv(f"{out}/summary.tsv", sep="\t", index=False)
    log(Sm.to_string(index=False))
    log(f"map: done in {time.time() - t0:.0f} s -> {out}")
    return Sm
