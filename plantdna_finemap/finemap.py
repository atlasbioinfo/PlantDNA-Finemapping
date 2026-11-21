"""expGWAS fine-mapping: weighted Bonferroni eGene detection and SuSiE, under the uniform and the attention prior.

Input: one row per association test (gene, chrom, pos, P or -log10 P) from any GWAS tool, e.g. GEMMA run per
gene on the variants in gene body +/- 10 kb, as in Voichek et al.; every row counts towards M.
1. Weights. q = within-gene rank-normal attention score (variants without a score: q = 0, neutral);
   w = exp(lambda * q), mean 1 within each gene. lambda by two-fold cross-fitting over chromosome halves
   (xfit.py), or fixed.
2. eGenes. A variant is significant if -log10 P > -log10(alpha * w / M) (uniform prior: w = 1). Because w has
   mean 1 within each gene, the family-wise error rate stays at alpha. An eGene has >= 1 significant variant.
3. Fine-mapping. SuSiE (L effects) on the significant variants of each eGene (mean-imputed genotypes,
   expression as is), prior inclusion probabilities pi = 0.8 * w / sum(w) + 0.2 / m over these m variants
   (uniform prior: 1 / m). An expGWAS signal is a 95% credible set with purity >= 0.5.
Single-gene mode (candidate-gene analysis): only the tests of one gene, M = its number of tests, fixed lambda.
"""
import os
import time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from .genotypes import open_genotypes, resolve_contig
from .prior import gene_weights, rank_q, susie_pi, thresholds
from .susie import susie
from .xfit import default_folds, gene_codes, xfit

PRIOR = "PlantDNA-FM"
_G = {}


def _table(path):
    return pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, sep=None if path.endswith(".txt") else "\t", engine="python")


def read_tests(path, logp=True):
    """Association table -> gene, chrom, pos (, logp). Accepted columns: gene; chrom or chr; pos or ps;
    logp (-log10 P) or p / p_wald / p_lrt / p_score (the first found)."""
    T = _table(path)
    ren = {c: c.lower() for c in T.columns}
    T = T.rename(columns=ren).rename(columns={"chr": "chrom", "ps": "pos", "-log10p": "logp"})
    for c in ("gene", "chrom", "pos"):
        if c not in T.columns:
            raise ValueError(f"{path}: column '{c}' is missing")
    T["gene"], T["chrom"], T["pos"] = T.gene.astype(str), T.chrom.astype(str), T.pos.astype(np.int64)
    if not logp:
        return T[["gene", "chrom", "pos"]]
    if "logp" not in T.columns:
        pc = next((c for c in ("p", "p_wald", "p_lrt", "p_score", "pvalue", "p_value") if c in T.columns), None)
        if pc is None:
            raise ValueError(f"{path}: needs a column logp (-log10 P) or p (P value)")
        T["logp"] = -np.log10(np.maximum(T[pc].astype(float), 1e-300))
    n0 = len(T)
    T = T.dropna(subset=["logp"]).drop_duplicates(["gene", "pos"]).reset_index(drop=True)
    if len(T) < n0:
        print(f"tests: {n0 - len(T)} rows without P or duplicated (gene, pos) dropped", flush=True)
    keep = ["gene", "chrom", "pos", "logp"] + [c for c in ("maf", "af", "ref", "alt") if c in T.columns]
    return T[keep]


def read_expression(path):
    E = pd.read_csv(path, sep="\t", index_col=0)
    E.index = E.index.astype(str)
    E.columns = E.columns.astype(str)
    return E.apply(pd.to_numeric, errors="coerce")


def _init(genotypes, Y, L):
    _G.update(geno=open_genotypes(genotypes), Y=Y, L=L)


def _finemap(job):
    """SuSiE for one gene under one or more priors. job: gene, chrom, {prior: (positions, pi or None)}."""
    gene, chrom, fits = job
    y = _G["Y"].loc[gene]
    y = y[y.notna()]
    out = []
    for name, (pos, pi) in fits.items():
        found, X = _G["geno"].fetch(chrom, pos, list(y.index))
        if not len(found):
            continue
        m = np.isnan(X)
        if m.any():                                   # mean imputation
            mu = np.nan_to_num(np.nansum(X, 0) / np.maximum((~m).sum(0), 1))
            X = np.where(m, mu[None, :], X)
        pw = None if pi is None else pd.Series(pi, index=pos).loc[found].values
        res = susie(X, y.values, L=_G["L"], prior_weights=pw)
        cs, purity = np.full(len(found), -1), np.full(len(found), np.nan)
        for c in res["sets"]:
            cs[c["idx"]] = c["effect"]
            purity[c["idx"]] = c["purity"]
        out.append(pd.DataFrame({"gene": gene, "prior": name, "pos": found, "pip": res["pip"].astype(np.float32), "cs": cs,
                                 "purity": purity.astype(np.float32), "n_samples": len(y), "converged": res["converged"]}))
    return pd.concat(out, ignore_index=True) if out else None


def _signals(P, T):
    """One row per expGWAS signal (95% credible set)."""
    rows = []
    lp = T.set_index(["gene", "pos"]).logp
    for (g, prior, c), d in P[P.cs >= 0].groupby(["gene", "prior", "cs"], sort=False):
        d = d.sort_values("pip", ascending=False)
        lead = d.iloc[0]
        rows.append(dict(gene=g, prior=prior, signal=f"{g}:{int(c)}", n_variants=len(d), purity=round(float(lead.purity), 4),
                         lead_pos=int(lead.pos), lead_pip=round(float(lead.pip), 4), lead_logp=round(float(lp[(g, int(lead.pos))]), 3),
                         variants=";".join(str(p) for p in sorted(d.pos.astype(int)))))
    cols = ["gene", "prior", "signal", "n_variants", "purity", "lead_pos", "lead_pip", "lead_logp", "variants"]
    return pd.DataFrame(rows, columns=cols)


def result_table(T, P, egenes, geno, genes=None, prior=None):
    """One row per variant that is significant with the uniform prior or with the attention prior (the table of the paper)."""
    sig = T.sig_uniform if prior is None else (T.sig_uniform | T[f"sig_{prior}"])
    F = T[sig].copy()
    F["n_tested_in_gene"] = T.groupby("gene").pos.transform("size")[sig]
    F["rs"] = F.chrom.astype(str) + "_" + F.pos.astype(str)
    al = {}
    if "ref" not in F.columns or "alt" not in F.columns:
        for c, d in F.groupby("chrom", sort=False):
            al.update({(c, p): a for p, a in geno.alleles(c, d.pos.unique()).items()})
    if "ref" not in F.columns or "alt" not in F.columns:      # alleles from the association table if given, else from the genotypes
        F["ref"] = [al.get((c, p), ("", ""))[0] for c, p in zip(F.chrom, F.pos)]
        F["alt"] = [al.get((c, p), ("", ""))[1] for c, p in zip(F.chrom, F.pos)]
    if "maf" not in F.columns:                                # GEMMA writes the allele frequency (af)
        F["maf"] = np.minimum(F.af, 1 - F.af).round(3) if "af" in F.columns else np.nan
    pip = lambda name: P[P.prior == name].drop_duplicates(["gene", "pos"]).set_index(["gene", "pos"]).pip.astype(float)
    key = pd.MultiIndex.from_arrays([F.gene, F.pos])
    F["sig_raw"], F["eGene_raw"] = F.sig_uniform, F.gene.isin(egenes["uniform"])
    F["PIP_raw"] = pip("uniform").reindex(key).round(5).values
    F["raw_pval (-log10)"] = F.logp.round(4)
    cols = ["gene", "rs", "ref", "alt", "maf", "raw_pval (-log10)", "n_tested_in_gene"]
    if prior is not None:
        F["FM_prior_weight"], F["threshold_FM (-log10)"] = F.weight.round(4), F[f"threshold_{prior}"].round(4)
        F["sig_FM"], F["eGene_FM"] = F[f"sig_{prior}"], F.gene.isin(egenes[prior])
        F["FM_new_sig"] = F.sig_FM & ~F.sig_raw
        F["PIP_FM"] = pip(prior).reindex(key).round(5).values
        cols += ["FM_prior_weight", "threshold_FM (-log10)", "sig_raw", "eGene_raw", "sig_FM", "eGene_FM", "FM_new_sig", "PIP_raw", "PIP_FM"]
    else:
        cols += ["sig_raw", "eGene_raw", "PIP_raw"]
    if genes is not None:
        G = pd.read_csv(genes, sep="\t", dtype=str, keep_default_na=False).drop_duplicates("gene").set_index("gene")
        for c in ("gene_name", "strand", "gene_start", "gene_end"):
            F[c] = F.gene.map(G[c]) if c in G.columns else ""
        F.loc[F.gene_name == F.gene, "gene_name"] = ""          # no name in the annotation
        cols = ["gene", "gene_name", "strand", "gene_start", "gene_end"] + cols[1:]
    return F.sort_values(["gene", "rs"])[cols].reset_index(drop=True)


def _lambdas(T, lam, folds, M, alpha, out, log):
    """Weights of the attention prior: fixed lambda (one value, or one per fold) or two-fold cross-fitting."""
    if lam is not None and len(lam) == 1:
        return gene_weights(T.attention.values, T.gene.values, lam[0]), {0: lam[0], 1: lam[0]}
    fmap = folds if folds is not None else default_folds(T.chrom.unique())
    fold = T.chrom.map({c: fmap.get(resolve_contig(c, list(fmap))) for c in T.chrom.unique()})
    if fold.isna().any():
        raise ValueError(f"chromosomes without a fold: {sorted(T.chrom[fold.isna()].unique())}")
    fold = fold.astype(int).values
    T["fold"] = fold
    log(f"finemap: folds {' / '.join(','.join(str(c) for c in fmap if fmap[c] == f) for f in (0, 1))}")
    if lam is not None:
        lam = {0: lam[0], 1: lam[1]}
        w = np.empty(len(T))
        for f in (0, 1):
            wf = gene_weights(T.attention.values, T.gene.values, lam[f])
            w[fold == f] = wf[fold == f]
        return w, lam
    if len(set(fold)) < 2:
        raise ValueError("cross-fitting needs tests on both folds; give --lambda for a single fold or a single gene")
    w, lam, curve = xfit(T.logp.values, T.q.values, gene_codes(T.gene.values), fold, M, alpha)
    curve.to_csv(f"{out}/lambda_curve.tsv", sep="\t", index=False)
    log(f"finemap: cross-fitted lambda: fold 0 {lam[0]}, fold 1 {lam[1]} (each learned on the other fold)")
    return w, lam


def run(tests, expression, genotypes, out, attention=None, lam=None, folds=None, alpha=0.05, n_tests=None, L=10,
        single_gene=None, susie_all=False, workers=4, log=print, genes=None, check=True):
    t0 = time.time()
    if check:                                              # stop before any computation if the inputs do not fit together
        from .check import Report, check_finemap, check_tests
        log("finemap: checking the input files")
        R = Report(log)
        Tc = check_tests(R, tests)
        if Tc is not None:
            check_finemap(R, Tc, expression, genotypes, attention, genes, single_gene, lam, folds)
        R.finish("finemap inputs")
    os.makedirs(out, exist_ok=True)
    T = read_tests(tests)
    if single_gene:
        T = T[T.gene == single_gene].reset_index(drop=True)
        if not len(T):
            raise ValueError(f"no tests for gene {single_gene}")
        if attention is not None and lam is None:
            raise ValueError("single-gene mode needs a fixed --lambda (for example the cross-fitted value of a genome-wide run)")
    M = len(T) if n_tests is None or single_gene else int(n_tests)
    log(f"finemap: {T.gene.nunique():,} genes, {len(T):,} tests, M = {M:,}" + (f" (single gene {single_gene})" if single_gene else ""))

    priors = {"uniform": np.ones(len(T))}
    lams = {}
    if attention is not None:
        S = _table(attention)
        S = S.assign(gene=S.gene.astype(str), pos=S.pos.astype(np.int64))[["gene", "pos", "attention"]].drop_duplicates(["gene", "pos"])
        T = T.merge(S, on=["gene", "pos"], how="left")
        log(f"finemap: attention for {T.attention.notna().mean():.2%} of tests (missing = neutral weight)")
        T["q"] = rank_q(T.attention.values, T.gene.values)
        w, lams = _lambdas(T, lam, folds, M, alpha, out, log)
        priors[PRIOR] = w
        T["weight"] = w
        pd.DataFrame({"fold": list(lams), "lambda": list(lams.values())}).to_csv(f"{out}/lambda.tsv", sep="\t", index=False)
    egenes, fits = {}, {}
    for name, w in priors.items():
        thr = thresholds(w, M, alpha)
        T[f"threshold_{name}"] = thr
        sig = T.logp.values > thr
        T[f"sig_{name}"] = sig
        egenes[name] = set(T.gene[sig])
        log(f"finemap: {name}: {len(egenes[name]):,} eGenes, {int(sig.sum()):,} significant variants")
        sub = T[T.gene.isin(egenes[name])] if susie_all else T[sig]
        for g, d in sub.groupby("gene", sort=False):
            ww = w[d.index.values]
            fits.setdefault(g, (d.chrom.iloc[0], {}))[1][name] = (d.pos.values, None if name == "uniform" else susie_pi(ww))

    E = read_expression(expression)
    geno = open_genotypes(genotypes)
    acc = [s for s in E.columns if s in set(geno.samples)]
    log(f"finemap: {len(acc)} accessions with expression and genotypes; SuSiE on {len(fits):,} eGenes")
    missing = sorted(set(fits) - set(E.index))
    if missing:
        log(f"finemap: {len(missing)} eGenes have no expression and are not fine-mapped")
    jobs = [(g, c, f) for g, (c, f) in fits.items() if g in E.index]
    Y = E.loc[[g for g, _, _ in jobs], acc]
    with Pool(workers, _init, (genotypes, Y, L)) as pool:
        res = [d for d in pool.imap(_finemap, jobs, chunksize=4) if d is not None]
    cols = ["gene", "prior", "pos", "pip", "cs", "purity", "n_samples", "converged"]
    P = pd.concat(res, ignore_index=True) if res else pd.DataFrame(columns=cols)
    Q = _signals(P, T)

    G = T.groupby("gene", sort=False).agg(chrom=("chrom", "first"), n_tests=("pos", "size"))
    lead = T.loc[T.groupby("gene", sort=False).logp.idxmax()].set_index("gene")
    G["lead_pos"], G["lead_logp"] = lead.pos, lead.logp.round(3)
    for name in priors:
        G[f"egene_{name}"] = G.index.isin(egenes[name])
        G[f"signals_{name}"] = Q[Q.prior == name].groupby("gene").size().reindex(G.index).fillna(0).astype(int)
    rows = []
    for name in priors:
        q, p = Q[Q.prior == name], P[P.prior == name]
        rows.append({"prior": name, "lambda": "-" if name == "uniform" else " / ".join(dict.fromkeys(str(lams[f]) for f in (0, 1))),
                     "genes_tested": T.gene.nunique(), "tests": len(T), "M": M, "eGenes": len(egenes[name]),
                     "eGenes_with_signal": q.gene.nunique(), "expGWAS_signals": len(q), "single_variant_signals": int((q.n_variants == 1).sum()),
                     "variants_PIP_gt_0_9": int((p.pip > 0.9).sum()), "mean_signal_size": round(float(q.n_variants.mean()), 3) if len(q) else np.nan})
    Sm = pd.DataFrame(rows)

    T.drop(columns=[c for c in ("q",) if c in T.columns]).to_csv(f"{out}/tests.tsv.gz", sep="\t", index=False, float_format="%.6g")
    G.reset_index().to_csv(f"{out}/genes.tsv", sep="\t", index=False)
    P.to_csv(f"{out}/pip.tsv.gz", sep="\t", index=False, float_format="%.6g")
    Q.to_csv(f"{out}/signals.tsv", sep="\t", index=False)
    Sm.to_csv(f"{out}/summary.tsv", sep="\t", index=False)
    F = result_table(T, P, egenes, geno, genes, PRIOR if attention is not None else None)
    F.to_csv(f"{out}/finemapping.csv", index=False)
    log(f"finemap: finemapping.csv: {len(F):,} significant variants in {F.gene.nunique():,} eGenes")
    log("\n" + Sm.to_string(index=False))
    log(f"finemap: done in {time.time() - t0:.0f} s -> {out}")
    return Sm
