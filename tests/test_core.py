import numpy as np
import pandas as pd

from plantdna_finemap.genotypes import decode_gt, resolve_contig
from plantdna_finemap.lmm import lmm_transform, marginal_logp
from plantdna_finemap.prior import gene_weights, rank_q, susie_pi, thresholds
from plantdna_finemap.susie import susie
from plantdna_finemap.windows import make_windows, parse_exons


def test_decode_gt():
    d = decode_gt(["0|0:1:2", "0/1", "1|1:3", ".:.", "./.", "1", "0"])
    assert np.array_equal(np.nan_to_num(d, nan=-1), [0, 1, 2, -1, -1, 2, 0])
    assert np.array_equal(decode_gt(["9:0/1", "9:1/1"], gt_index=1), [1, 2])


def test_resolve_contig():
    assert resolve_contig("5", ["Chr5", "Chr4"]) == "Chr5"
    assert resolve_contig("chr10", ["10"]) == "10"
    assert resolve_contig("3", ["1"]) is None


def test_weights_keep_gene_budget():
    g = np.repeat(["a", "b"], [5, 3])
    s = np.array([5, 1, 3, np.nan, 2, 0.1, 0.3, 0.2])
    w = gene_weights(s, g, lam=1.0)
    assert np.allclose(pd.Series(w).groupby(g).mean(), 1)
    assert rank_q(s, g)[3] == 0                      # missing score -> neutral
    assert w[0] == w[:5].max()
    assert np.allclose(thresholds(np.ones(4), 4, 0.05), -np.log10(0.05 / 4))
    pi = susie_pi(w[:5])
    assert np.isclose(pi.sum(), 1) and pi.min() >= 0.2 / 5 - 1e-12


def test_susie_and_lmm_recover_causal_variant():
    rng = np.random.default_rng(1)
    n, p = 300, 60
    Z = rng.standard_normal((n, 6))
    X = np.stack([(0.9 * Z[:, j // 10] + 0.45 * rng.standard_normal(n) > 0).astype(float) for j in range(p)], 1)
    X += rng.random((n, p)) < 0.2
    K = np.corrcoef(X)
    y = X[:, 17] * 0.8 + rng.standard_normal(n)
    Xt, yt, h2 = lmm_transform(X, y, K)
    assert 0 < h2 < 1
    lp = marginal_logp(Xt, yt)
    assert np.argmax(lp) in range(10, 20)
    res = susie(Xt, yt, L=5, intercept=False)
    assert any(17 in c["idx"] for c in res["sets"])
    w = np.ones(p)
    w[17] = 20
    res2 = susie(Xt, yt, L=5, intercept=False, prior_weights=susie_pi(w))
    assert res2["pip"][17] >= res["pip"][17]


def test_windows(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "1\tx\tgene\t5000\t7000\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n"
        "1\tx\tmRNA\t5000\t7000\t.\t+\t.\tID=transcript:G1.2;Parent=gene:G1\n"
        "1\tx\texon\t5000\t7000\t.\t+\t.\tParent=transcript:G1.2\n"
        "1\tx\tmRNA\t5000\t6500\t.\t+\t.\tID=transcript:G1.1;Parent=gene:G1\n"
        "1\tx\texon\t5000\t5100\t.\t+\t.\tParent=transcript:G1.1\n"
        "1\tx\texon\t6000\t6500\t.\t+\t.\tParent=transcript:G1.1\n"
        "1\tx\tgene\t9000\t9900\t.\t-\t.\tID=gene:G2;biotype=protein_coding\n"
        "1\tx\tmRNA\t9000\t9900\t.\t-\t.\tID=transcript:G2_T02;Parent=gene:G2\n"
        "1\tx\texon\t9000\t9900\t.\t-\t.\tParent=transcript:G2_T02\n"
        "1\tx\tgene\t100\t200\t.\t+\t.\tID=gene:R1;biotype=ncRNA\n")
    W = make_windows(str(gff)).set_index("gene")
    assert list(W.index) == ["G1", "G2"]
    assert (W.at["G1", "win_start"], W.at["G1", "win_end"]) == (2000, 9999)
    assert (W.at["G2", "tss"], W.at["G2", "win_start"], W.at["G2", "win_end"]) == (9900, 4901, 12900)
    assert parse_exons(W.at["G1", "exons"]) == [(5000, 5100), (6000, 6500)]      # <gene>.1 preferred
    assert parse_exons(W.at["G2", "exons"]) == [(9000, 9900)]                    # otherwise longest
