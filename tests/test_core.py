import numpy as np
import pandas as pd

from plantdna_finemap.attention import aggregate, tile
from plantdna_finemap.genes import make_genes, parse_exons
from plantdna_finemap.genotypes import Plink, decode_gt, resolve_contig
from plantdna_finemap.prior import gene_weights, rank_q, susie_pi, thresholds
from plantdna_finemap.susie import susie
from plantdna_finemap.xfit import default_folds, gene_codes, weights, xfit


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
    assert np.allclose(w, weights(rank_q(s, g), gene_codes(g), 1.0))      # same weights as the cross-fitting code
    assert np.allclose(thresholds(np.ones(4), 4, 0.05), -np.log10(0.05 / 4))
    pi = susie_pi(w[:5])
    assert np.isclose(pi.sum(), 1) and pi.min() >= 0.2 / 5 - 1e-12


def test_xfit_uses_the_other_fold():
    rng = np.random.default_rng(0)
    n_genes, m = 400, 50
    g = np.repeat(np.arange(n_genes), m)
    score = rng.standard_normal(len(g))
    q = rank_q(score, g)
    logp = -np.log10(rng.uniform(size=len(g)))
    top = q > 1.5                                    # informative prior: signals sit at high-score variants
    hit = (np.arange(len(g)) // m % 2 == 0) & top
    logp[hit] += 4.8
    fold = (g >= n_genes // 2).astype(int)
    w, lam, curve = xfit(logp, q, g, fold)
    assert lam[0] > 0 and lam[1] > 0
    for f in (0, 1):                                 # each fold's lambda maximises eGenes on the other fold
        other = curve[f"eGenes_fold{1 - f}"]
        best = curve["lambda"][other == other.max()].values
        assert lam[f] == best[len(best) // 2]                 # ties: median (upper middle) of the tied values
    assert np.allclose(pd.Series(w).groupby(g).mean(), 1)
    # no eGenes for any lambda: all 81 grid values tie, the median is 2.0
    w0, lam0, _ = xfit(np.zeros(len(g)), q, g, fold)
    assert lam0 == {0: 2.0, 1: 2.0}


def test_default_folds():
    assert default_folds(["5", "1", "3", "2", "4"]) == {"1": 0, "2": 0, "3": 0, "4": 1, "5": 1}
    f = default_folds([f"chr{i}" for i in range(10, 0, -1)])
    assert [c for c in f if f[c] == 0] == [f"chr{i}" for i in range(1, 6)]


def test_susie_recovers_causal_variant():
    rng = np.random.default_rng(1)
    n, p = 300, 60
    Z = rng.standard_normal((n, 6))
    X = np.stack([(0.9 * Z[:, j // 10] + 0.45 * rng.standard_normal(n) > 0).astype(float) for j in range(p)], 1)
    X += rng.random((n, p)) < 0.2
    y = X[:, 17] * 0.8 + rng.standard_normal(n)
    res = susie(X, y, L=5)
    assert any(17 in c["idx"] for c in res["sets"])
    w = np.ones(p)
    w[17] = 20
    res2 = susie(X, y, L=5, prior_weights=susie_pi(w))
    assert res2["pip"][17] >= res["pip"][17]


def test_tile():
    s = tile(20001, 22000, np.array([15000, 31000]), window=16000, step=8000, flank=10000)
    assert s == [10001, 16001]                       # 10001..32000 covered, last window ends at the end of the region
    s = tile(20001, 22000, np.array([5000, 40000]), window=16000, step=8000, flank=10000)
    assert s[0] == 5000 and s[-1] + 16000 - 1 == 40000
    assert tile(100, 200, np.array([150]), window=16000, step=8000, flank=1000) == [-900]    # region shorter than a window


def test_aggregate_prefers_exon_windows():
    W = pd.DataFrame({"gene": "g", "pos": [1, 1, 1, 2, 2], "window": [0, 1, 2, 1, 2],
                      "has_exon": [True, True, False, False, False], "attn_win": [1.0, 3.0, 10.0, 4.0, 6.0]})
    S = aggregate(W).set_index("pos")
    assert S.at[1, "attention"] == 2.0 and S.at[1, "attention_source"] == "exon windows"
    assert S.at[2, "attention"] == 5.0 and S.at[2, "attention_source"] == "windows without exon"


def test_plink_reader(tmp_path):
    X = np.array([[0, 1, 2, np.nan, 2], [2, 2, 0, 1, np.nan]])          # 2 variants x 5 samples, dosage of A1
    code = np.where(np.isnan(X), 1, np.where(X == 2, 0, np.where(X == 1, 2, 3))).astype(np.uint8)
    code = np.concatenate([code, np.zeros((2, 3), np.uint8)], 1).reshape(2, -1, 4)
    b = (code[:, :, 0] | code[:, :, 1] << 2 | code[:, :, 2] << 4 | code[:, :, 3] << 6).astype(np.uint8)
    (tmp_path / "g.bed").write_bytes(bytes([0x6C, 0x1B, 0x01]) + b.tobytes())
    (tmp_path / "g.bim").write_text("1\ta\t0\t100\tA\tG\n1\tb\t0\t200\tC\tT\n")
    (tmp_path / "g.fam").write_text("".join(f"s{i} s{i} 0 0 0 -9\n" for i in range(5)))
    G = Plink(str(tmp_path / "g"))
    pos, D = G.fetch("chr1", [200, 100, 300], ["s4", "s0", "s3"])
    assert list(pos) == [200, 100]
    assert np.array_equal(np.nan_to_num(D, nan=-1), [[-1, 2], [2, 0], [1, -1]])


def test_genes(tmp_path):
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
    G = make_genes(str(gff)).set_index("gene")
    assert list(G.index) == ["G1", "G2"]
    assert (G.at["G2", "strand"], G.at["G2", "gene_start"], G.at["G2", "gene_end"]) == ("-", 9000, 9900)
    assert parse_exons(G.at["G1", "exons"]) == [(5000, 5100), (6000, 6500)]      # <gene>.1 preferred
    assert G.at["G2", "exons"] == ""                     # no <gene>.1 / <gene>_T01 transcript: all positions are queries
    assert "R1" in set(make_genes(str(gff), biotype=None).gene)


def test_example_gene_without_dot1_uses_all_positions():
    """AT4G00520 has only .3/.4 transcripts: no exons, so every window queries all positions, as in the paper."""
    from pathlib import Path
    from plantdna_finemap.attention import window_sequence
    gff = Path(__file__).resolve().parents[1] / "examples/arabidopsis/TAIR10_example.gff3.gz"
    G = make_genes(str(gff), biotype=None).set_index("gene")
    assert G.at["AT4G00520", "exons"] == "" and G.at["AT5G10140", "exons"] != ""
    assert 0 < (G.exons == "").sum() < 10                # the few genes without a .1 transcript

    class Fa:
        def fetch(self, c, a, b):
            return "A" * (b - a)
    _, ex = window_sequence(Fa(), {"4": 10 ** 6}, "4", 220000, 16000, "+", parse_exons(G.at["AT4G00520", "exons"]))
    assert not ex.any()


def test_plink_gzipped(tmp_path):
    import gzip
    from plantdna_finemap.genotypes import open_genotypes
    with gzip.open(tmp_path / "g.bed.gz", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01, 0b11100100]))          # 4 samples: hom A1, missing, het, hom A2
    with gzip.open(tmp_path / "g.bim.gz", "wt") as f:
        f.write("1\t1_50\t0\t50\tA\tC\n")
    (tmp_path / "g.fam").write_text("".join(f"s{i} s{i} 0 0 0 -9\n" for i in range(4)))
    g = open_genotypes(str(tmp_path / "g"))
    P, X = g.fetch("1", [50], ["s0", "s1", "s2", "s3"])
    assert list(P) == [50] and np.array_equal(np.nan_to_num(X[:, 0], nan=-1), [2, -1, 1, 0])
    assert g.alleles("Chr1", [50, 60]) == {50: ("C", "A")}             # REF = A2, ALT = A1
