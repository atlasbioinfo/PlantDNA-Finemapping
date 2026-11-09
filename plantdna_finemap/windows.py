"""Gene windows and representative-transcript exons from a GFF3 file.

Each gene gets one 8-kb window, TSS -3 kb ... +5 kb in the direction of transcription
(plus strand [tss - 3000, tss + 4999], minus strand [tss - 4999, tss + 3000]); the TSS is the
5' end of the gene feature. The exons of one representative transcript per gene are the query
rows of the attention readout: <gene>.1 or <gene>_T01 if present, otherwise the transcript with
the longest total exon length.
"""
import gzip
import pandas as pd

UPSTREAM, DOWNSTREAM = 3000, 5000


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _attrs(s):
    return dict(x.split("=", 1) for x in s.split(";") if "=" in x)


def _strip(x):
    return x.split(":", 1)[1] if x.split(":", 1)[0] in ("gene", "transcript") else x


def read_gff(path, biotype="protein_coding"):
    """Return genes (DataFrame) and {gene: [(start, end), ...]} exons of the representative transcript.
    biotype: keep genes with this `biotype` attribute; genes without the attribute are kept. None keeps all."""
    genes, tx, ex = [], {}, {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            a = _attrs(p[8])
            if p[2] == "gene":
                if biotype and a.get("biotype", biotype) != biotype:
                    continue
                genes.append((_strip(a["ID"]), p[0], p[6], int(p[3]), int(p[4])))
            elif p[2] in ("mRNA", "transcript") and "ID" in a and "Parent" in a:
                tx[_strip(a["ID"])] = _strip(a["Parent"].split(",")[0])
            elif p[2] == "exon" and "Parent" in a:
                for t in a["Parent"].split(","):
                    ex.setdefault(_strip(t), []).append((int(p[3]), int(p[4])))
    G = pd.DataFrame(genes, columns=["gene", "chrom", "strand", "gene_start", "gene_end"]).drop_duplicates("gene")
    by_gene = {}
    for t, g in tx.items():
        if t in ex:
            by_gene.setdefault(g, []).append(t)
    exons = {}
    for g in G.gene:
        ts = by_gene.get(g, [])
        if not ts:
            continue
        t = next((c for c in (f"{g}.1", f"{g}_T01") if c in ts), None)
        if t is None:
            t = max(ts, key=lambda x: sum(e - s + 1 for s, e in ex[x]))
        exons[g] = sorted(set(ex[t]))
    return G, exons


def make_windows(gff, biotype="protein_coding", chroms=None):
    """Windows table: gene chrom strand gene_start gene_end tss win_start win_end exons ("start-end,..." genomic)."""
    G, exons = read_gff(gff, biotype)
    if chroms:
        G = G[G.chrom.isin(chroms)]
    plus = G.strand == "+"
    G = G.assign(tss=G.gene_start.where(plus, G.gene_end))
    G["win_start"] = (G.tss - UPSTREAM).where(plus, G.tss - DOWNSTREAM + 1)
    G["win_end"] = (G.tss + DOWNSTREAM - 1).where(plus, G.tss + UPSTREAM)
    G["exons"] = [",".join(f"{s}-{e}" for s, e in exons.get(g, [])) for g in G.gene]
    return G.reset_index(drop=True)


def parse_exons(s):
    return [tuple(int(v) for v in x.split("-")) for x in s.split(",")] if isinstance(s, str) and s else []


def read_windows(path):
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "gene": str}, keep_default_na=False)
