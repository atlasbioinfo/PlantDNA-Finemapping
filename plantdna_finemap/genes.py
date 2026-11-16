"""Gene table and representative-transcript exons from a GFF3 file.

One row per gene: chrom, strand, gene_start, gene_end (1-based, inclusive, from the gene feature) and
the exons of the representative transcript <gene>.1 (A. thaliana) or <gene>_T01 (maize), which are the query
positions of the attention readout. Genes without such a transcript get no exons, so every window of the gene
uses all positions as queries (as in the paper).
"""
import gzip
import pandas as pd


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
        if t is not None:
            exons[g] = sorted(set(ex[t]))
    return G, exons


def make_genes(gff, biotype="protein_coding", chroms=None):
    """Gene table: gene chrom strand gene_start gene_end exons ("start-end,..." genomic, 1-based)."""
    G, exons = read_gff(gff, biotype)
    if chroms:
        G = G[G.chrom.isin(chroms)]
    G = G.assign(exons=[",".join(f"{s}-{e}" for s, e in exons.get(g, [])) for g in G.gene])
    return G.reset_index(drop=True)


def parse_exons(s):
    return [tuple(int(v) for v in x.split("-")) for x in s.split(",")] if isinstance(s, str) and s else []


def read_genes(path):
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "gene": str}, keep_default_na=False)
