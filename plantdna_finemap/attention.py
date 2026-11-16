"""Per-variant attention scores from PlantDNA16K-FM (GPU recommended).

For each gene, the region gene body +/- 10 kb, extended to cover every tested variant of the gene, is tiled
with 16-kb windows at a step of 8 kb (the last window ends at the end of the region). Each window is read
in the direction of transcription (minus-strand genes as reverse complement; positions beyond the
chromosome ends as N) and run through the model once. In every layer, attention from the exon positions
of the gene's representative transcript (all positions if the window contains no exon of the gene) is
averaged over heads and query positions; the window score of a position is this attention to it,
averaged over all layers (36 layers x 12 heads = 432 attention matrices).

Score of a variant: the mean of its window scores over the windows that contain an exon of the gene,
otherwise over the windows that contain none (attention_source). Output: attention.parquet with columns
gene, pos, attention, attention_source; per-window values are kept in chunks/ (rerunning resumes).
"""
import glob
import os
import re
import time
import numpy as np
import pandas as pd
import pysam
from .genes import parse_exons, read_genes
from .genotypes import resolve_contig

_RC = str.maketrans("ACGTN", "TGCAN")
_NON = re.compile("[^ACGT]")


class AttentionReader:
    def __init__(self, model, device="cuda", query_chunk=2048):
        import torch
        from transformers import AttentionInterface, AutoModel, AutoTokenizer
        from transformers.utils import logging as hf_logging
        hf_logging.set_verbosity_error()          # the unused masked-LM head is reported as unexpected keys
        from transformers.integrations.sdpa_attention import sdpa_attention_forward
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        self.torch, self.device, self.cur = torch, device, {}
        cur = self.cur

        def query_rows(module, query, key, value, attention_mask, scaling=1.0, dropout=0.0, **kw):
            # query / key already carry the scaling and rotary embedding: (batch, heads, length, dim)
            qi, k = cur["q"], key[0].float().transpose(-1, -2)
            s = torch.zeros(key.shape[2], dtype=torch.float64, device=key.device)
            for i in range(0, len(qi), query_chunk):          # in chunks of query rows to bound GPU memory
                A = torch.softmax((query[0][:, qi[i:i + query_chunk]].float() @ k) * scaling, dim=-1)
                s += A.sum((0, 1)).double()
            cur["sum"] += s / (query.shape[1] * len(qi))
            return sdpa_attention_forward(module, query, key, value, attention_mask, scaling=scaling, dropout=dropout, **kw)

        AttentionInterface.register("pdfm_query_rows", query_rows)
        tok = AutoTokenizer.from_pretrained(model)
        v = tok.get_vocab()
        self.tok = {b: v[b] for b in "ACGTN"}
        self.cls, self.eos = tok.cls_token_id, tok.eos_token_id
        self.net = AutoModel.from_pretrained(model, dtype=torch.float32, attn_implementation="pdfm_query_rows",
                                             add_pooling_layer=False).to(device).eval()
        self.n_layers = self.net.config.num_hidden_layers

    def window(self, seq, exon_mask):
        """Layer-mean query attention to every position of seq (transcription orientation)."""
        torch = self.torch
        q = np.where(exon_mask)[0] if exon_mask.any() else np.arange(len(seq))
        self.cur.update(q=torch.tensor(q + 1, device=self.device), sum=torch.zeros(len(seq) + 2, dtype=torch.float64, device=self.device))
        ids = torch.tensor([[self.cls] + [self.tok[b] for b in seq] + [self.eos]], device=self.device)
        with torch.no_grad():
            self.net(input_ids=ids)
        return (self.cur["sum"] / self.n_layers).cpu().numpy()[1:-1]


def tile(gene_start, gene_end, pos, window=16000, step=8000, flank=10000):
    """Window starts covering min(gene_start - flank, first variant) .. max(gene_end + flank, last variant)."""
    lo, hi = min(int(gene_start) - flank, int(np.min(pos))), max(int(gene_end) + flank, int(np.max(pos)))
    starts = list(range(lo, max(lo, hi - window + 1) + 1, step))
    if starts[-1] + window - 1 < hi:
        starts.append(hi - window + 1)
    return starts


def window_sequence(fa, lengths, chrom, ws, window, strand, exons):
    """Sequence and exon mask of the window [ws, ws + window - 1] in the direction of transcription."""
    we, L = ws + window - 1, lengths[chrom]
    s = _NON.sub("N", fa.fetch(chrom, max(ws, 1) - 1, min(we, L)).upper()) if we >= 1 and ws <= L else ""
    s = ("N" * max(0, 1 - ws) + s + "N" * max(0, we - L))[:window].ljust(window, "N")
    ex = np.zeros(window, bool)
    for a, b in exons:
        lo, hi = max(a, ws) - ws, min(b, we) - ws
        if hi >= lo:
            ex[lo:hi + 1] = True
    if strand == "-":
        return s.translate(_RC)[::-1], ex[::-1]
    return s, ex


def aggregate(W):
    """Per-window rows (gene, pos, has_exon, attn_win) -> one score per variant."""
    k = ["gene", "pos"]
    ex = W[W.has_exon].groupby(k).attn_win.mean().rename("ex")
    al = W[~W.has_exon].groupby(k).attn_win.mean().rename("al")
    S = pd.concat([ex, al], axis=1).reset_index()
    S["attention"] = S.ex.fillna(S.al).astype(np.float32)
    S["attention_source"] = np.where(S.ex.notna(), "exon windows", "windows without exon")
    return S[["gene", "pos", "attention", "attention_source"]]


def run(model, genes, fasta, tests, out, window=16000, step=8000, flank=10000, only=None, device="cuda", log=print):
    from .finemap import read_tests
    G = read_genes(genes).drop_duplicates("gene").set_index("gene")
    T = read_tests(tests, logp=False)
    T = T[T.gene.isin(G.index)]
    if only:
        T = T[T.gene.isin(set(only))]
    P = {g: np.unique(d.values) for g, d in T.groupby("gene").pos}
    os.makedirs(f"{out}/chunks", exist_ok=True)
    done = set()
    for f in glob.glob(f"{out}/chunks/*.parquet"):
        done |= set(pd.read_parquet(f, columns=["gene"]).gene)
    todo = [g for g in P if g not in done]
    fa = pysam.FastaFile(fasta)
    lengths = {c: fa.get_reference_length(c) for c in fa.references}
    log(f"attention: {len(todo)} genes to score ({len(done)} done); window {window} nt, step {step} nt, flank {flank} nt")
    reader = AttentionReader(model, device) if todo else None
    rows, k, nw, t0 = [], len(glob.glob(f"{out}/chunks/*.parquet")), 0, time.time()
    for i, g in enumerate(todo):
        r, pos = G.loc[g], P[g]
        c = resolve_contig(r.chrom, lengths)
        if c is None:
            log(f"attention: chromosome {r.chrom} of {g} is not in the FASTA; skipped")
            continue
        exons = parse_exons(r.exons)
        for wi, ws in enumerate(tile(r.gene_start, r.gene_end, pos, window, step, flank)):
            we = ws + window - 1
            sel = pos[(pos >= ws) & (pos <= we)]
            if not len(sel):
                continue
            seq, ex = window_sequence(fa, lengths, c, ws, window, r.strand, exons)
            a = reader.window(seq, ex)
            idx = (we - sel) if r.strand == "-" else (sel - ws)
            rows.append(pd.DataFrame({"gene": g, "pos": sel, "window": wi, "has_exon": bool(ex.any()), "attn_win": a[idx].astype(np.float32)}))
            nw += 1
        if rows and ((i + 1) % 50 == 0 or i + 1 == len(todo)):
            pd.concat(rows).to_parquet(f"{out}/chunks/{k:05d}.parquet", index=False)
            k, rows = k + 1, []
            log(f"attention: {i + 1}/{len(todo)} genes, {nw} windows, {time.time() - t0:.0f} s")
    fs = sorted(glob.glob(f"{out}/chunks/*.parquet"))
    W = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame(columns=["gene", "pos", "has_exon", "attn_win"])
    S = aggregate(W) if len(W) else pd.DataFrame(columns=["gene", "pos", "attention", "attention_source"])
    S.to_parquet(f"{out}/attention.parquet", index=False)
    log(f"attention: {S.gene.nunique()} genes, {len(S):,} variants -> {out}/attention.parquet")
