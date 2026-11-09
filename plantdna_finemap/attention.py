"""Per-variant attention scores from PlantDNA-FM (GPU recommended).

For each gene, the 8-kb window is read in the direction of transcription (minus-strand genes as
reverse complement; positions beyond the chromosome ends as N) and run through the model once.
In every layer, attention from the exon positions of the representative transcript (all positions
if the window has no exon) is averaged over heads and query positions. The score of a variant is
this attention to its position, averaged over all layers (36 layers x 12 heads = 432 attention
matrices for PlantDNA-FM). Output: parquet with columns gene, pos, attention.
"""
import glob
import os
import re
import time
import numpy as np
import pandas as pd
import pysam
from .genotypes import VCF, resolve_contig
from .windows import parse_exons, read_windows

_RC = str.maketrans("ACGTN", "TGCAN")
_NON = re.compile("[^ACGT]")


class AttentionReader:
    def __init__(self, model, device="cuda"):
        import torch
        from transformers import AttentionInterface, AutoModel, AutoTokenizer
        from transformers.utils import logging as hf_logging
        hf_logging.set_verbosity_error()          # the unused masked-LM head is reported as unexpected keys
        from transformers.integrations.sdpa_attention import sdpa_attention_forward
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        self.torch, self.device, self.cur = torch, device, {}
        cur = self.cur

        def exon_rows(module, query, key, value, attention_mask, scaling=1.0, dropout=0.0, **kw):
            # query / key already carry the scaling and rotary embedding: (batch, heads, length, dim)
            qi = cur["q"]
            A = torch.softmax((query[0][:, qi].float() @ key[0].float().transpose(-1, -2)) * scaling, dim=-1)
            cur["layers"][module.layer_idx] = A.mean((0, 1)).cpu().numpy()
            return sdpa_attention_forward(module, query, key, value, attention_mask, scaling=scaling, dropout=dropout, **kw)

        AttentionInterface.register("pdfm_exon_rows", exon_rows)
        tok = AutoTokenizer.from_pretrained(model)
        v = tok.get_vocab()
        self.tok = {b: v[b] for b in "ACGTN"}
        self.cls, self.eos = tok.cls_token_id, tok.eos_token_id
        self.net = AutoModel.from_pretrained(model, dtype=torch.float32, attn_implementation="pdfm_exon_rows", add_pooling_layer=False).to(device).eval()
        for i, layer in enumerate(self.net.encoder.layer):
            layer.attention.self.layer_idx = i
        self.n_layers = self.net.config.num_hidden_layers

    def window(self, seq, exon_mask):
        """Layer-mean exon-query attention to every position of seq (transcription orientation)."""
        torch = self.torch
        q = np.where(exon_mask)[0] if exon_mask.any() else np.arange(len(seq))
        self.cur.update(q=torch.tensor(q + 1, device=self.device), layers=[None] * self.n_layers)
        ids = torch.tensor([[self.cls] + [self.tok[b] for b in seq] + [self.eos]], device=self.device)
        with torch.no_grad():
            self.net(input_ids=ids)
        return np.mean(self.cur["layers"], 0)[1:-1]


def window_sequence(fa, lengths, r):
    """Window sequence and exon mask in the direction of transcription."""
    c = resolve_contig(r.chrom, lengths)
    ws, we = int(r.win_start), int(r.win_end)
    L = lengths[c]
    s = _NON.sub("N", fa.fetch(c, max(ws, 1) - 1, min(we, L)).upper())
    s = "N" * max(0, 1 - ws) + s + "N" * max(0, we - L)
    ex = np.zeros(len(s), bool)
    for a, b in parse_exons(r.exons):
        lo, hi = max(a, ws) - ws, min(b, we) - ws
        if hi >= lo:
            ex[lo:hi + 1] = True
    if r.strand == "-":
        return s.translate(_RC)[::-1], ex[::-1]
    return s, ex


def run(model, windows, fasta, vcf, out, genes=None, device="cuda", log=print):
    W = read_windows(windows)
    if genes:
        W = W[W.gene.isin(set(genes))]
    os.makedirs(f"{out}/chunks", exist_ok=True)
    done = set()
    for f in glob.glob(f"{out}/chunks/*.parquet"):
        done |= set(pd.read_parquet(f, columns=["gene"]).gene)
    W = W[~W.gene.isin(done)]
    fa = pysam.FastaFile(fasta)
    lengths = {c: fa.get_reference_length(c) for c in fa.references}
    v = VCF(vcf)
    log(f"attention: {len(W)} genes to score ({len(done)} done)")
    reader = AttentionReader(model, device) if len(W) else None
    rows, k, t0 = [], len(glob.glob(f"{out}/chunks/*.parquet")), time.time()
    for i, r in enumerate(W.itertuples(index=False)):
        if resolve_contig(r.chrom, lengths) is None:
            continue
        pos = np.array([p for p, *_ in v.records(r.chrom, int(r.win_start), int(r.win_end), cols=[])], np.int64)
        if len(pos):
            seq, ex = window_sequence(fa, lengths, r)
            a = reader.window(seq, ex)
            idx = (int(r.win_end) - pos) if r.strand == "-" else (pos - int(r.win_start))
            rows.append(pd.DataFrame({"gene": r.gene, "pos": pos, "attention": a[idx].astype(np.float32)}))
        if rows and ((i + 1) % 100 == 0 or i + 1 == len(W)):
            pd.concat(rows).to_parquet(f"{out}/chunks/{k:05d}.parquet", index=False)
            k, rows = k + 1, []
            log(f"attention: {i + 1}/{len(W)} genes, {time.time() - t0:.0f} s")
    fs = sorted(glob.glob(f"{out}/chunks/*.parquet"))
    A = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True) if fs else pd.DataFrame(columns=["gene", "pos", "attention"])
    A.drop_duplicates(["gene", "pos"]).to_parquet(f"{out}/attention.parquet", index=False)
    log(f"attention: {A.gene.nunique()} genes, {len(A):,} variants -> {out}/attention.parquet")
