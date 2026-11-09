"""Biallelic SNV dosages from a bgzipped, tabix-indexed VCF.

Dosage = number of ALT alleles (0, 1, 2); missing calls are NaN. A haploid call counts as
homozygous. Records that are not single-base biallelic SNVs (REF and ALT in A/C/G/T) are skipped;
for duplicated positions the first record is kept.
"""
import numpy as np
import pysam

_ACGT = frozenset("ACGT")
_ONE, _DOT = ord("1"), ord(".")


def resolve_contig(name, contigs):
    """Match chromosome names that differ only by a chr / Chr prefix."""
    name = str(name)
    if name in contigs:
        return name
    base = name[3:] if name.lower().startswith("chr") else name
    for c in contigs:
        cb = c[3:] if c.lower().startswith("chr") else c
        if cb == base:
            return c
    return None


def decode_gt(fields, gt_index=0):
    """Dosages from the sample columns of one VCF record."""
    if gt_index:
        fields = [s.split(":")[gt_index] for s in fields]
    g = np.array(fields, dtype="S3").view(np.uint8).reshape(-1, 3)
    diploid = (g[:, 1] == ord("/")) | (g[:, 1] == ord("|"))
    d = np.where(diploid, (g[:, 0] == _ONE).astype(np.float32) + (g[:, 2] == _ONE), 2.0 * (g[:, 0] == _ONE))
    d = d.astype(np.float32)
    d[(g[:, 0] == _DOT) | (diploid & (g[:, 2] == _DOT))] = np.nan
    return d


class VCF:
    def __init__(self, path):
        self.path = path
        self.tbx = pysam.TabixFile(path)
        self.samples = [str(s) for s in self.tbx.header[-1].rstrip("\n").split("\t")[9:]]
        self.contigs = list(self.tbx.contigs)

    def records(self, chrom, start=None, end=None, cols=None):
        """Yield (pos, ref, alt, dosage) for SNVs in [start, end] (1-based, inclusive); cols selects samples."""
        c = resolve_contig(chrom, self.contigs)
        if c is None:
            return
        it = self.tbx.fetch(c) if start is None else self.tbx.fetch(c, max(0, start - 1), end)
        last = -1
        for line in it:
            f = line.split("\t", 9)
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1 or ref not in _ACGT or alt not in _ACGT:
                continue
            pos = int(f[1])
            if pos == last:
                continue
            last = pos
            fmt = f[8].split(":")
            s = f[9].rstrip("\n").split("\t")
            if cols is not None:
                s = [s[i] for i in cols]
            yield pos, ref, alt, decode_gt(s, fmt.index("GT"))

    def region(self, chrom, start, end, cols=None):
        """Positions (int64) and dosages (samples x SNVs, float32) in a region."""
        P, X = [], []
        for pos, _, _, d in self.records(chrom, start, end, cols):
            P.append(pos)
            X.append(d)
        n = len(self.samples) if cols is None else len(cols)
        return np.array(P, np.int64), (np.stack(X, 1) if X else np.zeros((n, 0), np.float32))

    def sample_index(self, samples):
        ix = {s: i for i, s in enumerate(self.samples)}
        return [ix[s] for s in samples]
