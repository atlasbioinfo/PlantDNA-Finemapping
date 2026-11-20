"""Genotype dosages (samples x variants, missing = NaN) from a VCF, a PLINK 1 fileset or a numpy .npz.

VCF: bgzipped and tabix-indexed; dosage = number of ALT alleles (0, 1, 2); a haploid call counts as
homozygous; records that are not single-base biallelic SNVs (REF and ALT in A/C/G/T) are skipped.
PLINK: dosage = number of A1 alleles. The sign of the dosage does not change SuSiE PIPs or credible sets.
For duplicated positions the first record is kept. Chromosome names may differ by a chr / Chr prefix.
"""
import os
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

    def chunks(self, cols=None, size=20000):
        """Yield (chrom, positions, dosages) for all SNVs of the file, `size` variants at a time."""
        n = len(self.samples) if cols is None else len(cols)
        for c in self.contigs:
            P, X = [], []
            for pos, _, _, d in self.records(c, cols=cols):
                P.append(pos)
                X.append(d)
                if len(P) == size:
                    yield c, np.array(P, np.int64), np.stack(X, 1)
                    P, X = [], []
            if P:
                yield c, np.array(P, np.int64), np.stack(X, 1)

    def sample_index(self, samples):
        ix = {s: i for i, s in enumerate(self.samples)}
        return [ix[s] for s in samples]

    def fetch(self, chrom, positions, samples):
        """Dosages (samples x variants) at the requested positions that are biallelic SNVs in the VCF."""
        want = set(int(p) for p in positions)
        cols = self.sample_index(samples)
        P, X = [], []
        for pos, _, _, d in self.records(chrom, min(want), max(want), cols):
            if pos in want:
                P.append(pos)
                X.append(d)
        return np.array(P, np.int64), (np.stack(X, 1) if X else np.zeros((len(cols), 0), np.float32))


class _Indexed:
    """Genotypes held as a variant table (chrom, pos) and a sample list; subclasses implement _dosage(idx)."""

    def _index(self, chrom, pos):
        self.contigs = list(dict.fromkeys(str(c) for c in chrom))
        self._ix = {}
        for i, (c, p) in enumerate(zip(chrom, pos)):
            self._ix.setdefault((str(c), int(p)), i)      # first record of a duplicated position
        first = np.array(sorted(self._ix.values()), int)
        cc, pp = np.array([str(c) for c in chrom])[first], np.asarray(pos, np.int64)[first]
        self._by = {}
        for c in self.contigs:                            # per chromosome: positions in ascending order, record index
            k = np.where(cc == c)[0]
            o = np.argsort(pp[k], kind="stable")
            self._by[c] = (pp[k][o], first[k][o])

    def region(self, chrom, start, end, cols=None):
        """Positions (int64) and dosages (samples x variants, float32) in [start, end] (1-based, inclusive)."""
        c = resolve_contig(chrom, self.contigs)
        n = len(self.samples) if cols is None else len(cols)
        if c is None:
            return np.zeros(0, np.int64), np.zeros((n, 0), np.float32)
        P, I = self._by[c]
        a, b = np.searchsorted(P, start, "left"), np.searchsorted(P, end, "right")
        if a == b:
            return np.zeros(0, np.int64), np.zeros((n, 0), np.float32)
        X = self._dosage(I[a:b])
        return P[a:b], (X if cols is None else X[np.asarray(cols, int)])

    def chunks(self, cols=None, size=20000):
        """Yield (chrom, positions, dosages) for all variants of the file, `size` variants at a time."""
        for c in self.contigs:
            P, I = self._by[c]
            for a in range(0, len(P), size):
                X = self._dosage(I[a:a + size])
                yield c, P[a:a + size], (X if cols is None else X[np.asarray(cols, int)])

    def sample_index(self, samples):
        ix = {s: i for i, s in enumerate(self.samples)}
        return [ix[s] for s in samples]

    def fetch(self, chrom, positions, samples):
        c = resolve_contig(chrom, self.contigs)
        hits = [(int(p), self._ix[(c, int(p))]) for p in positions if (c, int(p)) in self._ix] if c is not None else []
        cols = np.array(self.sample_index(samples), int)
        if not hits:
            return np.zeros(0, np.int64), np.zeros((len(cols), 0), np.float32)
        return np.array([p for p, _ in hits], np.int64), self._dosage(np.array([i for _, i in hits]))[cols]


class Plink(_Indexed):
    """PLINK 1 binary fileset (prefix.bed/.bim/.fam, SNP-major). Dosage = number of A1 alleles; samples = IID."""

    def __init__(self, prefix):
        prefix = prefix[:-4] if prefix.endswith((".bed", ".bim", ".fam")) else prefix
        bim = np.loadtxt(f"{prefix}.bim", dtype=str, ndmin=2)
        self.samples = [l.split()[1] for l in open(f"{prefix}.fam") if l.strip()]
        self._index(bim[:, 0], bim[:, 3].astype(np.int64))
        self._nb = (len(self.samples) + 3) // 4
        self._bed = np.memmap(f"{prefix}.bed", np.uint8, "r")
        if self._bed[:3].tolist() != [0x6C, 0x1B, 0x01]:
            raise ValueError(f"{prefix}.bed is not a SNP-major PLINK 1 .bed file")

    def _dosage(self, idx):
        b = np.stack([self._bed[3 + i * self._nb: 3 + (i + 1) * self._nb] for i in idx])
        codes = np.unpackbits(b[:, :, None], axis=2, bitorder="little").reshape(len(idx), -1, 2)
        codes = (codes[:, :, 0] + 2 * codes[:, :, 1])[:, :len(self.samples)]       # 0 hom A1, 1 missing, 2 het, 3 hom A2
        return np.array([2, np.nan, 1, 0], np.float32)[codes].T


class NumpyGenotypes(_Indexed):
    """.npz with arrays samples (n), chrom (p), pos (p) and dosage (n x p; NaN or negative = missing)."""

    def __init__(self, path):
        z = np.load(path, allow_pickle=False)
        self.samples = [str(s) for s in z["samples"]]
        self._index(z["chrom"], z["pos"])
        self._G = z["dosage"]

    def _dosage(self, idx):
        X = self._G[:, idx].astype(np.float32)
        X[X < 0] = np.nan
        return X.T


def open_genotypes(path):
    """VCF (.vcf.gz with .tbi), PLINK 1 fileset (prefix or .bed) or .npz."""
    if path.endswith(".npz"):
        return NumpyGenotypes(path)
    if path.endswith(".bed") or os.path.exists(f"{path}.bed"):
        return Plink(path)
    return VCF(path)
