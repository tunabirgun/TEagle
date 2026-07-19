"""Deterministic sequence builders + fixture loaders for the test suite."""
import os, random

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
_B = "ACGT"
_STOP = {"TAA", "TAG", "TGA"}


def fixture_seq(acc):
    from teagle_core import sequtil
    fa = open(os.path.join(FIX, acc + ".fasta"), encoding="utf-8").read()
    return sequtil.parse_fasta(fa)[0][1]


def make_ltr_element(seed=1, ltr_len=160, divergence=0.02, flank=250):
    """A synthetic LTR retrotransposon: flank + TSD + LTR + internal(ORF) + LTR' + TSD + flank."""
    r = random.Random(seed)
    rb = lambda n: "".join(r.choice(_B) for _ in range(n))

    def rc():
        while True:
            c = rb(3)
            if c not in _STOP:
                return c

    def orf(ncodons):
        return "ATG" + "".join(rc() for _ in range(ncodons)) + "TAA"

    def mut(s, rate):
        a = list(s)
        for i in range(len(a)):
            if r.random() < rate:
                a[i] = r.choice(_B)
        return "".join(a)

    ltr = rb(ltr_len)
    tsd = rb(5)
    internal = rb(400) + orf(210) + rb(300)
    return rb(flank) + tsd + ltr + internal + mut(ltr, divergence) + tsd + rb(flank)


def make_tir_element(seed=7, tir_len=24, flank=200):
    """A synthetic DNA transposon with terminal inverted repeats."""
    from teagle_core.sequtil import reverse_complement
    r = random.Random(seed)
    rb = lambda n: "".join(r.choice(_B) for _ in range(n))
    tir = rb(tir_len)
    internal = rb(800)
    return rb(flank) + tir + internal + reverse_complement(tir) + rb(flank)
