"""Deterministic demo LTR element (seeded), ported from the web UI so 'load sample' is reproducible
and fires the structural/domain detectors. Same mulberry32 seed as app.js -> same sequence."""
from __future__ import annotations

_B = "ACGT"
_STOP = {"TAA", "TAG", "TGA"}


def _imul(x, y):                      # Math.imul: 32-bit integer multiply, unsigned result
    return (x * y) & 0xFFFFFFFF


def _mulberry32(a):
    state = {"a": a & 0xFFFFFFFF}
    def rng():
        state["a"] = (state["a"] + 0x6D2B79F5) & 0xFFFFFFFF
        a2 = state["a"]
        t = _imul(a2 ^ (a2 >> 15), 1 | a2)                     # imul(a^a>>>15, 1|a)
        t = ((t + _imul(t ^ (t >> 7), 61 | t)) & 0xFFFFFFFF) ^ t   # (t + imul(t^t>>>7, 61|t)) ^ t
        return ((t ^ (t >> 14)) & 0xFFFFFFFF) / 4294967296
    return rng


def make_sample() -> str:
    r = _mulberry32(20260717)

    def rb(n):
        return "".join(_B[int(r() * 4)] for _ in range(n))

    def rc():
        while True:
            c = _B[int(r() * 4)] + _B[int(r() * 4)] + _B[int(r() * 4)]
            if c not in _STOP:
                return c

    def orf(ncodons):
        return "ATG" + "".join(rc() for _ in range(ncodons)) + "TAA"

    def mut(s, rate):
        a = list(s)
        for i in range(len(a)):
            if r() < rate:
                a[i] = _B[int(r() * 4)]
        return "".join(a)

    tsd = rb(5)
    ltr = rb(160)
    internal = rb(300) + orf(210) + rb(280)
    seq = rb(250) + tsd + ltr + internal + mut(ltr, 0.02) + tsd + rb(250)
    wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
    return ">sample_LTR_element  demo construct (illustrative)\n" + wrapped
