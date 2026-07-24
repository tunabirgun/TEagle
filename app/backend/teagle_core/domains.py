"""Real TE protein-domain detection (Layer C) via native HMMER (pyhmmer) against a
bundled CC0 Pfam TE-domain profile set. Translates ORFs, runs hmmsearch, maps hits
back to nucleotide coordinates. No WSL, no external binaries, fully offline."""
from __future__ import annotations
import os
from .sequtil import reverse_complement, translate, find_orfs
from . import appdirs

try:                                                # a broken/missing pyhmmer must not crash the engine —
    import pyhmmer                                  # domain detection degrades to "unavailable", everything else runs
    PYHMMER_ERROR = None
except Exception as _e:
    pyhmmer = None
    PYHMMER_ERROR = f"{type(_e).__name__}: {_e}"

HMM_PATH = appdirs.resource("data", "te_domains.hmm") or \
    os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "te_domains.hmm"))

PYHMMER_VERSION = getattr(pyhmmer, "__version__", "unavailable") if pyhmmer is not None else "unavailable"


def _hmm_sha256():                                  # pin the bundled profile set into provenance
    import hashlib
    try:
        with open(HMM_PATH, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except Exception:
        return None


HMM_SHA256 = _hmm_sha256()

# hmm profile name -> (short domain code, human label, functional class, Pfam accession)
DOMAIN_INFO = {
    "RVT_1": ("RT", "reverse transcriptase", "retro", "PF00078"),
    "RVT_2": ("RT", "reverse transcriptase", "retro", "PF07727"),
    "RVT_3": ("RT", "reverse transcriptase", "retro", "PF13456"),
    "rve": ("INT", "integrase", "retro", "PF00665"),
    "RNase_H": ("RNaseH", "RNase H", "retro", "PF00075"),
    "RVP": ("PR", "aspartic protease", "retro", "PF00077"),
    "PEG10_N-capsid": ("GAG", "gag capsid (retrotransposon/PEG10-type)", "retro", "PF03732"),
    # retroviral / ERV gag (matrix, capsid, nucleocapsid) and env (glycoprotein, TM, surface) — the models that
    # annotate the HERV-K(HML-2) Gag/Env polyproteins (UniProt P62684, HERV-K env entries). All Pfam-A (CC0).
    "Gag_p24": ("GAG", "gag capsid (CA)", "retro", "PF00607"),
    "Gag_p24_C": ("GAG", "gag capsid, C-terminal (CA)", "retro", "PF19317"),
    "Gag_p10": ("GAG", "gag matrix (MA)", "retro", "PF02337"),
    "zf-CCHC_5": ("GAG", "gag nucleocapsid zinc-finger (NC)", "retro", "PF14787"),
    "HERV-K_env_2": ("ENV", "envelope glycoprotein", "retro", "PF13804"),
    "GP41": ("ENV", "envelope, transmembrane (TM)", "retro", "PF00517"),
    "TLV_coat": ("ENV", "envelope, surface (SU)", "retro", "PF00429"),
    "Chromo": ("CHR", "chromodomain", "retro", "PF00385"),
    "HTH_Tnp_Tc3_2": ("TPase", "Tc1/mariner transposase", "dna:Tc1-Mariner", "PF01498"),
    "DDE_1": ("TPase", "DDE transposase", "dna:DDE", "PF03184"),
    "DDE_3": ("TPase", "DDE transposase", "dna:Tc1-Mariner", "PF13358"),
    "Transposase_1": ("TPase", "mariner-type transposase", "dna:Tc1-Mariner", "PF01359"),
    "Dimer_Tnp_hAT": ("TPase", "hAT transposase", "dna:hAT", "PF05699"),
    "hAT-like_RNase-H": ("TPase", "hAT-like transposase", "dna:hAT", "PF14372"),
}

_ABC = None
_HMMS = None


def _abc():                                         # lazy so a missing pyhmmer never fails at import
    global _ABC
    if _ABC is None:
        _ABC = pyhmmer.easel.Alphabet.amino()
    return _ABC


def _hmms():
    global _HMMS
    if _HMMS is None:
        with pyhmmer.plan7.HMMFile(HMM_PATH) as f:
            _HMMS = list(f)
    return _HMMS


def scan_domains(seq: str, max_orfs: int = 12, evalue: float = 1e-3):
    """Detect TE protein domains in the sequence's ORFs. Returns hits ordered along
    the element by genomic position, each with nucleotide coordinates."""
    if pyhmmer is None:                             # domain detection unavailable in this environment
        return []
    orfs = find_orfs(seq)[:max_orfs]
    seqs, meta = [], {}
    for n, o in enumerate(orfs):
        sub = seq[o["start"]:o["end"]] if o["strand"] == "+" else reverse_complement(seq[o["start"]:o["end"]])
        prot = translate(sub).rstrip("*")
        if len(prot) >= 40:
            seqs.append(pyhmmer.easel.TextSequence(name=f"orf{n}".encode(), sequence=prot).digitize(_abc()))
            meta[n] = o
    if not seqs:
        return []
    block = pyhmmer.easel.DigitalSequenceBlock(_abc(), seqs)
    hits = []
    for top in pyhmmer.hmmsearch(_hmms(), block, E=evalue):
        hmm_name = str(top.query.name)
        code, label, dclass, pfam = DOMAIN_INFO.get(hmm_name, (hmm_name, hmm_name, "other", ""))
        for h in top:
            n = int(str(h.name)[3:])
            o = meta[n]
            for d in h.domains:
                if d.i_evalue >= evalue:
                    continue
                aa0, aa1 = d.env_from - 1, d.env_to           # 1-based -> 0-based half-open
                if o["strand"] == "+":
                    nt = [o["start"] + aa0 * 3, o["start"] + aa1 * 3]
                    coding = seq[nt[0]:nt[1]]
                else:
                    nt = [o["end"] - aa1 * 3, o["end"] - aa0 * 3]
                    coding = reverse_complement(seq[nt[0]:nt[1]])
                iev = float(d.i_evalue)
                hits.append({
                    "domain": code, "label": label, "class": dclass, "hmm": hmm_name, "pfam": pfam,
                    "score": round(d.score, 1), "evalue": iev,
                    # per-domain call confidence = the HMMER i-Evalue (Eddy 2011); high when strongly significant
                    "confidence": "high" if iev <= 1e-10 else "moderate",
                    "orf": n, "strand": o["strand"], "aa": [d.env_from, d.env_to], "nt": nt,
                    "dna": coding, "protein": translate(coding).rstrip("*"),
                })
    return _dedup_domains(hits)


def _dedup_domains(hits):
    """Keep the best-scoring hit per (domain-code, STRAND, overlapping-nt-region), then order along the
    element by genomic position. The strand check keeps a genuine minus-strand hit that merely overlaps a
    higher-scoring plus-strand hit of the same code — they are different biological features."""
    kept = []
    for hh in sorted(hits, key=lambda x: -x["score"]):
        if any(o["domain"] == hh["domain"] and o["strand"] == hh["strand"]
               and not (hh["nt"][1] <= o["nt"][0] or hh["nt"][0] >= o["nt"][1]) for o in kept):
            continue                                          # overlapping same-domain same-strand, lower score -> drop
        kept.append(hh)
    return sorted(kept, key=lambda x: x["nt"][0])
