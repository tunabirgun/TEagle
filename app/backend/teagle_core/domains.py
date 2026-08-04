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
    # --- LINE / non-LTR coding modules. Until these, the panel modelled a LINE by its RT alone, so a
    # full-length L1 and a dead 5'-truncated fragment scored identically. Validated on M80343 (L1.2):
    # ORF1p E=1.6e-68 / 2.1e-35 / 1.5e-51, endonuclease E=1.0e-23, with no hit on copia, gypsy, Tc1 or Ac.
    "Transposase_22": ("ORF1", "LINE ORF1p, RNA-binding domain", "retro:LINE", "PF02994"),
    "Tnp_22_trimer": ("ORF1", "LINE ORF1p, trimerisation domain", "retro:LINE", "PF17489"),
    "Tnp_22_dsRBD": ("ORF1", "LINE ORF1p, dsRBD-like domain", "retro:LINE", "PF17490"),
    # BROAD family: Pfam's Exo_endo_phos also covers host DNase I and AP endonucleases, so an EN hit is
    # only element evidence when it sits N-terminal to an RT in the same ORF (classify enforces this).
    "Exo_endo_phos": ("EN", "apurinic-like endonuclease (ORF2p)", "retro:LINE?", "PF03372"),
    # BROAD family: tyrosine recombinases are widespread in hosts and phage. Reported as a YR signal that
    # makes a DIRS-group element plausible; never asserted as one on its own. Validated on M11340 (DIRS-1).
    "Phage_integrase": ("YR", "tyrosine recombinase (DIRS/Crypton-type)", "retro:YR?", "PF00589"),
    "Helitron_like_N": ("HEL", "Helitron helicase-like domain", "dna:Helitron", "PF14214"),
    "Transposase_24": ("TPase", "CACTA/En-Spm transposase", "dna:CACTA", "PF03004"),
    "MULE": ("TPase", "MULE/Mutator transposase", "dna:MULE", "PF10551"),
    "DDE_Tnp_1_7": ("TPase", "IS4-like DDE transposase (piggyBac-like)", "dna:IS4", "PF13843"),
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


# How many ORFs the profile search actually reads. find_orfs sorts descending by length, so this keeps the
# longest — where a TE's coding modules live — but on a sequence with more ORFs than this the remainder are
# never searched, and an unsearched ORF must not contribute to a "not detected" claim. Named so callers can
# compare it against the real ORF count and say so.
MAX_ORFS_SCANNED = 12


def scan_domains(seq: str, max_orfs: int = MAX_ORFS_SCANNED, evalue: float = 1e-3):
    """Detect TE protein domains in the sequence's ORFs. Returns hits ordered along
    the element by genomic position, each with nucleotide coordinates.

    Only the `max_orfs` longest ORFs are searched — see MAX_ORFS_SCANNED."""
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
