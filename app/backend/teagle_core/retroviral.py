"""Retroviral transcript architecture (Layer E) — the correct coding-organisation model for an
ENDOGENOUS RETROVIRUS, in place of a host exon-intron gene model.

A retrovirus does not express its genes the way a host gene does. gag-pro-pol is one continuous
coding unit translated from the full-length, UNSPLICED genomic RNA as a fused polyprotein (via
ribosomal frameshifting OR stop-codon readthrough — the exact mechanism is genus-dependent);
env is the only subgenomic mRNA spliced DIRECTLY from that full-length RNA — a splice donor in the
5' leader joins to an acceptor 5' of env, removing the whole gag-pro-pol span as ONE large intron.
Some lineages (notably HML-2) additionally sub-splice the env mRNA to accessory transcripts
(HERV-K rec / np9); that second splice is noted, not drawn.

This module reports the PREDICTED transcript architecture anchored on the LTR + protein-domain
positions (which are reliable), with the splice JUNCTIONS labelled approximate — the exact junction
bases cannot be recovered from proviral DNA alone and need a supplied transcript / RNA-seq. It never
fabricates a single-base donor/acceptor from motif guessing, and it is built only when the domain
geometry is clean and internally consistent (both strands handled; degenerate layouts return None).

Refs: Löwer et al. 1995 J Virol 69:141 (spliced env/rec transcripts); Schmitt et al. 2015
Mobile DNA 6:4 (HML-2 rec/np9, type-1/type-2)."""
from __future__ import annotations

_BODY_CODES = ("GAG", "PR", "RT", "RNaseH", "INT")   # gag-pro-pol polyprotein modules (the env intron)


def _strand(domains) -> str:
    # dominant strand of the gag/pol/env modules; determines whether the element reads 5'->3' with
    # ascending genomic coordinates (+) or descending (-, a reverse-complement-oriented provirus).
    votes = [d.get("strand", "+") for d in domains if d.get("domain") in _BODY_CODES + ("ENV",)]
    return "-" if votes.count("-") > votes.count("+") else "+"


def transcript_architecture(structural_ev, domains, classification, seq_len: int):
    """Build the ERV subgenomic-env transcript model in genomic coordinates. Returns None unless the
    element is an ERV (env + paired LTRs + RT) whose gag-pro-pol/env domains lay out consistently on
    one strand — so a non-ERV, or a garbled/degenerate domain layout, never yields a spurious model."""
    if not classification or not classification.get("is_erv"):
        return None
    ltr = next((e for e in structural_ev if e.get("type", "").startswith("LTR")), None)
    env = [d for d in domains if d.get("domain") == "ENV"]
    if not ltr or not env:
        return None
    lo_ltr, hi_ltr = ltr["five_prime"][1], ltr["three_prime"][0]   # inner edges of the 5'/3' genomic LTRs
    strand = _strand(domains)
    env_lo = min(d["nt"][0] for d in env)
    env_hi = max(d["nt"][1] for d in env)

    if strand == "+":
        # 5'->3' ascends: [5'LTR] leader gag-pro-pol env [3'LTR]. Body = gag/pol modules BEFORE env.
        body = [d for d in domains if d.get("domain") in _BODY_CODES and d["nt"][1] <= env_lo]
        if not body:
            return None
        body_start = min(d["nt"][0] for d in body)
        body_end = max(d["nt"][1] for d in body)
        leader_exon = [lo_ltr, body_start]
        intron = [body_start, env_lo]
        env_exon = [env_lo, min(env_hi, hi_ltr)]
    else:
        # reverse-complement-oriented: 5'->3' descends in genomic coords: [3'LTR] env gag-pro-pol leader [5'LTR].
        # Body = gag/pol modules ABOVE env; the leader sits between the top of gag and the (genomic-high) 5' LTR.
        body = [d for d in domains if d.get("domain") in _BODY_CODES and d["nt"][0] >= env_hi]
        if not body:
            return None
        body_end = max(d["nt"][1] for d in body)
        env_exon = [max(env_lo, lo_ltr), env_hi]
        intron = [env_hi, body_end]
        leader_exon = [body_end, hi_ltr]

    # strict geometry: every span positive-length, ordered leader/intron/env with no overlap, inside the LTRs.
    spans = [leader_exon, intron, env_exon]
    if any(a >= b for a, b in spans):
        return None
    ordered = (leader_exon, intron, env_exon) if strand == "+" else (env_exon, intron, leader_exon)
    xs = [c for span in ordered for c in span]
    if xs != sorted(xs) or xs[0] < lo_ltr or xs[-1] > hi_ltr:
        return None

    model = {
        "kind": "ERV subgenomic env mRNA (predicted)",
        "strand": strand,
        "leader_exon": leader_exon,
        "env_exon": env_exon,
        "intron": intron,                                 # gag-pro-pol span = the single env intron + the polyprotein
        "polyprotein": intron[:],
        "approximate": True,
        "note": ("env is expressed from a subgenomic mRNA spliced directly from the full-length RNA: the "
                 "gag-pro-pol region is removed as a single intron and is itself a fused polyprotein "
                 "(translated by ribosomal frameshift or stop-codon readthrough, the mechanism is "
                 "genus-dependent), not a set of host exons. The junctions are anchored on the LTR and domain "
                 "positions and are approximate — the true splice donor lies within the drawn leader (upstream "
                 "of gag) and the acceptor upstream of the env domain, so the leader exon is an upper bound; "
                 "exact splice bases need a supplied env transcript or RNA-seq."),
    }
    # HML-2-specific accessory sub-splicing (rec/np9) is only asserted CONDITIONALLY — TEagle has no
    # HML-2-vs-other-Gypsy signal, so it must not attach HERV-K biology to a non-HML-2 (e.g. insect) ERV.
    model["subsplice_note"] = ("Some ERV lineages additionally sub-splice the env mRNA to accessory transcripts "
                               "(e.g. HERV-K/HML-2 rec and np9); such a second splice is noted, not drawn, as it "
                               "needs a transcript to place.")
    return model
