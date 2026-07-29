"""TE classification from structural evidence + protein-domain architecture (Layer C).
Superfamily calls follow domain content, structural context, and — for LTR elements —
the diagnostic integrase-vs-RT order (Copia: INT before RT; Gypsy: INT after RT).
Transparent, evidence-derived; never a family/name call beyond what the evidence supports."""
from __future__ import annotations

from . import structural as structural_mod        # tsd_congruence: the expected-TSD-length table


def _pos(d):
    # translation-order position along the pol polyprotein (strand-aware): smaller = earlier.
    # domains.py maps +strand aa->nt ascending and -strand aa->nt descending in nt[0]. .get for robustness.
    nt = d.get("nt") or [0, 0]
    return nt[0] if d.get("strand", "+") == "+" else -nt[0]


def _rep(domains, code):
    # representative hit for a domain code = highest-scoring occurrence (not ORF length-rank)
    ds = [d for d in domains if d["domain"] == code]
    return max(ds, key=lambda d: d.get("score", 0.0)) if ds else None


def _brackets(term, d):
    # do a terminal-repeat pair's two arms enclose this domain hit? an autonomous TIR transposon's
    # transposase lies BETWEEN the ends its own product binds. None when either side is unavailable.
    if not term or not d:
        return None
    s = (term.get("five_prime") or [0, 0])[0]
    e = (term.get("three_prime") or [0, 0])[1]
    nt = d.get("nt") or [0, 0]
    return s <= nt[0] and nt[1] <= e


# structural core of gag = capsid or matrix (incl. PEG10-type capsid); the CCHC nucleocapsid zinc-knuckle
# (zf-CCHC_5) is promiscuous (shared with host nucleic-acid-binding proteins) and is NOT core gag evidence on its own.
_GAG_CORE = {"Gag_p24", "Gag_p24_C", "Gag_p10", "PEG10_N-capsid"}

_CODE_ORDER = ["ORF1", "GAG", "PR", "EN", "RT", "RNaseH", "INT", "YR", "CHR", "ENV", "TPase", "HEL"]


def _ordered(cset):
    # detected domain codes in canonical retroviral order, unknown codes appended alphabetically
    return [c for c in _CODE_ORDER if c in cset] + sorted(c for c in cset if c not in _CODE_ORDER)


def classify(structural, domains):
    has_ltr = any(e["type"].startswith("LTR") for e in structural)
    tir_ev = next((e for e in structural if e["type"].startswith("TIR")), None)
    has_tir = tir_ev is not None
    has_tsd = any(e["type"].startswith("TSD") for e in structural)
    has_polya = any(e["type"].startswith("poly") for e in structural)
    codes = [d["domain"] for d in domains]
    cset = set(codes)
    rt_d, int_d, tp_d = _rep(domains, "RT"), _rep(domains, "INT"), _rep(domains, "TPase")
    # ORF2p is a single reading frame with the endonuclease N-terminal to the RT. Pfam's Exo_endo_phos
    # also matches host DNase I and AP endonucleases, so an EN hit only counts as element evidence when
    # it holds that arrangement — same strand, upstream of the RT in translation order.
    en_d = _rep(domains, "EN")
    en_ok = bool(en_d and rt_d and en_d.get("strand", "+") == rt_d.get("strand", "+")
                 and _pos(en_d) < _pos(rt_d))
    rt, intg = rt_d is not None, int_d is not None
    tpase = "TPase" in cset
    # the element's own ends: a TIR pair that actually encloses the transposase. The scan resolves the two
    # arms together, so they stand or fall as a pair; an inverted repeat elsewhere in the record is not
    # evidence that THIS coding module has its termini.
    tir_encloses_tpase = _brackets(tir_ev, tp_d)
    tir_ok = has_tir and tir_encloses_tpase is not False
    ev, superfamily, te_class, order = [], None, None, None
    order_resolvable = False
    tpase_conflict = False

    # DIRS-group elements carry a tyrosine recombinase in place of the DDE integrase, so they must be
    # tested BEFORE the generic RT branch — otherwise every one of them falls through to "LINE (non-LTR)"
    # on the strength of an absent integrase. Wicker 2007 order DIRS; PASTEC's YR-vs-DDE logic
    # (Hoede 2014 PLoS ONE 9:e91929) makes the same distinction.
    yr = "YR" in cset
    if rt and yr and not intg:
        klass = "Class I · retrotransposon"
        superfamily, te_class = "DIRS-group (tyrosine-recombinase retroelement)", "DIRS"
        ev.append("reverse transcriptase with a tyrosine recombinase and no DDE integrase → DIRS-group "
                  "architecture, not a LINE")
        if has_ltr:
            ev.append("terminal repeats present, consistent with the split direct / inverted repeats of DIRS")
        ev.append("tyrosine recombinases also occur in host and phage proteins — this call rests on the "
                  "combination with reverse transcriptase, never on the recombinase alone")
        order = "–".join(_ordered(cset))
    elif rt:
        klass = "Class I · retrotransposon"
        ev.append("reverse-transcriptase domain present")
        if intg and has_ltr:
            # decide Copia vs Gypsy by strand-aware TRANSLATION order of INT vs RT, not ORF length-rank;
            # the call is only unambiguous when both are on the same strand with non-overlapping spans
            i_nt, r_nt = (int_d.get("nt") or [0, 0]), (rt_d.get("nt") or [0, 0])
            order_resolvable = (int_d.get("strand", "+") == rt_d.get("strand", "+") and
                                (i_nt[1] <= r_nt[0] or r_nt[1] <= i_nt[0]))
            if _pos(int_d) < _pos(rt_d):
                superfamily = "Copia (Ty1)"
                ev.append("integrase N-terminal to RT + paired LTRs → Copia/Ty1 order")
            else:
                superfamily = "Gypsy (Ty3)"
                if "CHR" in cset:
                    superfamily += " · chromovirus"
                ev.append("integrase C-terminal to RT + paired LTRs → Gypsy/Ty3 order")
            if not order_resolvable:
                ev.append("integrase/RT order not cleanly resolvable (different strands or overlapping spans) — superfamily call is tentative")
            te_class = "LTR/" + superfamily.split(" ")[0]
            if "PR" in cset:                                  # cset holds emitted domain codes; domains.py maps RVP -> code "PR"
                ev.append("aspartic-protease domain present")
            if "RNaseH" in cset:
                ev.append("RNase H domain present")
        elif not has_ltr and not intg:
            superfamily, te_class = "LINE (non-LTR)", "LINE"
            ev.append("RT without a DDE integrase and without LTRs → non-LTR retrotransposon (LINE)")
            if has_polya:                                     # name the tail that was actually detected, not always poly-A
                _pa = any(e["type"].startswith("poly-A") for e in structural)
                ev.append(("3′ poly-A tail" if _pa else "5′ poly-T tract") +
                          " consistent with target-primed reverse transcription")
            if "ORF1" in cset:
                ev.append("ORF1p domain present — the LINE-specific coding module, not shared with LTR elements")
            if en_ok:
                ev.append("apurinic-like endonuclease N-terminal to the RT in the same reading frame → "
                          "the ORF2p EN–RT architecture of a LINE")
            elif "EN" in cset:
                ev.append("an endonuclease domain was detected but does not sit N-terminal to the RT in the "
                          "same frame; Pfam's endonuclease family also covers host enzymes, so it is not "
                          "credited as ORF2p here")
            # A DIRS-group element is now tested for directly (above) rather than merely hedged against.
            # Penelope-like elements remain outside the panel: their GIY-YIG endonuclease is not modelled.
            ev.append("a Penelope-like element would also show RT without a DDE integrase; its GIY-YIG "
                      "endonuclease is not in the tested panel, so it is not excluded here")
        elif has_ltr:
            superfamily, te_class = "LTR retrotransposon (superfamily undetermined)", "LTR/unclassified"
            ev.append("LTRs + RT present but integrase/order not resolved")
        else:
            superfamily, te_class = "retrotransposon (partial)", "retro/partial"
    elif tpase:
        klass = "Class II · DNA transposon"
        tp_hits = [d for d in domains if d["domain"] == "TPase"]
        best = max(tp_hits, key=lambda d: d.get("score", 0.0))       # decide by the strongest hit, not a fixed hAT-first precedence
        bcl = best.get("class", "")
        if "hAT" in bcl:
            superfamily = "hAT"
        elif "Tc1-Mariner" in bcl:
            superfamily = "Tc1/Mariner"
        elif "CACTA" in bcl:
            superfamily = "CACTA (En/Spm)"
        elif "MULE" in bcl:
            superfamily = "MULE (Mutator)"
        elif "IS4" in bcl:
            # PF13843 is Pfam's generic IS4-like DDE family; piggyBac sits inside it but so do others,
            # so the call names the family that was actually matched, not the best-known member.
            superfamily = "IS4-like DDE (piggyBac group)"
        elif "DDE" in bcl:
            superfamily = "DDE transposon"
        else:
            superfamily = "DNA transposon"
        # take the leading token, then the part before any "/" — mirrors the LTR branch. Splitting only on
        # "/" garbled a name whose slash sits inside a parenthetical ("CACTA (En/Spm)" -> "CACTA (En").
        _sf_token = superfamily.split(" ")[0].split("/")[0]
        # the generic "DNA transposon" fallback would collapse to a redundant "DNA/DNA"; name it honestly
        te_class = "DNA/" + (_sf_token if _sf_token != "DNA" else "unclassified")
        ev.append("transposase domain present → Class II DNA transposon")
        tpase_conflict = len({d["class"] for d in tp_hits}) > 1
        if tpase_conflict:
            ev.append("multiple transposase classes detected — superfamily assigned from the strongest-scoring hit (ambiguous)")
        if has_tir:
            ev.append("terminal inverted repeats consistent with a cut-and-paste transposon")
        # The ends are half of what makes a DNA transposon autonomous, so state what was actually recovered.
        # The "insertion site captured" claim is a completeness statement built on a short exact direct
        # repeat, which arises by chance at ~4^-L — so it is GATED on the TSD length not contradicting the
        # superfamily's expected one, and the incongruent case is reported instead of suppressed.
        if tir_ok and has_tsd:
            _tsd = next((e for e in structural if e["type"].startswith("TSD")), None)
            _cong = structural_mod.tsd_congruence((_tsd or {}).get("length"), superfamily)
            if _cong["verdict"] == "incongruent":
                ev.append(f"a target-site duplication flanks the inverted repeats, but it is "
                          f"{_cong['observed']} bp where {superfamily} elements duplicate "
                          f"{_cong['expected']} bp ({_cong['basis']}) — the repeat may be coincidental, or the "
                          f"boundary or superfamily call may be wrong, so the ends are not credited from it")
            else:
                ev.append("a target-site duplication flanks the inverted repeats — the insertion site itself is "
                          "captured, so both element termini are present in the record"
                          + (f" (length congruent with {superfamily}, {_cong['basis']})"
                             if _cong["verdict"] == "congruent" else ""))
        elif has_tir and tir_encloses_tpase is False:
            ev.append("the detected inverted-repeat pair does not enclose the transposase — the termini and the "
                      "coding module may not belong to the same element, so the ends are not credited")
        elif not has_tir:
            ev.append("no terminal inverted-repeat pair recovered — neither element end is confirmed, which is "
                      "consistent with a 5′- or 3′-truncated (or internally deleted) copy, with a superfamily "
                      "that does not carry TIRs, or with termini too diverged for the scan")
        order = "–".join(_ordered(cset))     # coding architecture only — the TIR ends are structural, not domains
    else:
        klass = "unclassified"
        if cset:                                          # coding domain(s) recovered but no RT and no transposase
            _cs = "–".join(_ordered(cset))                # e.g. an ERV relic that kept env/capsid but lost pol
            if has_ltr:
                superfamily, te_class = f"LTR retroelement fragment ({_cs}; RT/pol not detected)", "LTR/partial"
            else:
                superfamily, te_class = f"coding fragment ({_cs}; RT/transposase not detected)", "retro/partial"
            ev.append(f"coding domain(s) recovered ({_cs}) but neither reverse transcriptase nor transposase detected")
        elif has_ltr:
            superfamily, te_class = "LTR retrotransposon (no coding domains detected)", "LTR/structural-only"
            ev.append("paired LTRs but no coding domain recovered")
        elif has_tir:
            # An inverted terminal repeat WITHOUT a transposase is not Class II evidence. Verified on the
            # real copia 5' LTR (276 bp): fed alone it trips the terminal-inverted-repeat scan and used to
            # be filed as a DNA transposon — a Class I fragment assigned to Class II. Solo LTRs outnumber
            # full-length elements in most genomes, so this is the common case, not an edge case.
            superfamily, te_class = "terminal inverted repeat, class unassigned", "repeat/structural-only"
            ev.append("an inverted terminal repeat was detected but no transposase — this does not establish "
                      "a DNA transposon. A solo LTR, a fragment of a larger element, or an unrelated "
                      "inverted repeat all produce the same signal, so the class is left unassigned")
        else:
            superfamily, te_class = "no clear TE signature", "none"

    # a transposase domain co-occurring with RT is never inspected by the RT branch (the `elif tpase` above is
    # unreachable once rt is True), so surface it explicitly: it flags a nested/composite locus rather than a
    # confident single element. Mirrors the tpase_conflict downgrade below.
    composite = rt and tpase
    if composite:
        ev.append("transposase domain also present alongside reverse transcriptase — possible nested / composite "
                  "element (e.g. a DNA transposon within a retroelement, or two overlapping TE fragments)")

    # envelope (env) domain + paired LTRs mark an ENDOGENOUS RETROVIRUS (ERV) / errantivirus — an env-bearing
    # LTR retroelement (potentially infection-capable), distinct from a plain LTR retrotransposon.
    has_env = "ENV" in cset
    gag_core = any(d.get("hmm") in _GAG_CORE for d in domains)   # capsid/matrix, not nucleocapsid zf-CCHC alone
    # a DIRS-group call (RT + tyrosine recombinase, no DDE integrase) is never an ERV, even with a spurious
    # env hit — routing it through the retroviral transcript-architecture model would contradict its own class.
    is_erv = bool(rt and has_ltr and has_env and not (yr and not intg))
    if is_erv:
        ev.append("envelope (env) domain present with paired LTRs → endogenous retrovirus (ERV) / errantivirus "
                  "lineage (an env-bearing LTR retroelement)")

    # domain-architecture order shows ONLY the domains actually detected, in the superfamily's
    # canonical order — never a full template presented as if it were observed. ENV closes the retroviral order.
    _CANON = {"Copia": ["GAG", "PR", "INT", "RT", "RNaseH", "CHR", "ENV"],
              "Gypsy": ["GAG", "PR", "RT", "RNaseH", "INT", "CHR", "ENV"]}
    _sf = superfamily.split(" ")[0] if superfamily else ""
    if _sf in _CANON:
        order = "–".join(c for c in _CANON[_sf] if c in cset) or None

    ndom = len(cset)
    is_ltr_super = bool(rt and superfamily and superfamily.split(" ")[0] in ("Copia", "Gypsy") and has_ltr)
    if is_ltr_super and ndom >= 2 and order_resolvable:
        confidence = "High"
    elif is_ltr_super:                              # Copia/Gypsy called but INT/RT order indeterminate -> not High
        confidence = "Moderate"
    elif te_class == "LINE":
        # the 3' tail is the only POSITIVE structural evidence for a LINE here; without it the call rests on
        # an absence (no integrase, no LTRs), which DIRS and Penelope-like elements share. `ndom >= 1` was
        # always true on this branch (RT is a domain), so it never discriminated anything.
        confidence = "Moderate" if has_polya else "Candidate"
    elif (rt or tpase) and (has_ltr or has_tir):
        confidence = "Moderate"
    elif rt or tpase:
        confidence = "Moderate" if ndom >= 2 else "Candidate"
    else:
        confidence = "Candidate"
    if tpase_conflict:                             # conflicting transposase classes -> don't overstate
        confidence = "Candidate"
    if composite and confidence == "High":         # a co-occurring transposase makes a confident single-element call unsafe
        confidence = "Moderate"

    completeness = _completeness(cset, rt, intg, tpase, has_ltr, has_tir, has_polya, is_erv, order_resolvable,
                                 gag_core, tir_ok, en_ok=en_ok, yr=yr)

    dom_str = " · ".join(codes) or "none"
    struct_str = ", ".join(e["type"].split(" (")[0] for e in structural) or "none"
    explanation = (f"Classified as {te_class} ({confidence.lower()} confidence). "
                   f"Structural: {struct_str}. Domains: {dom_str}." + (" " + "; ".join(ev) + "." if ev else ""))
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "confidence": confidence, "evidence": ev, "explanation": explanation, "n_domains": ndom,
            "is_erv": is_erv, "completeness": completeness}


# Domain families TEagle can test (its bundled Pfam profile panel). A module absent from a result is only
# meaningfully "missing" relative to THIS panel — a divergent or unmodelled domain reads as not-detected, not decay.
DOMAINS_TESTED = ("gag (matrix/capsid/nucleocapsid), protease, RT, RNase H, DDE integrase, envelope, "
                  "chromodomain, LINE ORF1p and apurinic-like endonuclease, tyrosine recombinase, "
                  "Helitron helicase, and transposases of the Tc1/Mariner, hAT, CACTA, MULE and IS4 groups")


def _completeness(cset, rt, intg, tpase, has_ltr, has_tir, has_polya, is_erv, order_resolvable, gag_core=False,
                  tir_ok=False, en_ok=False, yr=False):
    """A CATEGORICAL structural-completeness call (never a fabricated numeric score), scoped to the models tested.
    Tiers map to established terms: an element with its expected coding architecture + intact structural context is
    'intact / autonomous-consistent' (Wicker 2007 autonomous; TEsorter Complete; LTR_retriever intact); a core
    module missing is 'partial'; terminal repeats with no coding is 'structural-only'. Every branch derives its tier
    from an explicit expected/present/missing ledger, so the tier can be contradicted by its own bookkeeping — for a
    DNA transposon the ledger holds the transposase AND both terminal inverted repeats, not the transposase alone.
    The tier describes how much of the expected ARCHITECTURE is present at the domain level IN ONE SEQUENCE. It is
    not a claim that the ORFs are functional, that the element is transcribed, that it retains transposition or
    infection competence, or that any individual carries the insertion — the three conflations Lanciano & Cristofari
    2020 (Nat Rev Genet 21:721-736) identify as the field's central interpretive trap."""
    if rt and yr and not intg:                            # DIRS-group: tyrosine recombinase in place of DDE integrase
        # A DIRS-group element encodes Gag, an RT-RNaseH pol module and a tyrosine recombinase (YR); it lacks
        # both the DDE integrase and the aspartic protease of an LTR retrotransposon (Wicker 2007; Poulter &
        # Goodwin 2005 Cytogenet Genome Res 110:575-588). Scoring it against the LTR ledger would flag INT — a
        # domain the class never has — as permanently missing and omit YR entirely. Ledger scoped to the
        # diagnostic modules the panel tests; the split/inverted terminal repeats are structural, not domains.
        expected = ["GAG", "RT", "RNaseH", "YR"]
        present = [m for m in expected if m in cset]
        missing = [m for m in expected if m not in cset]
        # RT + YR are guaranteed on this branch; the top tier additionally needs a CAPSID/matrix gag
        # (gag_core), not a promiscuous zf-CCHC hit — the same discipline the LTR/ERV branch applies.
        tier = ("intact / autonomous-consistent" if (gag_core and not missing)
                else "partial (DIRS core present: RT + tyrosine recombinase)")
        kind = "DIRS-group retroelement (tyrosine recombinase)"
    elif rt and (has_ltr or intg):                        # LTR retroelement / ERV
        expected = ["GAG", "PR", "RT", "RNaseH", "INT"] + (["ENV"] if "ENV" in cset else [])
        present = [m for m in expected if m in cset]
        missing = [m for m in expected if m not in cset]
        # autonomous signature = a CORE gag (capsid/matrix) + RT + integrase; a nucleocapsid-only (zf-CCHC) gag hit
        # is detected and shown in `present` but does not by itself earn the intact/near-complete tier.
        core_ok = gag_core and rt and intg
        if core_ok and not missing and has_ltr and order_resolvable:
            tier = "intact / autonomous-consistent"
        elif core_ok and has_ltr:
            tier = "near-complete"
        elif core_ok:                                     # gag core + pol present; LTRs / order just not confirmed
            tier = "coding core present (gag + pol" + ("; env" if "ENV" in cset else "") + "); LTRs not confirmed"
        elif rt and intg:
            # keep the tier a short label — the results banner renders the 'not detected' list itself, so
            # embedding it here too would double-render it. Only flag the nucleocapsid-only case, which the
            # missing list cannot express (gag IS present, just not a capsid/matrix).
            tier = "partial (pol core present)" if missing else "partial (pol core present; gag capsid/matrix not confirmed)"
        else:
            tier = "fragment"
        if is_erv:
            kind = "ERV (env-bearing LTR retroelement)"
        elif has_ltr:
            kind = "LTR retroelement"
        else:
            kind = "retroelement (LTR not confirmed)"
    elif rt:                                              # non-LTR (LINE-like)
        # An autonomous LINE carries ORF1 (nucleic-acid chaperone) plus ORF2 (apurinic-like endonuclease
        # then RT in one frame), and target-primed reverse transcription leaves a 3' tail
        # (Wicker 2007; Ostertag & Kazazian 2001 Annu Rev Genet 35:501-538). All four are now testable,
        # so the ledger carries all four: a full-length L1 and a 5'-truncated fragment no longer return
        # the same verdict. EN is credited only in the ORF2p arrangement (see en_ok).
        _TAIL = "3′ tail (poly-A / poly-T)"
        expected = ["ORF1", "EN", "RT", _TAIL]
        present = [m for m in expected if (m == "RT"
                                           or (m == "ORF1" and "ORF1" in cset)
                                           or (m == "EN" and en_ok)
                                           or (m == _TAIL and has_polya))]
        missing = [m for m in expected if m not in present]
        if not missing:
            tier = "intact / autonomous-consistent"
        elif "ORF1" in present and "EN" in present:
            tier = "near-complete (ORF1 + ORF2 endonuclease and RT recovered)"
        elif "ORF1" in present or "EN" in present:
            tier = "partial (LINE coding core incomplete)"
        else:
            tier = "partial (RT only)"
        kind = "non-LTR retroelement (LINE-like)"
    elif tpase:                                           # DNA transposon (TIR-type, cut-and-paste)
        # an AUTONOMOUS cut-and-paste transposon needs its transposase ORF *and* both terminal inverted repeats —
        # the cis ends its own product binds and excises (Wicker 2007). Both arms are resolved together by the TIR
        # scan and are credited only when they enclose the transposase, so a locus whose ends were never recovered
        # (5'/3'-truncated, internally deleted, or too diverged to detect) now fails the ledger instead of being
        # asserted intact from the transposase alone.
        expected = ["TPase", "TIR (5′)", "TIR (3′)"]
        present = ["TPase"] + (["TIR (5′)", "TIR (3′)"] if tir_ok else [])
        missing = [m for m in expected if m not in present]
        # keep the tier a short label — the results banner renders the 'not detected' list itself
        tier = "intact / autonomous-consistent" if not missing else "partial (transposase present)"
        kind = "DNA transposon" if tir_ok else "DNA transposon (terminal repeats not confirmed)"
    elif has_ltr or has_tir:
        kind = "LTR retroelement" if has_ltr else "DNA transposon"
        if cset:                                          # coding domain(s) recovered but no RT/transposase — name them
            present = _ordered(cset)
            expected, missing = list(present), []
            tier = f"coding fragment ({'–'.join(present)} detected; RT/pol not confirmed)"
        else:
            expected, present, missing = [], [], []
            tier = "structural-only (terminal repeats, no coding domain detected)"
    elif cset:                                            # coding domain(s) only — no RT/tpase, no terminal repeats
        present = _ordered(cset)
        expected, missing = list(present), []
        tier = f"coding fragment ({'–'.join(present)} detected; no structural context)"
        kind = "coding fragment"
    else:
        return None
    return {"tier": tier, "kind": kind, "present": present, "missing": missing,
            "expected": expected, "scope": DOMAINS_TESTED}
