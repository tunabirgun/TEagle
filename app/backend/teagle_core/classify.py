"""TE classification from structural evidence + protein-domain architecture (Layer C).
Superfamily calls follow domain content, structural context, and — for LTR elements —
the diagnostic integrase-vs-RT order (Copia: INT before RT; Gypsy: INT after RT).
Transparent, evidence-derived; never a family/name call beyond what the evidence supports."""
from __future__ import annotations


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

_CODE_ORDER = ["GAG", "PR", "RT", "RNaseH", "INT", "CHR", "ENV", "TPase"]


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

    if rt:
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
            ev.append("RT without integrase and without LTRs → non-LTR retrotransposon (LINE)")
            if has_polya:                                     # name the tail that was actually detected, not always poly-A
                _pa = any(e["type"].startswith("poly-A") for e in structural)
                ev.append(("3′ poly-A tail" if _pa else "5′ poly-T tract") + " consistent with LINE")
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
        elif "DDE" in bcl:
            superfamily = "DDE transposon"
        else:
            superfamily = "DNA transposon"
        te_class = "DNA/" + superfamily.split("/")[0]
        ev.append("transposase domain present → Class II DNA transposon")
        tpase_conflict = len({d["class"] for d in tp_hits}) > 1
        if tpase_conflict:
            ev.append("multiple transposase classes detected — superfamily assigned from the strongest-scoring hit (ambiguous)")
        if has_tir:
            ev.append("terminal inverted repeats consistent with a cut-and-paste transposon")
        # the ends are half of what makes a DNA transposon autonomous, so state what was actually recovered
        if tir_ok and has_tsd:
            ev.append("a target-site duplication flanks the inverted repeats — the insertion site itself is "
                      "captured, so both element termini are present in the record")
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
            superfamily, te_class = "DNA transposon (TIR, no transposase detected)", "DNA/structural-only"
            ev.append("terminal inverted repeats but no transposase recovered")
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
    is_erv = bool(rt and has_ltr and has_env)
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
    elif te_class == "LINE" and (has_polya or ndom >= 1):
        confidence = "Moderate"
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
                                 gag_core, tir_ok)

    dom_str = " · ".join(codes) or "none"
    struct_str = ", ".join(e["type"].split(" (")[0] for e in structural) or "none"
    explanation = (f"Classified as {te_class} ({confidence.lower()} confidence). "
                   f"Structural: {struct_str}. Domains: {dom_str}." + (" " + "; ".join(ev) + "." if ev else ""))
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "confidence": confidence, "evidence": ev, "explanation": explanation, "n_domains": ndom,
            "is_erv": is_erv, "completeness": completeness}


# Domain families TEagle can test (its bundled Pfam profile panel). A module absent from a result is only
# meaningfully "missing" relative to THIS panel — a divergent or unmodelled domain reads as not-detected, not decay.
DOMAINS_TESTED = "gag (matrix/capsid/nucleocapsid), protease, RT, RNase H, integrase, envelope, chromodomain, transposases"


def _completeness(cset, rt, intg, tpase, has_ltr, has_tir, has_polya, is_erv, order_resolvable, gag_core=False,
                  tir_ok=False):
    """A CATEGORICAL structural-completeness call (never a fabricated numeric score), scoped to the models tested.
    Tiers map to established terms: an element with its expected coding architecture + intact structural context is
    'intact / autonomous-consistent' (Wicker 2007 autonomous; TEsorter Complete; LTR_retriever intact); a core
    module missing is 'partial'; terminal repeats with no coding is 'structural-only'. Every branch derives its tier
    from an explicit expected/present/missing ledger, so the tier can be contradicted by its own bookkeeping — for a
    DNA transposon the ledger holds the transposase AND both terminal inverted repeats, not the transposase alone.
    The tier describes how much of the expected ARCHITECTURE is present at the domain level — it is not a claim that
    the ORFs are functional."""
    if rt and (has_ltr or intg):                          # LTR retroelement / ERV
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
        expected, present = ["RT"], ["RT"]
        missing = []
        tier = "coding present (RT); non-LTR element" if has_polya else "coding present (RT)"
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
