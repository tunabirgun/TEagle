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


def classify(structural, domains):
    has_ltr = any(e["type"].startswith("LTR") for e in structural)
    has_tir = any(e["type"].startswith("TIR") for e in structural)
    has_polya = any(e["type"].startswith("poly") for e in structural)
    codes = [d["domain"] for d in domains]
    cset = set(codes)
    rt_d, int_d = _rep(domains, "RT"), _rep(domains, "INT")
    rt, intg = rt_d is not None, int_d is not None
    tpase = "TPase" in cset
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
    else:
        klass = "unclassified"
        if has_ltr:
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

    # domain-architecture order shows ONLY the domains actually detected, in the superfamily's
    # canonical order — never a full 5-domain template presented as if it were observed
    _CANON = {"Copia": ["GAG", "PR", "INT", "RT", "RNaseH", "CHR"],
              "Gypsy": ["GAG", "PR", "RT", "RNaseH", "INT", "CHR"]}
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

    dom_str = " · ".join(codes) or "none"
    struct_str = ", ".join(e["type"].split(" (")[0] for e in structural) or "none"
    explanation = (f"Classified as {te_class} ({confidence.lower()} confidence). "
                   f"Structural: {struct_str}. Domains: {dom_str}." + (" " + "; ".join(ev) + "." if ev else ""))
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "confidence": confidence, "evidence": ev, "explanation": explanation, "n_domains": ndom}
