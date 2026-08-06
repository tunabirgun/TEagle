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

# Retained for reference: the shortest record that could plausibly hold an autonomous non-LTR element.
# Autonomous LINEs run 3-6 kb (L1 ~6 kb, R2Bm 3.6 kb). Length is NOT used to license a LINE call - a
# multi-kilobase record whose sole detected domain is a reverse transcriptase is equally consistent with a
# host telomerase, and the bundled panel carries one RT model rather than a clade-resolved set, so the two
# cannot be separated on that evidence. The order is withheld instead.
_MIN_LINE_RECORD_BP = 2000


def _ordered(cset):
    # detected domain codes in canonical retroviral order, unknown codes appended alphabetically
    return [c for c in _CODE_ORDER if c in cset] + sorted(c for c in cset if c not in _CODE_ORDER)


def _refine_tsd(tsd, superfamily, seq):
    """Re-detect the TSD now that the superfamily is known, preferring its literature target-site length.
    detect_all runs before classification, so it picked the LONGEST exact flanking repeat — which a
    coincidental longer direct repeat in an AT-rich flank can win over the real short TSD (e.g. a Tc1/Mariner
    2 bp TA). Re-running find_tsd with `expect` returns the expected length only when it genuinely flanks
    (and it is always <= the longest pick), so this can only shorten a coincidental TSD to the diagnostic
    one, never lengthen or fabricate. No-op without a sequence, without an expected length, or when unchanged."""
    if not seq or not tsd:
        return
    up, dn = tsd.get("upstream"), tsd.get("downstream")
    if not up or not dn:                                  # a real TSD always carries flank coordinates
        return
    exp = structural_mod.tsd_congruence(tsd.get("length"), superfamily).get("expected")
    if not exp:                                           # no attributable expected length for this superfamily
        return
    refined = structural_mod.find_tsd(seq, up[1], dn[0], expect=exp)
    if refined and refined["length"] == exp and refined["length"] != tsd.get("length"):
        tsd.update(refined)                              # correct length/motif/coords/matched_expected in place


def classify(structural, domains, seq=None, domains_ok=True, orfs_unscanned=0, orfs=None):
    """`domains_ok=False` means the protein-domain scan RAISED, so `domains` is empty because nothing ran —
    not because nothing is there. Without that distinction every domain-derived negative ("not detected:
    GAG, PR, RNaseH") is asserted from a scan that never happened, which is absence-of-evidence sold as
    evidence-of-absence. The result is marked so the card can say "not assessed" instead of "not detected"."""
    # Advisory rows are REPORTED evidence, never CREDITED evidence. A sub-threshold terminal repeat is
    # shown to the user precisely because it did not qualify, so it must not satisfy any has_* test —
    # otherwise the classifier treats a rejection as a confirmation. Belt and braces with the naming rule
    # in structural.py: a future advisory type that happens to start with a credited prefix is caught here.
    credited = [e for e in structural if not e.get("advisory")]
    has_ltr = any(e["type"].startswith("LTR") for e in credited)
    tir_ev = next((e for e in credited if e["type"].startswith("TIR")), None)
    has_tir = tir_ev is not None
    has_tsd = any(e["type"].startswith("TSD") for e in credited)
    has_polya = any(e["type"].startswith("poly") for e in credited)
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
    hel = "HEL" in cset          # Helitron Rep/Helicase (PF14214) — bundled and advertised, so it must be consumed
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
            if not order_resolvable:
                # Do NOT name a superfamily from a comparison that cannot carry the answer. _pos negates the
                # minus-strand coordinate to get translation order, which is coherent WITHIN one strand and
                # meaningless ACROSS two: a minus-strand hit always yields a negative value, so the test
                # reduces to "which of the two is on the minus strand" and returns the SAME answer whether
                # INT sits upstream or downstream of RT (measured: INT-on-minus -> always Copia,
                # RT-on-minus -> always Gypsy, at every coordinate). Naming a tentative superfamily from that
                # puts a family in the card's headline on evidence containing no information about it, which
                # AGENTS.md rule 1 forbids. The overlapping-span case is equally unreadable. Fall through to
                # the same honest label the no-integrase path already uses.
                superfamily, te_class = "LTR retrotransposon (superfamily undetermined)", "LTR/unclassified"
                ev.append("integrase and RT are not in a readable arrangement (different strands, or "
                          "overlapping spans), so the diagnostic integrase-vs-RT order cannot be determined — "
                          "Copia/Gypsy is not called")
            elif _pos(int_d) < _pos(rt_d):
                superfamily, te_class = "Copia (Ty1)", "LTR/Copia"
                ev.append("integrase N-terminal to RT + paired LTRs → Copia/Ty1 order")
            else:
                superfamily = "Gypsy (Ty3)"
                if "CHR" in cset:
                    superfamily += " · chromovirus"
                te_class = "LTR/Gypsy"
                ev.append("integrase C-terminal to RT + paired LTRs → Gypsy/Ty3 order")
            if "PR" in cset:                                  # cset holds emitted domain codes; domains.py maps RVP -> code "PR"
                ev.append("aspartic-protease domain present")
            if "RNaseH" in cset:
                ev.append("RNase H domain present")
            # TSD-length congruence, same discipline as the DNA branch: Copia/Gypsy duplicate a 5 bp target
            # site (Ou 2019), so a flanking repeat of another length is likely coincidental and its ends are
            # not credited. TSD_EXPECT already carries the Copia/Gypsy lengths; this wires them for LTR too.
            if has_tsd:
                _ltsd = next((e for e in structural if e["type"].startswith("TSD")), None)
                _refine_tsd(_ltsd, superfamily, seq)      # prefer the superfamily's target-site length over a coincidental longer flank
                _lcong = structural_mod.tsd_congruence((_ltsd or {}).get("length"), superfamily)
                if _ltsd is not None:
                    _ltsd["tsd_congruence"] = _lcong["verdict"]
                if _lcong["verdict"] == "incongruent":
                    ev.append(f"a flanking target-site duplication is present but {_lcong['observed']} bp where "
                              f"{superfamily.split(' ')[0]} elements duplicate {_lcong['expected']} bp "
                              f"({_lcong['basis']}) — likely coincidental, so the ends are not credited from it")
                elif _lcong["verdict"] == "congruent":
                    ev.append(f"a {_lcong['observed']} bp target-site duplication flanks the element — congruent "
                              f"with the {superfamily.split(' ')[0]} target-site length ({_lcong['basis']})")
        elif not has_ltr and not intg and "CHR" in cset:
            # A chromodomain is carried by the chromoviruses, a Ty3/Gypsy lineage, and is not carried by
            # non-LTR elements. Reaching the LINE branch on an absent integrase alone therefore filed
            # chromovirus fragments as LINEs whenever the integrase went undetected — a Class I order error
            # from evidence that positively excluded it. The integrase is what resolves Copia from Gypsy, so
            # without it the superfamily is still withheld; only the order is named here.
            klass = "Class I · retrotransposon"
            superfamily = "LTR retrotransposon (superfamily undetermined)"
            te_class = "LTR/unclassified"
            ev.append("chromodomain present with reverse transcriptase — the chromoviruses are a Ty3/Gypsy "
                      "lineage and non-LTR elements do not carry a chromodomain, so this is an LTR "
                      "retrotransposon even though no integrase was detected and no terminal repeat was "
                      "accepted")
            ev.append("the superfamily is not called: Copia versus Gypsy needs the integrase-versus-RT "
                      "translation order, and no integrase was detected")
        elif (not has_ltr and not intg
              and not (has_polya or "ORF1" in cset or "EN" in cset or "RNaseH" in cset)):
            # A LONE reverse transcriptase is not a LINE. The branch below reads an absent integrase as
            # positive evidence for a non-LTR element, which holds only when something else about the record
            # is LINE-like: a poly-A tail from target-primed reverse transcription, an ORF1p, an
            # endonuclease, or an RNase H. With none of those the only fact in evidence is "a reverse
            # transcriptase is present", which every retroelement order satisfies — and which host genes
            # satisfy too. A telomerase reverse transcriptase presents exactly this way, as a multi-kilobase
            # record whose sole detected domain is RT, and cannot be separated from an R2-clade element on
            # that evidence: both are reverse transcriptases and the panel carries one RT model, not a
            # clade-resolved set. The order is therefore withheld rather than defaulted to LINE.
            superfamily, te_class = "retrotransposon (partial)", "retro/partial"
            ev.append("a reverse transcriptase was detected with no integrase, no terminal repeat, and none "
                      "of the features that distinguish a non-LTR element — no poly-A tail, no ORF1p, no "
                      "ORF2p endonuclease–RT arrangement, no RNase H. Every retrotransposon order carries a "
                      "reverse transcriptase, so the order is left unassigned rather than defaulted to LINE")
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
            # RT + a DDE integrase, but no terminal repeat was accepted. The order is still determined:
            # a DDE integrase is diagnostic of the LTR order and is not carried by non-LTR elements, which
            # is the same discriminator the LINE branch above relies on in its negative form. Refusing to
            # name the order here treated "the repeat was not recovered" as "the order is unknown", which
            # are different statements — a repeat goes unrecovered on a long deposit, on a truncated one,
            # or on an element old enough for its copies to have diverged past the identity floor.
            #
            # Measured: on the benchmark corpus this branch fired on 16 elements whose literature label is
            # LTR, all of them carrying RT + INT, and a domain-architecture classifier (TEsorter) named
            # every one of them correctly from the same sequence. The superfamily is still withheld —
            # that needs the INT-vs-RT translation order, which is a separate question from the order.
            superfamily, te_class = "LTR retrotransposon (superfamily undetermined)", "LTR/unclassified"
            ev.append("RT + DDE integrase present → LTR order, from domain architecture; no terminal "
                      "repeat was accepted, so the assignment does not rest on structural evidence")
            if "ENV" in cset:                     # has_env is not bound until later in this function
                ev.append("envelope glycoprotein present — carried by LTR elements and never by non-LTR "
                          "retrotransposons, corroborating the order")
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
        # take the leading token before any parenthetical ("CACTA (En/Spm)" -> "CACTA"), but KEEP a
        # compound superfamily name whose "/" joins two members ("Tc1/Mariner"): rewrite that "/" to "-"
        # so the full superfamily survives without a second slash (the "/" in te_class separates class
        # from subclass). Truncating at "/" dropped "Mariner" and asserted a narrower call than the
        # domain evidence supports (Tc1/Mariner is one Wicker superfamily).
        _sf_token = superfamily.split(" ")[0].replace("/", "-")
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
            _refine_tsd(_tsd, superfamily, seq)               # prefer the superfamily's target-site length over a coincidental longer flank
            _cong = structural_mod.tsd_congruence((_tsd or {}).get("length"), superfamily)
            if _tsd is not None:
                _tsd["tsd_congruence"] = _cong["verdict"]     # carry the verdict to the record (GFF3 export reads it)
            if _cong["verdict"] == "incongruent":
                ev.append(f"a target-site duplication flanks the inverted repeats, but it is "
                          f"{_cong['observed']} bp where {superfamily} elements duplicate "
                          f"{_cong['expected']} bp ({_cong['basis']}) — the repeat may be coincidental, or the "
                          f"boundary or superfamily call may be wrong, so the ends are not credited from it")
            elif _cong["verdict"] == "congruent":
                ev.append("a target-site duplication flanks the inverted repeats — the insertion site itself is "
                          f"captured, so both element termini are present in the record (length congruent with "
                          f"{superfamily}, {_cong['basis']})")
            else:
                # a TSD flanks the repeats but this superfamily has no literature-attributed expected length,
                # so the length cannot corroborate the termini — say that rather than imply congruence.
                ev.append("a target-site duplication flanks the inverted repeats — the insertion site itself is "
                          "captured, so both element termini are present in the record; no literature-attributed "
                          f"expected target-site length exists for {superfamily} to corroborate its length")
        elif has_tir and tir_encloses_tpase is False:
            ev.append("the detected inverted-repeat pair does not enclose the transposase — the termini and the "
                      "coding module may not belong to the same element, so the ends are not credited")
        elif not has_tir:
            ev.append("no terminal inverted-repeat pair recovered — neither element end is confirmed, which is "
                      "consistent with a 5′- or 3′-truncated (or internally deleted) copy, with a superfamily "
                      "that does not carry TIRs, or with termini too diverged for the scan")
        order = "–".join(_ordered(cset))     # coding architecture only — the TIR ends are structural, not domains
    elif hel:
        # Helitrons are Wicker Class II, SUBCLASS 2 — they move by rolling-circle replication, not by a
        # cut-and-paste DDE transposase, so they legitimately have no TPase, no TIR and no TSD (they insert
        # between an A and a T and duplicate nothing). Without a branch here the HEL hit fell through to the
        # generic fragment arm and came out as te_class "retro/partial" — a Class II element filed under
        # Class I retro — with the label "RT/transposase not detected", which reads as "nothing diagnostic
        # found" on a record where the diagnostic domain WAS found. The profile is bundled (domains.py:70)
        # and DOMAINS_TESTED advertises it, so the panel promised a call it could not make.
        # DELIBERATELY NOT a Helitron superfamily call. tests/test_docs_track_build.py
        # (test_doc_does_not_claim_a_helitron_superfamily_call) and both capability docs record a standing
        # decision that naming a Helitron superfamily from the helicase domain alone would over-state the
        # tool, and that decision is not this loop's to reverse. What IS fixed here is the misfiling: with
        # no branch at all, a HEL hit fell through to the generic arm and came out as te_class
        # "retro/partial" — a Class II rolling-circle element filed under Class I retro — labelled
        # "RT/transposase not detected" on a record where the diagnostic domain WAS found. Declining to
        # name the superfamily is honest; calling it a retroelement is simply wrong.
        klass = "unclassified"
        superfamily = "Helitron-family helicase present — superfamily not called"
        te_class = "unclassified"
        ev.append("Helitron-family Rep/Helicase domain (PF14214) present. Helitrons move by rolling-circle "
                  "replication, so they carry no terminal inverted repeat and duplicate no target site — "
                  "the absence of TIR/TSD here is expected and is not evidence against a Helitron")
        # The structural signature that would confirm it — a 5' TC / 3' CTRR terminus and the ~11-20 nt
        # hairpin upstream of the 3' end — is NOT among TEagle's detectors, so say so rather than let the
        # domain hit alone read as a complete Helitron call.
        ev.append("the diagnostic Helitron termini (5′ TC … CTRR 3′ and the subterminal hairpin) are not "
                  "among the tested structural detectors, so TEagle does not assign a Helitron superfamily "
                  "from the helicase domain alone — but this is NOT a retroelement, and the class is left "
                  "unassigned rather than guessed")
    else:
        klass = "unclassified"
        if cset:                                          # coding domain(s) recovered but no RT and no transposase
            _cs = "–".join(_ordered(cset))                # e.g. an ERV relic that kept env/capsid but lost pol
            # A domain shared with host enzymes is not, on its own, evidence of a transposable element.
            # The endonuclease model (Pfam Exo_endo_phos) also matches host AP endonucleases and DNase I;
            # the tyrosine-recombinase model matches phage and plasmid integrases and host site-specific
            # recombinases; RNase H and the chromodomain are likewise carried by host proteins. Asserting
            # Class I from any of those alone reported host DNA-repair and recombinase genes as
            # retroelement fragments. A retro call now requires a domain that non-mobile host genes do not
            # carry: a capsid, an envelope glycoprotein, or a DDE integrase.
            _RETRO_SPECIFIC = {"GAG", "ENV", "INT", "ORF1"}
            if has_ltr:
                superfamily, te_class = f"LTR retroelement fragment ({_cs}; RT/pol not detected)", "LTR/partial"
                ev.append(f"coding domain(s) recovered ({_cs}) with paired terminal repeats, but neither "
                          f"reverse transcriptase nor transposase detected")
            elif cset & _RETRO_SPECIFIC:
                superfamily, te_class = f"coding fragment ({_cs}; RT/transposase not detected)", "retro/partial"
                ev.append(f"coding domain(s) recovered ({_cs}) but neither reverse transcriptase nor "
                          f"transposase detected")
            else:
                klass = "unclassified"
                superfamily = f"{_cs} domain(s) detected; no transposable-element assignment"
                te_class = "unclassified"
                ev.append(f"the only coding evidence is {_cs}, which host genes also carry — an "
                          f"endonuclease model matches apurinic endonucleases and DNase I, and a tyrosine "
                          f"recombinase matches phage, plasmid and host site-specific recombinases. Without "
                          f"a reverse transcriptase, a transposase, a capsid, an envelope or an integrase "
                          f"there is no element-specific evidence, so no class is assigned")
        elif has_ltr:
            superfamily, te_class = "LTR retrotransposon (no coding domains detected)", "LTR/structural-only"
            ev.append("paired LTRs but no coding domain recovered")
        elif has_tir:
            # A SHORT inverted terminal repeat without a transposase is not Class II evidence. Verified on
            # the real copia 5' LTR (276 bp): fed alone it trips the terminal-inverted-repeat scan and used
            # to be filed as a DNA transposon — a Class I fragment assigned to Class II. Solo LTRs
            # outnumber full-length elements in most genomes, so this is the common case, not an edge case.
            #
            # A LONG one is different, and the distinction is measured rather than assumed. Chance inverted
            # repeats reached at most 15 bp at the recorded seed and 17 bp across five seeds over 4,500
            # random sequences (benchmarks/chance_tir.py; see structural.MIN_TIR_STANDALONE),
            # and the solo copia LTR that motivated this branch yields 12 bp. Above that ceiling the repeat
            # is not attributable to chance, and a terminal inverted repeat is the defining structural
            # feature of a Class II element. Refusing to call it left every NON-AUTONOMOUS DNA transposon
            # unassigned — a MITE has no transposase by definition, which is what makes it non-autonomous,
            # so requiring one to name the class excluded the whole category by construction.
            _tl = (tir_ev or {}).get("tir_len") or 0
            if _tl >= structural_mod.MIN_TIR_STANDALONE:
                klass = "Class II · DNA transposon"
                superfamily = "non-autonomous TIR element (superfamily undetermined)"
                te_class = "DNA/non-autonomous"
                ev.append(f"a {_tl} bp terminal inverted repeat was detected with no transposase. Chance "
                          f"inverted repeats reached at most {structural_mod.MIN_TIR_STANDALONE - 3} bp "
                          f"over 4,500 random sequences, so a repeat this "
                          f"long is structural evidence of a Class II element; the absence of a transposase "
                          f"is consistent with a non-autonomous derivative such as a MITE or an internally "
                          f"deleted copy, which carry terminal repeats and no coding capacity")
                ev.append("the superfamily is not called: that needs the transposase family, which was not "
                          "detected here")
            else:
                superfamily, te_class = "terminal inverted repeat, class unassigned", "repeat/structural-only"
                ev.append(f"an inverted terminal repeat was detected ({_tl} bp) but no transposase, and a "
                          f"repeat this short arises by chance on random sequence — this does not establish "
                          f"a DNA transposon. A solo LTR, a fragment of a larger element, or an unrelated "
                          f"inverted repeat all produce the same signal, so the class is left unassigned")
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
    # No CHR in the Copia canon: a chromodomain-bearing integrase is a chromovirus feature and
    # chromoviruses are a Ty3/Gypsy lineage, so listing it here would render a spurious PF00385 hit on a
    # Copia call as though it were canonical Copia architecture.
    _CANON = {"Copia": ["GAG", "PR", "INT", "RT", "RNaseH", "ENV"],
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
    # The transposase arm gates on tir_ok, not raw has_tir. This function already refuses to credit a
    # terminal repeat that does not enclose the transposase — the ledger drops the arms and the evidence
    # says "the ends are not credited" — because the termini and the coding module may belong to different
    # elements. Reading that same repeat as enough to lift the badge out of Candidate contradicted the
    # refusal one screen above it: two records with an identical completeness ledger got different badges.
    # A solo LTR fed alone trips the inverted-repeat scan (see the note below), so this is the common case.
    elif (rt and (has_ltr or has_tir)) or (tpase and (has_ltr or tir_ok)):
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
    # The coding axis rides ALONGSIDE the tier, never inside it. The tier value is benchmarked (moving it
    # would move published rows), and the defect was never that the tier miscounted domains — it was that a
    # domain-only ledger was the ONLY thing on screen, so "partial" was read as "decayed". Adding the axis
    # answers the question the reader was actually asking without perturbing a benchmarked number.
    coding = _coding_axis(domains, orfs)
    if coding:
        completeness["coding"] = coding

    dom_str = " · ".join(codes) or "none"
    struct_str = ", ".join(e["type"].split(" (")[0] for e in structural) or "none"
    explanation = (f"Classified as {te_class} ({confidence.lower()} confidence). "
                   f"Structural: {struct_str}. Domains: {dom_str}." + (" " + "; ".join(ev) + "." if ev else ""))
    if orfs_unscanned and domains_ok:
        # A partial version of the domains_ok problem: the search ran, but not over every ORF, so the
        # not-detected list is scoped to what was searched rather than to the sequence.
        ev.append(f"{orfs_unscanned} shorter ORF(s) were outside the profile search's ORF budget and were "
                  "not searched — a domain listed as not detected was not tested on those")
    if not domains_ok:
        # The scan raised. Every domain-derived statement below is therefore uninformative, and the
        # confidence tier — which counts domain evidence — cannot stand. Say so first, in the evidence the
        # card renders, and cap the tier at the lowest rung rather than letting a structural-only record
        # inherit a tier that domain hits were supposed to earn.
        ev.insert(0, "the protein-domain scan did not run for this record, so no domain was tested — "
                     "any domain named or not named below is unassessed, not absent")
        confidence = "Candidate"
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "confidence": confidence, "evidence": ev, "explanation": explanation, "n_domains": ndom,
            "is_erv": is_erv, "completeness": completeness, "domains_unavailable": not domains_ok,
            "orfs_unscanned": int(orfs_unscanned or 0)}


# Domain families TEagle can test (its bundled Pfam profile panel). A module absent from a result is only
# meaningfully "missing" relative to THIS panel — a divergent or unmodelled domain reads as not-detected, not decay.
def _tpase_groups():
    """The DNA-transposon groups the bundled panel can actually name, DERIVED from domains.DOMAIN_INFO.

    Hand-listing them let the sentence fall behind the panel: it named five groups while the panel also
    bundles Pfam's generic DDE family (PF03184, class dna:DDE), which the superfamily logic below names
    outright as "DDE transposon". A user reading the scope line was told a family that CAN be called is not
    tested, which breaks the not-detected-vs-not-tested distinction the line exists to preserve. Helitron
    is excluded here because it is a rolling-circle helicase, not a transposase, and is named separately.
    Imported inside the function: domains imports sequtil/appdirs, and this keeps the dependency one-way."""
    from . import domains as _dom
    rename = {"Tc1-Mariner": "Tc1/Mariner", "DDE": "unclassified DDE"}
    groups = []
    for code, _label, cls, _pfam in _dom.DOMAIN_INFO.values():
        if code != "TPase" or not cls.startswith("dna:"):
            continue
        g = cls.split(":", 1)[1]
        g = rename.get(g, g)
        if g not in groups:
            groups.append(g)
    # the catch-all reads last, however DOMAIN_INFO happens to be ordered
    groups.sort(key=lambda g: (g == "unclassified DDE", g.lower()))
    return groups


def _tpase_group_phrase():
    g = _tpase_groups()
    return (", ".join(g[:-1]) + " and " + g[-1]) if len(g) > 1 else (g[0] if g else "")


DOMAINS_TESTED = ("gag (matrix/capsid/nucleocapsid), protease, RT, RNase H, DDE integrase, envelope, "
                  "chromodomain, LINE ORF1p and apurinic-like endonuclease, tyrosine recombinase, "
                  "Helitron helicase, and transposases of the " + _tpase_group_phrase() + " groups")


def _coding_axis(domains, orfs):
    """Do the DETECTED domains sit in one uninterrupted reading frame?

    A SECOND, independent axis, reported beside the domain ledger rather than folded into it. The two
    answer different questions and the tier word conflated them: "partial" is emitted whenever an expected
    Pfam module was not detected, and a reader takes that as "this copy is decayed". But a module can be
    absent because it is divergent, unmodelled, or outside the bundled panel, in a copy whose reading frame
    is completely intact — which is why the field's other tools ("intact", "Complete") disagree with this
    one on the same element. Neither axis can substitute for the other, so both are shown.

    Returns None when there is nothing to assess (no domains, or no ORFs), rather than a reassuring default."""
    if not domains or not orfs:
        return None
    placed, unplaced, truncated = {}, 0, False
    for d in domains:
        # Attribute a domain to the ORF it was actually FOUND in, not to one that merely contains it.
        # scan_domains records that as `orf` — the index into the same length-sorted find_orfs list. Geometry
        # gets this wrong in exactly the case the axis exists to detect: ORFs overlap, find_orfs sorts
        # longest-first, so a containment search returns the LONGEST enclosing ORF for every domain inside it.
        # Two domains in genuinely different, frame-shifted ORFs both resolved to that one ORF and the axis
        # reported "single reading frame" — reproduced, and the exact opposite of the truth.
        idx = d.get("orf")
        host = orfs[idx] if isinstance(idx, int) and 0 <= idx < len(orfs) else None
        if host is None:
            # no provenance (a synthetic hit, or a record from an older run): fall back to geometry, which is
            # right whenever the ORFs do not overlap and is the best available guess when they do
            nt = d.get("nt") or []
            if len(nt) < 2:
                unplaced += 1
                continue
            ds, de = min(nt[0], nt[1]), max(nt[0], nt[1])
            host = next((o for o in orfs
                         if o.get("strand") == d.get("strand") and o["start"] <= ds and de <= o["end"]), None)
        if host is None:
            unplaced += 1                       # a hit straddling an ORF boundary is itself informative
            continue
        placed.setdefault((host["start"], host["end"], host["strand"], host.get("frame")), []).append(d["domain"])
        truncated = truncated or bool(host.get("open_end"))
    if not placed:
        return {"frames": 0, "unplaced": unplaced, "truncated_frame": False,
                "state": "not assessed",
                "note": "no detected domain could be placed inside a single called ORF"}
    n = len(placed)
    state = "single reading frame" if n == 1 and not unplaced else "split across reading frames"
    if n == 1 and not unplaced:
        note = (f"all {sum(len(v) for v in placed.values())} detected domain(s) lie in ONE uninterrupted ORF — "
                "the coding structure of what was found is intact, whatever the domain ledger above lists as "
                "not detected")
    else:
        note = (f"the detected domains fall in {n} separate reading frame(s)"
                + (f", and {unplaced} could not be placed in any single ORF" if unplaced else "")
                + " — consistent with an interrupted, rearranged or degenerate copy, though a genuine "
                  "multi-ORF architecture (LINE ORF1/ORF2, some ERVs) also looks like this")
    if truncated:
        note += "; at least one hosting ORF runs to the end of the supplied sequence, so it may be truncated"
    return {"frames": n, "unplaced": unplaced, "truncated_frame": truncated, "state": state, "note": note}


def _completeness(cset, rt, intg, tpase, has_ltr, has_tir, has_polya, is_erv, order_resolvable, gag_core=False,
                  tir_ok=False, en_ok=False, yr=False):
    """A CATEGORICAL domain-completeness call (never a fabricated numeric score), scoped to the models tested.
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
        elif gag_core and rt and has_ltr and not intg:
            # capsid/matrix gag + RT + paired LTRs, but the DDE integrase specifically was not detected — a
            # common state in aged/degenerate copies whose integrase ORF diverged first. Without this branch the
            # record collapsed to the bare "fragment" tier, ranking BELOW a 2-domain RT+INT "partial" despite
            # carrying more recovered core architecture. Integrase is a core module, so this is a 'partial' tier.
            tier = "partial (gag + RT core present; integrase not confirmed)"
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
