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


def _same_strand(domains, a, b):
    # is there an `a` hit and a `b` hit on one strand? (RT with a tyrosine recombinase is DIRS only when co-oriented)
    return any(x.get("domain") == a and y.get("domain") == b and x.get("strand", "+") == y.get("strand", "+")
               for x in domains for y in domains)


def _wicker_order(te_class, outcome):
    # the Wicker order the call sits in (or its withheld/unassigned state), as distinct from the domain architecture
    if outcome == "withheld":
        return "withheld"
    tc = te_class or ""
    if tc.startswith("LTR"):
        return "LTR"
    if tc.startswith("DNA/"):
        return "TIR"
    return tc if tc in ("LINE", "DIRS") else "unassigned"


def _tpase_family(tp_hits):
    """Family of the strongest transposase hit (not a fixed hAT-first precedence): (superfamily, te_class,
    classes_conflict). Shared by the structure-aware and structure-free paths."""
    best = max(tp_hits, key=lambda d: d.get("score", 0.0))
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
    token = superfamily.split(" ")[0].replace("/", "-")
    # the generic "DNA transposon" fallback would collapse to a redundant "DNA/DNA"; name it honestly
    return (superfamily, "DNA/" + (token if token != "DNA" else "unclassified"),
            len({d.get("class") for d in tp_hits}) > 1)


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

# Profiles of retrovirus-specific capsid/matrix and envelope proteins. Retroviruses and ERVs share the
# Gypsy pol order (PR-RT-RH-INT), so on a hit from one of these the integrase-after-RT test cannot separate
# Gypsy from Retrovirus/ERV (Wicker 2007 sibling superfamilies) and the superfamily is withheld. TLV_coat is
# excluded on purpose: it also annotates the envelope of insect Ty3/gypsy elements (errantiviruses).
_RETROVIRAL_PROFILES = frozenset({"Gag_p24", "Gag_p24_C", "Gag_p10", "HERV-K_env_2", "GP41"})
RETROVIRAL_UNDETERMINED = ("LTR retrotransposon (superfamily undetermined; retroviral-type gag/env — Gypsy, "
                           "Bel-Pao and Retrovirus/ERV are not separated by the bundled panel)")

# Longest chance inverted repeat over five seeds and the number of random sequences behind it
# (benchmarks/chance_tir.py). Constants because the benchmark output is not shipped; a claim about the live
# detector is checked by a seeded test against the detector itself, not against the recorded file.
CHANCE_TIR_MAX_BP = 17
CHANCE_TIR_SEQUENCES = 22500

# What has to accompany a reverse transcriptase before a Class I call is licensed.
#
# A reverse transcriptase on its own does not identify a transposable element. Telomerase carries one,
# group II intron maturases carry one, and the bundled panel holds a single RT model rather than a
# clade-resolved set, so it cannot say which kind of reverse transcriptase it found. Each code below is a
# structural or coding module that a retrotransposon carries and a host reverse-transcribing gene does
# not. LINE order requires LINE-specific evidence: ORF1p or an accepted ORF2p endonuclease--RT
# arrangement. Core gag, protease, envelope and RNase H can support a retroelement call but cannot choose
# LTR versus non-LTR order when terminal repeats and integrase are absent.
#
# Terminal repeats, a DDE integrase, a chromodomain and a tyrosine recombinase are all handled by
# branches that run before this test, so they are deliberately absent here rather than forgotten.
#
# A poly-A tail is deliberately NOT in this set. Every mRNA and cDNA deposit has one, so it distinguishes
# nothing; admitting it is what let a telomerase mRNA be called a LINE.
_LINE_CORROBORATION = {"ORF1"}
_RETRO_CORROBORATION = {"RNaseH", "PR", "ENV"}

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


def _recover_target_site(structural, tir_ev, superfamily, seq):
    """A Tc1/Mariner TIR on a flanked record runs through the palindromic TA target site, so detect_all reports
    the arms 2 bp too long and no TSD. Only a caller that knows the superfamily may trim them (structural.
    palindromic_target_site); the rows are updated in place so the caller's list and the verdict stay in step.
    Returns True when the TIR row was trimmed or a TSD row was set."""
    motif = structural_mod.TSD_EXPECT["Tc1/Mariner"][3]
    found = structural_mod.palindromic_target_site(seq, tir_ev, motif=motif) if seq and tir_ev else None
    if not found or not found[1]:
        return False
    tir, tsd = found
    tir_ev.clear()
    tir_ev.update(tir)
    old = next((e for e in structural if e["type"].startswith("TSD") and not e.get("advisory")), None)
    if old is None:
        structural.insert(next(i for i, e in enumerate(structural) if e is tir_ev) + 1, tsd)
    else:
        old.clear()
        old.update(tsd)
    return True


def _tsd_expectation(c):
    lo, hi = c["expected_range"]
    return (f"{lo} bp" if lo == hi else f"{lo}-{hi} bp") + (f" ({c['expected_motif']})" if c.get("expected_motif") else "")


def _tsd_observed(c):
    lo, hi = c["expected_range"]
    return f"{c['observed']} bp" + (f" ({c['observed_motif']})" if c.get("observed_motif") and lo <= c["observed"] <= hi else "")


def _tsd_chance_clause(tsd):
    """Why a congruent or undiagnostic TSD does not establish the insertion site, or None when it does."""
    flag = (tsd or {}).get("chance_level")
    if flag is False:
        return None
    if flag is None:
        return "its chance probability was not assessed"
    p = tsd.get("chance_probability")
    return (f"a flanking direct repeat of {tsd.get('length')} bp is short enough to arise by chance" +
            (f" (p = {p:.2g})" if p is not None else ""))


def _scan_gap(seq, kind):
    """Parenthetical saying which part of a long record the terminal-repeat scan could not see, or ''.
    Appended where a terminal repeat is reported absent: absence is only meaningful within the scanned windows."""
    if not seq:
        return ""
    sc = structural_mod.scan_scope(seq)
    n = sc["record_bp"]
    if kind == "LTR":
        wide = narrow = sc["ltr_seed_window_bp"]
    else:
        wide, narrow = sc["tir_seed_window_bp"], sc["tir_joint_window_bp"]
    if 2 * wide >= n and 2 * narrow >= n:
        return ""
    if 2 * wide >= n:
        return (f" (the record is {n:,} bp and short or imperfect terminal repeats were searched only within "
                f"{narrow:,} bp of each end)")
    return f" (the record is {n:,} bp and terminal repeats were searched only within {wide:,} bp of each end)"


def _classify_without_structure(domains, domains_ok=True, orfs_unscanned=0):
    """Classify only positive coding evidence when structural detection did not run."""
    codes = [d.get("domain") for d in domains if d.get("domain")]
    cset = set(codes)
    gag_core = any(d.get("hmm") in _GAG_CORE for d in domains)
    ledger = [{"dependency": "structural", "status": "unavailable",
               "statement": "the structural scan did not run; repeats, tails and target-site duplications are unassessed"}]
    outcome = "complete"
    klass, superfamily, te_class, order, completeness = "unclassified", None, "unclassified", None, None
    if not domains_ok:
        outcome = "withheld"
        superfamily = "classification withheld (structural and domain scans unavailable)"
        ledger.append({"dependency": "domains", "status": "unavailable",
                       "statement": "the protein-domain scan did not run; coding evidence is unassessed"})
    elif not cset:
        outcome = "withheld"
        superfamily = "classification withheld (structural scan unavailable)"
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "no tested protein domain was detected; structural evidence is unavailable"})
    elif "INT" not in cset and _same_strand(domains, "RT", "YR"):
        klass, superfamily, te_class = ("Class I · retrotransposon",
                                        "DIRS-group (coding architecture; structure unassessed)", "DIRS")
        order = "–".join(_ordered(cset))
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "reverse transcriptase plus tyrosine recombinase supports DIRS-group coding architecture"})
    elif "RT" in cset and "INT" in cset:
        klass, superfamily, te_class = ("Class I · retrotransposon",
                                        "LTR retroelement (coding architecture; termini unassessed)", "LTR/unclassified")
        order = "–".join(_ordered(cset))
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "reverse transcriptase plus DDE integrase supports the LTR order from coding architecture"})
    elif "TPase" in cset:
        _fam, te_class, _conflict = _tpase_family([d for d in domains if d.get("domain") == "TPase"])
        klass = "Class II · DNA transposon"
        superfamily = f"{_fam} (advisory: strongest transposase profile; termini unassessed)"
        order = "–".join(_ordered(cset))
        completeness = {"tier": "coding evidence present (transposase); terminal structure unassessed",
                        "kind": "DNA transposon (terminal structure unassessed)", "present": ["TPase"],
                        "missing": [], "expected": ["TPase", "TIR (5′)", "TIR (3′)"],
                        "unassessed": ["TIR (5′)", "TIR (3′)"], "scope": DOMAINS_TESTED}
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "a transposase domain supports a Class II coding assignment; termini are unassessed"})
        if _conflict:
            ledger.append({"dependency": "domains", "status": "complete",
                           "statement": "multiple transposase classes detected — family assigned from the strongest-scoring hit (ambiguous)"})
    elif "HEL" in cset:
        superfamily = "Helitron-family helicase present — superfamily not called"
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "a Helitron-family helicase was detected; the untested terminal architecture prevents a superfamily call"})
    elif "RT" in cset and (gag_core or cset & (_RETRO_CORROBORATION | _LINE_CORROBORATION)):
        klass, superfamily, te_class = ("Class I · retrotransposon",
                                        "retrotransposon order unassigned (structure unavailable)", "retro/unassigned")
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "reverse transcriptase plus corroborating coding modules supports a retroelement; order is unassigned"})
    elif cset == {"RT"}:
        outcome = "withheld"
        superfamily = "reverse transcriptase, class unassigned"
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "reverse transcriptase alone cannot distinguish a transposable element from another reverse-transcribing gene"})
    else:
        superfamily = f"{'–'.join(_ordered(cset))} domain(s) detected; transposable-element class unassigned"
        ledger.append({"dependency": "domains", "status": "complete",
                       "statement": "coding domains were detected but do not establish a transposable-element class"})
    evidence = [item["statement"] for item in ledger]
    prefix = (f"Classification withheld; reported as {te_class} (candidate confidence). "
              if outcome == "withheld" else f"Classified as {te_class} (candidate confidence). ")
    explanation = (prefix + "Structural: unavailable. "
                   f"Domains: {' · '.join(codes) or 'none'}. " + "; ".join(evidence) + ".")
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "domain_architecture": ("unavailable" if not domains_ok else "–".join(_ordered(cset)) or "none detected"),
            "wicker_order": _wicker_order(te_class, outcome),
            "confidence": "Candidate", "evidence": evidence, "evidenceLedger": ledger,
            "outcome": outcome, "explanation": explanation, "n_domains": len(cset), "is_erv": False,
            "completeness": completeness, "domains_unavailable": not domains_ok,
            "structural_unavailable": True, "orfs_unscanned": int(orfs_unscanned or 0)}


def classify(structural, domains, seq=None, domains_ok=True, orfs_unscanned=0, orfs=None,
             structural_ok=True):
    """`domains_ok=False` means the protein-domain scan RAISED, so `domains` is empty because nothing ran —
    not because nothing is there. Without that distinction every domain-derived negative ("not detected:
    GAG, PR, RNaseH") is asserted from a scan that never happened, which is absence-of-evidence sold as
    evidence-of-absence. ``structural_ok=False`` applies the same rule to structural negatives: an empty
    list then means the detector failed, not that terminal evidence was absent."""
    if not structural_ok:
        return _classify_without_structure(domains, domains_ok=domains_ok, orfs_unscanned=orfs_unscanned)

    # Advisory rows are REPORTED evidence, never CREDITED evidence. A sub-threshold terminal repeat is
    # shown to the user precisely because it did not qualify, so it must not satisfy any has_* test —
    # otherwise the classifier treats a rejection as a confirmation. Belt and braces with the naming rule
    # in structural.py: a future advisory type that happens to start with a credited prefix is caught here.
    credited = [e for e in structural if not e.get("advisory")]
    has_ltr = any(e["type"].startswith("LTR") for e in credited)
    tir_ev = next((e for e in credited if e["type"].startswith("TIR")), None)
    has_tir = tir_ev is not None
    has_tsd = any(e["type"].startswith("TSD") for e in credited)
    codes = [d["domain"] for d in domains]
    cset = set(codes)
    gag_core = any(d.get("hmm") in _GAG_CORE for d in domains)   # capsid/matrix, not nucleocapsid zf-CCHC alone
    rt_d, int_d, tp_d = _rep(domains, "RT"), _rep(domains, "INT"), _rep(domains, "TPase")
    # ORF2p is a single reading frame with the endonuclease N-terminal to the RT. Pfam's Exo_endo_phos
    # also matches host DNase I and AP endonucleases, so an EN hit only counts as element evidence when
    # it holds that arrangement — same strand, upstream of the RT in translation order.
    en_d = _rep(domains, "EN")
    en_rt_pairs = [
        (en_hit, rt_hit)
        for en_hit in domains if en_hit.get("domain") == "EN"
        for rt_hit in domains if rt_hit.get("domain") == "RT"
        if en_hit.get("orf") is not None
        and en_hit.get("orf") == rt_hit.get("orf")
        and en_hit.get("strand", "+") == rt_hit.get("strand", "+")
        and _pos(en_hit) < _pos(rt_hit)
    ]
    en_ok = bool(en_rt_pairs)
    line_corroborated = en_ok or bool(cset & _LINE_CORROBORATION)
    retro_corroborated = line_corroborated or gag_core or bool(cset & _RETRO_CORROBORATION)
    rt, intg = rt_d is not None, int_d is not None
    # a 3' poly-A tail belongs to a plus-strand ORF2 and a 5' poly-T to a minus-strand one; the opposite
    # tract is an unrelated T/A run (or the other strand's copy) and is reported but not credited. The
    # advisory polyA-signal motif row is neither: tail types are matched exactly.
    _tails = {e["type"].split(" ")[0] for e in credited if e["type"].startswith(("poly-A", "poly-T"))}
    has_polya = bool(_tails)
    tail_ok = ("poly-A" if (rt_d or {}).get("strand", "+") == "+" else "poly-T") in _tails
    tpase = "TPase" in cset
    # the element's own ends: a TIR pair that actually encloses the transposase. The scan resolves the two
    # arms together, so they stand or fall as a pair; an inverted repeat elsewhere in the record is not
    # evidence that THIS coding module has its termini.
    tir_encloses_tpase = _brackets(tir_ev, tp_d)
    tir_ok = has_tir and tir_encloses_tpase is not False
    ev, superfamily, te_class, order = [], None, None, None
    order_resolvable = False
    tpase_conflict = False
    classification_withheld = False
    rt_class_withheld = False
    order_conflict = scan_unassessed = False

    # DIRS-group elements carry a tyrosine recombinase in place of the DDE integrase, so they must be
    # tested BEFORE the generic RT branch — otherwise every one of them falls through to "LINE (non-LTR)"
    # on the strength of an absent integrase. Wicker 2007 order DIRS; PASTEC's YR-vs-DDE logic
    # (Hoede 2014 PLoS ONE 9:e91929) makes the same distinction.
    yr = rt and _same_strand(domains, "RT", "YR")          # a recombinase on the other strand is not the DIRS architecture
    hel = "HEL" in cset          # Helitron Rep/Helicase (PF14214) — bundled and advertised, so it must be consumed
    if yr and not intg:
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
            r_nt = rt_d.get("nt") or [0, 0]
            # every integrase hit that can be read against the RT (same strand, non-overlapping). Several can
            # exist in a nested or adjacent locus, so the call must not follow the highest bit score alone:
            # an integrase in the RT's own ORF is preferred, and integrases on both sides of the RT with no
            # such tie-break are a conflict.
            ints = [d for d in domains if d["domain"] == "INT" and d.get("strand", "+") == rt_d.get("strand", "+")
                    and ((d.get("nt") or [0, 0])[1] <= r_nt[0] or r_nt[1] <= (d.get("nt") or [0, 0])[0])]
            own_orf = [d for d in ints if rt_d.get("orf") is not None and d.get("orf") == rt_d.get("orf")]
            ints = own_orf or ints
            sides = {_pos(d) < _pos(rt_d) for d in ints}
            order_resolvable = len(sides) == 1
            if len(sides) == 2:
                superfamily, te_class = "LTR retrotransposon (superfamily undetermined)", "LTR/unclassified"
                ev.append("integrase hits lie on both sides of the RT on its strand (possible nested or adjacent "
                          "elements), so the integrase-vs-RT order is conflicting — Copia/Gypsy is not called")
            elif not order_resolvable:
                # Do NOT name a superfamily from a comparison that cannot carry the answer. _pos negates the
                # minus-strand coordinate to get translation order, which is coherent WITHIN one strand and
                # meaningless ACROSS two: a minus-strand hit always yields a negative value, so the test
                # reduces to "which of the two is on the minus strand" and returns the SAME answer whether
                # INT sits upstream or downstream of RT (measured: INT-on-minus -> always Copia,
                # RT-on-minus -> always Gypsy, at every coordinate). Naming a tentative superfamily from that
                # puts a family in the card's headline on evidence containing no information about it, which
                # would exceed the evidence available. The overlapping-span case is equally unreadable. Fall through to
                # the same honest label the no-integrase path already uses.
                superfamily, te_class = "LTR retrotransposon (superfamily undetermined)", "LTR/unclassified"
                ev.append("integrase and RT are not in a readable arrangement (different strands, or "
                          "overlapping spans), so the diagnostic integrase-vs-RT order cannot be determined — "
                          "Copia/Gypsy is not called")
            elif sides == {True}:
                superfamily, te_class = "Copia (Ty1)", "LTR/Copia"
                ev.append("integrase N-terminal to RT (INT before RT, the Copia/Ty1-type domain order) + paired LTRs → "
                          "Copia superfamily")
            elif _RETROVIRAL_PROFILES & {d.get("hmm") for d in domains}:
                # integrase C-terminal to RT is the pol order of Gypsy, Bel-Pao, Retrovirus and ERV alike
                superfamily, te_class = RETROVIRAL_UNDETERMINED, "LTR/unclassified"
                ev.append("integrase C-terminal to RT + paired LTRs fits the Gypsy-type pol order, but a "
                          "retrovirus-specific capsid/matrix or envelope profile is present: that order is "
                          "shared by Gypsy, Bel-Pao and Retrovirus/ERV, so the superfamily is not called")
            else:
                superfamily = "Gypsy (Ty3)"
                if "CHR" in cset:
                    superfamily += " · chromovirus"
                te_class = "LTR/Gypsy"
                ev.append("integrase C-terminal to RT (RT before INT, the Gypsy/Ty3-type domain order) + paired LTRs → "
                          "Gypsy superfamily")
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
                _lcong = structural_mod.tsd_congruence((_ltsd or {}).get("length"), superfamily, motif=(_ltsd or {}).get("motif"))
                if _ltsd is not None:
                    _ltsd["tsd_congruence"] = _lcong["verdict"]
                if _lcong["verdict"] == "incongruent":
                    ev.append(f"a flanking target-site duplication is present but {_tsd_observed(_lcong)} where "
                              f"{superfamily.split(' ')[0]} elements duplicate {_tsd_expectation(_lcong)} "
                              f"({_lcong['basis']}) — the repeat may be coincidental, or the boundary or superfamily "
                              f"call may be wrong, so the ends are not credited from it")
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
                      "accepted" + _scan_gap(seq, "LTR"))
            ev.append("the superfamily is not called: Copia versus Gypsy needs the integrase-versus-RT "
                      "translation order, and no integrase was detected")
        elif not has_ltr and not intg and not retro_corroborated:
            # A LONE reverse transcriptase establishes neither the order NOR the class.
            #
            # An earlier form of this guard withheld only the order, and admitted a poly-A tail as evidence
            # that the record was LINE-like. Both were wrong, and the same record exposed both. Telomerase
            # reverse transcriptase is a host enzyme, not a transposable element, and the human catalytic
            # subunit is deposited as an mRNA — so it carries a poly-A tail that belongs to the transcript,
            # not to target-primed reverse transcription. The tail defeated the guard, the record fell
            # through to the LINE branch, and a telomerase gene was called a non-LTR retrotransposon.
            #
            # A poly-A tail cannot license this call at all. Every mRNA and every cDNA deposit carries one,
            # so it separates nothing; it is retained as reported evidence once a call is licensed on coding
            # grounds, and never as the thing that licenses it.
            #
            # The class is withheld for the same reason the order is. A reverse transcriptase is carried by
            # every retrotransposon order, by telomerase, and by group II intron maturases — the panel holds
            # one RT model, not a clade-resolved set, so on this evidence a Class I element cannot be
            # separated from a host gene that reverse-transcribes. Asserting Class I here is the error that
            # a negative-control panel exists to catch.
            klass = "unclassified"
            superfamily, te_class = "reverse transcriptase, class unassigned", "unclassified"
            rt_class_withheld = True
            ev.append("a reverse transcriptase was detected with no integrase, no terminal repeat and no "
                      "credited coding architecture that distinguishes a retrotransposon from another "
                      "reverse-transcribing gene — no ORF1p, no ORF2p endonuclease–RT arrangement, no "
                      "RNase H, no core capsid/matrix gag, no protease, and no envelope. Telomerase and "
                      "group II intron maturases also present as a "
                      "lone reverse transcriptase, so neither the class nor the order is assigned")
            if has_polya:
                ev.append("a poly-A tail is present but is not credited: every mRNA and cDNA deposit carries "
                          "one, so it does not distinguish a non-LTR element from a transcript")
        elif not has_ltr and not intg and not line_corroborated:
            superfamily, te_class = "retrotransposon order unassigned", "retro/unassigned"
            observed = (["core GAG"] if gag_core else []) + [
                label for code, label in (("PR", "protease"), ("RNaseH", "RNase H"), ("ENV", "envelope"))
                if code in cset
            ]
            ev.append("reverse transcriptase plus " + ", ".join(observed) +
                      " supports a retroelement, but this evidence does not establish the non-LTR order; "
                      "no LINE-specific ORF1p or accepted ORF2p endonuclease–RT architecture was detected")
        elif not has_ltr and not intg:
            superfamily, te_class = "LINE (non-LTR)", "LINE"
            ev.append("RT without a DDE integrase and without LTRs → non-LTR retrotransposon (LINE)")
            if tail_ok:                                       # name the tail that was actually detected, not always poly-A
                ev.append(("3′ poly-A tail" if "poly-A" in _tails else "5′ poly-T tract") +
                          " consistent with target-primed reverse transcription")
            elif has_polya:
                ev.append(("a 5′ poly-T tract" if "poly-T" in _tails else "a 3′ poly-A tail") +
                          f" is present but the RT is on the {rt_d.get('strand', '+')} strand, where a tail of that "
                          "kind does not belong, so it is not credited as this element's tail")
            if "ORF1" in cset:
                ev.append("ORF1p domain present — the LINE-specific coding module, not shared with LTR elements")
            if en_ok:
                ev.append("apurinic/apyrimidinic endonuclease-like (APE) domain N-terminal to the RT in the same called ORF → "
                          "the ORF2p EN–RT architecture of a LINE")
            elif "EN" in cset:
                ev.append("an endonuclease domain was detected but does not sit N-terminal to an RT in the "
                          "same called ORF; Pfam's endonuclease family also covers host enzymes, so it is not "
                          "credited as ORF2p here")
            # A DIRS-group element is now tested for directly (above) rather than merely hedged against.
            # Penelope-like elements remain outside the panel: their GIY-YIG endonuclease is not modelled.
            ev.append("a Penelope-like element would also show RT without a DDE integrase; its GIY-YIG "
                      "endonuclease is not in the tested panel, so it is not excluded here")
        elif has_ltr and (en_ok or "ORF1" in cset):
            # LINE-specific coding modules and terminal direct repeats disagree (e.g. a LINE inside a
            # segmental or tandem duplication). Neither outranks the other, so the order is withheld.
            superfamily, te_class = "retrotransposon order unassigned", "retro/unassigned"
            order_conflict = True
            ev.append("LINE-specific coding evidence (" +
                      " and ".join(x for x, ok in (("ORF1p", "ORF1" in cset), ("an ORF2p endonuclease–RT arrangement", en_ok)) if ok) +
                      ") conflicts with paired terminal direct repeats, which suggest an LTR element, and no "
                      "integrase was detected, so the order is withheld")
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
                      "repeat was accepted" + _scan_gap(seq, "LTR") + ", so the assignment does not rest on structural evidence")
            if "ENV" in cset:                     # has_env is not bound until later in this function
                ev.append("envelope glycoprotein present — carried by LTR elements and never by non-LTR "
                          "retrotransposons, corroborating the order")
    elif tpase:
        klass = "Class II · DNA transposon"
        tp_hits = [d for d in domains if d["domain"] == "TPase"]
        superfamily, te_class, tpase_conflict = _tpase_family(tp_hits)
        ev.append("transposase domain present → Class II DNA transposon")
        if tpase_conflict:
            ev.append("multiple transposase classes detected — superfamily assigned from the strongest-scoring hit (ambiguous)")
        elif superfamily.startswith("Tc1/Mariner") and _recover_target_site(structural, tir_ev, superfamily, seq):
            has_tsd = True                                   # the trimmed arms may enclose the transposase differently
            tir_encloses_tpase = _brackets(tir_ev, tp_d)
            tir_ok = tir_encloses_tpase is not False
        if has_tir:
            ev.append("terminal inverted repeats consistent with a cut-and-paste transposon")
        # The ends are half of what makes a DNA transposon autonomous, so state what was actually recovered.
        # The "insertion site captured" claim is a completeness statement built on a short exact direct
        # repeat, which arises by chance at ~4^-L — so it is GATED on the TSD length not contradicting the
        # superfamily's expected one, and the incongruent case is reported instead of suppressed.
        if tir_ok and has_tsd:
            _tsd = next((e for e in structural if e["type"].startswith("TSD")), None)
            _refine_tsd(_tsd, superfamily, seq)               # prefer the superfamily's target-site length over a coincidental longer flank
            _cong = structural_mod.tsd_congruence((_tsd or {}).get("length"), superfamily, motif=(_tsd or {}).get("motif"))
            if _tsd is not None:
                _tsd["tsd_congruence"] = _cong["verdict"]     # carry the verdict to the record (GFF3 export reads it)
            _chance = _tsd_chance_clause(_tsd)
            _fit = (f"length congruent with {superfamily}, {_cong['basis']}" if _cong["verdict"] == "congruent" else
                    f"no literature-attributed expected target-site length exists for {superfamily} to corroborate its length")
            if _cong["verdict"] == "incongruent":
                ev.append(f"a target-site duplication flanks the inverted repeats, but it is "
                          f"{_tsd_observed(_cong)} where {superfamily} elements duplicate "
                          f"{_tsd_expectation(_cong)} ({_cong['basis']}) — the repeat may be coincidental, or the "
                          f"boundary or superfamily call may be wrong, so the ends are not credited from it")
            elif _chance:
                # a short exact direct repeat arises by chance at ~4^-L, so neither a congruent length nor its
                # absence turns it into a captured insertion site
                ev.append(f"a target-site duplication flanks the inverted repeats ({_fit}), but {_chance}, so it "
                          "is not credited as the insertion site and does not confirm the element ends")
            else:
                ev.append("a target-site duplication flanks the inverted repeats — the insertion site itself is "
                          f"captured, so both element termini are present in the record ({_fit})")
        elif has_tir and tir_encloses_tpase is False:
            ev.append("the detected inverted-repeat pair does not enclose the transposase — the termini and the "
                      "coding module may not belong to the same element, so the ends are not credited")
        elif not has_tir:
            ev.append("no terminal inverted-repeat pair recovered" + _scan_gap(seq, "TIR") +
                      " — neither element end is confirmed, which is "
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
                  "replication, so they carry no terminal inverted repeat and duplicate no target site" +
                  (" — the absence of TIR/TSD here is expected and is not evidence against a Helitron"
                   if not has_tir else
                   "; a terminal inverted repeat was detected in this record, which is not a Helitron feature "
                   "(it may belong to another element or be coincidental)"))
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
            # Direct terminal repeats are not, by themselves, specific to an LTR retroelement. They can also
            # be a tandem/satellite pair, and an isolated domain does not resolve that ambiguity. Preserve the
            # empirically supported H0080-like fragment boundary: a capsid/matrix GAG plus another retroviral
            # module is a coherent coding architecture even when RT/pol was not recovered.
            _trusted_ltr_fragment = has_ltr and gag_core and bool(cset & {"PR", "INT", "ENV"})
            if _trusted_ltr_fragment:
                superfamily, te_class = f"LTR retroelement fragment ({_cs}; RT/pol not detected)", "LTR/partial"
                ev.append(f"coding domain(s) recovered ({_cs}) with paired terminal repeats, but neither "
                          f"reverse transcriptase nor transposase detected")
            elif has_ltr:
                superfamily = "terminal direct-repeat candidate (LTR origin unconfirmed)"
                te_class = "repeat/structural-only"
                classification_withheld = True
                ev.append(f"paired terminal direct repeats and coding domain(s) ({_cs}) were detected, but this "
                          "does not establish an LTR retroelement: the detected domains do not form the trusted "
                          "capsid-plus-retroviral-module architecture and no reverse transcriptase was detected")
            elif cset & _RETRO_SPECIFIC:
                superfamily, te_class = f"coding fragment ({_cs}; RT/transposase not detected)", "retro/partial"
                ev.append(f"coding domain(s) recovered ({_cs}) but neither reverse transcriptase nor "
                          f"transposase detected")
            elif has_tir and (tir_ev.get("tir_len") or 0) >= structural_mod.MIN_TIR_STANDALONE:
                klass = "Class II · DNA transposon"
                superfamily = "non-autonomous TIR element (superfamily undetermined)"
                te_class = "DNA/non-autonomous"
                ev.append(f"a {tir_ev.get('tir_len')} bp terminal inverted repeat was detected, beyond what chance "
                          f"inverted repeats reach ({CHANCE_TIR_MAX_BP} bp across {CHANCE_TIR_SEQUENCES:,} random "
                          f"sequences), so it is structural evidence of a Class II element; the only coding "
                          f"evidence is {_cs}, which host genes also carry and which may be a captured gene "
                          f"fragment (as in Pack-MULEs), so no transposase is credited and the superfamily is "
                          f"not called")
            else:
                klass = "unclassified"
                superfamily = f"{_cs} domain(s) detected; no transposable-element assignment"
                te_class = "unclassified"
                ev.append(f"the only coding evidence is {_cs}, which host genes also carry — the "
                          f"APE-family endonuclease model matches host AP endonucleases and DNase I, and a tyrosine "
                          f"recombinase matches phage, plasmid and host site-specific recombinases. Without "
                          f"a reverse transcriptase, a transposase, a capsid, an envelope or an integrase "
                          f"there is no element-specific evidence, so no class is assigned")
                if has_tir:
                    ev.append(f"an inverted terminal repeat ({tir_ev.get('tir_len') or 0} bp) was also detected but "
                              f"is too short to attribute to a Class II element")
        elif has_ltr:
            superfamily, te_class = "terminal direct-repeat candidate (LTR origin unconfirmed)", "repeat/structural-only"
            classification_withheld = True
            if domains_ok:
                ev.append("paired terminal direct repeats were detected, but no trusted LTR coding architecture "
                          "was recovered; direct repeats alone do not distinguish an LTR retroelement from a "
                          "tandem or satellite repeat, so the order is withheld")
            else:
                ev.append("paired terminal direct repeats were detected, but the protein-domain scan was "
                          "unavailable; direct repeats alone do not distinguish an LTR retroelement from a "
                          "tandem or satellite repeat, so the order is withheld")
        elif not domains_ok:
            # Every domain-absence the TIR and no-signal arms below would assert (no transposase, no coding
            # capacity) comes from a scan that never ran, so nothing is concluded from structure alone.
            klass = "unclassified"
            classification_withheld = scan_unassessed = True
            if has_tir:
                superfamily, te_class = "terminal inverted repeat, class unassigned", "repeat/structural-only"
                ev.append(f"an inverted terminal repeat was detected ({(tir_ev or {}).get('tir_len') or 0} bp), but "
                          "the protein-domain scan was unavailable, so whether a transposase is present is "
                          "unassessed and no class is assigned")
            else:
                superfamily, te_class = "classification withheld (protein-domain scan unavailable)", "unclassified"
                ev.append("the protein-domain scan was unavailable, so the absence of a transposable-element "
                          "signature cannot be concluded from the structural evidence alone")
        elif has_tir:
            # A SHORT inverted terminal repeat without a transposase is not Class II evidence. Verified on
            # the real copia 5' LTR (276 bp): fed alone it trips the terminal-inverted-repeat scan and used
            # to be filed as a DNA transposon — a Class I fragment assigned to Class II. Solo LTRs
            # outnumber full-length elements in most genomes, so this is the common case, not an edge case.
            #
            # A LONG one is different, and the distinction is measured rather than assumed. Chance inverted
            # repeats reached at most 15 bp at the recorded seed (4,500 sequences) and 17 bp across five
            # seeds (22,500 sequences; benchmarks/chance_tir.py, see structural.MIN_TIR_STANDALONE),
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
                          f"inverted repeats reached at most {CHANCE_TIR_MAX_BP} bp across "
                          f"{CHANCE_TIR_SEQUENCES:,} random sequences (five seeds), so a repeat this "
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

    # envelope (env) domain + paired LTRs mark an env-bearing LTR retroelement (potentially infection-capable),
    # distinct from a plain LTR retrotransposon. "Errantivirus" is the insect Ty3/gypsy lineage, so it is named
    # only when the element is Gypsy-called; a vertebrate ERV is not one.
    has_env = "ENV" in cset
    # a DIRS-group call (RT + tyrosine recombinase, no DDE integrase) is never an ERV, even with a spurious
    # env hit — routing it through the retroviral transcript-architecture model would contradict its own class.
    is_erv = bool(rt and has_ltr and has_env and not (yr and not intg))
    if is_erv:
        ev.append("envelope (env) domain present with paired LTRs → env-bearing LTR retroelement (ERV architecture" +
                  ("; errantivirus lineage for a Gypsy element)" if (superfamily or "").startswith("Gypsy") else ")"))

    # domain-architecture order shows ONLY the domains actually detected, in the superfamily's
    # canonical order — never a full template presented as if it were observed. ENV closes the retroviral order.
    # No CHR in the Copia canon: a chromodomain-bearing integrase is a chromovirus feature and
    # chromoviruses are a Ty3/Gypsy lineage, so listing it here would render a spurious PF00385 hit on a
    # Copia call as though it were canonical Copia architecture.
    _CANON = {"Copia": ["GAG", "PR", "INT", "RT", "RNaseH", "ENV"],
              "Gypsy": ["GAG", "PR", "RT", "RNaseH", "INT", "CHR", "ENV"]}
    _sf = superfamily.split(" ")[0] if superfamily else ""
    if superfamily == RETROVIRAL_UNDETERMINED:                  # the pol order is Gypsy-type, only the superfamily is withheld
        _sf = "Gypsy"
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
        confidence = "Moderate" if tail_ok else "Candidate"
    # The transposase arm gates on tir_ok, not raw has_tir. This function already refuses to credit a
    # terminal repeat that does not enclose the transposase — the ledger drops the arms and the evidence
    # says "the ends are not credited" — because the termini and the coding module may belong to different
    # elements. Reading that same repeat as enough to lift the badge out of Candidate contradicted the
    # refusal one screen above it: two records with an identical completeness ledger got different badges.
    # A solo LTR fed alone trips the inverted-repeat scan (see the note below), so this is the common case.
    elif (rt and (has_ltr or has_tir)) or (tpase and tir_ok):
        confidence = "Moderate"
    elif rt or tpase:
        confidence = "Moderate" if ndom >= 2 else "Candidate"
    else:
        confidence = "Candidate"
    if tpase_conflict or order_conflict or rt_class_withheld or classification_withheld:
        confidence = "Candidate"                   # conflicting or withheld evidence never earns a higher badge than the call it refuses
    if composite and confidence == "High":         # a co-occurring transposase makes a confident single-element call unsafe
        confidence = "Moderate"

    completeness = _completeness(cset, rt, intg, tpase, has_ltr, has_tir, tail_ok, is_erv, order_resolvable,
                                 gag_core, tir_ok, en_ok=en_ok, yr=yr, te_class=te_class)
    if scan_unassessed:
        completeness = ({"tier": "structural-only (protein-domain scan unavailable; coding content unassessed)",
                         "kind": "terminal inverted repeat (class unassigned)", "present": [], "missing": [],
                         "expected": [], "unassessed": ["protein domains"], "scope": DOMAINS_TESTED}
                        if has_tir else None)
    elif classification_withheld:
        completeness = {
            "tier": "structural candidate (no trusted LTR coding architecture)",
            "kind": "terminal direct-repeat candidate (LTR origin unconfirmed)",
            "present": _ordered(cset), "missing": [], "expected": [], "scope": DOMAINS_TESTED,
        }
        if not domains_ok:
            completeness["unassessed"] = ["protein domains"]
    # The coding axis rides ALONGSIDE the tier, never inside it. The tier value is benchmarked (moving it
    # would move published rows), and the defect was never that the tier miscounted domains — it was that a
    # domain-only ledger was the ONLY thing on screen, so "partial" was read as "decayed". Adding the axis
    # answers the question the reader was actually asking without perturbing a benchmarked number.
    coding = _coding_axis(domains, orfs)
    if coding and completeness:
        completeness["coding"] = coding

    dom_str = (" · ".join(codes) or "none detected in tested panel") if domains_ok else "unavailable (scan failed)"
    struct_str = ", ".join(e["type"].split(" (")[0] for e in structural) or "none"
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
    prefix = (f"Classification withheld; reported as {te_class} ({confidence.lower()} confidence). "
              if classification_withheld or rt_class_withheld
              else f"Classified as {te_class} ({confidence.lower()} confidence). ")
    explanation = (prefix + f"Structural: {struct_str}. Domains: {dom_str}." +
                   (" " + "; ".join(ev) + "." if ev else ""))
    if classification_withheld and not scan_unassessed:
        evidence_ledger = [{"dependency": "structural", "status": "complete",
                            "statement": "paired terminal direct repeats were detected"}]
        if domains_ok:
            evidence_ledger.append({"dependency": "domains", "status": "complete",
                                    "statement": (f"detected coding domains: {'–'.join(_ordered(cset))}"
                                                  if cset else "no coding domain was detected in the tested panel")})
        else:
            evidence_ledger.append({"dependency": "domains", "status": "unavailable",
                                    "statement": "the protein-domain scan did not run; coding evidence is unassessed"})
    elif not domains_ok:
        evidence_ledger = [
            {"dependency": "domains", "status": "unavailable",
             "statement": "the protein-domain scan did not run; coding evidence is unassessed"},
            *[{"dependency": "structural", "status": "complete", "statement": item}
              for item in ev if not item.startswith("the protein-domain scan did not run")],
        ]
    else:
        evidence_ledger = [{"dependency": "combined", "status": "complete", "statement": item} for item in ev]
    outcome = ("withheld" if classification_withheld or rt_class_withheld or (not domains_ok and not credited)
               else "complete")
    if not domains_ok:
        domain_architecture = "unavailable"
    elif not cset:
        domain_architecture = "none detected"
    else:
        domain_architecture = order if (_sf in _CANON and order) else "–".join(_ordered(cset))
    return {"class": klass, "superfamily": superfamily, "te_class": te_class, "order": order,
            "domain_architecture": domain_architecture, "wicker_order": _wicker_order(te_class, outcome),
            "confidence": confidence, "evidence": ev,
            "evidenceLedger": evidence_ledger,
            "outcome": outcome,
            "explanation": explanation, "n_domains": ndom,
            "is_erv": is_erv, "completeness": completeness, "domains_unavailable": not domains_ok,
            "structural_unavailable": not structural_ok, "orfs_unscanned": int(orfs_unscanned or 0)}


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
                  "chromodomain, LINE ORF1p and apurinic/apyrimidinic endonuclease-like (APE) domain, tyrosine recombinase, "
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
        # n counts host ORFs; a sequence has only six reading frames, so the frames are listed, not counted
        frames = sorted({f"{k[2]}{k[3]}" for k in placed if k[3] is not None})
        note = (f"the detected domains fall in {n} separate ORFs" + (f" (frames {', '.join(frames)})" if frames else "")
                + (f", and {unplaced} could not be placed in any single ORF" if unplaced else ""))
        line_layout = False
        if n == 2 and not unplaced:
            (ka, ca), (kb, cb) = sorted(placed.items(), key=lambda kv: kv[0][0])
            if ka[2] == kb[2]:                            # ORF1 upstream of the endonuclease-RT ORF in translation order
                up, down = (ca, cb) if ka[2] == "+" else (cb, ca)
                line_layout = set(up) == {"ORF1"} and "RT" in down and "ORF1" not in down
        note += (" — consistent with the canonical LINE ORF1/ORF2 organisation (ORF1p, then the endonuclease–RT ORF)"
                 if line_layout else
                 " — consistent with an interrupted, rearranged or degenerate copy, though a genuine "
                 "multi-ORF architecture (LINE ORF1/ORF2, some ERVs) also looks like this")
    if truncated:
        note += "; at least one hosting ORF runs to the end of the supplied sequence, so it may be truncated"
    return {"frames": n, "unplaced": unplaced, "truncated_frame": truncated, "state": state, "note": note}


def _completeness(cset, rt, intg, tpase, has_ltr, has_tir, has_polya, is_erv, order_resolvable, gag_core=False,
                  tir_ok=False, en_ok=False, yr=False, te_class=None):
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
    if rt and yr and not intg:                            # DIRS-group: tyrosine recombinase in place of DDE integrase (yr is the co-oriented RT+YR test)
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
    elif rt and (has_ltr or intg or str(te_class or "").startswith("LTR")) and te_class != "retro/unassigned":
        # LTR retroelement / ERV; an LTR-order call that rests on RT + chromodomain has no LTRs or integrase
        # yet still gets this ledger, not the LINE one
        expected = ["GAG", "PR", "RT", "RNaseH", "INT"] + [c for c in ("ENV", "CHR") if c in cset]
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
        elif rt and not has_ltr and not intg:
            tier = "partial (RT + chromodomain; integrase not confirmed)"
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
    elif rt and te_class == "unclassified":
        expected = []
        present = _ordered(cset)
        missing = []
        tier = "reverse-transcriptase evidence only; class unassigned"
        kind = "reverse transcriptase, class unassigned"
    elif rt and te_class == "retro/unassigned":
        expected = []
        present = _ordered(cset)
        missing = []
        tier = "coding evidence present; order unassigned"
        kind = "retroelement (order unassigned)"
    elif rt:                                              # non-LTR (LINE-like)
        # An autonomous LINE carries ORF1 (nucleic-acid chaperone) plus ORF2 (APE-family endonuclease
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
        # the kind follows the call: a short TIR, or a TIR beside only non-transposase domains, is not Class II
        kind = ("LTR retroelement" if has_ltr
                else "DNA transposon" if str(te_class or "").startswith("DNA/")
                else "retroelement fragment (order unassigned)" if str(te_class or "").startswith("retro")
                else "terminal inverted repeat (class unassigned)")
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
