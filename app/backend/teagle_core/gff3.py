"""GFF3 and BED export of a TEagle annotation.

TEagle is otherwise a terminal node: a result that cannot leave the application undercuts the one thing
it does that a hosted search does not — dissect a candidate, design the confirming assay, and seal the
run. This is the interoperability half of that.

Column 3 of GFF3 is an ASSERTION that downstream software trusts. It is therefore gated on the same
completeness tier the report shows: a structural-only call exports as the generic `repeat_region`, never
as `LTR_retrotransposon`, because a specific subclass claims coding evidence that a structural-only call
does not have. Wicker codes have no SO accession and go in LOWERCASE attributes — GFF3 reserves
capitalised tags for its own predefined set, so the `TEID=` convention seen in the wild is non-compliant.

Term names and accessions were verified on 2026-07-28 against the Sequence Ontology's own
`so-simple.obo` (The-Sequence-Ontology/SO-Ontologies, CC-BY-4.0). EDTA's mapping file is GPL-3.0 and is
deliberately not used. Format: GFF3 specification; Diesh et al. 2023 Genome Biol 24:74.
"""
from __future__ import annotations

# name -> accession. Verified against so-simple.obo; the accession travels in Ontology_term= so a
# consumer can resolve the exact class rather than string-matching column 3.
SO = {
    "repeat_region": "SO:0000657",
    "transposable_element": "SO:0000101",
    "LTR_retrotransposon": "SO:0000186",
    "non_LTR_retrotransposon": "SO:0000189",
    "LINE_element": "SO:0000194",
    "SINE_element": "SO:0000206",
    "YR_retrotransposon": "SO:0002286",
    "terminal_inverted_repeat_element": "SO:0000208",
    "DNA_transposon": "SO:0000182",
    "helitron": "SO:0000544",
    "MITE": "SO:0000338",
    "long_terminal_repeat": "SO:0000286",
    "terminal_inverted_repeat": "SO:0000481",
    "target_site_duplication": "SO:0000434",
    "primer_binding_site": "SO:0005850",
    "RR_tract": "SO:0000435",
    "ORF": "SO:0000236",
    "polypeptide_domain": "SO:0000417",
}

# te_class prefix -> the SO term the evidence supports when coding evidence IS present.
_CLASS_TERM = {
    "LTR/": "LTR_retrotransposon",
    "LINE": "LINE_element",
    "SINE": "SINE_element",
    "DIRS": "YR_retrotransposon",
    "DNA/": "terminal_inverted_repeat_element",
    "HEL": "helitron",
}

# structural evidence type -> SO term for the sub-feature rows
_STRUCT_TERM = {
    "LTR": "long_terminal_repeat",
    "TIR": "terminal_inverted_repeat",
    "TSD": "target_site_duplication",
    "PBS": "primer_binding_site",
    "PPT": "RR_tract",
}


def so_term_for(classification) -> str:
    """The SO term an annotation may claim, gated on the evidence its completeness ledger records.

    Column 3 is a machine-trusted assertion, so a subclass term is emitted only when the run actually
    established that subclass's defining evidence, and degrades to the nearest supported parent otherwise:
      - a 'structural-only' tier (terminal repeats, no coding domain) or an unresolved class -> repeat_region.
      - an LTR retrotransposon subclass claims reverse-transcriptase evidence; a fragment that kept a
        coding domain (e.g. env/capsid) but lost pol has no RT in the ledger -> repeat_region.
      - terminal_inverted_repeat_element claims the inverted repeats themselves; a transposase whose TIRs
        were not recovered -> the generic DNA_transposon, which asserts only DNA-mediated transposition
        (what the transposase supports), never the absent ends."""
    cl = classification or {}
    te_class = (cl.get("te_class") or "")
    comp = cl.get("completeness") or {}
    tier = (comp.get("tier") or "").lower()
    present = comp.get("present") or []
    missing = comp.get("missing") or []
    if not te_class or te_class == "none" or "structural-only" in tier or "structural-only" in te_class:
        return "repeat_region"
    for prefix, term in _CLASS_TERM.items():
        if te_class.startswith(prefix):
            if term == "terminal_inverted_repeat_element" and any("TIR" in m for m in missing):
                return "DNA_transposon"
            if term == "LTR_retrotransposon" and "RT" not in present:
                return "repeat_region"
            return term
    return "repeat_region"


def _esc(v) -> str:
    """GFF3 attribute escaping: the reserved characters must be percent-encoded, not stripped."""
    s = "" if v is None else str(v)
    for ch, rep in ((";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C"),
                    ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D")):
        s = s.replace(ch, rep)
    return s


def _attrs(pairs) -> str:
    return ";".join(f"{k}={_esc(v)}" for k, v in pairs if v not in (None, "", []))


def build_features(rec, seqid="locus"):
    """One annotation record -> a list of GFF3 feature tuples, parent first."""
    cl = rec.get("classification") or {}
    comp = cl.get("completeness") or {}
    length = (rec.get("composition") or {}).get("length") or len(rec.get("seq_preview") or "") or 1
    term = so_term_for(cl)
    rows = []
    gene_id = "TE_1"
    # parent: the element itself. Wicker code and the tier ride in lowercase custom attributes.
    rows.append((seqid, "TEagle", term, 1, length, ".", ".", ".", _attrs([
        ("ID", gene_id),
        ("Name", cl.get("superfamily") or term),
        ("Ontology_term", SO.get(term, "")),
        ("wicker_code", cl.get("te_class")),
        ("te_class", cl.get("te_class")),
        ("classification_confidence", cl.get("confidence")),
        ("completeness_tier", comp.get("tier")),
        ("domains_tested_scope", comp.get("scope")),
    ])))
    n = 0
    for ev in rec.get("structural") or []:
        kind = (ev.get("type") or "").split(" ")[0]
        sub = _STRUCT_TERM.get(kind)
        if not sub:
            continue
        spans = []
        if ev.get("five_prime"):
            spans.append(ev["five_prime"])
        if ev.get("three_prime"):
            spans.append(ev["three_prime"])
        if ev.get("pos"):                     # PBS, PPT: a single [start, end) span
            spans.append(ev["pos"])
        if ev.get("upstream"):                # TSD: two flanking direct repeats, both emitted as sub-features
            spans.append(ev["upstream"])
        if ev.get("downstream"):
            spans.append(ev["downstream"])
        if not spans and ev.get("start") is not None:
            spans.append([ev["start"], ev.get("end", ev["start"])])
        for span in spans:
            n += 1
            rows.append((seqid, "TEagle", sub, int(span[0]) + 1, int(span[1]), ".", ".", ".", _attrs([
                ("ID", f"{gene_id}.s{n}"), ("Parent", gene_id), ("Name", kind),
                ("Ontology_term", SO.get(sub, "")),
                # a length still limited by the record is declared here too, not silently exported
                ("length_is_lower_bound", "true" if ev.get("length_is_lower_bound") else None),
            ])))
    for i, d in enumerate(rec.get("domains") or [], 1):
        nt = d.get("nt") or [0, 0]
        rows.append((seqid, "TEagle", "polypeptide_domain", int(nt[0]) + 1, int(nt[1]),
                     f"{d.get('score', '.')}", d.get("strand", "."), ".", _attrs([
                         ("ID", f"{gene_id}.d{i}"), ("Parent", gene_id),
                         ("Name", d.get("label") or d.get("domain")),
                         ("Ontology_term", SO["polypeptide_domain"]),
                         ("domain_code", d.get("domain")), ("hmm_profile", d.get("hmm")),
                         ("pfam_accession", d.get("pfam")), ("evalue", d.get("evalue")),
                     ])))
    return rows


def to_gff3(rec, seqid="locus", sequence=None, source_note=None) -> str:
    """A complete, self-contained GFF3 document.

    Coordinates are LOCUS-RELATIVE unless the caller supplies a genome-anchored seqid, because exporting
    locus offsets under a chromosome name would silently misplace every feature. The sequence is embedded
    after ##FASTA so the file stands alone in a browser that has no matching reference."""
    rows = build_features(rec, seqid)
    length = (rec.get("composition") or {}).get("length") or 1
    out = ["##gff-version 3", f"##sequence-region {seqid} 1 {length}"]
    if source_note:
        out.append(f"#!source {source_note}")
    out.append("#!provenance column 3 is gated on the completeness tier; a structural-only call exports "
               "as repeat_region, never a specific subclass")
    out += ["\t".join(str(c) for c in r) for r in rows]
    if sequence:
        out.append("##FASTA")
        out.append(f">{seqid}")
        out += [sequence[i:i + 60] for i in range(0, len(sequence), 60)]
    return "\n".join(out) + "\n"


def to_bed(rec, seqid="locus") -> str:
    """BED12-less BED6: chrom, start (0-based), end, name, score, strand — the lowest common denominator
    every browser reads. GFF3 carries the ontology; BED carries the intervals."""
    lines = []
    for r in build_features(rec, seqid):
        chrom, _src, term, start, end, score, strand, _ph, attrs = r
        name = term
        for part in attrs.split(";"):
            if part.startswith("Name="):
                name = part[5:]
                break
        lines.append("\t".join([str(chrom), str(int(start) - 1), str(int(end)), name,
                                str(_bed_score(score)),
                                strand if strand in ("+", "-") else "."]))
    return "\n".join(lines) + "\n"


def _bed_score(score) -> int:
    """BED's score column is a 0-1000 integer that drives shading, not a free-form float. The GFF3 column 6
    carries the real HMMER bit score; here it is rounded and clamped into the valid range (monotone with
    confidence) so browsers accept the file instead of rejecting a raw, unbounded bit score."""
    if score in (".", None, ""):
        return 0
    try:
        return max(0, min(1000, int(round(float(score)))))
    except (TypeError, ValueError, OverflowError):        # OverflowError: int(round(inf))
        return 0
