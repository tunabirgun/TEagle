"""GFF3 and BED export of a TEagle annotation.

TEagle is otherwise a terminal node: a result that cannot leave the application undercuts the one thing
it does that a hosted search does not — dissect a candidate, design the confirming assay, and seal the
run. This is the interoperability half of that.

Column 3 of GFF3 is an ASSERTION that downstream software trusts. It is therefore gated on the same
domain completeness the report shows: a structural-only call exports as the generic `repeat_region`, never
as `LTR_retrotransposon`, because a specific subclass claims coding evidence that a structural-only call
does not have. Wicker codes have no SO accession and go in LOWERCASE attributes — GFF3 reserves
capitalised tags for its own predefined set, so the `TEID=` convention seen in the wild is non-compliant.

Term names and accessions were verified on 2026-07-28 against the Sequence Ontology's own
`so-simple.obo` (The-Sequence-Ontology/SO-Ontologies, CC-BY-4.0). EDTA's mapping file is GPL-3.0 and is
deliberately not used. Format: GFF3 specification; Diesh et al. 2023 Genome Biol 24:74.
"""
from __future__ import annotations

import hashlib
import re

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
    # verified against the Sequence Ontology (OLS4, 2026-10-04): "A repeat where the same sequence is repeated
    # in the same direction." The orientation-neutral term for an unconfirmed terminal direct-repeat pair.
    "direct_repeat": "SO:0000314",
    "polypeptide_domain": "SO:0000417",
    # A LINE's genomic 3' poly-A / 5' poly-T tail. SO's polyA_sequence (SO:0000610) is specifically the
    # post-transcriptional poly-A ADDED TO AN mRNA, so it would misdescribe this genomic tract; the honest
    # generic term is sequence_feature. Verified on OLS4: "Any extent of continuous biological sequence."
    "sequence_feature": "SO:0000110",
    # verified against the Sequence Ontology (OLS4): "The recognition sequence necessary for endonuclease
    # cleavage of an RNA transcript that is followed by polyadenylation; consensus=AATAAA."
    "polyA_signal_sequence": "SO:0000551",
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

# structural evidence type (first token of ev["type"]) -> SO term for the sub-feature rows
_STRUCT_TERM = {
    "LTR": "long_terminal_repeat",
    "TIR": "terminal_inverted_repeat",
    "TSD": "target_site_duplication",
    "PBS": "primer_binding_site",
    "PPT": "RR_tract",
    # LINE tail: no genomic-specific SO subclass exists, so emit the generic sequence_feature (Name=
    # "poly-A"/"poly-T", the same short-token convention as the other rows) rather than drop this
    # on-screen track from the file.
    "poly-A": "sequence_feature",
    "poly-T": "sequence_feature",
    # the advisory poly(A)-SIGNAL motif inside an LTR (first token of "poly(A) signal (motif)"). SO's
    # polyA_signal_sequence is exactly this recognition hexamer ("consensus=AATAAA"), NOT the tail above,
    # so the two never share a term. The row still carries the motif's own hedge in `note`.
    "polyA-signal": "polyA_signal_sequence",
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
    unassessed = comp.get("unassessed") or []
    if not te_class or te_class == "none" or "structural-only" in tier or "structural-only" in te_class:
        return "repeat_region"
    for prefix, term in _CLASS_TERM.items():
        if te_class.startswith(prefix):
            if term == "terminal_inverted_repeat_element" and any("TIR" in m for m in list(missing) + list(unassessed)):
                return "DNA_transposon"
            if term == "LTR_retrotransposon" and "RT" not in present:
                return "repeat_region"
            return term
    return "repeat_region"


def _esc(v) -> str:
    """GFF3 attribute escaping: the reserved characters must be percent-encoded, not stripped.
    '%' is encoded FIRST so the '%' we introduce for the others (e.g. '%3B') is not re-escaped, and a
    literal '%' in a value (e.g. a PBS note's '61.1%') is preserved rather than left as a bare '%'."""
    s = "" if v is None else str(v)
    for ch, rep in (("%", "%25"), (";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C"),
                    ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D")):
        s = s.replace(ch, rep)
    return s


def _attrs(pairs) -> str:
    return ";".join(f"{k}={_esc(v)}" for k, v in pairs if v not in (None, "", []))


def domain_scan_state(cl) -> tuple:
    """(domain_evidence_status, ORFs outside the profile search's ORF budget) for one classification.
    The status keeps its three values; a partial search is declared beside it so a negative is never
    read as a complete one."""
    cl = cl or {}
    if cl.get("domains_unavailable"):
        return "unavailable_failed", 0
    status = "tested_empty" if not cl.get("n_domains") else "tested_with_hits"
    return status, max(0, int(cl.get("orfs_unscanned") or 0))


def _orders(cl) -> tuple:
    """(wicker_order, domain_architecture) from a classification; a hand-built one without them falls back."""
    cl = cl or {}
    withheld = cl.get("outcome") == "withheld"
    return (cl.get("wicker_order") or ("withheld" if withheld else "unassigned"),
            cl.get("domain_architecture") or cl.get("order"))


def is_neutral_tdr(cl) -> bool:
    """A withheld terminal direct-repeat candidate: the repeats are exported, but nothing in the file may
    assert LTR (retroviral) origin that the classification itself declines to claim."""
    cl = cl or {}
    kind = str((cl.get("completeness") or {}).get("kind") or "")
    return cl.get("te_class") == "repeat/structural-only" and kind.startswith("terminal direct-repeat candidate")


def _element_span(cl, structural, length):
    """0-based half-open [start, end) of the element from the non-advisory terminal-repeat evidence behind
    the call, or None when no evidence bounds it (domain-only and poly-A-only calls keep the record extent)."""
    te_class = (cl.get("te_class") or "")
    if te_class.startswith("LTR/"):
        kind = "LTR"
    elif te_class.startswith("DNA/"):
        comp = cl.get("completeness") or {}
        if any("TIR" in m for m in list(comp.get("missing") or []) + list(comp.get("unassessed") or [])):
            return None                                # TIRs not credited to this call do not bound the element
        kind = "TIR"
    elif te_class == "repeat/structural-only":
        kind = "TIR" if "inverted" in (cl.get("superfamily") or "") else "LTR"
    else:
        return None
    for ev in structural or []:
        span = ev.get("element_span")
        if ((ev.get("type") or "").split(" ")[0] == kind and not ev.get("advisory") and
                isinstance(span, (list, tuple)) and len(span) == 2 and 0 <= span[0] < span[1] <= length):
            return int(span[0]), int(span[1])
    return None


def build_features(rec, seqid="locus", *, identity_trust="local_input", record_id_sha256=None):
    """One annotation record -> a list of GFF3 feature tuples, parent first.

    The parent is the element itself: its terminal-repeat span when the evidence bounds it, else the
    analysed record. Flanking target-site duplications lie outside the element, so a flanked record also
    gets a repeat_region wrapper (TE_1.region) that parents them."""
    original_seqid = str(seqid or "")
    identity = _writer_identity(original_seqid, identity_trust, record_id_sha256)
    seqid = _safe_seqid(identity["id"])
    cl = rec.get("classification") or {}
    comp = cl.get("completeness") or {}
    length = (rec.get("composition") or {}).get("length") or len(rec.get("seq_preview") or "") or 1
    term = so_term_for(cl)
    outcome = cl.get("outcome") or "complete"
    wicker_order, architecture = _orders(cl)
    domain_status, orfs_unscanned = domain_scan_state(cl)
    neutral = is_neutral_tdr(cl)
    gene_id = "TE_1"
    structural = rec.get("structural") or []
    span = _element_span(cl, structural, length)
    children = []                                  # (start0, end0, is_tsd, row builder taking the Parent)
    n = 0
    for ev in structural:
        kind = (ev.get("type") or "").split(" ")[0]
        sub = _STRUCT_TERM.get(kind)
        if not sub:
            continue
        name = kind
        if neutral and kind == "LTR":
            sub, name = "direct_repeat", "terminal_direct_repeat"
        elif neutral and kind == "PPT":
            sub = "sequence_feature"
        elif kind == "PBS" and ev.get("confident") is False:
            sub, name = "sequence_feature", "PBS_candidate"      # below the confident threshold: advisory, not a PBS call
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
        for sp in spans:
            n += 1

            def build(parent, sp=sp, ev=ev, sub=sub, name=name, kind=kind, n=n):
                return (seqid, "TEagle", sub, int(sp[0]) + 1, int(sp[1]), ".", ".", ".", _attrs([
                    ("ID", f"{gene_id}.s{n}"), ("Parent", parent), ("Name", name),
                    ("Ontology_term", SO.get(sub, "")),
                    # the repeats of an unconfirmed-origin candidate must not read as established LTR cis-elements
                    ("context", "ltr_origin_unconfirmed" if neutral and kind in ("LTR", "PBS", "PPT", "polyA-signal")
                     else None),
                    # a length still limited by the record is declared here too, not silently exported
                    ("length_is_lower_bound", "true" if ev.get("length_is_lower_bound") else None),
                    # carry the evidence's OWN confidence/hedge into the file — dropping it would let a browser
                    # read an unhedged call the on-screen panel qualified (e.g. a diverged PBS whose priming tRNA
                    # is "undetermined"). _attrs omits any field that is None/""/[] , so these appear only when set.
                    ("identity_pct", ev.get("identity")),
                    ("priming_trna", ev.get("priming_trna")),
                    ("priming_trna_best_match", ev.get("best_match")),
                    ("confident", ("true" if ev.get("confident") else "false") if "confident" in ev else None),
                    ("note", ev.get("note")),
                    ("purine_frac", ev.get("purine_frac")),   # PPT tract quality (fraction of purines)
                    # the TSD-length congruence verdict classify() computed against the superfamily's expected
                    # length — the hedge shown on screen ("ends not credited from it") must reach the file too.
                    ("tsd_congruence", ev.get("tsd_congruence")),
                    ("tsd_state", ev.get("tsd_state")),
                    # poly(A)-signal motif: which hexamer, how common that variant is, whether the gating
                    # downstream element was found (or could not be assessed), and whether both LTR copies
                    # agree. A browser must not read this motif as a located cleavage site.
                    ("polya_variant", ev.get("variant")),
                    ("downstream_element", ev.get("dse_state")),
                    ("in_both_ltr_copies", ("true" if ev.get("in_both_ltr_copies") else "false")
                     if "in_both_ltr_copies" in ev else None),
                ]))
            children.append((int(sp[0]), int(sp[1]), kind == "TSD", build))
    for i, d in enumerate(rec.get("domains") or [], 1):
        nt = d.get("nt") or [0, 0]

        def build(parent, d=d, nt=nt, i=i):
            return (seqid, "TEagle", "polypeptide_domain", int(nt[0]) + 1, int(nt[1]),
                    f"{d.get('score', '.')}", d.get("strand", "."), ".", _attrs([
                        ("ID", f"{gene_id}.d{i}"), ("Parent", parent),
                        ("Name", d.get("label") or d.get("domain")),
                        ("Ontology_term", SO["polypeptide_domain"]),
                        ("domain_code", d.get("domain")), ("hmm_profile", d.get("hmm")),
                        ("pfam_accession", d.get("pfam")), ("evalue", d.get("evalue")),
                    ]))
        children.append((int(nt[0]), int(nt[1]), False, build))
    te_lo, te_hi = span if span else (0, length)
    wrapper = None                                 # element plus its TSDs only; never stretched over other flank features
    if span:
        tsds = [c for c in children if c[2]]
        if tsds:
            wrapper = (min([te_lo] + [c[0] for c in tsds]), max([te_hi] + [c[1] for c in tsds]))
    rows = []
    wrapper_id = f"{gene_id}.region"
    name = cl.get("superfamily") or term
    if wrapper:
        rows.append((seqid, "TEagle", "repeat_region", wrapper[0] + 1, wrapper[1], ".", ".", ".", _attrs([
            ("ID", wrapper_id), ("Name", f"{name} region"), ("Ontology_term", SO["repeat_region"]),
            ("note", "the element plus the target-site duplications that flank it")])))
    # parent: the element itself. Wicker code and the tier ride in lowercase custom attributes.
    rows.append((seqid, "TEagle", term, te_lo + 1, te_hi, ".", ".", ".", _attrs([
        ("ID", gene_id),
        ("Parent", wrapper_id if wrapper else None),
        ("Name", name),
        ("Ontology_term", SO.get(term, "")),
        ("wicker_code", cl.get("te_class")),
        ("te_class", cl.get("te_class")),
        ("classification_confidence", cl.get("confidence")),
        ("classification_outcome", outcome),
        ("wicker_order", wicker_order),
        ("domain_architecture", architecture),
        ("classification_order", architecture or wicker_order),      # deprecated alias of domain_architecture
        ("domain_evidence_status", domain_status),
        ("domain_scan_scope", "partial" if orfs_unscanned else None),
        ("orfs_not_scanned", orfs_unscanned or None),
        ("structural_evidence_status", "unavailable_failed" if cl.get("structural_unavailable") else None),
        ("classification_explanation", cl.get("explanation")),
        ("completeness_tier", comp.get("tier")),
        ("domains_tested_scope", comp.get("scope")),
        ("original_record_sha256", identity.get("idSha256")),
    ])))
    for lo, hi, is_tsd, build in children:
        if is_tsd and wrapper:
            parent = wrapper_id
        else:                                      # a feature outside the element has no parent rather than a wrong one
            parent = gene_id if te_lo <= lo and hi <= te_hi else None
        rows.append(build(parent))
    return rows


def _shift(row, offset: int):
    """Move one GFF3 row into genome space. Columns 4/5 are 1-based inclusive, so a locus feature at 1-based
    p sits at p + offset, where offset is the 0-based genome coordinate of the locus's first base."""
    chrom, src, term, start, end, score, strand, phase, attrs = row
    return (chrom, src, term, int(start) + offset, int(end) + offset, score, strand, phase, attrs)


def to_gff3(rec, seqid="locus", sequence=None, source_note=None, offset: int = 0, region_span=None, *,
            identity_trust="local_input", record_id_sha256=None) -> str:
    """A complete, self-contained GFF3 document.

    Coordinates are LOCUS-RELATIVE unless the caller supplies a genome-anchored seqid, because exporting
    locus offsets under a chromosome name would silently misplace every feature. The sequence is embedded
    after ##FASTA so the file stands alone in a browser that has no matching reference."""
    original_seqid = str(seqid or "")
    identity = _writer_identity(original_seqid, identity_trust, record_id_sha256)
    seqid = _safe_seqid(identity["id"])
    rows = build_features(rec, original_seqid, identity_trust=identity_trust,
                          record_id_sha256=record_id_sha256)
    if offset:
        rows = [_shift(r, offset) for r in rows]
    length = (rec.get("composition") or {}).get("length") or 1
    lo, hi = region_span if region_span else (1, length)
    out = ["##gff-version 3", f"##sequence-region {seqid} {lo} {hi}"]
    if offset:
        out.append(f"#!coordinates genome-anchored on {seqid}; a +{offset} bp locus offset is applied to "
                   "every feature. The sequence is deliberately NOT embedded below: this file describes a "
                   "slice of a chromosome, and a ##FASTA block under the chromosome's name would supply the "
                   "slice as if it were the whole sequence.")
    if source_note:
        note = str(source_note).replace("\t", " ").replace("\r", " ").replace("\n", " ")
        out.append(f"#!source {note}")
    out.append("#!provenance column 3 is gated on the domain completeness; a structural-only call exports "
               "as repeat_region, never a specific subclass")
    out += ["\t".join(str(c) for c in r) for r in rows]
    if sequence and not offset:
        out.append("##FASTA")
        out.append(f">{seqid}")
        out += [sequence[i:i + 60] for i in range(0, len(sequence), 60)]
    return "\n".join(out) + "\n"


def to_bed(rec, seqid="locus", offset: int = 0, *, identity_trust="local_input", record_id_sha256=None) -> str:
    """BED12-less BED6: chrom, start (0-based), end, name, score, strand — the lowest common denominator
    every browser reads. GFF3 carries the ontology; BED carries the intervals.

    `offset` shifts every interval into genome space (the 0-based genome coordinate of the locus's first
    base), so the file can be intersected against a chromosome-space track. Without it a locus-relative BED
    intersected against a real genome track returns zero overlaps — silently, which is the dangerous part."""
    lines = []
    original_seqid = str(seqid or "")
    identity = _writer_identity(original_seqid, identity_trust, record_id_sha256)
    seqid = _safe_seqid(identity["id"])
    cl = rec.get("classification") or {}
    outcome = cl.get("outcome") or "complete"
    order = _orders(cl)[0]
    domain_status, orfs_unscanned = domain_scan_state(cl)
    for r in build_features(rec, original_seqid, identity_trust=identity_trust,
                            record_id_sha256=record_id_sha256):
        if offset:
            r = _shift(r, offset)
        chrom, _src, term, start, end, score, strand, _ph, attrs = r
        name = term
        for part in attrs.split(";"):
            if part.startswith("Name="):
                name = part[5:]
                break
        if "ID=TE_1" in attrs.split(";"):
            name += f"|outcome={outcome}|order={order}|domains={domain_status}"
            if orfs_unscanned:
                name += f"|domain_scan=partial|orfs_not_scanned={orfs_unscanned}"
        if identity.get("idSha256"):
            name += f"|record_sha256={identity['idSha256']}"
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


_SAFE_SEQID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
_SAFE_FILENAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_PUBLIC_ACCESSION = re.compile(
    r"(?:[A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{3}\d{5}|[A-Z]{4}\d{8,10}|[A-Z]{6}\d{9,11}|"
    r"[A-Z]{2}_\d+|NZ_[A-Z]{2,6}\d+|GC[AF]_\d+)(?:\.\d+)?"
)
_PUBLIC_LOCUS = re.compile(
    rf"(?:{_PUBLIC_ACCESSION.pattern}|chr[A-Za-z0-9][A-Za-z0-9_.-]*)(?::\d+(?:-\d+)?)?"
)
_GENERATED_PUBLIC_ID = re.compile(r"(?:record_[0-9a-f]{12,64}|P[1-9]\d*)")
_GENERATED_NEUTRAL_IDS = {"record", "input_sequence", "multiple_records", "batch", "splice_locus", "locus"}
_MODEL_FEATURES = (("exons", "exon", "SO:0000147"),
                   ("introns", "intron", "SO:0000188"),
                   ("cds", "CDS", "SO:0000316"))
_DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
_WINDOWS_RESERVED = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
                     *(f"LPT{i}" for i in range(1, 10))}


def _identity_digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]


def _identity_sha256(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def public_record_identity(record_id, *, trust: str = "local_input") -> dict:
    """Return a public record ID and, when rewritten, its collision-resistant identity digest.

    Only recognized public sequence accessions/loci and application-generated neutral IDs remain readable.
    User headers and sample-like labels become neutral content-derived IDs; the full digest keeps scientific
    bindings collision-resistant even when a compact display prefix collides.
    """
    original = str(record_id or "")
    if trust not in {"local_input", "validated_fetch", "generated"}:
        raise ValueError("record identity trust must be local_input, validated_fetch, or generated")
    trusted = ((trust == "validated_fetch" and bool(_PUBLIC_LOCUS.fullmatch(original))) or
               (trust == "generated" and
                (original in _GENERATED_NEUTRAL_IDS or bool(_GENERATED_PUBLIC_ID.fullmatch(original)))))
    if trusted:
        return {"id": original}
    digest = _identity_sha256(original)
    return {"id": f"record_{_identity_digest(original)}", "idSha256": digest}


def public_record_id(record_id, *, trust: str = "local_input") -> str:
    """Public string form of :func:`public_record_identity`."""
    return public_record_identity(record_id, trust=trust)["id"]


def public_record_label(value, *, trust: str = "local_input") -> str:
    """Public record ID with an application-owned PCR role suffix preserved."""
    raw = str(value or "")
    for suffix in (" (input)", " (custom background)"):
        if raw.endswith(suffix):
            return public_record_id(raw[:-len(suffix)], trust=trust) + suffix
    return public_record_id(raw, trust=trust)


def bound_public_record_label(value, id_sha256=None) -> str:
    """Preserve a rewritten public label only when its full digest proves the compact-ID binding."""
    raw = str(value or "")
    base, suffix = raw, ""
    for candidate in (" (input)", " (custom background)"):
        if base.endswith(candidate):
            base, suffix = base[:-len(candidate)], candidate
            break
    if id_sha256 is not None:
        digest = str(id_sha256).lower()
        if re.fullmatch(r"[0-9a-f]{64}", digest) and base == f"record_{digest[:12]}":
            return base + suffix
    return public_record_label(raw)


def is_public_accession(value) -> bool:
    """Return whether ``value`` is a recognized public sequence/assembly accession."""
    return bool(_PUBLIC_ACCESSION.fullmatch(str(value or "")))


def is_public_record_id(value) -> bool:
    """Return whether ``value`` may be emitted without rewriting."""
    raw = str(value or "")
    return (raw in _GENERATED_NEUTRAL_IDS or bool(_GENERATED_PUBLIC_ID.fullmatch(raw)) or
            bool(_PUBLIC_LOCUS.fullmatch(raw)))


def _safe_seqid(value: str) -> str:
    return value if _SAFE_SEQID.fullmatch(value) else f"record_{_identity_digest(value)}"


def _writer_identity(record_id, trust="local_input", record_id_sha256=None) -> dict:
    original = str(record_id or "")
    if record_id_sha256 is None:
        return public_record_identity(original, trust=trust)
    digest = str(record_id_sha256).lower()
    if not re.fullmatch(r"[0-9a-f]{64}", digest) or original != f"record_{digest[:12]}":
        raise ValueError("record identity SHA-256 does not match the compact public ID")
    return {"id": original, "idSha256": digest}


def sanitize_seqid(record_id, *, trust="local_input") -> str:
    """Return a deterministic GFF3/BED-safe sequence identifier.

    Path separators, control characters and GFF3 delimiters cause a content-derived replacement rather
    than a lossy transliteration. Writers preserve the original identifier in an attribute when this
    function rewrites it.
    """
    return _safe_seqid(public_record_id(record_id, trust=trust))


def safe_filename_stem(record_id, default: str = "record", *, trust="local_input") -> str:
    """Return a Windows-safe deterministic stem; never interpret a record ID as a path."""
    original = str(record_id or "")
    public = public_record_id(original, trust=trust)
    fallback = str(default or "record")
    if (not _SAFE_FILENAME.fullmatch(fallback) or fallback in {".", ".."} or
            fallback.split(".", 1)[0].upper() in _WINDOWS_RESERVED):
        fallback = "record"
    if (_SAFE_FILENAME.fullmatch(public) and
            public.split(".", 1)[0].upper() not in _WINDOWS_RESERVED):
        return public.rstrip(". ") or fallback
    return f"{fallback}_{_identity_digest(original)}"


def _model_rows(model: dict, sequence_length: int | None = None):
    if not isinstance(model, dict):
        raise ValueError("gene model must be an object")
    if sequence_length is not None and (isinstance(sequence_length, bool) or
                                        not isinstance(sequence_length, int) or sequence_length < 1):
        raise ValueError("sequence length must be a positive integer")
    rows = []
    for key, feature_type, ontology in _MODEL_FEATURES:
        features = model.get(key, [])
        if features is None:
            features = []
        if not isinstance(features, list):
            raise ValueError(f"{key} must be a list")
        for source_index, feature in enumerate(features):
            if not isinstance(feature, dict):
                raise ValueError(f"{key}[{source_index}] must be an object")
            start, end = feature.get("start"), feature.get("end")
            if (isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) or
                    not isinstance(end, int) or start < 0 or end <= start):
                raise ValueError(f"{key}[{source_index}] is not a positive 0-based half-open interval")
            if sequence_length is not None and end > sequence_length:
                raise ValueError(f"{key}[{source_index}] exceeds the source sequence")
            strand = feature.get("strand", model.get("strand", ".")) or "."
            if strand not in {"+", "-", "."}:
                raise ValueError(f"{key}[{source_index}] has invalid strand {strand!r}")
            rows.append({"key": key, "type": feature_type, "ontology": ontology,
                         "start": start, "end": end, "strand": strand,
                         "source_index": source_index, "feature": feature})
    order = {feature_type: i for i, (_key, feature_type, _ontology) in enumerate(_MODEL_FEATURES)}
    rows.sort(key=lambda row: (row["start"], row["end"], order[row["type"]], row["source_index"]))
    return rows


def _clean_sequence(sequence: str) -> str:
    if not isinstance(sequence, str) or ">" in sequence:
        raise ValueError("source sequence must be one raw nucleotide sequence, not FASTA")
    clean = "".join(sequence.split()).upper().replace("U", "T")
    if not clean or not _DNA.fullmatch(clean):
        raise ValueError("source sequence must contain only IUPAC nucleotide symbols")
    return clean


def gene_model_to_gff3(model: dict, seqid, *, sequence: str | None = None, source: str = "TEagle",
                       model_kind: str = "fetched_gene_model", identity_trust="local_input",
                       record_id_sha256=None) -> str:
    """Export fetched or splice-derived features whose input coordinates are 0-based half-open."""
    original = str(seqid or "")
    identity = _writer_identity(original, identity_trust, record_id_sha256)
    safe = _safe_seqid(identity["id"])
    clean = _clean_sequence(sequence) if sequence is not None else None
    rows = _model_rows(model, len(clean) if clean is not None else None)
    if not rows:
        raise ValueError("gene model has no exon, intron or CDS intervals")
    maximum = len(clean) if clean is not None else max(row["end"] for row in rows)
    source_name = str(source or "TEagle").replace("\t", " ").replace("\r", " ").replace("\n", " ")
    kind = str(model_kind or "gene_model")
    out = ["##gff-version 3", f"##sequence-region {safe} 1 {maximum}",
           "#!coordinate-system input=0-based-half-open; gff3=1-based-inclusive"]
    counters = {"exon": 0, "intron": 0, "CDS": 0}
    for row in rows:
        counters[row["type"]] += 1
        feature = row["feature"]
        phase = "."
        if row["type"] == "CDS":
            phase = feature.get("phase")
            if isinstance(phase, bool) or not isinstance(phase, int) or phase not in {0, 1, 2}:
                raise ValueError(f"cds[{row['source_index']}] has no valid GFF3 phase; BED and FASTA remain available")
            phase = str(phase)
        group = feature.get("cds_group") if row["type"] == "CDS" else None
        # the segments of one spliced CDS are ONE discontinuous feature, which GFF3 writes as lines sharing an ID
        # (and one Name, so validators do not see conflicting attributes on the shared ID)
        attrs = [("ID", f"annotation.{group}" if group else f"annotation.{row['type'].lower()}.{counters[row['type']]}"),
                 ("Name", group or f"{row['type']}_{counters[row['type']]}"),
                 ("Ontology_term", row["ontology"]), ("model_kind", kind),
                 ("original_record_sha256", identity.get("idSha256")),
                 ("derived", "true" if feature.get("derived") else None),
                 ("cds_group", group),
                 ("note", feature.get("note")), ("donor", feature.get("donor")),
                 ("acceptor", feature.get("acceptor")),
                 ("canonical", ("true" if feature.get("canonical") else "false")
                  if "canonical" in feature else None)]
        out.append("\t".join((safe, source_name, row["type"], str(row["start"] + 1),
                              str(row["end"]), ".", row["strand"], phase, _attrs(attrs))))
    if clean is not None:
        out.extend(("##FASTA", f">{safe}"))
        out.extend(clean[i:i + 60] for i in range(0, len(clean), 60))
    return "\n".join(out) + "\n"


def gene_model_to_bed(model: dict, seqid, *, sequence_length: int | None = None,
                      source: str = "TEagle", model_kind: str = "fetched_gene_model",
                      identity_trust="local_input", record_id_sha256=None) -> str:
    """Export BED6 without changing the model's 0-based half-open coordinates."""
    del source, model_kind  # BED6 has no source/model columns; GFF3 retains these fields.
    original = str(seqid or "")
    identity = _writer_identity(original, identity_trust, record_id_sha256)
    safe = _safe_seqid(identity["id"])
    rows = _model_rows(model, sequence_length)
    if not rows:
        raise ValueError("gene model has no exon, intron or CDS intervals")
    counters = {"exon": 0, "intron": 0, "CDS": 0}
    out = []
    for row in rows:
        counters[row["type"]] += 1
        name = f"{row['type']}_{counters[row['type']]}"
        if row["feature"].get("derived"):
            name += "|derived=true"
        if identity.get("idSha256"):
            name += f"|record_sha256={identity['idSha256']}"
        out.append("\t".join((safe, str(row["start"]), str(row["end"]), name, "0", row["strand"])))
    return "\n".join(out) + "\n"


def gene_model_to_fasta(model: dict, sequence: str, seqid, *, feature_types=("exon", "CDS"),
                        line_width: int = 60, identity_trust="local_input", record_id_sha256=None) -> str:
    """Export exon/CDS intervals from the exact source sequence, in deterministic genomic order."""
    clean = _clean_sequence(sequence)
    if isinstance(line_width, bool) or not isinstance(line_width, int) or line_width < 1:
        raise ValueError("line width must be a positive integer")
    requested = tuple(feature_types)
    if not requested or any(t not in {"exon", "intron", "CDS"} for t in requested):
        raise ValueError("feature_types must contain exon, intron and/or CDS")
    rows = [row for row in _model_rows(model, len(clean)) if row["type"] in requested]
    if not rows:
        raise ValueError("gene model has no requested intervals")
    original = str(seqid or "")
    identity = _writer_identity(original, identity_trust, record_id_sha256)
    safe = _safe_seqid(identity["id"])
    complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    counters = {"exon": 0, "intron": 0, "CDS": 0}
    out = []
    for row in rows:
        counters[row["type"]] += 1
        fragment = clean[row["start"]:row["end"]]
        if row["strand"] == "-":
            fragment = fragment.translate(complement)[::-1]
        digest_note = f" record_sha256={identity['idSha256']}" if identity.get("idSha256") else ""
        digest_note += " derived=true" if row["feature"].get("derived") else ""
        desc = (f">{safe}.{row['type'].lower()}{counters[row['type']]}{digest_note} "
                f"coordinates={row['start']}-{row['end']} coordinate_system=0-based-half-open "
                f"strand={row['strand']}")
        out.append(desc)
        out.extend(fragment[i:i + line_width] for i in range(0, len(fragment), line_width))
    return "\n".join(out) + "\n"
