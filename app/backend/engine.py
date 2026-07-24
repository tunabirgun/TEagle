"""TEagle engine adapter — the single source of truth for every operation.

Each `run_*` function takes an already-parsed request dict and returns a result dict,
or raises `BadRequest` for a malformed client request (never an unhandled 500 for user error).
Both the stdlib HTTP server (server.py) and the native PySide6 app (app/native) call these
functions, so request validation and scientific behaviour cannot drift between the two front ends.

Error taxonomy (consumed by the Qt worker, obs 4702):
  - returns dict              -> success
  - raises BadRequest         -> user-correctable input error (HTTP 400 / inline message)
  - raises anything else      -> unexpected fault (HTTP 500 / error banner + traceback to log)
"""
from __future__ import annotations
import hashlib, re, time

from teagle_core import (sequtil, structural, primers, provenance, fetch,
                         domains, classify, refs, wsl, timing, genomepcr, oligoqc, retroviral)
import envcheck


class BadRequest(Exception):
    """A malformed client request — surfaced as HTTP 400 / an inline user message, never a 500."""


# ---------- request-value coercion (encodes hard-won fixes: NaN max_mm, float tp, non-string species) ----------
def _num(v, default):
    """Coerce a scalar to a number: None -> default, a non-numeric value -> BadRequest (not a 500 crash)."""
    if v is None:
        return default
    if isinstance(v, bool):
        raise BadRequest("expected a number, got a boolean")
    if isinstance(v, (int, float)):
        if isinstance(v, float) and v != v or v in (float("inf"), float("-inf")):
            raise BadRequest("expected a finite number")
        return v
    try:
        f = float(v)
        if f != f or f in (float("inf"), float("-inf")):     # NaN / Infinity literals -> 400, not a 500
            raise BadRequest("expected a finite number")
        return int(f) if f.is_integer() and "." not in str(v) and "e" not in str(v).lower() else f
    except (TypeError, ValueError):
        raise BadRequest(f"expected a number, got {v!r}")


def _clean_params(p):
    """Coerce every value of a params object through _num so a mistyped OR non-finite (NaN/Infinity)
    parameter is a clean 400, never a 500 or a non-strict-JSON token in the sealed manifest."""
    if p is None:
        return {}
    if not isinstance(p, dict):
        raise BadRequest("params must be a JSON object")
    return {k: (_num(v, v) if isinstance(v, (str, int, float)) and not isinstance(v, bool) else v)
            for k, v in p.items()}


def _require_nt(seq, label="sequence"):
    """Reject input that is clearly not nucleotide (e.g. an accession pasted into the sequence box),
    so those operations fail with a clear message instead of silently succeeding with no result."""
    s = (seq or "").upper()
    if not s:
        raise BadRequest(f"empty {label}")
    if sum(1 for c in s if c in "ACGTUN") / len(s) < 0.5:
        raise BadRequest(f"{label} does not look like nucleotide sequence — if you meant an accession, fetch it first")


_AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZUO*")
def _require_aa(seq, label="protein"):
    """Reject a reference protein that is actually nucleotide (miniprot needs amino acids), so a
    pasted DNA sequence fails clearly instead of mis-aligning as a degenerate protein."""
    s = (seq or "").upper()
    if not s:
        raise BadRequest(f"empty {label}")
    if sum(1 for c in s if c in "ACGTUN") / len(s) > 0.9:
        raise BadRequest(f"{label} looks like nucleotide, not amino acids — miniprot needs a protein sequence")
    if sum(1 for c in s if c in _AA_ALPHABET) / len(s) < 0.9:
        raise BadRequest(f"{label} does not look like a valid amino-acid sequence")


# ---------- health / environment ----------
def run_health():
    return {"ok": True, "primer3": primers.PRIMER3_VERSION, "core": provenance.TEAGLE_VERSION}


def run_env():
    try:
        env = envcheck.check()
        env["wsl"] = wsl.available()          # fast availability probe
        return env
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


def run_wsl_status():
    try:
        return wsl.env_status()
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}


def run_wsl_install():
    return wsl.start_install()


def run_wsl_install_log():
    return {"log": wsl.install_log(60)}


def run_wsl_install_wsl2():
    return wsl.install_wsl2()


def run_wsl2_install_log():
    return {"log": wsl.wsl2_install_log(200)}


def run_wsl_components():
    try:
        return wsl.components_status()
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}", "components": []}


def run_wsl_repair(body):
    return wsl.repair_component((body or {}).get("component", ""))


def run_wsl_integrity():
    try:
        return wsl.integrity_check()
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}", "checks": []}


# ---------- analysis ----------
def run_analyze(body):
    return analyze(body.get("sequence", ""), body.get("source"))


def analyze(seq_text: str, source: dict | None = None):
    recs = sequtil.parse_fasta(seq_text)
    # RNA input (normalized U->T by parse_fasta): detect U in the SEQUENCE body only, never in a header line
    _body = "".join(l for l in seq_text.splitlines() if not l.startswith(">")) if isinstance(seq_text, str) else ""
    rna = "U" in _body.upper()
    out = []
    for rid, seq in recs:
        ok, bad = sequtil.validate_iupac(seq)
        rec = {"id": rid, "valid": ok,
               "invalid": [{"pos": p, "char": c} for p, c in bad],
               "composition": sequtil.composition(seq),
               "structural": [], "domains": [], "classification": None, "orfs": [],
               "seq_preview": seq[:120], "notes": []}
        if ok:
            try:
                rec["structural"] = structural.detect_all(seq)      # isolate: a detector fault must not 500 the whole analyze
            except Exception as e:
                rec["notes"].append(f"structural detection unavailable: {type(e).__name__}")
            try:
                rec["domains"] = domains.scan_domains(seq)
            except Exception as e:
                rec["notes"].append(f"domain scan unavailable: {type(e).__name__}")
            rec["classification"] = classify.classify(rec["structural"], rec["domains"])
            rec["orfs"] = sequtil.find_orfs(seq)
            # retroviral transcript architecture (ERV only) — the correct coding-organisation model in place
            # of a host exon-intron gene model; None for non-ERVs, so it never overrides a normal record.
            rec["retroviral"] = None
            try:
                rec["retroviral"] = retroviral.transcript_architecture(
                    rec["structural"], rec["domains"], rec["classification"], len(seq))
            except Exception as e:                          # a retroviral-model fault must not fail the analysis
                rec["notes"].append(f"retroviral architecture unavailable: {type(e).__name__}")
            if rna:
                rec["notes"].append("RNA input — U was read as T (DNA equivalent) for all steps")
            if any(o.get("open_end") for o in rec["orfs"]):
                rec["notes"].append("sequence appears 3′-truncated — some ORFs are open-ended; domain calls on them may be partial")
            if len(seq) > 200_000:
                rec["notes"].append("large input — structural/ORF heuristics are tuned for single loci; treat as approximate")
        out.append(rec)
    joined = "".join(s for _, s in recs)
    any_dom = any(r["domains"] for r in out)
    references = refs.for_run("analysis", domains=any_dom, fetched=bool(source))
    dbs = ([{"name": "Pfam TE-domain profiles (bundled, CC0)", "file": "data/te_domains.hmm",
             "sha256": domains.HMM_SHA256}] if domains.HMM_SHA256 else [])
    result = {"records": out, "references": references,
              "provenance": provenance.build_manifest("analysis", joined,
                            recs[0][0] if recs else "empty", {"orf_min_aa": 40},
                            source=source, references=references, databases=dbs)}
    if not recs:
        result["warning"] = "No sequence provided — paste, upload, or fetch a sequence first."
    elif len(recs) > 1:
        result["warning"] = f"{len(recs)} records found; downstream primer / family steps use the first record only."
    return result


# ---------- accession fetch ----------
def run_fetch(body):
    acc = body.get("accession", "")
    if not isinstance(acc, str):
        raise BadRequest("accession must be a string")
    try:
        meta = fetch.retrieve(acc, refresh=bool(body.get("refresh")))
        return {"ok": True, **meta}
    except fetch.FetchError as fe:
        return {"ok": False, "error": str(fe)}


# ---------- coordinate-based fetch (organism + chr:start-end) ----------
def run_fetch_coords(body):
    regions = body.get("regions")
    if not isinstance(regions, str) or not regions.strip():
        raise BadRequest("no regions provided (e.g. chr13:33,016,423-33,066,143)")
    strand = body.get("strand", "+")
    if strand not in ("+", "-"):
        raise BadRequest("strand must be '+' or '-'")
    org = body.get("organism", "")
    custom = body.get("customQuery", "")
    if not isinstance(org, str) or not isinstance(custom, str):
        raise BadRequest("organism / customQuery must be strings")
    try:
        if org in fetch.COORD_ASSEMBLIES:
            a = fetch.COORD_ASSEMBLIES[org]
            meta = fetch.retrieve_coords(regions, a["assemblyAccession"], a["assemblyName"], org,
                                         taxid=a.get("taxid", ""), strand=strand, refresh=bool(body.get("refresh")))
        elif custom.strip():
            asm = fetch.resolve_assembly(custom)          # organism/taxon name or GCF/GCA accession
            meta = fetch.retrieve_coords(regions, asm["assemblyAccession"], asm["assemblyName"],
                                         asm["organism"], taxid=asm["taxid"], strand=strand,
                                         refresh=bool(body.get("refresh")))
        else:
            raise fetch.CoordError("select an organism or enter a custom organism / assembly accession")
        return {"ok": True, **meta}
    except fetch.FetchError as fe:
        return {"ok": False, "error": str(fe)}


# ---------- calibrated ETA ----------
def run_eta(body):
    return timing.estimate(str(body.get("job", "")), int(_num(body.get("size"), 0)))


# ---------- Dfam / RepeatMasker family annotation (WSL) ----------
def run_annotate(body):
    recs = sequtil.parse_fasta(body.get("sequence", ""))
    if not recs:
        raise BadRequest("no sequence")
    seq = recs[0][1]
    _require_nt(seq)
    sp = body.get("species")
    if sp is not None and not isinstance(sp, str):
        raise BadRequest("species must be a string")
    t0 = time.time()
    r = wsl.annotate(f">{recs[0][0]}\n{seq}", species=(sp or None),
                     timeout=int(_num(body.get("timeout"), 600)))
    if r.get("ok"):
        el = time.time() - t0
        timing.record("annotate", len(seq), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("annotate", fetched=bool(body.get("source")))
        r["references"] = references
        r["provenance"] = provenance.build_manifest(
            "annotate", seq, recs[0][0],
            {"engine": "RMBLAST", "mode": "-qq", "species": r.get("species")},
            not_run=["External NCBI Primer-BLAST", "De novo family discovery"],
            databases=[{"name": "Dfam (curated)", "version": r.get("dfam_version") or "4.0"},
                       {"name": "RepeatMasker", "version": r.get("repeatmasker_version")}],
            source=body.get("source"), references=references)
    return r


# ---------- de-novo splice detection (WSL / minimap2) ----------
def run_splice(body):
    grecs = sequtil.parse_fasta(body.get("sequence", ""))
    trecs = sequtil.parse_fasta(body.get("transcript", ""))
    if not grecs:
        raise BadRequest("no genomic sequence")
    if not trecs:
        raise BadRequest("no transcript sequence to align")
    _require_nt(grecs[0][1], "genomic sequence")
    _require_nt(trecs[0][1], "transcript")
    t0 = time.time()
    r = wsl.splice_align(f">{grecs[0][0]}\n{grecs[0][1]}", f">{trecs[0][0]}\n{trecs[0][1]}",
                         timeout=int(_num(body.get("timeout"), 180)))
    if r.get("ok"):
        el = time.time() - t0
        timing.record("splice", len(grecs[0][1]), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("splice", fetched=bool(body.get("source")))
        r["references"] = references
        r["provenance"] = provenance.build_manifest(
            "splice", grecs[0][1], grecs[0][0],
            {"tool": "minimap2", "preset": "splice", "minimap2": r.get("minimap2_version"),
             "transcript_sha256": hashlib.sha256(trecs[0][1].encode()).hexdigest()},   # seal the transcript too
            source=body.get("source"), references=references)
    return r


# ---------- homology-based coding/intron recovery (WSL / miniprot) ----------
def run_miniprot(body):
    grecs = sequtil.parse_fasta(body.get("sequence", ""))
    if not grecs:
        raise BadRequest("no genomic sequence")
    _require_nt(grecs[0][1], "genomic sequence")
    prot_text = body.get("protein", "")
    if not isinstance(prot_text, str):
        raise BadRequest("protein must be a string of FASTA / raw amino-acid sequence")
    precs = sequtil.parse_protein(prot_text) if prot_text.strip() else []
    if not precs:                                    # library-default source is wired in a later step
        raise BadRequest("no reference protein supplied — paste or fetch a related TE protein "
                         "(Gag-Pol, transposase, reverse transcriptase, integrase)")
    for pid, pseq in precs:
        _require_aa(pseq, f"reference protein '{pid}'")
    protein_fasta = "\n".join(f">{pid}\n{pseq}" for pid, pseq in precs)
    # hash the FASTA text (headers + record boundaries) so IDs and multi-record splits seal distinctly
    protein_sha = hashlib.sha256(protein_fasta.encode()).hexdigest()
    t0 = time.time()
    r = wsl.protein_align(f">{grecs[0][0]}\n{grecs[0][1]}", protein_fasta,
                          timeout=int(_num(body.get("timeout"), 180)))
    if r.get("ok"):
        el = time.time() - t0
        timing.record("homology", len(grecs[0][1]), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("homology", fetched=bool(body.get("source")))
        r["references"] = references
        # the reference protein IS the evidence: its sha256 must seal the run (different protein -> different result)
        r["provenance"] = provenance.build_manifest(
            "homology", grecs[0][1], grecs[0][0],
            {"tool": "miniprot", "mode": "--gff", "miniprot": r.get("miniprot_version"),
             "protein_source": "user", "protein_sha256": protein_sha, "n_reference_proteins": len(precs)},
            not_run=["De novo ab-initio gene prediction (no transcript, no protein)",
                     "Curated TE-protein library (bundled default)"],
            source=body.get("source"), references=references)
    return r


# ---------- primer design (Primer3) ----------
def run_primers(body):
    recs = sequtil.parse_fasta(body.get("sequence", ""))
    if not recs:
        raise BadRequest("no sequence")
    seq = recs[0][1]
    _require_nt(seq)
    params = _clean_params(body.get("params", {}))
    inc = body.get("included")
    if inc is not None:                              # clamp the domain region to the template (never past the end)
        if not (isinstance(inc, (list, tuple)) and len(inc) == 2):
            raise BadRequest("included must be a [start, length] pair")
        s0 = max(0, min(int(_num(inc[0], 0)), len(seq) - 1))
        inc = [s0, max(1, min(int(_num(inc[1], 1)), len(seq) - s0))]
    try:
        res = primers.design_primers(seq, params, body.get("target"), inc)
    except (ImportError, RuntimeError):                # environment fault (Primer3 unavailable) -> 500, not a user 400
        raise
    except Exception as e:
        msg = str(e)
        if "PRIMER_PRODUCT_SIZE_RANGE" in msg or "SEQUENCE_INCLUDED_REGION" in msg:
            raise BadRequest("requested product size is larger than the sequence — lower the product-size range")
        raise BadRequest(f"primer parameters rejected: {type(e).__name__}: {e}")
    for c in res.get("candidates", []):               # secondary-structure QC (hairpin/dimer/3'-end, dual-engine) per pair
        try:
            c["qc"] = oligoqc.qc_pair(c["left_seq"], c["right_seq"])
        except Exception as e:                        # a QC fault must never drop the designed pair
            c["qc"] = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    res["oligoqc_engines"] = oligoqc.available()
    res["oligoqc_references"] = refs.oligoqc_refs()    # advisory-method citations — reported, NOT sealed
    references = refs.for_run("primer")
    res["references"] = references
    seal_params = dict(params)                        # the target/included region drives the design -> must be in the seal
    if body.get("target") is not None:
        seal_params["target"] = body.get("target")
    if inc is not None:
        seal_params["included"] = inc
    res["provenance"] = provenance.build_manifest("primer", seq, recs[0][0], seal_params, references=references)
    return res


# ---------- in-silico PCR ----------
def run_pcr(body):
    recs = sequtil.parse_fasta(body.get("sequence", ""))
    if not recs:
        raise BadRequest("no sequence")
    seq = recs[0][1]
    fwd, rev = body.get("fwd"), body.get("rev")
    if not fwd or not rev:
        raise BadRequest("fwd and rev primer sequences are required")
    if not isinstance(fwd, str) or not isinstance(rev, str):
        raise BadRequest("fwd and rev primers must be strings")
    _require_nt(seq)
    _require_nt(fwd, "forward primer")     # a non-nucleotide primer must error, not scan as empty
    _require_nt(rev, "reverse primer")
    p = _clean_params(body.get("params", {}))
    backgrounds = [{"id": recs[0][0] + " (input)", "seq": seq}]
    bg_raw = body.get("background", "")
    if not isinstance(bg_raw, str):
        raise BadRequest("background must be a string of FASTA / raw sequence")
    bg_text = bg_raw.strip()
    if bg_text:
        for bid, bseq in sequtil.parse_fasta(bg_text):
            _require_nt(bseq, f"background '{bid}'")
            backgrounds.append({"id": bid + " (custom background)", "seq": bseq})
    ts = body.get("target_span")
    if ts is not None and not (isinstance(ts, (list, tuple)) and len(ts) == 2
            and all(isinstance(x, (int, float)) and not isinstance(x, bool)
                    and x == x and abs(x) != float("inf") for x in ts)):    # reject NaN/Inf (would seal into the manifest)
        raise BadRequest("target_span must be a [start, end] pair of finite numbers")
    tp = max(0, int(_num(p.get("tp", 5), 5)))            # non-negative int; float/NaN already rejected by _clean_params
    mm = max(0, int(_num(p.get("max_mm", 2), 2)))
    try:
        amps = []
        for k, bg in enumerate(backgrounds):
            amps += primers.in_silico_pcr(
                fwd, rev, bg["seq"], bg["id"],
                max_mm=mm, tp=tp,
                prod_min=int(_num(p.get("prod_min", 70), 70)), prod_max=int(_num(p.get("prod_max", 1000), 1000)),
                target_span=(list(ts) if (ts and k == 0) else None))  # on-target only on the input sequence
    except (ValueError, TypeError) as e:
        raise BadRequest(f"in-silico PCR parameters rejected: {type(e).__name__}: {e}")
    seal_p = dict(p)                                     # primers/target/background drive the result -> must be in the seal
    seal_p.update({"fwd": fwd, "rev": rev, "target_span": (list(ts) if ts else None),
                   "background_sha256": (hashlib.sha256(bg_text.encode()).hexdigest() if bg_text else None)})
    return {
        "amplicons": amps,
        "backgroundsSearched": [bg["id"] for bg in backgrounds],
        "notSearched": ["Reference assembly (needs WSL/NCBI)", "External NCBI Primer-BLAST",
                        "IUPAC-ambiguous template bases"],
        "criteria": {"max_mismatch": mm, "three_prime_strict": tp,
                     "product_size": [int(_num(p.get("prod_min", 70), 70)), int(_num(p.get("prod_max", 1000), 1000))]},
        "references": refs.for_run("in-silico-pcr"),
        "provenance": provenance.build_manifest("in-silico-pcr", seq, recs[0][0], seal_p,
                      references=refs.for_run("in-silico-pcr"))}


_PRIMER_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")


def _resolve_assembly(body):
    """Resolve the organism (or a directly supplied accession) to a pinned RefSeq assembly. Returns
    (organism, accession, assembly_name, taxid) or raises BadRequest."""
    org = body.get("organism", "")
    if not isinstance(org, str):
        raise BadRequest("organism must be a string")
    acc = str(body.get("assemblyAccession", "") or (fetch.COORD_ASSEMBLIES.get(org, {}) or {}).get("assemblyAccession", "")).strip()
    if not re.match(r"^GC[AF]_\d+\.\d+$", acc):
        raise BadRequest("select an organism (or provide a RefSeq assembly accession GCF_…) for the whole-genome scan")
    info = fetch.COORD_ASSEMBLIES.get(org, {}) or {}
    name = str(body.get("assemblyName", "") or info.get("assemblyName", "")).strip()
    taxid = str(body.get("taxid", "") or info.get("taxid", "")).strip()
    return org, acc, name, taxid


def run_genome_prepare(body):
    """One-time download + cache of an organism's RefSeq assembly for LOCAL whole-genome scanning (no
    remote query at scan time). Slow for large genomes; runs off the UI thread. Idempotent."""
    org, acc, name, taxid = _resolve_assembly(body)
    r = wsl.genome_prepare(acc, name, timeout=int(_num(body.get("timeout"), 3600)))
    if not r.get("ok"):
        return r
    return {"ok": True, "organism": org, "taxid": taxid, "assemblyAccession": acc, "assemblyName": name, **r}


def run_genome_prepare_log(body):
    """Latest genome-download milestone (for the UI liveness indicator during a long prepare)."""
    return {"log": wsl.genome_prepare_log(int(_num((body or {}).get("tail"), 1)))}


def run_genome_list(body):
    """List locally cached genomes (accession, size, contig count) for the genome manager."""
    return wsl.genome_list()


def run_genome_remove(body):
    """Delete a cached genome to reclaim disk."""
    acc = str(body.get("assemblyAccession", "")).strip()
    return wsl.genome_remove(acc)


def run_genome_pcr(body):
    """Whole-genome off-target scan of one primer pair against a LOCAL, downloaded RefSeq assembly (isPcr).
    No remote query; a local safety timeout applies (body.timeout, default 600s — a very large genome may need
    it raised). The assembly is a fixed checksummed file, so the result IS reproducible and sealed (accession +
    source-FASTA sha256 + isPcr version + params). Every genome-wide product is an off-target candidate — a
    specificity screen, not validated bands."""
    fwd, rev = body.get("fwd"), body.get("rev")
    if not isinstance(fwd, str) or not isinstance(rev, str):
        raise BadRequest("a forward and reverse primer are required for the genome scan")
    fwd, rev = fwd.strip().upper(), rev.strip().upper()
    if not _PRIMER_RE.match(fwd) or not _PRIMER_RE.match(rev):
        raise BadRequest("primers must be DNA (A/C/G/T plus IUPAC ambiguity codes)")
    org, acc, name, taxid = _resolve_assembly(body)
    p = _clean_params(body.get("params", {}))
    min_perfect = int(_num(p.get("min_perfect"), 15))
    min_good = int(_num(p.get("min_good"), 15))
    max_size = int(_num(p.get("prod_max"), 4000))
    min_size = int(_num(p.get("prod_min"), 0))
    r = wsl.genome_scan(acc, genomepcr.query_rows(fwd, rev), max_size=max_size, min_size=min_size,
                        min_perfect=min_perfect, min_good=min_good, timeout=int(_num(body.get("timeout"), 600)))
    if not r.get("ok"):
        return r                                     # {ok:False, error, [need_prepare]} — surfaced to the UI
    amps = genomepcr.parse_ispcr(r["raw"])
    if min_size > 0:                                 # isPcr has no -minSize; apply the lower product bound here
        amps = [a for a in amps if a["length"] >= min_size]
    seal_params = {"fwd": fwd, "rev": rev, "engine": "isPcr", "isPcr_version": r.get("isPcr_version"),
                   "min_perfect": min_perfect, "min_good": min_good, "max_size": max_size, "min_size": min_size}
    # Seal the genome by its accession (with version) + source-FASTA sha256 ONLY. The human-readable
    # assemblyName is a label (like input.assemblyName in the seal-exclude list) — folding it into the sealed
    # database entry would make the same genome seal differently on the dropdown path (name populated) vs a
    # bare-accession call (name empty), so the db name is a stable, accession-independent constant.
    manifest = provenance.build_manifest(
        "genome-scan", f"{fwd}/{rev}", f"{org or acc} whole-genome scan", seal_params,
        databases=[{"name": "RefSeq genome assembly", "version": acc, "sha256": r.get("sha256")}],
        references=refs.for_run("genome-scan"),
        not_run=["Wet-lab validation of predicted products", "Off-target scan of unlisted organisms"])
    # on-target from the DESIGN LOCUS: if the specimen sits at a known position in this assembly, the product
    # overlapping it is the on-target and the rest are off-target paralogs; with no locus there is no single
    # intended target, so the products are neutral 'genomic priming sites'. This labeling is DERIVED, not sealed
    # (like the summary): the isPcr product set — the sealed content — is identical with or without a design locus.
    raw_loci = body.get("design_locus")
    raw_loci = raw_loci if isinstance(raw_loci, list) else ([raw_loci] if isinstance(raw_loci, dict) else [])
    loci = [{"acc": str(L["accession"]),
             "lo": min(int(_num(L.get("start"), 0)), int(_num(L.get("stop"), 0))),
             "hi": max(int(_num(L.get("start"), 0)), int(_num(L.get("stop"), 0)))}
            for L in raw_loci if isinstance(L, dict) and L.get("accession")]
    has_locus = bool(loci)
    n_on = 0
    if has_locus:
        def _overlap(a):                                 # best 1-based-inclusive overlap of product a with any region
            best = 0
            for L in loci:
                if a["source"] == L["acc"]:
                    best = max(best, min(a["end"], L["hi"]) - max(a["start"], L["lo"]) + 1)
            return best
        overlapping = [(a, _overlap(a)) for a in amps if not a.get("single_primer")]
        overlapping = [(a, ov) for a, ov in overlapping if ov > 0]
        if overlapping:                                  # exactly ONE on-target = the best-overlapping product (never
            best = max(overlapping, key=lambda t: (t[1], -t[0]["start"]))[0]   # every overlapper, or a clustered paralog
            best["on_target"] = True                     # inside the window would falsely read 'copy-specific')
            n_on = 1
    # interpretation (on/off split, per-chromosome spread, verdict). DERIVED, NOT sealed — never passed to
    # build_manifest, so the same primers+genome seal identically regardless of the design locus or the products.
    summary = genomepcr.summarize(amps, has_locus=has_locus)
    return {"ok": True, "organism": org, "taxid": taxid, "scope": "whole-genome", "assemblyAccession": acc,
            "assemblyName": name, "amplicons": amps, "n_amplicons": len(amps), "n_seqs": r.get("n_seqs"),
            "summary": summary, "has_locus": has_locus, "n_on_target": n_on,
            "advisory_note": ("candidate priming sites (isPcr, ≥15 bp 3'-perfect match) — not wet-lab-validated "
                              "amplicons; products priming off a shorter 3' match are not flagged"),
            "references": refs.for_run("genome-scan"), "provenance": manifest}
