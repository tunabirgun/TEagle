"""TEagle engine adapter — the single source of truth for every operation.

Each `run_*` function takes an already-parsed request dict and returns a result dict,
or raises `BadRequest` for a malformed client request (never an unhandled fault for user error).
The native PySide6 app (app/native) calls these functions through its engine worker, so request
validation and scientific behaviour live in one place.

Error taxonomy (consumed by the Qt worker):
  - returns dict              -> success
  - raises BadRequest         -> user-correctable input error (inline message)
  - raises anything else      -> unexpected fault (error banner + traceback to log)
"""
from __future__ import annotations
import copy, hashlib, re, time

from teagle_core import (sequtil, structural, primers, provenance, fetch,
                         domains, classify, refs, wsl, timing, genomepcr, oligoqc, retroviral, gff3)
import envcheck


class BadRequest(Exception):
    """A malformed client request, surfaced as an inline user message and never as an unhandled fault."""


ANALYSIS_RESULT_SCHEMA = "1.1"
# 1.0 documents are validated under the 4.0.0 rules: no domains "not_run" state and no invalidReason field.
ANALYSIS_LEGACY_SCHEMA = "1.0"
ANALYSIS_ACCEPTED_SCHEMAS = frozenset({ANALYSIS_LEGACY_SCHEMA, ANALYSIS_RESULT_SCHEMA})
ANALYSIS_RESULT_FORMAT = "TEagle analysis result"
ANALYSIS_STATUS_COMPLETE = "complete"
ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS = "complete_with_invalid_records"
ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES = "complete_with_stage_failures"
ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES_AND_INVALID_RECORDS = "complete_with_stage_failures_and_invalid_records"
ANALYSIS_STATUS_NO_INPUT = "no_input"
ANALYSIS_STATUS_NO_VALID_RECORDS = "no_valid_records"
ANALYSIS_STATUSES = frozenset({
    ANALYSIS_STATUS_COMPLETE,
    ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS,
    ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES,
    ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES_AND_INVALID_RECORDS,
    ANALYSIS_STATUS_NO_INPUT,
    ANALYSIS_STATUS_NO_VALID_RECORDS,
})
ANALYSIS_EXPORTABLE_STATUSES = frozenset({
    ANALYSIS_STATUS_COMPLETE,
    ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS,
    ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES,
    ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES_AND_INVALID_RECORDS,
})
ANALYSIS_STAGE_NAMES = ("structural", "orfs", "domains", "classification", "retroviral")


# Documented single-locus scope: above this length the structural/ORF heuristics are flagged as approximate.
LARGE_INPUT_BP = 200_000


def derive_analysis_status(records: list) -> str:
    """Derive the only valid top-level status from record validity and structured stage outcomes."""
    if not records:
        return ANALYSIS_STATUS_NO_INPUT
    valid = [rec for rec in records if rec.get("valid") is True]
    if not valid:
        return ANALYSIS_STATUS_NO_VALID_RECORDS
    invalid_present = len(valid) != len(records)
    stage_failure = any(rec.get("processing", {}).get(stage, {}).get("status") == "failed"
                        for rec in valid for stage in ANALYSIS_STAGE_NAMES)
    if stage_failure and invalid_present:
        return ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES_AND_INVALID_RECORDS
    if stage_failure:
        return ANALYSIS_STATUS_COMPLETE_WITH_STAGE_FAILURES
    if invalid_present:
        return ANALYSIS_STATUS_COMPLETE_WITH_INVALID_RECORDS
    return ANALYSIS_STATUS_COMPLETE


def stage_summary(records: list) -> dict:
    """Records per stage state across the valid records, so a run with partial or withheld stages is not
    mistaken for a fully searched one by a reader who filters on the top-level status alone."""
    summary = {stage: {} for stage in ANALYSIS_STAGE_NAMES}
    for rec in records:
        if rec.get("valid") is True:
            for stage in ANALYSIS_STAGE_NAMES:
                state = rec.get("processing", {}).get(stage, {}).get("status")
                summary[stage][state] = summary[stage].get(state, 0) + 1
    return summary


def _nonnegative_int(value) -> bool:
    return not isinstance(value, bool) and isinstance(value, int) and value >= 0


def _analysis_processing_valid(rec: dict, schema: str = ANALYSIS_RESULT_SCHEMA) -> bool:
    """Validate the complete per-record stage universe and each outcome's exact schema."""
    legacy = schema == ANALYSIS_LEGACY_SCHEMA
    required = {"id", "valid", "invalid", "composition", "structural", "domains", "classification",
                "orfs", "seq", "seq_preview", "notes", "processing", "identitySource"}
    if (not isinstance(rec, dict) or not required.issubset(rec) or
            not isinstance(rec["id"], str) or not isinstance(rec["seq"], str) or
            not isinstance(rec["seq_preview"], str) or
            not isinstance(rec["invalid"], list) or not isinstance(rec["composition"], dict) or
            not isinstance(rec["structural"], list) or not isinstance(rec["domains"], list) or
            not isinstance(rec["orfs"], list) or not isinstance(rec["notes"], list) or
            not isinstance(rec["processing"], dict) or
            set(rec["processing"]) != set(ANALYSIS_STAGE_NAMES)):
        return False
    identity = {"id": rec["id"]}
    if "idSha256" in rec:
        identity["idSha256"] = rec["idSha256"]
    source_kind = rec.get("identitySource")
    if source_kind == "local_input":
        digest = identity.get("idSha256")
        if (not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest) or
                identity["id"] != f"record_{digest[:12]}"):
            return False
    elif source_kind == "validated_fetch":
        if "idSha256" in identity or not gff3.is_public_record_id(identity["id"]):
            return False
    elif source_kind == "generated":
        if identity != gff3.public_record_identity(identity["id"], trust="generated"):
            return False
    else:
        return False
    processing = rec["processing"]
    reason = rec.get("invalidReason")
    if legacy and "invalidReason" in rec:
        return False
    if rec.get("valid") is False:
        return ((bool(rec["invalid"]) or (isinstance(reason, str) and bool(reason))) and not rec["structural"] and not rec["domains"] and not rec["orfs"] and
                rec["classification"] is None and rec.get("retroviral") is None and
                all(processing[stage] == {"status": "not_run"} for stage in ANALYSIS_STAGE_NAMES))
    if rec.get("valid") is not True or rec["invalid"] or reason is not None:
        return False

    structural_state = processing["structural"]
    if structural_state.get("status") == "complete":
        if set(structural_state) != {"status", "count"} or not _nonnegative_int(structural_state["count"]):
            return False
        if structural_state["count"] != len(rec["structural"]):
            return False
    elif (structural_state.get("status") != "failed" or
          set(structural_state) != {"status", "errorType"} or
          not isinstance(structural_state["errorType"], str) or not structural_state["errorType"] or
          rec["structural"]):
        return False

    orf_state = processing["orfs"]
    if orf_state.get("status") == "complete":
        if (set(orf_state) != {"status", "count"} or not _nonnegative_int(orf_state["count"]) or
                orf_state["count"] != len(rec["orfs"])):
            return False
    elif (orf_state.get("status") != "failed" or set(orf_state) != {"status", "errorType"} or
          not isinstance(orf_state["errorType"], str) or not orf_state["errorType"] or rec["orfs"]):
        return False

    domain_state = processing["domains"]
    if domain_state.get("status") in ({"complete", "partial"} if legacy else {"complete", "partial", "not_run"}):
        keys = {"status", "count", "orfsScanned", "orfsNotScanned"}
        if set(domain_state) != keys or not all(_nonnegative_int(domain_state[k]) for k in keys - {"status"}):
            return False
        if domain_state["count"] != len(rec["domains"]):
            return False
        if domain_state["status"] == "not_run":
            if domain_state["orfsScanned"] or domain_state["count"]:
                return False
        elif (domain_state["status"] == "partial") != (domain_state["orfsNotScanned"] > 0):
            return False
        if orf_state["status"] == "complete":
            expected_scanned = min(orf_state["count"], domains.MAX_ORFS_SCANNED)
            if (domain_state["orfsScanned"] != expected_scanned or
                    domain_state["orfsNotScanned"] != orf_state["count"] - expected_scanned):
                return False
    elif (domain_state.get("status") != "failed" or set(domain_state) != {"status", "errorType"} or
          not isinstance(domain_state["errorType"], str) or not domain_state["errorType"] or rec["domains"]):
        return False

    classification_state = processing["classification"]
    classification_status = classification_state.get("status")
    if classification_status in {"complete", "withheld"}:
        if (set(classification_state) != {"status"} or not isinstance(rec["classification"], dict) or
                rec["classification"].get("outcome") != classification_status):
            return False
    elif (classification_status != "failed" or set(classification_state) != {"status", "errorType"} or
          not isinstance(classification_state["errorType"], str) or not classification_state["errorType"] or
          rec["classification"] is not None):
        return False

    retroviral_state = processing["retroviral"]
    retroviral_status = retroviral_state.get("status")
    if retroviral_status == "complete":
        if (set(retroviral_state) != {"status", "detected"} or
                not isinstance(retroviral_state["detected"], bool) or
                retroviral_state["detected"] != (rec.get("retroviral") is not None) or
                classification_status != "complete"):
            return False
    elif retroviral_status == "not_run":
        if retroviral_state != {"status": "not_run"} or rec.get("retroviral") is not None or classification_status == "complete":
            return False
    elif (retroviral_status != "failed" or set(retroviral_state) != {"status", "errorType"} or
          not isinstance(retroviral_state["errorType"], str) or not retroviral_state["errorType"] or
          rec.get("retroviral") is not None or classification_status != "complete"):
        return False
    return True


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
    # 200 lines, not 60: the Dfam download reports progress every 15 s, and at 60 lines the step header
    # that names what is downloading scrolled out of the panel within a quarter of an hour.
    return {"log": wsl.install_log(200)}


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


# Seal key -> (module, function name, parameter) whose signature default is the value that function applies.
# Every numeric default of a domain or ORF function is sealed here, and of a structural function either here or
# under a derived key (_other_structural_thresholds); tests/test_every_applied_threshold_is_sealed.py checks it.
_SEALED_DEFAULTS = {
    "orf_min_aa": (sequtil, "find_orfs", "min_aa"),
    "domain_evalue_max": (domains, "scan_domains", "evalue"),
    "domain_orfs_scanned_max": (domains, "scan_domains", "max_orfs"),
    "ltr_seed_k": (structural, "find_ltr", "k"),
    "ltr_min_len": (structural, "find_ltr", "min_ltr"),
    "ltr_min_anchors": (structural, "find_ltr", "min_anchors"),
    "ltr_extend_x_drop": (structural, "_extend", "x_drop"),
    "tir_seed_k": (structural, "find_tir", "k"),
    "tir_min_len": (structural, "find_tir", "min_tir"),
    "tir_max_len": (structural, "find_tir", "max_tir"),
    "tir_min_anchors": (structural, "find_tir", "min_anchors"),
    "polya_min_run": (structural, "find_polya", "min_run"),
    "polya_signal_dse_lo": (structural, "find_polya_signal", "lo"),
    "polya_signal_dse_hi": (structural, "find_polya_signal", "hi"),
    "polya_signal_dse_min_gt_frac": (structural, "find_polya_signal", "min_gt_frac"),
    "polya_signal_dse_min_t": (structural, "find_polya_signal", "min_t"),
    "pbs_search_window": (structural, "find_pbs", "search"),
    "pbs_min_identity": (structural, "find_pbs", "min_ident"),
    "ppt_window": (structural, "find_ppt", "window"),
    "ppt_min_len": (structural, "find_ppt", "min_len"),
    "ppt_min_purine": (structural, "find_ppt", "min_purine"),
    "ppt_max_defects": (structural, "find_ppt", "max_defects"),
}
# Applied from module scope rather than a signature default: seal key -> structural attribute. These keys are
# stable names; every other numeric constant and default in structural.py is sealed below under a derived key.
_SEALED_STRUCTURAL_CONSTANTS = {
    "ltr_min_identity_pct": "MIN_LTR_IDENTITY",
    "tir_min_identity_pct": "MIN_TIR_IDENTITY",
    # MIN_TIR_STANDALONE is the length at or above which a terminal inverted repeat licenses a Class II call
    # with no transposase, and the bound the benchmark measures a chance ceiling against.
    "tir_min_standalone": "MIN_TIR_STANDALONE",
    # _LTR_SEED_WINDOW caps the region searched for terminal-repeat seeds, so on a long record it decides what
    # can be found at all.
    "ltr_seed_window": "_LTR_SEED_WINDOW",
    "tsd_min_len": "MIN_TSD",
    "tsd_max_len": "MAX_TSD",
}
# The random-sequence TIR ceiling the classification explanations quote (17 bp across 22,500 sequences); it is
# reported with the call, so it is sealed with the run. Other numeric names in classify.py are not listed here
# because none of them is applied or reported (_MIN_LINE_RECORD_BP is read by nothing).
_SEALED_CLASSIFY_CONSTANTS = {
    "classify_chance_tir_max_bp": "CHANCE_TIR_MAX_BP",
    "classify_chance_tir_sequences": "CHANCE_TIR_SEQUENCES",
}
# Numeric names in structural.py that decide nothing: an overflow-safe minimum sentinel and accumulator starts.
_NOT_THRESHOLDS = {"_NEG", ("_ltr_stats", "events"), ("_ltr_stats", "gap_cols")}


def _threshold_value(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, range):
        return [value.start, value.stop]
    if isinstance(value, tuple) and value and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value):
        return list(value)
    return None


def _other_structural_thresholds():
    """Every numeric module constant and signature default in structural.py that has no stable key above, under
    a key derived from its name. A threshold added to the detectors is therefore sealed without anyone having to
    remember to list it here."""
    import inspect
    explicit_attrs = set(_SEALED_STRUCTURAL_CONSTANTS.values())
    explicit_defaults = {(fname, name) for _module, fname, name in _SEALED_DEFAULTS.values()}
    out = {}
    for name, value in sorted(vars(structural).items()):
        if (re.fullmatch(r"_?[A-Z][A-Z0-9_]*", name) and name not in explicit_attrs and name not in _NOT_THRESHOLDS
                and _threshold_value(value) is not None):
            out["structural_" + name.lstrip("_").lower()] = _threshold_value(value)
    for fname, fn in sorted(inspect.getmembers(structural, inspect.isfunction)):
        if fn.__module__ != structural.__name__:
            continue
        for pname, param in inspect.signature(fn).parameters.items():
            if ((fname, pname) not in explicit_defaults and (fname, pname) not in _NOT_THRESHOLDS
                    and _threshold_value(param.default) is not None):
                out[f"structural_{fname.lstrip('_')}_{pname}"] = _threshold_value(param.default)
    return out


def _detector_parameters():
    """Every threshold that decides a structural or domain call, READ FROM THE CODE THAT APPLIES IT.

    The analysis manifest once sealed a single hand-written parameter while the numbers that determine what is
    found (the HMMER E-value, the LTR seed/length/anchor/identity floors, the TIR, PPT and poly-A thresholds,
    the TSD window, the ORF budget) lived only as function defaults. A manifest that omits them cannot
    reproduce the run it claims to seal.

    Values come from `inspect.signature` and from the defining module's constants, never restated here: a
    default changed in a detector changes the seal automatically."""
    import inspect

    params = {key: inspect.signature(getattr(module, fname)).parameters[name].default
              for key, (module, fname, name) in _SEALED_DEFAULTS.items()}
    params.update({key: getattr(structural, attr) for key, attr in _SEALED_STRUCTURAL_CONSTANTS.items()})
    params.update(_other_structural_thresholds())
    params.update({key: getattr(classify, attr) for key, attr in _SEALED_CLASSIFY_CONSTANTS.items()})
    # the expected TSD length (and accepted range and motif, where the table defines them) per superfamily
    # decides whether a TSD is called congruent or incongruent
    params["tsd_expectation_by_superfamily"] = {
        name: {"length": entry[0],
               **({"range": list(entry[1])} if len(entry) > 1 and isinstance(entry[1], (tuple, list)) else {}),
               **({"motif": entry[3]} if len(entry) > 3 and entry[3] else {})}
        for name, entry in structural.TSD_EXPECT.items()}
    return params


def _validated_operation_source(source, records):
    return provenance.validate_source_context(source, records)


def _operation_record_identity(source_context, record_id, index=0):
    if source_context.get("inputType") in {"fetched_accession", "fetched_coordinate"}:
        bindings = source_context.get("records") or []
        if index < len(bindings):
            return {**gff3.public_record_identity(bindings[index]["id"], trust="validated_fetch"),
                    "identitySource": "validated_fetch"}
    if record_id == "input_sequence":
        # parse_fasta's label for header-less text, which analyze() seals as a generated identity; hashing
        # the constant would give every bare specimen the same hash-derived record id
        return {**gff3.public_record_identity(record_id, trust="generated"), "identitySource": "generated"}
    identity = gff3.public_record_identity(record_id)
    return {**identity, "identitySource": "local_input"}


def _rna_flags(seq_text) -> list:
    """Per record, whether its sequence lines contain U. parse_fasta normalizes U to T, so this reads the raw
    text, splitting records as parse_fasta does: header lines never count and text before the first header
    belongs to no record."""
    if not isinstance(seq_text, str):
        return []
    if ">" not in seq_text:
        return ["U" in seq_text.upper()]
    flags = []
    for line in seq_text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        if line.startswith(">"):
            flags.append(False)
        elif flags and "U" in line.upper():
            flags[-1] = True
    return flags


_ORIENTING_DOMAINS = frozenset({"RT", "INT", "TPase"})


def _domain_strand_vote(hits) -> str | None:
    """The strand every RT, integrase and transposase hit reads on, or None when there is no such hit or
    they disagree (nested or adjacent elements, or a cross-scoring INT/TPase fold on the other strand)."""
    strands = {h.get("strand", "+") for h in hits if h.get("domain") in _ORIENTING_DOMAINS}
    return strands.pop() if len(strands) == 1 and strands <= {"+", "-"} else None


def analyze(seq_text: str, source: dict | None = None):
    recs = sequtil.parse_fasta(seq_text)
    validated_source = provenance.validate_source_context(source, recs)
    input_type = validated_source["inputType"]
    identities = []
    if input_type == "fetched_accession":
        identities = [gff3.public_record_identity(validated_source["accession"], trust="validated_fetch")]
    elif input_type == "fetched_coordinate":
        identities = [gff3.public_record_identity(
            f"{region['chrAccession']}:{region['start']}-{region['end']}", trust="validated_fetch")
            for region in validated_source["regions"]]
    else:
        generated_bare = isinstance(seq_text, str) and ">" not in seq_text and len(recs) == 1
        identities = [gff3.public_record_identity(rid, trust="generated" if generated_bare else "local_input")
                      for rid, _seq in recs]
    out = []
    for (rid, seq), identity, rna in zip(recs, identities, _rna_flags(seq_text)):
        ok, bad = sequtil.validate_iupac(seq)
        # A record with no A/C/G/T has no evidence to search; every detector would return nothing and the
        # result would read as a complete negative classification.
        reason = None
        if ok and not any(base in seq for base in "ACGT"):
            ok = False
            reason = ("empty record (a header with no sequence)" if not seq
                      else "no informative bases (only N and ambiguity codes)")
        not_run = {"status": "not_run"}
        identity_source = ("validated_fetch" if input_type != "local_input" else
                           ("generated" if "idSha256" not in identity else "local_input"))
        rec = {**identity, "identitySource": identity_source, "valid": ok,
               "invalid": [{"pos": p, "char": c} for p, c in bad],
               "composition": sequtil.composition(seq),
               "structural": [], "domains": [], "classification": None, "orfs": [],
               # each record carries its OWN full sequence (the exact coordinate system its structural /
               # domain / ORF calls index), so a per-record export or dot plot never has to re-derive it
               # from the whole-input textbox — which concatenated every record for a multi-FASTA paste.
               "seq": seq, "seq_preview": seq[:120], "notes": [],
               "processing": {stage: dict(not_run) for stage in ANALYSIS_STAGE_NAMES}}
        if reason:
            rec["invalidReason"] = reason
        if ok:
            structural_ok = True
            try:
                rec["structural"] = structural.detect_all(seq)      # isolate: a detector fault must not 500 the whole analyze
                rec["processing"]["structural"] = {"status": "complete", "count": len(rec["structural"])}
            except Exception as e:
                structural_ok = False
                rec["processing"]["structural"] = {"status": "failed", "errorType": type(e).__name__}
                rec["notes"].append(f"structural detection unavailable: {type(e).__name__}")
            # The classification needs the ORF count before it can say anything honest about what was
            # searched, so the list is built here. scan_domains still calls find_orfs itself — find_orfs is
            # a deterministic pure function of seq, so both calls return the identical list and the domain
            # hits' `orf` indices address THIS list correctly. (An earlier version of this comment claimed
            # the computation was shared; it is not, and a comment that misstates a code property is the
            # kind of thing that later gets relied on.)
            try:
                rec["orfs"] = sequtil.find_orfs(seq)
                rec["processing"]["orfs"] = {"status": "complete", "count": len(rec["orfs"])}
                unscanned = max(0, len(rec["orfs"]) - domains.MAX_ORFS_SCANNED)
            except Exception as e:
                rec["processing"]["orfs"] = {"status": "failed", "errorType": type(e).__name__}
                rec["notes"].append(f"ORF detection unavailable: {type(e).__name__}")
                unscanned = 0
            domains_ok = True
            try:
                rec["domains"] = domains.scan_domains(seq)
                scanned = min(len(rec["orfs"]), domains.MAX_ORFS_SCANNED)
                rec["processing"]["domains"] = {
                    # with no ORF to read the profile search had nothing to search: not run, not an empty result
                    "status": "partial" if unscanned else ("complete" if scanned else "not_run"),
                    "count": len(rec["domains"]),
                    "orfsScanned": scanned,
                    "orfsNotScanned": unscanned,
                }
            except Exception as e:
                domains_ok = False        # empty because nothing RAN — classify must not read that as "absent"
                rec["processing"]["domains"] = {"status": "failed", "errorType": type(e).__name__}
                rec["notes"].append(f"domain scan unavailable: {type(e).__name__}")
            if unscanned and domains_ok:
                # The profile search reads only the longest MAX_ORFS_SCANNED ORFs. Anything beyond that was
                # never looked at, so a domain "not detected" on this record is partly a statement about the
                # cap rather than about the sequence. Same not-tested-vs-not-found rule as domains_ok.
                rec["notes"].append(
                    f"protein-domain search read the {domains.MAX_ORFS_SCANNED} longest ORFs of "
                    f"{len(rec['orfs'])}; {unscanned} shorter ORF(s) were not searched, so a domain listed as "
                    "not detected was not necessarily absent from the whole sequence")
            vote = _domain_strand_vote(rec["domains"])
            ltr_row = next((e for e in rec["structural"] if e["type"].startswith("LTR (terminal")), None)
            if vote and structural_ok and (ltr_row or {}).get("cis_orientation") == "undetermined":
                try:
                    # the record does not state its orientation; the coding strand of RT/INT/TPase does
                    redone = structural.detect_all(seq, orientation=vote)
                    new_ltr = next(e for e in redone if e["type"].startswith("LTR (terminal"))
                    new_ltr["cis_orientation_source"] = "domain-strand"
                    rec["structural"] = redone
                    rec["processing"]["structural"] = {"status": "complete", "count": len(redone)}
                except Exception as e:
                    rec["notes"].append(f"orientation from the domain strand unavailable: {type(e).__name__}")
            if (rec["processing"]["structural"]["status"] == "complete" and
                    not any(e["type"].startswith(("LTR (terminal", "TIR (terminal")) for e in rec["structural"])):
                sc = structural.scan_scope(seq)
                n = sc["record_bp"]
                partial = [f"{name} seeds within the outer {sc[key]} bp of each end"
                           for name, key in (("LTR", "ltr_seed_window_bp"), ("TIR", "tir_seed_window_bp")) if 2 * sc[key] < n]
                if partial:
                    rec["notes"].append(
                        f"the record is {n} bp and no terminal repeat was found, but the terminal-repeat scans read only "
                        f"{' and '.join(partial)} and poly(A)/poly(T) tails within {sc['tail_window_bp']} bp of an end, "
                        "so a repeat or tail further in was not searched for"
                        + (f" — in a genomic window, more than {sc['tir_seed_window_bp']} bp of flank on a side hides "
                           "the element's ends from the TIR scan" if 2 * sc["tir_seed_window_bp"] < n else ""))
            try:
                rec["classification"] = classify.classify(rec["structural"], rec["domains"], seq=seq,
                                                          domains_ok=domains_ok, orfs_unscanned=unscanned,
                                                          orfs=rec["orfs"], structural_ok=structural_ok)
                outcome = ((rec["classification"] or {}).get("outcome") or
                           ("withheld" if rec["classification"] is None else "complete"))
                rec["processing"]["classification"] = {"status": outcome}
            except Exception as e:
                rec["processing"]["classification"] = {"status": "failed", "errorType": type(e).__name__}
                rec["notes"].append(f"classification unavailable: {type(e).__name__}")
            # retroviral transcript architecture (ERV only) — the correct coding-organisation model in place
            # of a host exon-intron gene model; None for non-ERVs, so it never overrides a normal record.
            rec["retroviral"] = None
            if rec["processing"]["classification"]["status"] == "complete":
                try:
                    rec["retroviral"] = retroviral.transcript_architecture(
                        rec["structural"], rec["domains"], rec["classification"], len(seq))
                    rec["processing"]["retroviral"] = {
                        "status": "complete", "detected": rec["retroviral"] is not None}
                except Exception as e:                          # a retroviral-model fault must not fail the analysis
                    rec["processing"]["retroviral"] = {"status": "failed", "errorType": type(e).__name__}
                    rec["notes"].append(f"retroviral architecture unavailable: {type(e).__name__}")
            if rna:
                rec["notes"].append("RNA input — U was read as T (DNA equivalent) for all steps")
            # Only an open-ended ORF that hosts a detected domain says anything about truncation: a 40-aa ATG..end
            # stretch in any of six frames is common by chance. find_orfs leaves a forward ORF open at the
            # sequence's 3′ end and a reverse-strand ORF open at its 5′ end.
            hosts = {d.get("orf") for d in rec["domains"]}
            ends = sorted({"3′" if o["strand"] == "+" else "5′"
                           for i, o in enumerate(rec["orfs"]) if i in hosts and o.get("open_end")})
            if ends:
                rec["notes"].append(
                    f"an ORF carrying a detected domain runs off the {' and '.join(ends)} end"
                    f"{'s' if len(ends) > 1 else ''} of the supplied sequence — the element may be truncated "
                    "there, and domain calls on that ORF may be partial")
            if len(seq) > LARGE_INPUT_BP:
                rec["notes"].append("large input — structural/ORF heuristics are tuned for single loci; treat as approximate")
        out.append(rec)
    joined = "".join(s for _, s in recs)
    # HMMER and Pfam produced a result wherever the profile search read at least one ORF, including a negative
    any_searched = any(r["processing"]["domains"].get("orfsScanned") for r in out)
    references = refs.for_run("analysis", domains=any_searched, fetched=input_type != "local_input")
    dbs = ([{"name": "Pfam TE-domain profiles (bundled, CC0)", "file": "data/te_domains.hmm",
             "sha256": domains.HMM_SHA256}] if domains.HMM_SHA256 else [])
    result = {"resultSchemaVersion": ANALYSIS_RESULT_SCHEMA, "records": out, "references": references}
    result["status"] = derive_analysis_status(out)
    result["stageSummary"] = stage_summary(out)
    unusable = [f"record {i} ({r['invalidReason']})" for i, r in enumerate(out, 1) if r.get("invalidReason")]
    skipped = f" Not analysed: {'; '.join(unusable)}." if unusable else ""
    if not recs:
        result["warning"] = "No sequence provided — paste, upload, or fetch a sequence first."
    elif result["status"] == ANALYSIS_STATUS_NO_VALID_RECORDS:
        result["warning"] = "No valid nucleotide record was available for analysis." + skipped
    elif len(recs) > 1:
        valid_n = sum(rec.get("valid") is True for rec in out)
        result["warning"] = (f"{len(recs)} records found and {valid_n} were valid for analysis. The detail view and "
                             f"every downstream step (primers, in-silico PCR, family naming, splice) act on the "
                             f"ONE record selected in the summary table — the first valid record by default." + skipped)
    preamble = len(sequtil.fasta_preamble(seq_text))
    if preamble:
        lead = (f"{preamble} characters before the first FASTA header belong to no record and were not analysed."
                if preamble > 1 else
                "1 character before the first FASTA header belongs to no record and was not analysed.")
        result["warning"] = f"{result['warning']} {lead}" if result.get("warning") else lead
    result["resultPayloadSha256"] = scientific_result_sha256(result)
    manifest_identity = (gff3.public_record_identity("multiple_records", trust="generated") if len(recs) > 1 else
                         (identities[0] if identities else gff3.public_record_identity("record", trust="generated")))
    manifest_identity["identitySource"] = (
        "generated" if len(recs) != 1 else out[0]["identitySource"])
    result["provenance"] = provenance.build_manifest(
        "analysis", joined, manifest_identity["id"], _detector_parameters(),
        source=validated_source, references=references, databases=dbs, input_identity=manifest_identity,
        input_records=[{**identity, "identitySource": rec["identitySource"], "sequence": seq}
                       for identity, rec, (_rid, seq) in zip(identities, out, recs)],
        result_descriptor={"resultSchemaVersion": ANALYSIS_RESULT_SCHEMA,
                           "resultPayloadSha256": result["resultPayloadSha256"]})
    return result


def _scientific_result_payload(result: dict) -> dict:
    payload = copy.deepcopy(result)
    payload.pop("resultPayloadSha256", None)
    payload.pop("format", None)
    payload.pop("provenance", None)
    return payload


def scientific_result_sha256(result: dict) -> str:
    """Hash the canonical scientific core without provenance or this digest (non-circular by design)."""
    import json
    encoded = json.dumps(_scientific_result_payload(result), ensure_ascii=False, allow_nan=False,
                         sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def verify_scientific_result_payload(result: dict) -> bool:
    try:
        claimed = result.get("resultPayloadSha256")
        return (isinstance(claimed, str) and len(claimed) == 64 and
                scientific_result_sha256(result) == claimed.lower())
    except (AttributeError, TypeError, ValueError):
        return False


def analysis_result_document(result: dict) -> dict:
    """Build the versioned, complete JSON document exported for one analysis run."""
    if not isinstance(result, dict) or result.get("resultSchemaVersion") not in ANALYSIS_ACCEPTED_SCHEMAS:
        raise ValueError("unsupported analysis result")
    if result.get("status") not in ANALYSIS_EXPORTABLE_STATUSES:
        raise ValueError("analysis has no exportable valid result")
    document = copy.deepcopy(result)
    document["format"] = ANALYSIS_RESULT_FORMAT
    if not validate_analysis_result_document(document):
        raise ValueError("analysis result is not bound to its provenance manifest")
    return document


def analysis_result_json_bytes(result: dict) -> bytes:
    """Serialize an exportable result as strict UTF-8 JSON and reverify its bindings."""
    import json

    document = analysis_result_document(result)
    encoded = json.dumps(
        document, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2,
    ).encode("utf-8")
    decoded = json.loads(encoded)
    if not validate_analysis_result_document(decoded):
        raise ValueError("serialized analysis result failed seal verification")
    return encoded


def validate_analysis_result_document(document: dict) -> bool:
    """Check the sealed result descriptor, scientific core, stages, and ordered input binding."""
    try:
        manifest = document.get("provenance")
        if (document.get("format") != ANALYSIS_RESULT_FORMAT or
                document.get("resultSchemaVersion") not in ANALYSIS_ACCEPTED_SCHEMAS or
                document.get("status") not in ANALYSIS_EXPORTABLE_STATUSES or
                not provenance.verify_manifest(manifest)):
            return False
        descriptor = {"resultSchemaVersion": document["resultSchemaVersion"],
                      "resultPayloadSha256": document.get("resultPayloadSha256")}
        if manifest.get("resultDescriptor") != descriptor or not verify_scientific_result_payload(document):
            return False
        records = document.get("records")
        bound = manifest["input"].get("records")
        if (not isinstance(records, list) or not all(_analysis_processing_valid(rec, document["resultSchemaVersion"]) for rec in records) or
                not isinstance(bound, list) or len(records) != len(bound)):
            return False
        observed = [{**({"id": str(rec["id"]), "identitySource": rec["identitySource"]} |
                        ({"idSha256": rec["idSha256"]} if rec.get("idSha256") else {})),
                     "length": len(rec["seq"]),
                     "sha256": hashlib.sha256(rec["seq"].encode("utf-8")).hexdigest()}
                    for rec in records]
        return (observed == bound and document["status"] == derive_analysis_status(records) and
                document.get("stageSummary", stage_summary(records)) == stage_summary(records))
    except (AttributeError, KeyError, TypeError, ValueError):
        return False


# ---------- accession fetch ----------
def run_fetch(body):
    acc = body.get("accession", "")
    if not isinstance(acc, str):
        raise BadRequest("accession must be a string")
    try:
        meta = fetch.retrieve(acc, refresh=bool(body.get("refresh")))
        meta["sourceContext"] = provenance.normalize_source_context(
            meta, sequtil.parse_fasta(meta.get("fasta") or ""))
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
        asm = fetch.all_assemblies()                          # curated + user-added, so a custom organism resolves here too
        if org in asm:
            a = asm[org]
            meta = fetch.retrieve_coords(regions, a["assemblyAccession"], a["assemblyName"], org,
                                         taxid=a.get("taxid", ""), strand=strand, refresh=bool(body.get("refresh")))
        elif custom.strip():
            asm = fetch.resolve_assembly(custom)          # organism/taxon name or GCF/GCA accession
            meta = fetch.retrieve_coords(regions, asm["assemblyAccession"], asm["assemblyName"],
                                         asm["organism"], taxid=asm["taxid"], strand=strand,
                                         refresh=bool(body.get("refresh")))
        else:
            raise fetch.CoordError("select an organism or enter a custom organism / assembly accession")
        meta["sourceContext"] = provenance.normalize_source_context(
            meta.get("source"), sequtil.parse_fasta(meta.get("fasta") or ""))
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
    source_context = _validated_operation_source(body.get("source"), recs)
    identity = _operation_record_identity(source_context, recs[0][0])
    sp = body.get("species")
    if sp is not None and not isinstance(sp, str):
        raise BadRequest("species must be a string")
    unc = bool(body.get("include_uncurated"))
    t0 = time.time()
    r = wsl.annotate(f">{identity['id']}\n{seq}", species=(sp or None),
                     timeout=int(_num(body.get("timeout"), 600)), include_uncurated=unc)
    if r.get("ok"):
        el = time.time() - t0
        timing.record("annotate", len(seq), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("annotate", fetched=source_context["inputType"] != "local_input")
        r["references"] = references
        # The Dfam entry is whatever famdb reported for the library actually installed — never a literal.
        # A library that cannot be identified seals an explicit marker instead of a version, because a run
        # whose database is unknown must not carry a version claim. The installed PARTITIONS are sealed with
        # it: the same Dfam version with a different partition set searches a different family universe.
        lib = r.get("dfam_library") or {}
        # the name asserts no partition qualifier — the exact installed partition set is sealed below in
        # `partitions`, and a hardcoded "curated consensus" would contradict a library with uncurated
        # partitions installed (the same Dfam version searches a different family universe).
        dfam_db = {"name": "Dfam (via famdb)",
                   "version": lib.get("version") or "unresolved — famdb reported no library version"}
        for src_k, dst_k in (("date", "releaseDate"), ("famdbFormat", "famdbFormat"),
                             ("partitions", "partitions"), ("consensusSequences", "consensusSequences")):
            if lib.get(src_k) is not None:
                dfam_db[dst_k] = lib[src_k]
        if not lib.get("version"):
            r["warning"] = ("The Dfam library version could not be read from famdb, so this run's record "
                            "identifies the database as unresolved rather than asserting a version. The "
                            "family calls are unaffected; the provenance record is weaker.")
        r["provenance"] = provenance.build_manifest(
            "annotate", seq, recs[0][0],
            # the searched family universe is a PARAMETER of the run: the same library answers a
            # different question with and without the uncurated families, so it is sealed, not implied
            {"engine": "RMBLAST", "mode": "-qq", "species": r.get("species"),
             "library": r.get("library_kind")},
            not_run=["External NCBI Primer-BLAST", "De novo family discovery"],
            databases=[dfam_db, {"name": "RepeatMasker", "version": r.get("repeatmasker_version")}],
            source=source_context, references=references, input_identity=identity)
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
    source_context = _validated_operation_source(body.get("source"), grecs)
    identity = _operation_record_identity(source_context, grecs[0][0])
    t0 = time.time()
    r = wsl.splice_align(f">{identity['id']}\n{grecs[0][1]}", f">{trecs[0][0]}\n{trecs[0][1]}",
                         timeout=int(_num(body.get("timeout"), 180)))
    if r.get("ok"):
        el = time.time() - t0
        timing.record("splice", len(grecs[0][1]), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("splice", fetched=source_context["inputType"] != "local_input")
        r["references"] = references
        r["provenance"] = provenance.build_manifest(
            "splice", grecs[0][1], grecs[0][0],
            {"tool": "minimap2", "preset": "splice", "minimap2": r.get("minimap2_version"),
             "transcript_sha256": hashlib.sha256(trecs[0][1].encode()).hexdigest()},   # seal the transcript too
            source=source_context, references=references, input_identity=identity)
    return r


# ---------- homology-based coding/intron recovery (WSL / miniprot) ----------
def run_miniprot(body):
    grecs = sequtil.parse_fasta(body.get("sequence", ""))
    if not grecs:
        raise BadRequest("no genomic sequence")
    _require_nt(grecs[0][1], "genomic sequence")
    source_context = _validated_operation_source(body.get("source"), grecs)
    identity = _operation_record_identity(source_context, grecs[0][0])
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
    r = wsl.protein_align(f">{identity['id']}\n{grecs[0][1]}", protein_fasta,
                          timeout=int(_num(body.get("timeout"), 180)))
    if r.get("ok"):
        el = time.time() - t0
        timing.record("homology", len(grecs[0][1]), el)
        r["elapsed_s"] = round(el, 1)
        references = refs.for_run("homology", fetched=source_context["inputType"] != "local_input")
        r["references"] = references
        # the reference protein IS the evidence: its sha256 must seal the run (different protein -> different result)
        r["provenance"] = provenance.build_manifest(
            "homology", grecs[0][1], grecs[0][0],
            {"tool": "miniprot", "mode": "--gff", "miniprot": r.get("miniprot_version"),
             "protein_source": "user", "protein_sha256": protein_sha, "n_reference_proteins": len(precs)},
            not_run=["De novo ab-initio gene prediction (no transcript, no protein)",
                     "Curated TE-protein library (bundled default)"],
            source=source_context, references=references, input_identity=identity)
    return r


# ---------- primer design (Primer3) ----------
def _sequence_record_binding(record_id, sequence, identity=None):
    return {"index": 0, **(identity or gff3.public_record_identity(record_id)), "length": len(sequence),
            "sha256": hashlib.sha256(sequence.encode()).hexdigest()}


def _validate_sealed_primer_candidate(manifest, design_record, design_seq):
    """Recompute candidate geometry against the exact design record and reject internal contradictions."""
    if not provenance.verify_manifest(manifest) or manifest.get("runType") != "primer":
        raise BadRequest("design_provenance is invalid or has an unsupported seal")
    params = manifest.get("parameters") or {}
    cand = params.get("candidate")
    if not isinstance(cand, dict):
        raise BadRequest("design_provenance does not bind a primer candidate")
    top_input = manifest.get("input")
    if not isinstance(top_input, dict):
        raise BadRequest("design_provenance has no input binding")
    for field in ("id", "length", "sha256"):
        if top_input.get(field) != design_record[field]:
            raise BadRequest("design_provenance top-level input contradicts its selected record")
    if cand.get("selected_record") != design_record:
        raise BadRequest("primer candidate is bound to a different design record")
    try:
        left = primers.endpoint_span(cand.get("left_pos"), len(design_seq), "candidate left_pos")
        right = primers.endpoint_span(cand.get("right_pos"), len(design_seq), "candidate right_pos")
        intended = primers.endpoint_span(
            cand.get("intended_product_span"), len(design_seq), "candidate intended_product_span")
        target = primers.start_length_region(cand.get("primer3_target"), len(design_seq), "candidate primer3_target")
        included = primers.start_length_region(cand.get("included_region"), len(design_seq), "candidate included_region")
        tolerance = primers.endpoint_tolerance(cand.get("intended_endpoint_tolerance_bp"))
    except (TypeError, ValueError) as e:
        raise BadRequest(str(e))
    fwd, rev = cand.get("left_seq"), cand.get("right_seq")
    if not isinstance(fwd, str) or not isinstance(rev, str):
        raise BadRequest("design_provenance has invalid primer sequences")
    if (left[1] - left[0] != len(fwd) or right[1] - right[0] != len(rev) or
            left[1] > right[0]):
        raise BadRequest("sealed primer spans contradict oligo lengths or orientation")
    if design_seq[left[0]:left[1]].upper() != fwd.upper():
        raise BadRequest("sealed forward-primer position contradicts the design template")
    if sequtil.reverse_complement(design_seq[right[0]:right[1]]).upper() != rev.upper():
        raise BadRequest("sealed reverse-primer position contradicts the design template")
    derived_intended = [left[0], right[1]]
    if intended != derived_intended or cand.get("product_size") != intended[1] - intended[0]:
        raise BadRequest("sealed intended product or product size contradicts primer coordinates")
    if cand.get("intended_contributors") != {"left": "F", "right": "R", "orientation": "F+R"}:
        raise BadRequest("sealed intended primer contributors or orientation are invalid")
    if target is not None and (intended[0] > target[0] or intended[1] < target[0] + target[1]):
        raise BadRequest("sealed intended product does not span primer3_target")
    if included is not None:
        included_end = included[0] + included[1]
        if (intended[0] < included[0] or intended[1] > included_end or
                (target is not None and (target[0] < included[0] or target[0] + target[1] > included_end))):
            raise BadRequest("sealed candidate geometry leaves included_region")
    return cand, {"primer3_target": target, "included_region": included,
                  "intended_product_span": intended,
                  "intended_endpoint_tolerance_bp": tolerance}


def run_primers(body):
    if "target" in body or "included" in body:
        raise BadRequest("use primer3_target and included_region; target/included are not accepted")
    recs = sequtil.parse_fasta(body.get("sequence", ""))
    if not recs:
        raise BadRequest("no sequence")
    seq = recs[0][1]
    _require_nt(seq)
    source_context = _validated_operation_source(body.get("source"), recs)
    identity = _operation_record_identity(source_context, recs[0][0])
    params = _clean_params(body.get("params", {}))
    try:
        effective_params = primers.resolved_design_parameters(params)
        target = primers.start_length_region(body.get("primer3_target"), len(seq), "primer3_target")
        inc = primers.start_length_region(body.get("included_region"), len(seq), "included_region")
        if target is not None and inc is not None:
            if target[0] < inc[0] or target[0] + target[1] > inc[0] + inc[1]:
                raise ValueError("primer3_target must be contained in included_region")
    except (TypeError, ValueError) as e:
        raise BadRequest(str(e))
    try:
        res = primers.design_primers(seq, effective_params, target, inc)
    except (ImportError, RuntimeError):                # environment fault (Primer3 unavailable) -> 500, not a user 400
        raise
    except primers.TemplateError as e:
        raise BadRequest(str(e))
    except Exception as e:
        msg = str(e)
        if "PRIMER_PRODUCT_SIZE_RANGE" in msg or "SEQUENCE_INCLUDED_REGION" in msg:
            raise BadRequest("requested product size is larger than the sequence — lower the product-size range")
        if "invalid sequence character" in msg:
            raise BadRequest(f"the template contains a character Primer3 cannot use: {msg}")
        raise BadRequest(f"primer parameters rejected: {type(e).__name__}: {e}")
    _prime_vienna(res.get("candidates", []))          # optional out-of-process cross-check, one round trip
    for c in res.get("candidates", []):               # secondary-structure QC (hairpin/dimer/3'-end) per pair
        try:
            c["qc"] = oligoqc.qc_pair(c["left_seq"], c["right_seq"])
        except Exception as e:                        # a QC fault must never drop the designed pair
            c["qc"] = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    res["oligoqc_engines"] = oligoqc.available()
    res["design_tm_conditions"] = dict(primers.TM_CONDITIONS)     # conditions behind the Tm columns (sealed)
    res["oligoqc_conditions"] = dict(oligoqc.CONDITIONS)          # conditions behind QC Tm / ΔG (advisory, not sealed)
    ambiguity = res.setdefault("template_ambiguity", primers.template_ambiguity(0))
    res["cautions"] = ([f"{ambiguity['masked_to_n']} ambiguous template base(s) were treated as N for Primer3; "
                        "no primer covers them."] if ambiguity["masked_to_n"] else [])
    res["oligoqc_references"] = refs.oligoqc_refs()    # advisory-method citations — reported, NOT sealed
    references = refs.for_run("primer", fetched=source_context["inputType"] != "local_input")
    res["references"] = references
    seal_params = dict(effective_params)
    seal_params.update({"coordinateSystem": "0-based-half-open",
                        "primer3_target": target, "included_region": inc,
                        "tm_conditions": dict(primers.TM_CONDITIONS),
                        "max_ns_accepted": primers.MAX_NS_ACCEPTED,
                        "template_ambiguity_masked": ambiguity["masked_to_n"]})
    res["provenance"] = provenance.build_manifest(
        "primer", seq, identity["id"], seal_params, source=source_context, references=references,
        input_identity=identity)
    record_binding = _sequence_record_binding(identity["id"], seq, identity)
    for c in res.get("candidates", []):
        candidate = {
            "id": c["id"], "left_seq": c["left_seq"], "right_seq": c["right_seq"],
            "left_pos": c["left_pos"], "right_pos": c["right_pos"],
            "product_size": c["product_size"], "primer3_target": c["primer3_target"],
            "included_region": c["included_region"],
            "intended_product_span": c["intended_product_span"],
            "intended_endpoint_tolerance_bp": c["intended_endpoint_tolerance_bp"],
            "intended_contributors": {"left": "F", "right": "R", "orientation": "F+R"},
            "selected_record": record_binding,
        }
        c["design_provenance"] = provenance.build_manifest(
            "primer", seq, identity["id"], {**seal_params, "candidate": candidate},
            source=source_context, references=references, input_identity=identity)
    return res


# ---------- in-silico PCR ----------
def run_pcr(body):
    if "target_span" in body:
        raise BadRequest("target_span is ambiguous; supply primer3_target and intended_product_span separately")
    recs = sequtil.parse_fasta(body.get("sequence", ""))
    if not recs:
        raise BadRequest("no sequence")
    seq = recs[0][1]
    source_context = _validated_operation_source(body.get("source"), recs)
    fwd, rev = body.get("fwd"), body.get("rev")
    if not fwd or not rev:
        raise BadRequest("fwd and rev primer sequences are required")
    if not isinstance(fwd, str) or not isinstance(rev, str):
        raise BadRequest("fwd and rev primers must be strings")
    _require_nt(seq)
    try:                                   # whitespace/case/U normalised; a non-IUPAC primer must error, not scan as empty
        fwd = primers.normalize_primer(fwd, "forward primer")
        rev = primers.normalize_primer(rev, "reverse primer")
    except ValueError as e:
        raise BadRequest(str(e))
    p = _clean_params(body.get("params", {}))
    assay_identity = _operation_record_identity(source_context, recs[0][0])
    backgrounds = [{"id": assay_identity["id"], "source": assay_identity["id"] + " (input)", "seq": seq,
                    "idSha256": assay_identity.get("idSha256")}]
    bg_raw = body.get("background", "")
    if not isinstance(bg_raw, str):
        raise BadRequest("background must be a string of FASTA / raw sequence")
    bg_text = bg_raw.strip()
    if bg_text:
        for bid, bseq in sequtil.parse_fasta(bg_text):
            identity = gff3.public_record_identity(bid)
            _require_nt(bseq, f"background '{identity['id']}'")
            backgrounds.append({"id": identity["id"],
                                "source": identity["id"] + " (custom background)", "seq": bseq,
                                "idSha256": identity.get("idSha256")})
    design_manifest = body.get("design_provenance")
    assay_record = _sequence_record_binding(assay_identity["id"], seq, assay_identity)
    design_input = None
    if design_manifest is not None:
        design_text = body.get("design_sequence")
        if not isinstance(design_text, str):
            raise BadRequest("design_sequence is required to verify design_provenance")
        design_recs = sequtil.parse_fasta(design_text)
        if not design_recs:
            raise BadRequest("design_sequence has no sequence record")
        design_id, design_seq = design_recs[0]
        _require_nt(design_seq, "design sequence")
        if not provenance.verify_manifest(design_manifest):
            raise BadRequest("design_provenance is invalid or has an unsupported seal")
        bound_identity = design_manifest.get("input") or {}
        source_kind = bound_identity.get("identitySource")
        candidate_identity = gff3.public_record_identity(
            design_id, trust=source_kind if source_kind in {"validated_fetch", "generated"} else "local_input")
        candidate_identity["identitySource"] = source_kind or "local_input"
        design_record = _sequence_record_binding(candidate_identity["id"], design_seq, candidate_identity)
        cand, geometry = _validate_sealed_primer_candidate(
            design_manifest, design_record, design_seq)
        target, inc = geometry["primer3_target"], geometry["included_region"]
        intended = geometry["intended_product_span"]
        tolerance = geometry["intended_endpoint_tolerance_bp"]
        expected_candidate = {
            "left_seq": fwd, "right_seq": rev, "primer3_target": target,
            "included_region": inc, "intended_product_span": intended,
            "intended_endpoint_tolerance_bp": tolerance,
        }
        if any(cand.get(k) != v for k, v in expected_candidate.items()):
            raise BadRequest("primer candidate geometry or sequence disagrees with design_provenance")
        if any(body.get(k) != v for k, v in geometry.items()):
            raise BadRequest("submitted primer geometry disagrees with design_provenance")
        design_input = copy.deepcopy(design_manifest["input"])
        same_design_record = design_record == assay_record
        expected_assessment = "assess" if same_design_record else "withheld_template_mismatch"
    else:
        if body.get("design_sequence") is not None:
            raise BadRequest("design_sequence requires design_provenance")
        if any(body.get(k) is not None for k in (
                "primer3_target", "included_region", "intended_product_span")):
            raise BadRequest("design_provenance is required to assess designed-product geometry")
        target = inc = intended = None
        try:
            tolerance = primers.endpoint_tolerance(body.get(
                "intended_endpoint_tolerance_bp", primers.INTENDED_ENDPOINT_TOLERANCE_BP))
        except (TypeError, ValueError) as e:
            raise BadRequest(str(e))
        same_design_record = False
        expected_assessment = "not_run"
    design_assessment = body.get("design_assessment", expected_assessment)
    if design_assessment != expected_assessment:
        raise BadRequest("design_assessment contradicts the design and assay record bindings")

    tp = _num(p.get("tp", 5), 5)
    mm = _num(p.get("max_mm", 2), 2)
    # The shortest 3'-proximal core that may be called a binding site. It changes which products are
    # reported for a tailed primer, so it is resolved here and written into the seal below rather than
    # left as a bare default inside the matcher where no manifest would record it.
    ma = _num(p.get("min_anneal", primers.MIN_PRIMER_SIZE), primers.MIN_PRIMER_SIZE)
    pmin = _num(p.get("prod_min", 70), 70)
    pmax = _num(p.get("prod_max", 1000), 1000)
    for name, value in (("tp", tp), ("max_mm", mm), ("min_anneal", ma),
                        ("prod_min", pmin), ("prod_max", pmax)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise BadRequest(f"{name} must be an integer")
    if tp < 0 or mm < 0:
        raise BadRequest("tp and max_mm must be non-negative")
    if ma < primers.MIN_PRIMER_SIZE or len(fwd) < ma or len(rev) < ma:
        raise BadRequest(f"min_anneal must be at least {primers.MIN_PRIMER_SIZE} and no longer than either primer")
    shortest_primer = min(len(fwd), len(rev))
    if tp > shortest_primer or mm > shortest_primer:
        raise BadRequest("tp and max_mm cannot exceed the shorter primer length")
    if pmin < 1 or pmax < pmin:
        raise BadRequest("product-size bounds must be positive and prod_min <= prod_max")
    try:
        amps = []
        pcr_stats = {}
        truncated = []
        # The specimen alone fixes each primer's annealing core; a background is capped at it, so an identical
        # site gets the same call there, a background that does not reach it keeps its own lone-search core,
        # and the specimen result does not depend on the backgrounds.
        cores = (primers.specimen_cores(fwd, rev, backgrounds[0]["seq"], mm, tp, ma)
                 if len(backgrounds) > 1 else None)
        used_cores = {o: {"forward_strand": None, "reverse_strand": None} for o in "FR"}
        for k, bg in enumerate(backgrounds):
            st = {}
            found = primers.in_silico_pcr(
                fwd, rev, bg["seq"], bg["source"],
                max_mm=mm, tp=tp, min_anneal=ma, prod_min=pmin, prod_max=pmax,
                primer3_target=(target if k == 0 and same_design_record else None),
                intended_product_span=(intended if k == 0 and same_design_record else None),
                intended_endpoint_tolerance_bp=tolerance,
                assess_design_geometry=(k == 0 and same_design_record),
                stats=st, cores=(cores if k else None))
            for o, per_strand in (st.get("cores") or {}).items():
                for strand, n in per_strand.items():
                    if n is not None and n > (used_cores[o][strand] or 0):
                        used_cores[o][strand] = n
            if bg.get("idSha256"):
                for amplicon in found:
                    amplicon["sourceIdSha256"] = bg["idSha256"]
                    amplicon["product_key"] = f"{amplicon.get('product_key', '')}:idsha256={bg['idSha256']}"
            amps += found
            if st.get("amplicons_capped") or st.get("sites_capped"):
                truncated.append(bg["source"])
            pcr_stats = st or pcr_stats
    except (ValueError, TypeError) as e:
        raise BadRequest(f"in-silico PCR parameters rejected: {type(e).__name__}: {e}")
    # Seal the RESOLVED values, not the submitted ones. `dict(p)` alone records only what the caller
    # happened to type, so a run left entirely at defaults sealed no thresholds at all — the same defect
    # _detector_parameters() was written to remove from the analysis path, still present here. Every number
    # that selects which products are reported now goes in, whether the caller set it or not.
    assay_binding = dict(assay_record)
    if source_context["inputType"] != "local_input":
        assay_binding.update(provenance.public_source_fields(source_context, recs))
    search_inputs = [{"role": "assay", **assay_binding}]
    for bg in backgrounds[1:]:
        binding = {"role": "background", **_sequence_record_binding(bg["id"], bg["seq"],
                                                                       {"id": bg["id"], **({"idSha256": bg["idSha256"]}
                                                                        if bg.get("idSha256") else {})})}
        if bg.get("idSha256"):
            binding["idSha256"] = bg["idSha256"]
        search_inputs.append(binding)
    seal_p = dict(p)                                     # primers/target/background drive the result -> must be in the seal
    seal_p.update({"fwd": fwd, "rev": rev,
                   "coordinateSystem": "0-based-half-open",
                   "primer3_target": target, "included_region": inc,
                   "intended_product_span": intended,
                   "intended_endpoint_tolerance_bp": tolerance,
                   "design_assessment": design_assessment,
                   "design_manifest_sha256": (design_manifest.get("manifestSha256") if design_manifest else None),
                   "design_input": design_input,
                   "assay_input": assay_binding,
                   "search_inputs": search_inputs,
                   "max_mm": mm, "tp": tp, "min_anneal": ma,
                   "prod_min": pmin, "prod_max": pmax,
                   "background_sha256": (hashlib.sha256(bg_text.encode()).hexdigest() if bg_text else None)})
    # A capped search must say so. Both caps bite on exactly the repetitive templates where "no further
    # product" would otherwise read as a specificity result rather than as the search stopping early.
    warning = None
    if truncated:
        warning = (f"in-silico PCR stopped early on {', '.join(truncated)}: the search caps at "
                   f"{pcr_stats.get('max_amplicons', 4000)} products and {pcr_stats.get('max_sites', 4000)} "
                   "binding sites per primer. The product list below is therefore INCOMPLETE — treat it as "
                   "a lower bound on priming, not as an exhaustive product search.")
    status = "partial_capped" if truncated else ("complete" if amps else "complete_zero_products")
    unique_products = len({a["product_key"] for a in amps})
    references = refs.for_run("in-silico-pcr", fetched=source_context["inputType"] != "local_input")
    return {
        "status": status,
        "amplicons": amps,
        "mechanismCount": len(amps),
        "uniqueProductCount": unique_products,
        "truncated": bool(truncated),
        "truncatedOn": truncated,
        "warning": warning,
        "backgroundsSearched": [bg["source"] for bg in backgrounds],
        "bindings": {"design": design_input, "assay": assay_binding,
                     "backgrounds": search_inputs[1:], "designAssessment": design_assessment},
        "notSearched": ["Reference assembly (needs WSL/NCBI)", "External NCBI Primer-BLAST",
                        "IUPAC-ambiguous template bases"],
        "cautions": ([f"3′ exact bases ({tp}) exceed the shortest annealed core ({ma}): the exact window is "
                      "measured inside the annealed core, and primer bases 5′ of it are reported as a tail "
                      "that is not scored for mismatches."] if tp > ma else []),
        "criteria": {"max_mismatch": mm, "three_prime_strict": tp, "min_anneal": ma,
                     "three_prime_strict_effective": min(tp, ma),
                     "annealing_core": used_cores,
                     "product_size": [pmin, pmax], "coordinateSystem": "0-based-half-open",
                     "primer3_target": target, "included_region": inc,
                     "intended_product_span": intended,
                     "intended_endpoint_tolerance_bp": tolerance,
                     "intended_product_status": ("assessed" if same_design_record and intended is not None
                                                 else ("withheld" if design_assessment == "withheld_template_mismatch"
                                                       else "not_run")),
                     "target_coverage_status": ("assessed" if same_design_record and target is not None
                                                else ("withheld" if design_assessment == "withheld_template_mismatch"
                                                      else "not_run"))},
        "references": references,
        "provenance": provenance.build_manifest(
            "in-silico-pcr", seq, assay_identity["id"], seal_p, source=source_context,
            references=references, input_identity=assay_identity)}


def _prime_vienna(candidates):
    """Fetch every ViennaRNA fold a primer batch needs in ONE WSL round trip, when the in-process
    module is absent but the optional backend component is installed. Per-metric remote calls would be
    roughly nine round trips per pair, which is unusable for a full design. An installed in-process
    ViennaRNA always wins, and any failure here degrades to primer3-only — the cross-check is
    advisory and must never break a design."""
    if not candidates or oligoqc.RNA is not None:
        return
    try:
        from teagle_core import oligoqc_wsl
        if not oligoqc_wsl.available():
            return oligoqc.prime_remote({})
        seqs, pairs = set(), set()
        for c in candidates:
            left, right = c.get("left_seq", ""), c.get("right_seq", "")
            seqs.update(s for s in (left, right) if s)
            if left:
                pairs.add((left, left))               # self-dimer
            if right:
                pairs.add((right, right))
            if left and right:
                pairs.add((left, right))              # cross-dimer
        oligoqc.prime_remote(oligoqc_wsl.compute(seqs, pairs, oligoqc.CONDITIONS),
                             version=oligoqc_wsl.version())
    except Exception:
        oligoqc.prime_remote({})


_PRIMER_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$")


def _resolve_assembly(body):
    """Resolve the organism (or a directly supplied accession) to a pinned RefSeq assembly. Returns
    (organism, accession, assembly_name, taxid) or raises BadRequest."""
    org = body.get("organism", "")
    if not isinstance(org, str):
        raise BadRequest("organism must be a string")
    info = fetch.all_assemblies().get(org, {}) or {}          # curated + user-added, so a custom genome resolves here too
    acc = str(body.get("assemblyAccession", "") or info.get("assemblyAccession", "")).strip()
    if not re.match(r"^GC[AF]_\d+\.\d+$", acc):
        raise BadRequest("select an organism (or provide a RefSeq assembly accession GCF_…) for the whole-genome scan")
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


def run_annotate_budget(body):
    """CPU/RAM/disk available to a whole-genome annotation, plus the selected Dfam search-set count.

    A user-supplied FASTA replaces Dfam, so its budget must not query or report an unrelated Dfam lineage
    count. The request identity is echoed so the UI can discard a slower response from an earlier library
    choice instead of rendering it under the newly selected one.
    """
    b = body or {}
    size = {}
    if b.get("fasta_bytes") is not None:
        fasta_bytes = _num(b["fasta_bytes"], 0)
        cached = b.get("fasta_cached")
        if fasta_bytes < 0 or (cached is not None and not isinstance(cached, bool)):
            raise BadRequest("fasta_bytes must be non-negative and fasta_cached a boolean")
        size = {"fasta_bytes": int(fasta_bytes), "fasta_cached": cached}
    r = wsl.annotate_budget(int(_num(b.get("genome_bytes"), 0)), **size)
    mode = b.get("library_mode")
    if mode is None:
        mode = "dfam_uncurated" if bool(b.get("include_uncurated")) else "dfam_curated"
    if mode not in ("dfam_curated", "dfam_uncurated", "custom"):
        raise BadRequest("library_mode must be dfam_curated, dfam_uncurated or custom")
    r["library_mode"] = mode
    if "budget_request_id" in b:
        request_id = b.get("budget_request_id")
        if isinstance(request_id, bool) or not isinstance(request_id, int) or request_id < 1:
            raise BadRequest("budget_request_id must be a positive integer")
        r["budget_request_id"] = request_id
    sp = b.get("species")
    if mode == "dfam_uncurated" and not sp:
        r["library_ready"] = False
        r["library_error"] = "select an organism/lineage before including uncurated Dfam families"
        return r
    if mode != "custom" and sp:
        if not isinstance(sp, str):
            raise BadRequest("species must be a string")
        lib = wsl.dfam_lineage_families(sp, include_uncurated=(mode == "dfam_uncurated"))
        r["library_families_for_species"] = lib.get("families") if lib.get("ok") else None
        if mode == "dfam_uncurated":
            r["uncurated_available"] = bool(lib.get("ok") and lib.get("uncurated_available") is not False)
            if not r["uncurated_available"]:
                r["library_ready"] = False
                r["library_error"] = lib.get("error") or "the optional uncurated Dfam partitions are unavailable"
        r["species"] = sp
    return r


def run_genome_annotate(body):
    """Whole-genome TE annotation: RepeatMasker over a cached assembly, chunked and resumable.

    Homology-bound by construction: a Dfam-backed run can place only represented Dfam families, while a
    custom FASTA replaces Dfam for interspersed-repeat homology. The result identifies the selected library
    and states its coverage limit rather than treating an absent match as a fact about the genome."""
    acc = str(body.get("assemblyAccession", "")).strip()
    if not acc:
        raise BadRequest("choose a downloaded genome to annotate")
    sp = body.get("species")
    if sp is not None and not isinstance(sp, str):
        raise BadRequest("species must be a string")
    sens = str(body.get("sensitivity", "default"))
    if sens not in ("default", "quick", "slow"):
        raise BadRequest("sensitivity must be default, quick or slow")
    t0 = time.time()
    custom_path = body.get("custom_library") or None
    if custom_path is not None and not isinstance(custom_path, str):
        raise BadRequest("custom_library must be a path")
    r = wsl.genome_annotate(acc, species=(sp or None), threads=int(_num(body.get("threads"), 4)),
                            sensitivity=sens, chunk_mb=int(_num(body.get("chunk_mb"), 40)),
                            timeout=int(_num(body.get("timeout"), 86400)), custom_lib=custom_path,
                            include_uncurated=bool(body.get("include_uncurated")))
    if r.get("ok"):
        r["elapsed_s"] = round(time.time() - t0, 1)
        custom = r.get("custom_library") if isinstance(r.get("custom_library"), dict) else None
        if custom_path:
            custom_sha = (custom or {}).get("sha256")
            custom_n = (custom or {}).get("sequences")
            if (not isinstance(custom_sha, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", custom_sha)
                    or isinstance(custom_n, bool) or not isinstance(custom_n, int) or custom_n < 1):
                return {"ok": False, "error": "annotation completed but the requested custom library did not "
                        "return a valid checksum and positive sequence count; the result was not sealed",
                        "partial": True}
        references = refs.for_run("annotate", fetched=False, annotate_dfam=not bool(custom))
        r["references"] = references
        dfam_lib = r.get("dfam_library") or {}
        if custom:
            # File name is a display label, not scientific identity. The content hash and sequence count
            # seal the library, so renaming the same bytes does not create a different experiment.
            library_db = {"name": "User-supplied repeat library", "sha256": custom.get("sha256"),
                          "sequences": custom.get("sequences"), "bytes": custom.get("bytes")}
            library_db = {k: v for k, v in library_db.items() if v is not None}
            not_run = ["Dfam database search (replaced by the user-supplied FASTA library)",
                       "De novo family discovery (RepeatModeler/EDTA)"]
        else:
            library_db = {"name": "Dfam (via famdb)",
                          "version": dfam_lib.get("version") or "unresolved — famdb reported no library version"}
            for src_k, dst_k in (("date", "releaseDate"), ("famdbFormat", "famdbFormat"),
                                 ("partitions", "partitions"), ("consensusSequences", "consensusSequences")):
                if dfam_lib.get(src_k) is not None:
                    library_db[dst_k] = dfam_lib[src_k]
            not_run = ["De novo family discovery (RepeatModeler/EDTA)"]
        genome_sha = r.get("genome_sha256")
        genome_bp = r.get("genome_bp")
        if not isinstance(genome_sha, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", genome_sha):
            return {"ok": False, "error": "annotation completed but its verified genome checksum was missing; "
                    "the result was not sealed", "partial": True}
        if isinstance(genome_bp, bool) or not isinstance(genome_bp, int) or genome_bp < 1:
            return {"ok": False, "error": "annotation completed but its verified genome length was missing; "
                    "the result was not sealed", "partial": True}
        if r.get("accession", acc) != acc or not provenance.is_public_accession(acc):
            return {"ok": False, "error": "annotation completed but its assembly accession did not match the "
                    "request; the result was not sealed", "partial": True}
        # The seal identifies the GENOME by the checksum recorded at prepare time and records every
        # parameter that changes which families could be found. `complete` is sealed too: a run that
        # stopped early must never be indistinguishable from one that covered the whole assembly.
        r["provenance"] = provenance.build_manifest(
            "genome_annotate", "", acc,
            {"engine": "RMBLAST", "species": r.get("species"), "sensitivity": sens,
             "threads": int(_num(body.get("threads"), 4)), "chunk_mb": int(_num(body.get("chunk_mb"), 40)),
             "chunks": r.get("chunks"), "complete": bool(r.get("complete")),
             # which library was searched is part of what the result MEANS, so it is sealed with the run:
             # the same genome against a different library is a different experiment.
             "library_kind": r.get("library_kind"),
             "include_uncurated": bool(r.get("include_uncurated")),
             "custom_library_sha256": (custom or {}).get("sha256")},
            not_run=not_run,
            # the assembly is named by its accession beside the checksum, as in the genome scan, so a reader
            # can tell which genome was annotated; the input identity below hashes only the accession text
            databases=[library_db, {"name": "RepeatMasker", "version": r.get("repeatmasker_version")},
                       {"name": "RefSeq genome assembly", "version": acc, "sha256": genome_sha}],
            # the genome's sealed checksum comes from the PREPARE step's own meta.txt (wsl returns it),
            # not from whatever the caller passed in — a caller-supplied hash would let the manifest
            # claim an identity the analysed sequence was never checked against.
            references=references, input_sha256=genome_sha, input_length=genome_bp)
    return r


def run_genome_annotate_reset(body):
    """Discard a part-finished annotation so a run with different settings can start clean."""
    return wsl.genome_annotate_reset(str(body.get("assemblyAccession", "")).strip())


def run_genome_annotate_log(body):
    """Latest annotation milestone (drives the N-of-M progress line during a multi-hour run)."""
    return {"log": wsl.genome_annotate_log(int(_num((body or {}).get("tail"), 1)))}


def run_genome_annotate_status(body):
    """How much of an annotation already exists on disk (resume + honest cost estimate)."""
    return wsl.genome_annotate_status(str(body.get("assemblyAccession", "")).strip())


def run_genome_remove(body):
    """Delete a cached genome to reclaim disk."""
    acc = str(body.get("assemblyAccession", "")).strip()
    return wsl.genome_remove(acc)


def run_add_custom_assembly(body):
    """Resolve a typed organism / assembly accession against NCBI and pin it to the user store. Runs off the
    UI thread via the engine worker so the resolve never freezes the window; a bad query raises BadRequest
    (surfaced as a clean user error, not a traceback)."""
    q = (body.get("query") or "").strip()
    if not q:
        raise BadRequest("enter an organism name or an assembly accession")
    try:
        return fetch.add_custom_assembly(q)
    except fetch.CoordError as e:
        raise BadRequest(e.args[0] if getattr(e, "args", None) else str(e))


def _design_loci(raw):
    """Genome positions of the designed product from `design_locus` entries (dict or list of dicts).

    An entry is {accession, start, end (or stop), strand (1/2 or '+'/'-'), intended_product_span [a, b],
    intended_endpoint_tolerance_bp}: the fetched window in 1-based inclusive assembly coordinates plus the
    designed product's 0-based half-open span within that fetched record. The span is mapped through the
    window, strand-aware: a reverse-complement fetch counts from the window's end. An entry that lacks the
    accession, either window end, the strand or the span names no intended product and is ignored."""
    raw = raw if isinstance(raw, list) else ([raw] if isinstance(raw, dict) else [])
    loci = []
    for L in raw:
        if not isinstance(L, dict) or not L.get("accession") or L.get("intended_product_span") is None:
            continue
        end = L.get("end", L.get("stop"))
        if L.get("start") is None or end is None or L.get("strand") is None:
            continue
        start, end, strand = int(_num(L["start"], 0)), int(_num(end, 0)), L["strand"]
        if isinstance(strand, bool) or strand not in (1, 2, "+", "-"):
            raise BadRequest("design_locus strand must be 1 or + for the forward strand, 2 or - for reverse")
        if start < 1 or end < start:
            raise BadRequest("design_locus window must satisfy 1 <= start <= end")
        try:
            a, b = primers.endpoint_span(L["intended_product_span"], end - start + 1, "intended_product_span")
            tol = primers.endpoint_tolerance(L.get("intended_endpoint_tolerance_bp",
                                                   primers.INTENDED_ENDPOINT_TOLERANCE_BP))
        except (TypeError, ValueError) as e:
            raise BadRequest(str(e))
        lo, hi = (start + a, start + b - 1) if strand in (1, "+") else (end - b + 1, end - a)
        loci.append({"acc": str(L["accession"]), "lo": lo, "hi": hi, "tol": tol})
    return loci


def run_genome_pcr(body):
    """Whole-genome off-target scan of one primer pair against a LOCAL, downloaded RefSeq assembly (isPcr).
    No remote query; a local safety timeout applies (body.timeout, default 600s — a very large genome may need
    it raised). The assembly is a fixed checksummed file, so the result IS reproducible and sealed (accession +
    source-FASTA sha256 + isPcr version + params). Products are candidate priming sites, not validated bands;
    on/off-target labels are applied only when the request maps the designed product to a genome position."""
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
    loci = _design_loci(body.get("design_locus"))
    r = wsl.genome_scan(acc, genomepcr.query_rows(fwd, rev), max_size=max_size, min_size=min_size,
                        min_perfect=min_perfect, min_good=min_good, timeout=int(_num(body.get("timeout"), 600)))
    if not r.get("ok"):
        return r                                     # {ok:False, error, [need_prepare]} — surfaced to the UI
    amps = genomepcr.parse_ispcr(r["raw"])
    for amp in amps:
        # These identifiers come from the verified, checksummed assembly downloaded through the managed
        # reference-genome path. Recognized accessions remain readable; malformed labels are still
        # neutralized by the validated-fetch identity rules at the display/export boundary.
        amp["sourceIdentitySource"] = "validated_fetch"
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
    # on-target from the DESIGNED PRODUCT: when the request maps the designed product to a genome position
    # (see _design_loci), the pair product at that position is the on-target and the rest are off-target candidate
    # priming sites; with no designed position there is no single intended target, so the products are neutral
    # 'genomic priming sites'. This labeling is DERIVED, not sealed
    # (like the summary): the isPcr product set — the sealed content — is identical with or without a design locus.
    has_locus = bool(loci)
    n_on = 0
    if has_locus:
        best = None                                      # exactly ONE on-target: the product whose two endpoints match the
        for a in amps:                                   # designed product's genome position; never another product that
            if a.get("single_primer"):                   # merely overlaps it, and nothing when the designed product is absent
                continue
            for L in loci:
                err = max(abs(a["start"] - L["lo"]), abs(a["end"] - L["hi"]))
                if a["source"] == L["acc"] and err <= L["tol"] and (best is None or (err, a["start"]) < best[0]):
                    best = ((err, a["start"]), a)
        if best:
            best[1]["on_target"] = True
            n_on = 1
    # interpretation (on/off split, per-assembly-sequence spread, verdict). DERIVED, NOT sealed — never passed to
    # build_manifest, so the same primers+genome seal identically regardless of the design locus or the products.
    summary = genomepcr.summarize(amps, has_locus=has_locus)
    return {"ok": True, "organism": org, "taxid": taxid, "scope": "whole-genome", "assemblyAccession": acc,
            "assemblyName": name, "amplicons": amps, "n_amplicons": len(amps), "n_seqs": r.get("n_seqs"),
            "summary": summary, "has_locus": has_locus, "n_on_target": n_on,
            "advisory_note": ("candidate priming sites (isPcr default: ≥15 bp minimum perfect 3'-end match) — not wet-lab-validated "
                              "amplicons; products priming off a shorter 3' match are not flagged"),
            "references": refs.for_run("genome-scan"), "provenance": manifest}
