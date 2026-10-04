"""Live sequence retrieval from NCBI (E-utilities) — real web fetch, server-side.
Verifies metadata (esummary) before pulling the sequence (efetch). Errors are typed
and explicit; a failed or empty fetch is never returned as a sequence."""
from __future__ import annotations
import os, ssl, urllib.request, urllib.parse, urllib.error, json, re, datetime, hashlib, time, threading, bisect

# Verify TLS against certifi's CA bundle when available (some Windows Pythons lack the roots for UCSC/EBI/
# NCBI chains); falls back to the system default. Hardens every HTTPS call the app makes.
#
# Built on FIRST USE, not at import. Parsing the CA bundle costs ~140 ms, which was paid on every launch
# even though most sessions start with a pasted sequence and never open a socket. The context is created
# once and cached, so the first network call pays it and no later call does.
_SSL_CTX_CACHE = []


def _ssl_ctx():
    if not _SSL_CTX_CACHE:
        try:
            import certifi
            _SSL_CTX_CACHE.append(ssl.create_default_context(cafile=certifi.where()))
        except Exception:
            _SSL_CTX_CACHE.append(None)
    return _SSL_CTX_CACHE[0]



EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ENA = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
TOOL = "TEagle"
ACC_RE = re.compile(r"^(?:NZ_)?[A-Za-z]{1,6}_?\d+(\.\d+)?$")   # GenBank/RefSeq nucleotide accession
_PROTEIN_ACC = re.compile(r"^(?:NP|XP|YP|WP|AP|ZP)_\d", re.I)       # RefSeq protein prefixes
_ASSEMBLY_ACC = re.compile(r"^GC[AF]_\d", re.I)
GENE_MODEL_VERSION = 2      # bump when build_gene_model output changes, so cached models are rebuilt

# On-disk cache: an accession.version identifies a fixed sequence version, so a fetched sequence is reused instead of
# re-downloaded on every run (the cached FASTA + original retrieval time travel with results).
from . import appdirs
CACHE_DIR = os.path.join(appdirs.user_data_dir(), "cache", "fetch")


def _cache_path(acc: str) -> str:
    return os.path.join(CACHE_DIR, re.sub(r"[^A-Za-z0-9._-]", "_", acc) + ".json")


def _read_cache(acc: str):
    try:
        p = _cache_path(acc)
        if os.path.isfile(p):
            with open(p, encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return None


def _write_cache(acc: str, meta: dict):
    try:                                                   # best-effort; never fail a fetch on cache error
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(_cache_path(acc), "w", encoding="utf-8") as f:
            json.dump(meta, f)
    except Exception:
        pass


class FetchError(Exception):
    pass


_NCBI_GAP = 0.34            # NCBI allows 3 requests/s without an API key
_last_ncbi = [0.0]
_ncbi_lock = threading.Lock()


def _pace(url: str):
    if "ncbi.nlm.nih.gov" not in url:
        return
    with _ncbi_lock:
        wait = _last_ncbi[0] + _NCBI_GAP - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        _last_ncbi[0] = time.monotonic()


def _retry_delay(e, attempt: int) -> float:
    try:
        return min(float(e.headers.get("Retry-After")), 10.0)
    except (AttributeError, TypeError, ValueError):
        return float(2 ** attempt)


def _get(url: str, timeout: int = 25) -> str:
    from . import __version__
    req = urllib.request.Request(url, headers={"User-Agent": f"TEagle/{__version__} (+local)"})
    for attempt in range(3):                               # 429/5xx are retried twice; nothing else is
        _pace(url)
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx()) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if (e.code == 429 or 500 <= e.code < 600) and attempt < 2:
                time.sleep(_retry_delay(e, attempt))
                continue
            if e.code == 429:
                raise FetchError("HTTP 429 from source — the server is rate-limiting requests; "
                                 "wait a few seconds and try again")
            if e.code >= 500:
                raise FetchError(f"HTTP {e.code} from source — the service is temporarily unavailable; try again shortly")
            raise FetchError(f"HTTP {e.code} from source")
        except urllib.error.URLError as e:
            raise FetchError(f"network error: {getattr(e, 'reason', e)} — check the internet connection and try "
                             "again (previously fetched accessions still load from the local cache)")
        except Exception as e:                             # timeout etc.
            raise FetchError(str(e))


def validate_accession(acc: str) -> str:
    acc = (acc or "").strip()
    if not acc:
        raise FetchError("no accession provided")
    if _ASSEMBLY_ACC.match(acc):
        raise FetchError(f"'{acc}' is an assembly accession, not a sequence record — use 'Fetch by coordinate'")
    if _PROTEIN_ACC.match(acc):
        raise FetchError(f"'{acc}' is a protein accession — TEagle fetches nucleotide records (e.g. X04456)")
    if not ACC_RE.match(acc):
        raise FetchError(f"'{acc}' is not a valid nucleotide accession (e.g. X04456, AB087615)")
    return acc


def resolve(acc: str) -> dict:
    """Verify the accession and return metadata BEFORE downloading the sequence."""
    acc = validate_accession(acc)
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "retmode": "json", "tool": TOOL})
    try:                                                       # a 200 with a non-JSON/empty body (maintenance/rate-limit)
        data = json.loads(_get(EUTILS + "esummary.fcgi?" + q))
    except (json.JSONDecodeError, ValueError):
        raise FetchError("NCBI returned an unexpected (non-JSON) response — try again shortly")
    res = data.get("result", {})
    uids = res.get("uids", [])
    if not uids:
        raise FetchError(f"accession '{acc}' not found in NCBI nuccore")
    r = res[uids[0]]
    if isinstance(r, dict) and r.get("error"):
        raise FetchError(f"NCBI: {r['error']}")
    return {
        "accession": r.get("accessionversion", acc),
        "organism": r.get("organism", ""),
        "taxid": r.get("taxid"),
        "title": r.get("title", ""),
        "length": r.get("slen"),
        "moltype": r.get("moltype", ""),
        "source": "NCBI nuccore",
        "status": r.get("status") or "",                       # display-only: NCBI record status (live/replaced/...)
        "replacedBy": r.get("replacedby") or "",
    }


def fetch_fasta(acc: str, served: list | None = None) -> str:
    """Fetch the FASTA for an accession. If `served` is given, the (source-label, endpoint) that
    actually served the bytes is appended to it, so the caller records the real serving database
    (NCBI or the ENA fallback) rather than always claiming NCBI."""
    acc = validate_accession(acc)
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "fasta",
                                "retmode": "text", "tool": TOOL})
    src = ("NCBI nuccore", EUTILS + "efetch.fcgi")
    try:
        txt = _get(EUTILS + "efetch.fcgi?" + q)
    except FetchError:
        txt = ""                                              # NCBI HTTP/URL error -> still try the ENA fallback below
    if not txt.lstrip().startswith(">"):
        # fall back to ENA once before giving up (covers both a non-FASTA 200 body and an NCBI request error)
        try:
            txt = _get(ENA + urllib.parse.quote(acc))
            src = ("ENA (EMBL-EBI)", ENA + acc)
        except FetchError:
            pass
    if not txt.lstrip().startswith(">"):
        raise FetchError("no FASTA returned (check the accession)")
    if served is not None:
        served.append(src)
    return txt


_COORD_RE = re.compile(r"^[<>]?\d+$")


def _interval(a: str, b: str) -> dict:
    def n(x): return int(str(x).lstrip("<>").strip())
    s, e = n(a), n(b)
    if s <= e:
        return {"start": s - 1, "end": e, "strand": "+"}      # 1-based inclusive -> 0-based half-open
    return {"start": e - 1, "end": s, "strand": "-"}          # complement interval (start>end)


def parse_feature_table(text: str) -> list:
    """Parse an NCBI feature table (efetch rettype=ft). Spliced features list several
    coordinate lines before their qualifiers; qualifier lines are tab-indented."""
    feats, cur = [], None
    for line in text.splitlines():
        if not line or line.startswith(">Feature"):
            continue
        parts = line.split("\t")
        if parts[0] != "" and len(parts) >= 2 and parts[1] != "" and _COORD_RE.match(parts[0].strip()):
            key = parts[2] if len(parts) >= 3 and parts[2] != "" else None
            iv = _interval(parts[0], parts[1])
            if key:                                            # new feature
                cur = {"type": key, "intervals": [iv], "qualifiers": {}}
                feats.append(cur)
            elif cur is not None:                              # another interval of a spliced feature
                cur["intervals"].append(iv)
        elif parts[0] == "" and cur is not None:               # qualifier line
            q = [p for p in parts if p != ""]
            if q:
                cur["qualifiers"].setdefault(q[0], q[1] if len(q) > 1 else "")
    return feats


def _dedup(ivs: list, identity_fields=()) -> list:
    seen, out = set(), []
    for iv in ivs:
        k = (iv["start"], iv["end"], *(iv.get(field) for field in identity_fields))
        if k not in seen:
            seen.add(k)
            out.append(iv)
    return sorted(out, key=lambda x: x["start"])


def _transcript_accession(raw: str):
    """The nucleotide accession inside a /transcript_id value ('NM_1.2', 'gb|AB1.1|', 'ref|NM_1.2|'), else None.
    Local identifiers (gnl|, lcl|) and anything that is not a whole accession field are never guessed at."""
    raw = (raw or "").strip()
    if raw.lower().startswith(("gnl|", "lcl|")):
        return None
    for field in raw.split("|"):
        field = field.strip()
        if ACC_RE.match(field) and not _PROTEIN_ACC.match(field) and not _ASSEMBLY_ACC.match(field):
            return field
    return None


def build_gene_model(feats: list) -> dict:
    """Reduce a feature table to exon/intron/CDS segments. Introns and exons are taken from
    explicit features when present, else derived from spliced (joined) CDS/mRNA coordinates.
    Introns are only ever derived WITHIN one join (one CDS/mRNA feature) or one gene's explicit exons, never
    across separate features or genes, and every derived intron is flagged derived."""
    def collect(kind):
        out = []
        group_index = 0
        for fi, f in enumerate(feats):
            if f["type"] == kind:
                group_index += 1
                q = f["qualifiers"]
                label = q.get("number") or q.get("gene") or q.get("product") or ""
                intervals = [dict(iv) for iv in f["intervals"]]
                # `_unit` names what one chain of exons belongs to; it never leaves this function
                unit = f"{kind}_{fi}" if kind != "exon" else "exon:" + (q.get("gene") or q.get("locus_tag") or
                                                                       q.get("transcript_id") or "")
                if kind == "CDS":
                    strands = {iv.get("strand", "+") for iv in intervals}
                    strand = next(iter(strands)) if len(strands) == 1 else None
                    intervals.sort(key=lambda iv: iv["start"], reverse=(strand == "-"))
                    # INSDC: codon_start defaults to 1, and NCBI's feature table omits it in that case; only a
                    # qualifier that is present but not 1/2/3 leaves the phase unavailable
                    raw_start = f["qualifiers"].get("codon_start")
                    phase = 0 if raw_start is None else (int(raw_start) - 1 if str(raw_start) in {"1", "2", "3"} else None)
                    for iv in intervals:
                        row = {**iv, "note": label, "cds_group": f"cds_{group_index}", "phase": phase, "_unit": unit}
                        out.append(row)
                        length = iv["end"] - iv["start"]
                        phase = ((3 - ((length - phase) % 3)) % 3
                                 if phase is not None and length > phase else None)
                else:
                    out.extend({**iv, "note": label, "_unit": unit + (":" + iv["strand"] if kind == "exon" else "")}
                               for iv in intervals)
        return out
    exons, introns, cds, mrna = collect("exon"), collect("intron"), collect("CDS"), collect("mRNA")
    derived_exons = derived_introns = False
    if not exons and (mrna or cds):                            # derive exons from a spliced transcript/CDS
        # flag each so the viewer marks them exon* (honest); the CDS-only keys do not describe an exon
        exons = [{**{k: v for k, v in x.items() if k not in ("cds_group", "phase")}, "derived": True}
                 for x in (mrna or cds)]
        derived_exons = True
    chains = {}
    for e in exons:
        chains.setdefault(e["_unit"], []).append(e)
    if not introns:                                            # introns = gaps between consecutive exons OF ONE UNIT
        for chain in chains.values():
            chain = _dedup(chain)
            for a, b in zip(chain, chain[1:]):
                if b["start"] > a["end"]:
                    introns.append({"start": a["end"], "end": b["start"], "strand": a.get("strand", "+"),
                                    "note": "", "derived": True})
        derived_introns = bool(introns)
    def clean(rows):
        return [{k: v for k, v in r.items() if k != "_unit"} for r in rows]
    exons, introns = _dedup(clean(exons)), _dedup(clean(introns))
    cds = _dedup(clean(cds), ("strand", "phase", "cds_group"))
    # the record's OWN annotated transcripts (mRNA /transcript_id = NM_/XM_) — a splice picker offers these so a
    # user does not have to hunt for the matching transcript. Aligning one back to its own record is a CONSISTENCY
    # check (same annotation source), NOT independent confirmation — the UI relabels accordingly.
    transcripts, seen_tx = [], set()
    for f in feats:
        acc = _transcript_accession(f.get("qualifiers", {}).get("transcript_id"))
        if not acc or acc in seen_tx:
            continue
        seen_tx.add(acc)
        q = f["qualifiers"]
        transcripts.append({"accession": acc, "product": q.get("product", ""), "gene": q.get("gene", "")})
    return complete_gene_model({
        "exons": exons, "introns": introns, "cds": cds,
        "counts": {"exons": len(exons), "introns": len(introns), "cds": len(cds)},
        "derived_exons": derived_exons, "derived_introns": derived_introns,
        "transcripts": transcripts,
        "source": "NCBI feature table (efetch rettype=ft)",
    })


def _subtract(iv: dict, cuts: list) -> list:
    pieces = [(iv["start"], iv["end"])]
    for c in cuts:
        nxt = []
        for lo, hi in pieces:
            if c["end"] <= lo or c["start"] >= hi:
                nxt.append((lo, hi))
                continue
            if lo < c["start"]:
                nxt.append((lo, c["start"]))
            if c["end"] < hi:
                nxt.append((c["end"], hi))
        pieces = nxt
    return [{**iv, "start": lo, "end": hi} for lo, hi in pieces]


def complete_gene_model(gm: dict) -> dict:
    """Fill in exons implied by the CDS/mRNA when a record annotates some exons but omits others — e.g.
    J00265 annotates exon 1 and the last exon but NOT the middle exon that carries the first CDS segment.
    A CDS segment must lie inside an exon; where none does, synthesize a derived exon spanning the flanking
    intron boundaries (its true extent), or the CDS segment span if introns don't bracket it. Idempotent, so
    it is safe to re-apply to an already-complete or cached model."""
    exons = [dict(e) for e in gm.get("exons", [])]
    introns, cds = gm.get("introns", []), gm.get("cds", [])
    derived_introns = bool(gm.get("derived_introns"))
    bounds = [] if derived_introns else introns                # gaps derived from the exons say nothing about an absent exon
    ordered = sorted(exons, key=lambda e: e["start"])
    starts, reach, far = [e["start"] for e in ordered], [], -1
    for e in ordered:
        far = max(far, e["end"])
        reach.append(far)
    new = []
    def covered(iv):
        i = bisect.bisect_right(starts, iv["start"]) - 1       # farthest-reaching exon starting at or before iv
        return (i >= 0 and reach[i] >= iv["end"]) or any(e["start"] <= iv["start"] and iv["end"] <= e["end"] for e in new)
    for c in cds:
        if covered(c):
            continue
        left = max([i["end"] for i in bounds if i["end"] <= c["start"]], default=c["start"])
        right = min([i["start"] for i in bounds if i["start"] >= c["end"]], default=c["end"])
        new.append({"start": min(left, c["start"]), "end": max(right, c["end"]),
                    "strand": c.get("strand", "+"), "note": c.get("note", ""), "derived": True})
    if not new:
        return gm
    exons = _dedup(exons + new)
    out = dict(gm)
    out["exons"] = exons
    out["derived_exons"] = True
    # derived introns are the gaps between the OLD exons of one unit: a newly inserted exon splits the intron it
    # lands in. Explicit introns are authoritative and untouched.
    if derived_introns:
        introns = _dedup([piece for i in introns for piece in _subtract(i, new)])
        out["introns"] = introns
    out["counts"] = {"exons": len(exons), "introns": len(out.get("introns", introns)), "cds": len(cds)}
    return out


def cross_check_models(annotation_introns: list, aligned_introns: list, tol: int = 2) -> dict:
    """Compare the record's (completed) annotation introns against an INDEPENDENT alignment's introns
    (an external transcript aligned to the genome — NOT a transcript rebuilt from the same annotation,
    which would be circular). An intron matches when both boundaries agree within `tol` bp. Returns the
    match count plus the introns seen on only one side, so the UI can flag agreement vs discrepancies."""
    aln = list(aligned_introns)
    used, matched, ann_only = set(), 0, []
    for a in annotation_introns:
        hit = next((j for j, b in enumerate(aln) if j not in used
                    and abs(a["start"] - b["start"]) <= tol and abs(a["end"] - b["end"]) <= tol), None)
        if hit is None:
            ann_only.append(a)
        else:
            used.add(hit); matched += 1
    aln_only = [aln[j] for j in range(len(aln)) if j not in used]
    return {"matched": matched, "annotation_total": len(annotation_introns), "aligned_total": len(aln),
            "annotation_only": ann_only, "aligned_only": aln_only}


def fetch_features(acc: str) -> dict:
    """Fetch and parse the accession's annotation feature table into a gene model."""
    acc = validate_accession(acc)
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "ft", "retmode": "text", "tool": TOOL})
    return build_gene_model(parse_feature_table(_get(EUTILS + "efetch.fcgi?" + q)))


_ENA_HEADER = re.compile(r"^>ENA\|[^|\s]+\|(\S+)(.*)$")


def _ena_header(fasta: str, accession: str) -> str:
    """ENA names its record '>ENA|acc|acc.ver ...'. When that third field is exactly the accession NCBI
    resolved, rewrite the header to '>acc.ver ...' (as fetch_fasta_range does); anything else is left untouched,
    so a record that does not match the verified metadata is never presented as the fetched accession."""
    head, nl, rest = fasta.partition("\n")
    m = _ENA_HEADER.match(head.strip())
    return f">{accession}{m.group(2)}{nl}{rest}" if m and m.group(1) == accession else fasta


def _canonical_urls(meta: dict):
    """The endpoint and record URL of the resolved accession.version (a cache written from the typed text
    names the unversioned form, which never matched the verified record)."""
    acc = meta.get("accession")
    if not acc:
        return
    if str(meta.get("source", "")).startswith("ENA"):
        meta["endpoint"], meta["sourceUrl"] = ENA + acc, "https://www.ebi.ac.uk/ena/browser/view/" + acc
    else:
        meta["sourceUrl"] = "https://www.ncbi.nlm.nih.gov/nuccore/" + acc


def _attach_features(acc: str, meta: dict):
    """Best-effort gene annotation (exon/intron/CDS); never blocks the fetch. featuresStatus records whether the
    table was read ('ok'), held nothing ('none') or could not be fetched ('failed' — retried on the next use)."""
    meta.pop("features", None)
    try:
        gm = fetch_features(acc)
    except Exception:
        meta["featuresStatus"] = "failed"
    else:
        if gm["counts"]["exons"] or gm["counts"]["cds"]:
            meta["features"], meta["featuresStatus"] = gm, "ok"
        else:
            meta["featuresStatus"] = "none"
    meta["featuresModel"] = GENE_MODEL_VERSION


def retrieve(acc: str, refresh: bool = False) -> dict:
    """Full retrieval: metadata + sequence, with provenance fields. Cached on disk per accession
    so the same organism is not re-downloaded on every run; `refresh=True` forces a fresh pull."""
    acc = validate_accession(acc)
    if not refresh:
        cached = _read_cache(acc)
        if cached and cached.get("fasta"):
            _canonical_urls(cached)
            if cached.get("featuresStatus") == "failed" or cached.get("featuresModel") != GENE_MODEL_VERSION:
                _attach_features(acc, cached)                  # sequence and retrievedUtc stay as first retrieved
                _write_cache(acc, cached)
            cached["fromCache"] = True
            cached["servedUtc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
            return cached
    meta = resolve(acc)
    served = []
    fasta = fetch_fasta(meta["accession"], served)             # the verified accession.version, from NCBI or ENA alike
    if served and served[0][0].startswith("ENA"):
        fasta = _ena_header(fasta, meta["accession"])
    seq = "".join(l for l in fasta.splitlines() if not l.startswith(">"))
    meta["fasta"] = fasta
    meta["seq_length"] = len(seq)
    meta["retrievedUtc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    if served:                                                 # record the DB that actually served the sequence
        meta["source"], meta["endpoint"] = served[0]
    else:
        meta["endpoint"] = EUTILS + "efetch.fcgi"
    _canonical_urls(meta)
    meta["fromCache"] = False
    _attach_features(acc, meta)
    _write_cache(acc, meta)
    return meta


# ============================ coordinate-based fetch (UCSC-style) ============================
# organism -> pinned RefSeq reference assembly. The accession is the reproducibility anchor and taxid the
# organism pin (both verified live via NCBI Datasets, so a curated coordinate run seals the same organism
# identity the "Other organism" resolver would); regenerate the bundled chromosome maps with build_assembly_maps.py.
COORD_ASSEMBLIES = {
    "Homo sapiens":              {"assemblyName": "GRCh38.p14", "assemblyAccession": "GCF_000001405.40", "taxid": "9606"},
    "Mus musculus":              {"assemblyName": "GRCm39",     "assemblyAccession": "GCF_000001635.27", "taxid": "10090"},
    "Rattus norvegicus":         {"assemblyName": "GRCr8",      "assemblyAccession": "GCF_036323735.1", "taxid": "10116"},
    "Bos taurus":                {"assemblyName": "ARS-UCD2.0", "assemblyAccession": "GCF_002263795.3", "taxid": "9913"},
    # GCF_000003025.6 (Sscrofa11.1) was pinned here but no longer resolves at NCBI Datasets — it returned
    # zero reports on 2026-07-28, so selecting pig failed only once the user had waited for a download.
    "Sus scrofa":                {"assemblyName": "Sscrofa_HxYL_maternal_Yorkshire_Landrace_v2.0",
                                  "assemblyAccession": "GCF_054392235.1", "taxid": "9823"},
    "Canis lupus familiaris":    {"assemblyName": "UU_Cfam_GSD_1.0", "assemblyAccession": "GCF_011100685.1", "taxid": "9615"},
    "Macaca mulatta":            {"assemblyName": "Mmul_10",    "assemblyAccession": "GCF_003339765.1", "taxid": "9544"},
    "Gallus gallus":             {"assemblyName": "GRCg7b",     "assemblyAccession": "GCF_016699485.2", "taxid": "9031"},
    "Xenopus tropicalis":        {"assemblyName": "UCB_Xtro_10.0", "assemblyAccession": "GCF_000004195.4", "taxid": "8364"},
    "Danio rerio":               {"assemblyName": "GRCz11",     "assemblyAccession": "GCF_000002035.6", "taxid": "7955"},
    "Drosophila melanogaster":   {"assemblyName": "Release 6 plus ISO1 MT", "assemblyAccession": "GCF_000001215.4", "taxid": "7227"},
    "Caenorhabditis elegans":    {"assemblyName": "WBcel235",   "assemblyAccession": "GCF_000002985.6", "taxid": "6239"},
    "Saccharomyces cerevisiae":  {"assemblyName": "R64",        "assemblyAccession": "GCF_000146045.2", "taxid": "559292"},
    "Schizosaccharomyces pombe": {"assemblyName": "ASM294v3",   "assemblyAccession": "GCF_000002945.2", "taxid": "4896"},
    "Arabidopsis thaliana":      {"assemblyName": "TAIR10.1",   "assemblyAccession": "GCF_000001735.4", "taxid": "3702"},
    "Oryza sativa":              {"assemblyName": "IRGSP-1.0",  "assemblyAccession": "GCF_001433935.1", "taxid": "39947"},
    "Zea mays":                  {"assemblyName": "Zm-B73-REFERENCE-NAM-5.0", "assemblyAccession": "GCF_902167145.1", "taxid": "4577"},
    # Insect and vector genomics — the community that built the piggyBac, Minos and Hermes transgenesis
    # systems — previously had a single entry, and it was a model organism rather than a vector. Every
    # accession below was resolved through NCBI Datasets (fetch.resolve_assembly) on 2026-07-28, not
    # transcribed from memory.
    "Aedes aegypti":             {"assemblyName": "AaegL5.0",   "assemblyAccession": "GCF_002204515.2", "taxid": "7159"},
    "Anopheles gambiae":         {"assemblyName": "idAnoGambNW_F1_1", "assemblyAccession": "GCF_943734735.2", "taxid": "7165"},
    "Bombyx mori":               {"assemblyName": "ASM3026992v2", "assemblyAccession": "GCF_030269925.1", "taxid": "7091"},
    "Tribolium castaneum":       {"assemblyName": "icTriCast1.1", "assemblyAccession": "GCF_031307605.1", "taxid": "7070"},
    # aquaculture / fish models
    "Salmo salar":               {"assemblyName": "fSalSal1",   "assemblyAccession": "GCF_965601325.1", "taxid": "8030"},
    "Oryzias latipes":           {"assemblyName": "Olatipes_Hd-rR_3.1", "assemblyAccession": "GCF_053564925.1", "taxid": "8090"},
    # crop plants, where most TE families are uncurated in Dfam (see the family-naming coverage note)
    "Glycine max":               {"assemblyName": "Glycine_max_v4.0", "assemblyAccession": "GCF_000004515.6", "taxid": "3847"},
    "Solanum lycopersicum":      {"assemblyName": "SLM_r2.1",   "assemblyAccession": "GCF_036512215.1", "taxid": "4081"},
    "Triticum aestivum":         {"assemblyName": "IWGSC CS RefSeq v2.1", "assemblyAccession": "GCF_018294505.1", "taxid": "4565"},
}
_DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
_REGION_RE = re.compile(r"^(.+?)\s*:\s*([\d,\s]+)\s*-\s*([\d,\s]+)$")
_ASM_ACC_RE = re.compile(r"^GC[AF]_\d+\.\d+$", re.I)


class CoordError(FetchError):
    """A malformed / unresolvable coordinate request — surfaced to the user, never a 500."""


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def parse_regions(text: str) -> list:
    """Parse UCSC-style loci (1-based inclusive) into normalized regions. Split on newline/';' ONLY
    (comma is a thousands separator inside a coordinate). Organism-specific chrom names (2L, II, X, MT)
    are kept verbatim for the resolver to match against the assembly map."""
    if not isinstance(text, str) or not text.strip():
        raise CoordError("no regions provided")
    regions = []
    for tok in re.split(r"[\n;]+", text.strip()):
        tok = tok.strip()
        if not tok:
            continue
        m = _REGION_RE.match(tok)
        if not m:
            raise CoordError(f"could not parse '{tok}' — use e.g. chr13:33,016,423-33,066,143")
        chrom_raw = m.group(1).strip()
        chrom_key = re.sub(r"^chr", "", chrom_raw, flags=re.I).strip() or chrom_raw
        s_start, s_end = re.sub(r"[,\s]", "", m.group(2)), re.sub(r"[,\s]", "", m.group(3))
        if not s_start or not s_end:                          # a group of only commas/spaces -> int('') would raise
            raise CoordError(f"could not parse coordinates in '{tok}' — use e.g. chr13:33,016,423-33,066,143")
        start, end = int(s_start), int(s_end)
        if start < 1 or end < 1:
            raise CoordError(f"coordinates must be >= 1 in '{tok}'")
        if start > end:
            raise CoordError(f"start > end in '{tok}' — a coordinate range must be ascending")
        label = chrom_raw if chrom_raw.lower().startswith("chr") else "chr" + chrom_key
        regions.append({"chromKey": chrom_key, "chromLabel": label, "start": start, "end": end})
    if not regions:
        raise CoordError("no regions provided")
    return regions


def _datasets_json(path: str) -> dict:
    txt = _get(_DATASETS + path, timeout=40)                  # _get already maps HTTP/URL errors to FetchError
    try:
        return json.loads(txt)
    except (json.JSONDecodeError, ValueError):                # 200 with a non-JSON body (maintenance/rate-limit page)
        raise CoordError("NCBI Datasets returned an unexpected (non-JSON) response — try again shortly")


def _assembly_map_dir() -> str:
    return appdirs.resource("teagle_core", "data", "assemblies") or \
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "assemblies")


def _sequence_reports(assembly_accession: str) -> list:
    """Every sequence report of an assembly. Datasets pages them (organelles come after the unplaced
    scaffolds), so the next_page_token is followed until it runs out."""
    reps, token, seen = [], None, set()
    while True:
        path = f"genome/accession/{assembly_accession}/sequence_reports?page_size=1000"
        doc = _datasets_json(path + (f"&page_token={urllib.parse.quote(token)}" if token else ""))
        reps += doc.get("reports", [])
        token = doc.get("next_page_token")
        if not token or token in seen:
            return reps
        seen.add(token)


def _fetch_assembly_map_live(assembly_accession: str) -> dict:
    """Build a chromosome map for an assembly not bundled (custom organism) from NCBI Datasets."""
    reps = _sequence_reports(assembly_accession)
    mols = []
    for r in reps:
        if r.get("role") != "assembled-molecule":
            continue
        ref = r.get("refseq_accession") or r.get("genbank_accession")
        if ref:
            mols.append({"chrName": r.get("chr_name", ""), "ucscStyleName": r.get("ucsc_style_name", ""),
                         "refseqAccession": ref, "length": int(r.get("length", 0))})
    if not mols:
        raise CoordError(f"no assembled chromosomes found for assembly {assembly_accession}")
    return {"assemblyAccession": assembly_accession, "molecules": mols, "paginated": True}


def _load_assembly_map(assembly_accession: str, refresh: bool = False) -> dict:
    """Bundled JSON for a curated assembly, else the live Datasets map resolved-and-cached on disk
    (so a custom assembly stays reproducible after first resolve)."""
    p = os.path.join(_assembly_map_dir(), assembly_accession + ".json")
    if not refresh and os.path.isfile(p):
        try:                                                  # a truncated/corrupt shipped map falls through to the live fetch
            with open(p, encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, ValueError, OSError):
            pass
    cached = _read_cache("asm_" + assembly_accession)
    if cached and cached.get("paginated") and not refresh:     # a map built before pagination lacks the organelles
        return cached
    doc = _fetch_assembly_map_live(assembly_accession)
    _write_cache("asm_" + assembly_accession, doc)
    return doc


_MT_NAMES = {"M", "MT", "CHRM", "CHRMT"}
_PLASTID_NAMES = {"C", "CHRC", "PT", "PLTD", "CHLOROPLAST"}


def resolve_chrom(assembly_accession: str, chrom_key: str, refresh: bool = False) -> dict:
    """Resolve a chromosome name to its versioned RefSeq accession + length. Matches against both
    chr_name and ucsc_style_name (with/without 'chr'); M<->MT aliased."""
    amap = _load_assembly_map(assembly_accession, refresh)
    want = {chrom_key.upper(), re.sub(r"^chr", "", chrom_key, flags=re.I).upper()}

    def find(names):
        for mol in amap["molecules"]:
            keys = set()
            for v in (mol.get("chrName", ""), mol.get("ucscStyleName", "")):
                if v:
                    keys.add(v.upper()); keys.add(re.sub(r"^chr", "", v, flags=re.I).upper())
            if keys & names:
                return mol
    mol = find(want)                                           # an exact name is never shadowed by an alias
    for group in (_MT_NAMES, _PLASTID_NAMES):                  # organelle synonyms biologists paste (Ensembl, TAIR, UCSC)
        if mol is None and want & group:
            mol = find(group)
    if mol is not None:
        return {"refseqAccession": mol["refseqAccession"], "length": int(mol["length"]),
                "ucscStyleName": mol.get("ucscStyleName") or mol.get("chrName") or ("chr" + chrom_key),
                "chrName": mol.get("chrName")}
    names = [m.get("chrName") or m.get("ucscStyleName") or "" for m in amap["molecules"]]
    shown = ", ".join(n for n in names[:40] if n) + (f", … (+{len(names) - 40} more)" if len(names) > 40 else "")
    raise CoordError(f"chromosome '{chrom_key}' not found in assembly {assembly_accession} (available: {shown})")


def resolve_assembly(query: str) -> dict:
    """Resolve a free-text query (organism/taxon name OR a GCF/GCA accession) to a pinned assembly, for
    the 'other organism' path. An accession is fully reproducible; a name resolves the current RefSeq
    reference (the resolved versioned accession is what the seal records)."""
    query = (query or "").strip()
    if not query:
        raise CoordError("enter an organism name or an assembly accession (e.g. GCF_000001405.40)")
    if _ASM_ACC_RE.match(query):
        reps = _datasets_json(f"genome/accession/{query}/dataset_report").get("reports", [])
        if not reps:
            raise CoordError(f"assembly '{query}' not found on NCBI")
    else:
        reps = _datasets_json("genome/taxon/" + urllib.parse.quote(query) +
                              "/dataset_report?filters.reference_only=true&filters.assembly_source=RefSeq&page_size=1").get("reports", [])
        if not reps:
            raise CoordError(f"no RefSeq reference assembly found for '{query}' — try an assembly accession (GCF_...)")
    r = reps[0]
    ai = r.get("assembly_info", {})
    return {"organism": r.get("organism", {}).get("organism_name", query),
            "taxid": str(r.get("organism", {}).get("tax_id", "")),
            "assemblyName": ai.get("assembly_name", ""), "assemblyAccession": r.get("accession", "")}


# ---- user-added assemblies (organism -> pinned assembly), resolved ONCE at add-time -------------------
# A lookup store ONLY. It lives OUTSIDE every content-addressed seal key, so a user-pinned assembly seals
# byte-identically to a curated COORD_ASSEMBLIES entry — the versioned accession is the sole reproducibility
# anchor (a bare organism NAME can be promoted to a new RefSeq build over time, so it is frozen here).
_CUSTOM_ASM_PATH = os.path.join(appdirs.user_data_dir(), "custom_assemblies.json")
_CUSTOM_ASM_LOCK = threading.Lock()     # serialise the read-modify-write so concurrent adds cannot clobber an entry

def load_custom_assemblies() -> dict:
    try:
        with open(_CUSTOM_ASM_PATH, encoding="utf-8") as f:
            data = json.load(f)
        return {k: v for k, v in data.items() if isinstance(v, dict) and v.get("assemblyAccession")}
    except Exception:
        return {}

def all_assemblies() -> dict:
    """COORD_ASSEMBLIES unioned with the user's custom store. A curated name always wins a collision, so a user
    cannot shadow a pinned reference assembly. Every downstream lookup reads THIS, so custom + curated genomes
    share one identical code path (and therefore seal identically)."""
    merged = dict(COORD_ASSEMBLIES)
    for org, e in load_custom_assemblies().items():
        merged.setdefault(org, {"assemblyName": e.get("assemblyName", ""),
                                "assemblyAccession": e["assemblyAccession"],
                                "taxid": str(e.get("taxid", "")), "custom": True})
    return merged

def add_custom_assembly(query: str) -> dict:
    """Resolve `query` (organism name OR GCF/GCA accession) ONCE, pin the versioned accession, and persist.
    Returns {organism, assemblyName, assemblyAccession, taxid}. Raises CoordError on a bad name/accession."""
    a = resolve_assembly(query)
    org = a["organism"]
    entry = {"assemblyName": a["assemblyName"], "assemblyAccession": a["assemblyAccession"],
             "taxid": str(a.get("taxid", "")),
             "resolvedUtc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    with _CUSTOM_ASM_LOCK:                                # atomic read-modify-write: a second add cannot read a stale store and clobber this entry
        store = load_custom_assemblies()
        pinned = COORD_ASSEMBLIES.get(org) or store.get(org)
        if pinned:                                            # one assembly per organism: refuse, never a false success
            if pinned["assemblyAccession"] == entry["assemblyAccession"]:
                return {"organism": org, **entry}
            raise CoordError(f"{org} is already pinned to {pinned['assemblyName']} ({pinned['assemblyAccession']}); "
                             f"one assembly per organism is supported, so {entry['assemblyName']} "
                             f"({entry['assemblyAccession']}) was not added")
        store[org] = entry
        os.makedirs(os.path.dirname(_CUSTOM_ASM_PATH), exist_ok=True)
        tmp = _CUSTOM_ASM_PATH + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(store, f, indent=1)
        os.replace(tmp, _CUSTOM_ASM_PATH)
    return {"organism": org, **entry}


def fetch_fasta_range(acc_version: str, start: int, end: int, strand: str = "+") -> str:
    """efetch a subrange of a nuccore accession. seq_start/seq_stop are 1-based inclusive — identical to
    the UCSC browser display, so the pasted numbers pass through with NO conversion. No ENA fallback
    (it ignores the range and would return the whole chromosome). Header rewritten deterministically."""
    strand_num = "2" if strand == "-" else "1"
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc_version, "rettype": "fasta", "retmode": "text",
                                "seq_start": start, "seq_stop": end, "strand": strand_num, "tool": TOOL})
    txt = _get(EUTILS + "efetch.fcgi?" + q, timeout=60)
    if ">" not in txt:
        raise FetchError(f"no FASTA returned for {acc_version}:{start}-{end}")
    seq = "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))
    if len(seq) != end - start + 1:                        # efetch silently clamps an out-of-range stop
        raise FetchError(f"{acc_version}:{start}-{end} returned {len(seq)} bp, expected {end-start+1} "
                         "(coordinates past the end of the sequence?)")
    hdr = f">{acc_version}:{start}-{end}{'(-)' if strand == '-' else ''}"
    return hdr + "\n" + "\n".join(seq[i:i+70] for i in range(0, len(seq), 70)) + "\n"


def retrieve_coords(regions_text: str, assembly_accession: str, assembly_name: str, organism: str,
                    taxid: str = "", strand: str = "+", refresh: bool = False) -> dict:
    """Fetch one or more chromosomal regions of an assembly by coordinate. Returns a dict shaped like
    retrieve() (so the UI reuses _on_fetch) plus regions/assembly identity for the provenance seal."""
    regions = parse_regions(regions_text)
    strand_num = 2 if strand == "-" else 1
    resolved = []
    for reg in regions:                                    # resolve + bounds-check all regions before any fetch
        cm = resolve_chrom(assembly_accession, reg["chromKey"], refresh)
        if reg["end"] > cm["length"]:
            raise CoordError(f"{reg['chromLabel']}:{reg['start']}-{reg['end']} is past the chromosome "
                             f"length ({cm['length']:,} bp on {cm['ucscStyleName']})")
        resolved.append({"chrAccession": cm["refseqAccession"], "start": reg["start"], "stop": reg["end"],
                         "strand": strand_num, "chromLabel": cm["ucscStyleName"] or reg["chromLabel"]})
    key = "coord_" + hashlib.sha256("|".join(                # assembly + taxid + versioned accession + coords + strand
        [assembly_accession, taxid] +                        # so patch assemblies sharing a chr accession never collide
        [f"{r['chrAccession']}:{r['start']}-{r['stop']}:{r['strand']}" for r in resolved]).encode()).hexdigest()[:40]
    if not refresh:
        cached = _read_cache(key)
        if cached and cached.get("fasta"):
            cached["fromCache"] = True
            cached["servedUtc"] = _now()                    # source.retrievedUtc stays the original retrieval time
            return cached
    fastas = []
    for i, (reg, r) in enumerate(zip(regions, resolved)):
        if i:
            time.sleep(0.34)                                # polite pacing: <= 3 req/s without an API key
        fastas.append(fetch_fasta_range(r["chrAccession"], r["start"], r["stop"], strand))
    fasta = "".join(fastas)
    seq_len = sum(r["stop"] - r["start"] + 1 for r in resolved)
    r0 = resolved[0]
    display = f"{r0['chromLabel']}:{r0['start']:,}-{r0['stop']:,}" + (f" (+{len(resolved)-1} more)" if len(resolved) > 1 else "")
    seal_regions = [{"chrAccession": r["chrAccession"], "start": r["start"], "stop": r["stop"], "strand": r["strand"]} for r in resolved]
    source = {
        "accession": r0["chrAccession"],                    # first region's chr accession (reuses accession-path fields)
        "organism": organism, "taxid": taxid,
        "assemblyAccession": assembly_accession, "assemblyName": assembly_name,
        "regions": seal_regions, "coordSystem": "1-based-inclusive", "retrievalType": "coordinate",
        "displayLocus": display, "chromName": r0["chromLabel"],
        "source": "NCBI E-utilities (efetch seq_start/seq_stop)", "endpoint": EUTILS + "efetch.fcgi",
        "retrievedUtc": _now(),
    }
    meta = {"fasta": fasta, "seq_length": seq_len, "organism": organism, "assemblyName": assembly_name,
            "assemblyAccession": assembly_accession, "runType": "coordinate", "regions": resolved,
            "displayLocus": display, "fromCache": False, "source": source,
            "retrievedUtc": source["retrievedUtc"], "endpoint": source["endpoint"]}
    _write_cache(key, meta)
    return meta
