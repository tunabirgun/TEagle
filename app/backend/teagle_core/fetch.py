"""Live sequence retrieval from NCBI (E-utilities) — real web fetch, server-side.
Verifies metadata (esummary) before pulling the sequence (efetch). Errors are typed
and explicit; a failed or empty fetch is never returned as a sequence."""
from __future__ import annotations
import os, urllib.request, urllib.parse, urllib.error, json, re, datetime, hashlib, time

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ENA = "https://www.ebi.ac.uk/ena/browser/api/fasta/"
TOOL = "TEagle"
ACC_RE = re.compile(r"^[A-Za-z]{1,6}_?\d+(\.\d+)?$")   # GenBank/RefSeq nucleotide accession

# On-disk cache: an accession version is immutable, so a fetched sequence is reused instead of
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


def _get(url: str, timeout: int = 25) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "TEagle/0.1 (+local)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        raise FetchError(f"HTTP {e.code} from source")
    except urllib.error.URLError as e:
        raise FetchError(f"network error: {getattr(e, 'reason', e)}")
    except Exception as e:                                 # timeout etc.
        raise FetchError(str(e))


def validate_accession(acc: str) -> str:
    acc = (acc or "").strip()
    if not acc:
        raise FetchError("no accession provided")
    if not ACC_RE.match(acc):
        raise FetchError(f"'{acc}' is not a valid nucleotide accession (e.g. M11240, NC_003075.7)")
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


def _dedup(ivs: list) -> list:
    seen, out = set(), []
    for iv in ivs:
        k = (iv["start"], iv["end"])
        if k not in seen:
            seen.add(k)
            out.append(iv)
    return sorted(out, key=lambda x: x["start"])


def build_gene_model(feats: list) -> dict:
    """Reduce a feature table to exon/intron/CDS segments. Introns and exons are taken from
    explicit features when present, else derived from spliced (joined) CDS/mRNA coordinates."""
    def collect(kind):
        out = []
        for f in feats:
            if f["type"] == kind:
                label = f["qualifiers"].get("number") or f["qualifiers"].get("gene") or f["qualifiers"].get("product") or ""
                for iv in f["intervals"]:
                    out.append({**iv, "note": label})
        return out
    exons, introns, cds, mrna = collect("exon"), collect("intron"), collect("CDS"), collect("mRNA")
    derived_exons = derived_introns = False
    if not exons and (mrna or cds):                            # derive exons from a spliced transcript/CDS
        exons = [dict(x) for x in (mrna or cds)]
        derived_exons = True
    exons = _dedup(exons)
    if not introns and len(exons) > 1:                         # introns = gaps between consecutive exons
        for a, b in zip(exons, exons[1:]):
            if b["start"] > a["end"]:
                introns.append({"start": a["end"], "end": b["start"], "strand": a.get("strand", "+"), "note": ""})
        derived_introns = bool(introns)
    exons, introns, cds = _dedup(exons), _dedup(introns), _dedup(cds)
    return {
        "exons": exons, "introns": introns, "cds": cds,
        "counts": {"exons": len(exons), "introns": len(introns), "cds": len(cds)},
        "derived_exons": derived_exons, "derived_introns": derived_introns,
        "source": "NCBI feature table (efetch rettype=ft)",
    }


def fetch_features(acc: str) -> dict:
    """Fetch and parse the accession's annotation feature table into a gene model."""
    acc = validate_accession(acc)
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "ft", "retmode": "text", "tool": TOOL})
    return build_gene_model(parse_feature_table(_get(EUTILS + "efetch.fcgi?" + q)))


def retrieve(acc: str, refresh: bool = False) -> dict:
    """Full retrieval: metadata + sequence, with provenance fields. Cached on disk per accession
    so the same organism is not re-downloaded on every run; `refresh=True` forces a fresh pull."""
    acc = validate_accession(acc)
    if not refresh:
        cached = _read_cache(acc)
        if cached and cached.get("fasta"):
            cached["fromCache"] = True
            cached["servedUtc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
            return cached
    meta = resolve(acc)
    served = []
    fasta = fetch_fasta(acc, served)
    seq = "".join(l for l in fasta.splitlines() if not l.startswith(">"))
    meta["fasta"] = fasta
    meta["seq_length"] = len(seq)
    meta["retrievedUtc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    if served:                                                 # record the DB that actually served the sequence
        meta["source"], meta["endpoint"] = served[0]
    else:
        meta["endpoint"] = EUTILS + "efetch.fcgi"
    meta["sourceUrl"] = ("https://www.ebi.ac.uk/ena/browser/view/" + acc if meta.get("source", "").startswith("ENA")
                         else "https://www.ncbi.nlm.nih.gov/nuccore/" + acc)
    meta["fromCache"] = False
    try:                                                       # best-effort gene annotation (exon/intron/CDS); never blocks the fetch
        gm = fetch_features(acc)
        if gm["counts"]["exons"] or gm["counts"]["cds"]:
            meta["features"] = gm
    except Exception:
        pass
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
    "Sus scrofa":                {"assemblyName": "Sscrofa11.1","assemblyAccession": "GCF_000003025.6", "taxid": "9823"},
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


def _fetch_assembly_map_live(assembly_accession: str) -> dict:
    """Build a chromosome map for an assembly not bundled (custom organism) from NCBI Datasets."""
    reps = _datasets_json(f"genome/accession/{assembly_accession}/sequence_reports?page_size=1000").get("reports", [])
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
    return {"assemblyAccession": assembly_accession, "molecules": mols}


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
    if cached and not refresh:
        return cached
    doc = _fetch_assembly_map_live(assembly_accession)
    _write_cache("asm_" + assembly_accession, doc)
    return doc


def resolve_chrom(assembly_accession: str, chrom_key: str, refresh: bool = False) -> dict:
    """Resolve a chromosome name to its versioned RefSeq accession + length. Matches against both
    chr_name and ucsc_style_name (with/without 'chr'); M<->MT aliased."""
    amap = _load_assembly_map(assembly_accession, refresh)
    want = {chrom_key.upper(), re.sub(r"^chr", "", chrom_key, flags=re.I).upper()}
    if want & {"M", "MT", "CHRM", "CHRMT"}:
        want |= {"M", "MT", "CHRM"}
    for mol in amap["molecules"]:
        keys = set()
        for v in (mol.get("chrName", ""), mol.get("ucscStyleName", "")):
            if v:
                keys.add(v.upper()); keys.add(re.sub(r"^chr", "", v, flags=re.I).upper())
        if keys & want:
            return {"refseqAccession": mol["refseqAccession"], "length": int(mol["length"]),
                    "ucscStyleName": mol.get("ucscStyleName") or mol.get("chrName") or ("chr" + chrom_key),
                    "chrName": mol.get("chrName")}
    raise CoordError(f"chromosome '{chrom_key}' not found in assembly {assembly_accession}")


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
            cached["source"]["retrievedUtc"] = _now()       # display-only; excluded from the seal
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
