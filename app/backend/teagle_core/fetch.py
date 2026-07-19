"""Live sequence retrieval from NCBI (E-utilities) — real web fetch, server-side.
Verifies metadata (esummary) before pulling the sequence (efetch). Errors are typed
and explicit; a failed or empty fetch is never returned as a sequence."""
from __future__ import annotations
import os, urllib.request, urllib.parse, urllib.error, json, re, datetime

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
    data = json.loads(_get(EUTILS + "esummary.fcgi?" + q))
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


def fetch_fasta(acc: str) -> str:
    acc = validate_accession(acc)
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "fasta",
                                "retmode": "text", "tool": TOOL})
    txt = _get(EUTILS + "efetch.fcgi?" + q)
    if not txt.lstrip().startswith(">"):
        # fall back to ENA once before giving up
        try:
            txt = _get(ENA + urllib.parse.quote(acc))
        except FetchError:
            pass
    if not txt.lstrip().startswith(">"):
        raise FetchError("no FASTA returned (check the accession)")
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
    fasta = fetch_fasta(acc)
    seq = "".join(l for l in fasta.splitlines() if not l.startswith(">"))
    meta["fasta"] = fasta
    meta["seq_length"] = len(seq)
    meta["retrievedUtc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    meta["endpoint"] = EUTILS + "efetch.fcgi"
    meta["fromCache"] = False
    try:                                                       # best-effort gene annotation (exon/intron/CDS); never blocks the fetch
        gm = fetch_features(acc)
        if gm["counts"]["exons"] or gm["counts"]["cds"]:
            meta["features"] = gm
    except Exception:
        pass
    _write_cache(acc, meta)
    return meta
