"""Recover element coordinates from each accession's GenBank feature table.

The corpus records which element it means but not where that element sits inside the deposited record, so
every case analysed the whole clone or contig (see CORPUS_DEFECT.md). This reads the feature table and
pulls the coordinates of the annotated mobile elements, so a case can be narrowed to the element the
literature label actually refers to.

Matching is by name where the corpus names a family, and by size and rank otherwise. Every assignment
records HOW it was made, so a weak match can be filtered out at scoring time instead of being trusted
silently. A case with no usable feature is left without coordinates and reported, never guessed.

    python benchmarks/extract_coords.py            # writes benchmarks/corpus_coords.tsv
    python benchmarks/extract_coords.py --limit 5

Output: benchmarks/corpus_coords.tsv - the corpus with resolved start/end and a match_method column.
"""
from __future__ import annotations
import argparse, csv, json, os, re, sys, time, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join(ROOT, "benchmarks", "corpus.tsv")
OUT = os.path.join(ROOT, "benchmarks", "corpus_coords.tsv")
CACHE = os.path.join(ROOT, "benchmarks", "raw", "genbank")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Features that delimit a transposable element in a GenBank record.
FEATURE_KEYS = ("mobile_element", "repeat_region", "LTR", "misc_feature")
NAME_QUALS = ("mobile_element_type", "rpt_family", "note", "standard_name", "label", "gene", "product")


def fetch_gb(acc: str) -> str:
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, acc.replace("/", "_") + ".gb")
    if os.path.exists(path) and os.path.getsize(path) > 200:
        return open(path, encoding="utf-8", errors="replace").read()
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "gbwithparts", "retmode": "text"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=120) as r:
                txt = r.read().decode("utf-8", "replace")
            if txt.strip():
                open(path, "w", encoding="utf-8").write(txt)
                return txt
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 + 3 * attempt)                 # E-utilities throttles; back off rather than hammer
    return ""


LOC = re.compile(r"^\s{5}(\S+)\s+(\S.*)$")
QUAL = re.compile(r'^\s{21}/(\w+)=?"?([^"]*)"?')


def parse_features(gb: str):
    """Yield {key, location, quals} for every feature. A minimal parser: the full GenBank grammar is not
    needed, only the location line and the qualifiers that carry a family name."""
    feats, cur, in_ft = [], None, False
    for line in gb.splitlines():
        if line.startswith("FEATURES"):
            in_ft = True
            continue
        if in_ft and line[:1] not in (" ", ""):          # ORIGIN / CONTIG ends the table
            break
        if not in_ft:
            continue
        m = LOC.match(line)
        if m:
            if cur:
                feats.append(cur)
            cur = {"key": m.group(1), "location": m.group(2).strip(), "quals": {}}
            continue
        if cur is None:
            continue
        q = QUAL.match(line)
        if q:
            cur["quals"][q.group(1)] = q.group(2).strip()
            cur["_last"] = q.group(1)
        elif line.startswith(" " * 21) and cur.get("_last"):
            cur["quals"][cur["_last"]] += " " + line.strip().strip('"')
        elif line.startswith(" " * 21):
            cur["location"] += line.strip()
    if cur:
        feats.append(cur)
    return feats


SPAN = re.compile(r"(\d+)\s*\.\.\s*(\d+)")


def span_of(location: str):
    """Outermost span of a location, including join() and complement(). The element occupies everything
    between its first and last coordinate; intervening gaps are other insertions nested inside it, which
    belong to the element's span even though they are not its sequence."""
    nums = [int(a) for pair in SPAN.findall(location) for a in pair]
    if not nums:
        return None
    return min(nums), max(nums), location.strip().startswith("complement")


def name_of(f):
    return " ".join(str(f["quals"].get(k, "")) for k in NAME_QUALS).lower()


def tokens(s: str):
    return {t for t in re.split(r"[^a-z0-9]+", (s or "").lower()) if len(t) > 2 and
            t not in {"family", "element", "the", "and", "ty1", "ty3", "clade", "related", "lineage"}}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default=CORPUS)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.corpus, encoding="utf-8"), delimiter="\t"))
    if a.limit:
        rows = rows[:a.limit]

    gb_cache, stats = {}, {"named": 0, "sole": 0, "ranked": 0, "none": 0, "fetch_failed": 0}
    used_by_acc = {}
    for i, r in enumerate(rows, 1):
        acc = r["accession"].strip()
        r["resolved_start"] = r["resolved_end"] = r["resolved_strand"] = ""
        r["match_method"] = ""
        # A corpus coordinate that is already a plain span needs no lookup.
        m = re.match(r"^\s*(\d[\d,]*)\s*(?:\.\.|-|–)\s*(\d[\d,]*)", r.get("coords") or "")
        if m:
            r["resolved_start"], r["resolved_end"] = m.group(1).replace(",", ""), m.group(2).replace(",", "")
            r["match_method"] = "corpus span"
            stats["named"] += 1
            continue
        # A join()/complement() expression in the corpus was copied from this record's own feature table,
        # so it is authoritative and must not be re-derived by name matching - the family name appears on
        # several features and the matcher would pick a different one.
        gb_span = span_of(r.get("coords") or "") if ".." in (r.get("coords") or "") else None
        if gb_span:
            r["resolved_start"], r["resolved_end"] = str(gb_span[0]), str(gb_span[1])
            r["resolved_strand"] = "-" if gb_span[2] else "+"
            r["match_method"] = "corpus join() span, taken from the record's own feature table"
            stats["named"] += 1
            print(f"[{i}/{len(rows)}] {acc:14} {gb_span[0]}-{gb_span[1]} corpus join() span")
            continue
        try:
            if acc not in gb_cache:
                gb_cache[acc] = parse_features(fetch_gb(acc))
            feats = gb_cache[acc]
        except Exception as e:
            r["match_method"] = f"fetch failed: {type(e).__name__}"
            stats["fetch_failed"] += 1
            print(f"[{i}/{len(rows)}] {acc:14} FETCH FAILED {e}")
            continue

        cands = []
        for f in feats:
            if f["key"] not in FEATURE_KEYS:
                continue
            sp = span_of(f["location"])
            if not sp or sp[1] - sp[0] < 200:
                continue
            cands.append((sp, name_of(f), f["key"]))
        if not cands:
            r["match_method"] = "no delimiting feature in record"
            stats["none"] += 1
            print(f"[{i}/{len(rows)}] {acc:14} no mobile_element/repeat_region feature")
            continue

        want = tokens(r.get("expected_superfamily", "")) | tokens(r.get("coords", ""))
        scored = [(len(want & tokens(nm)), sp, nm, key) for sp, nm, key in cands]
        best = max(scored, key=lambda x: x[0])
        if best[0] > 0:
            method = f"name match ({best[0]} token(s))"
            stats["named"] += 1
        elif len(cands) == 1:
            best = scored[0]
            method = "sole delimiting feature in record"
            stats["sole"] += 1
        else:
            # Deterministic fallback: take the next unused feature, largest first, so several corpus rows
            # on one accession map to different elements instead of all collapsing onto the same one.
            used = used_by_acc.setdefault(acc, set())
            order = sorted(scored, key=lambda x: -(x[1][1] - x[1][0]))
            pick = next((c for c in order if (c[1][0], c[1][1]) not in used), None)
            if pick is None:
                r["match_method"] = f"all {len(cands)} features already assigned"
                stats["none"] += 1
                continue
            best, method = pick, f"rank by size, {len(cands)} candidates - WEAK"
            stats["ranked"] += 1
        (s, e, comp), nm = best[1], best[2]
        used_by_acc.setdefault(acc, set()).add((s, e))
        r["resolved_start"], r["resolved_end"] = str(s), str(e)
        r["resolved_strand"] = "-" if comp else "+"
        r["match_method"] = method
        print(f"[{i}/{len(rows)}] {acc:14} {s}-{e} ({e - s} bp) {method}")

    cols = list(rows[0].keys())
    with open(a.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    resolved = sum(1 for r in rows if r["resolved_start"])
    weak = sum(1 for r in rows if "WEAK" in r["match_method"])
    print(f"\nresolved {resolved}/{len(rows)} ({weak} weak) -> {a.out}")
    print("  " + json.dumps(stats))
    return 0


if __name__ == "__main__":
    sys.exit(main())
