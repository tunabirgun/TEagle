"""Dev-only build script (NOT on the fetch path). Regenerates the bundled per-assembly chromosome
maps under data/assemblies/ from NCBI Datasets v2 `sequence_reports`. Machine-generated so the
chr -> RefSeq-accession map that seals a coordinate run carries no hand-transcription error.
Run: python app/backend/teagle_core/build_assembly_maps.py   (updates the committed JSON)
Source: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/sequence_reports
"""
import json, os, sys, time, urllib.request, urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "assemblies")
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"

sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))   # .../app/backend
from teagle_core.fetch import COORD_ASSEMBLIES               # single source of truth


def sequence_report(acc):
    url = API + acc + "/sequence_reports?" + urllib.parse.urlencode({"page_size": "1000"})
    req = urllib.request.Request(url, headers={"User-Agent": "TEagle"})
    return json.loads(urllib.request.urlopen(req, timeout=60).read()).get("reports", [])


def build():
    os.makedirs(OUT, exist_ok=True)
    for org, a in COORD_ASSEMBLIES.items():
        acc = a["assemblyAccession"]
        reps = sequence_report(acc)
        mols = []
        for r in reps:
            if r.get("role") != "assembled-molecule":
                continue
            ref = r.get("refseq_accession") or r.get("genbank_accession")
            if not ref:
                continue
            mols.append({"chrName": r.get("chr_name", ""), "ucscStyleName": r.get("ucsc_style_name", ""),
                         "refseqAccession": ref, "length": int(r.get("length", 0))})
        doc = {"organism": org, "assemblyName": a["assemblyName"], "assemblyAccession": acc,
               "source": "NCBI Datasets v2 sequence_reports", "molecules": mols}
        with open(os.path.join(OUT, acc + ".json"), "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=1)
        print(f"{org:28s} {acc:18s} {a['assemblyName']:26s} {len(mols):3d} molecules")
        time.sleep(0.4)


if __name__ == "__main__":
    build()
