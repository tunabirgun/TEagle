"""Run provenance manifest — packs exact tool versions, checksums, params, environment
with every result (see Deliverables/09). Values are REAL (captured at runtime)."""
from __future__ import annotations
import sys, platform, json, hashlib, datetime
from . import __version__ as TEAGLE_VERSION
from .primers import PRIMER3_VERSION            # guarded: "unavailable" if primer3 fails to load
from .domains import PYHMMER_VERSION            # guarded likewise

# wall-clock / machine-specific fields that must NOT enter the content-addressed seal
_SEAL_EXCLUDE_TOP = ("createdUtc", "environment")
# volatile / label-only input fields: retrieval time, and the raw FASTA-header id (differs between the
# NCBI and ENA-fallback fetch paths for the SAME accession+sequence). The sequence sha256 + the resolved
# accession/organism/taxid still seal the scientific identity.
_SEAL_EXCLUDE_INPUT = ("retrievedUtc", "id", "assemblyName", "displayLocus", "chromName",
                       "source", "endpoint", "sourceUrl")   # serving-DB labels: recorded, never sealed (NCBI vs ENA fallback
                                                            # must not change the content-addressed seal for identical bytes)


def _software(run_type: str) -> list:
    """Only the tools that actually produced THIS run enter the software list (and thus the seal),
    so an unused tool's version/availability can never change a byte-identical run's checksum."""
    sw = [{"name": "python", "version": platform.python_version()},
          {"name": "TEagle core", "version": TEAGLE_VERSION}]
    if run_type == "primer":
        sw.insert(0, {"name": "primer3 (primer3-py)", "version": PRIMER3_VERSION})
    elif run_type == "analysis":
        sw.insert(0, {"name": "HMMER (pyhmmer)", "version": PYHMMER_VERSION})
    return sw


def build_manifest(run_type: str, input_seq: str, input_id: str, params: dict,
                   not_run: list | None = None, source: dict | None = None,
                   references: list | None = None, databases: list | None = None):
    inp = {
        "id": input_id,
        "length": len(input_seq),
        "sha256": hashlib.sha256(input_seq.encode()).hexdigest(),
    }
    if isinstance(source, dict):                # accession/coordinate retrieval provenance, if fetched
        inp.update({k: source[k] for k in ("accession", "organism", "taxid", "source", "endpoint", "retrievedUtc",
                    "assemblyAccession", "regions", "coordSystem", "retrievalType",     # coordinate: sealed identity
                    "assemblyName", "displayLocus", "chromName") if k in source})        # coordinate: recorded-but-excluded labels
    m = {
        "teagleVersion": TEAGLE_VERSION,
        "schemaVersion": "0.1",
        "runType": run_type,
        "createdUtc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "input": inp,
        "software": _software(run_type),
        "databases": databases or [],
        "parameters": params,
        "references": references or [],
        "environment": {
            "os": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
        },
        "notRun": not_run or [
            "DB-backed family naming (Dfam / RepeatMasker) — available via the WSL backend; not run in this step",
            "External NCBI Primer-BLAST",
            "De novo family discovery",
        ],
    }
    # Content-addressed seal over the scientific content only: exclude wall-clock time, the machine,
    # and the volatile retrieval timestamp, so an identical run (same input, versions, params, DBs)
    # yields an identical seal across machines and refetches.
    seal_src = {k: v for k, v in m.items() if k not in _SEAL_EXCLUDE_TOP}
    if isinstance(seal_src.get("input"), dict):
        seal_src["input"] = {k: v for k, v in seal_src["input"].items() if k not in _SEAL_EXCLUDE_INPUT}
    m["manifestSha256"] = hashlib.sha256(
        json.dumps(seal_src, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return m
