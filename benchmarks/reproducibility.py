"""Measure whether TEagle's sealed record actually delivers what a reproducibility claim requires.

Three properties are tested, because a manuscript sentence about reproducibility is worth only as much as
the weakest of them:

  1. STABILITY   - the same input, analysed repeatedly, yields a byte-identical seal once the wall-clock
                   timestamp is excluded. A seal that varies run to run identifies nothing.
  2. SENSITIVITY - a different input yields a different seal. A seal that is constant across inputs is
                   equally useless; stability alone is trivially satisfied by a constant.
  3. FIDELITY    - every threshold in the seal equals the value the code actually applies, read back from
                   the function signature that applies it. This is the property that fails silently: a
                   seal can be perfectly stable and perfectly sensitive while recording numbers the
                   detector no longer uses.

    python benchmarks/reproducibility.py

Output: benchmarks/raw/reproducibility.json
"""
from __future__ import annotations
import copy, hashlib, inspect, json, os, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "app", "backend")]

import engine                                              # noqa: E402
from teagle_core import structural, domains, sequtil, primers   # noqa: E402

OUT = os.path.join(ROOT, "benchmarks", "raw", "reproducibility.json")
REPLICATES = 10

# Deterministic synthetic inputs. Real accessions would make the test depend on a network service, which
# is exactly the kind of hidden variability the seal exists to expose.
def _elem(seed: str, n: int) -> str:
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


INPUTS = {
    "A": ">A\n" + _elem("A", 3000),
    "B": ">B\n" + _elem("B", 3000),
    "A_one_base_changed": None,        # filled below - A with a single substitution
}
_a = INPUTS["A"].split("\n", 1)[1]
INPUTS["A_one_base_changed"] = ">A\n" + ("T" if _a[0] != "T" else "G") + _a[1:]
# A primer pair taken straight out of the synthetic template, so it must amplify: a 22-nt forward at 100
# and the reverse complement of a 22-nt window at 400, giving a 322 bp product inside the default window.
_pcr_rev = sequtil.reverse_complement(_a[400:422])

VOLATILE = ("createdUtc",)             # wall-clock; excluded by design, recorded here so it is explicit


def seal_hash(prov: dict) -> str:
    stable = {k: v for k, v in prov.items() if k not in VOLATILE}
    return hashlib.sha256(json.dumps(stable, sort_keys=True).encode()).hexdigest()


def main():
    t0 = time.time()
    result = {"teagle_version": None, "replicates": REPLICATES, "volatile_fields_excluded": list(VOLATILE)}

    # ---- 1. stability ------------------------------------------------------------------------------
    hashes, seals = {}, {}
    for name, seq in INPUTS.items():
        hs = []
        for _ in range(REPLICATES):
            prov = engine.run_analyze({"sequence": seq})["provenance"]
            hs.append(seal_hash(prov))
        hashes[name] = hs
        seals[name] = prov
        result["teagle_version"] = prov["teagleVersion"]
    result["stability"] = {n: {"distinct_seal_hashes": len(set(h)), "stable": len(set(h)) == 1,
                               "seal_sha256": h[0]} for n, h in hashes.items()}

    # ---- 2. sensitivity ----------------------------------------------------------------------------
    distinct = {n: hashes[n][0] for n in INPUTS}
    result["sensitivity"] = {
        "inputs": len(INPUTS),
        "distinct_seals": len(set(distinct.values())),
        "all_distinct": len(set(distinct.values())) == len(INPUTS),
        "one_base_change_alters_seal": distinct["A"] != distinct["A_one_base_changed"],
        "input_sha256": {n: seals[n]["input"]["sha256"] for n in INPUTS},
    }

    # ---- 3. fidelity -------------------------------------------------------------------------------
    # Two distinct questions, and only the second has teeth.
    #
    # (a) Transport. engine._detector_parameters() derives every value with inspect.signature, so a
    #     re-derivation cannot disagree with the code by construction. What it CAN disagree with is the
    #     serialised manifest, if a value is rounded, coerced or cached on the way into JSON. That is what
    #     this compares: live derivation against what actually came back inside the seal.
    #
    # (b) Coverage. The failure this seal was built to prevent is a threshold that decides a call but is
    #     never recorded - the state the manifest was in before, sealing one hand-written parameter while
    #     nineteen others lived only as function defaults. So enumerate every numeric default in the
    #     detectors the seal draws from, and report any that no seal key accounts for.
    sealed = seals["A"]["parameters"]
    live = engine._detector_parameters()
    mismatches = [{"key": k, "sealed": sealed.get(k), "live": v}
                  for k, v in live.items() if sealed.get(k) != v]
    missing_keys = sorted(set(live) - set(sealed))
    extra_keys = sorted(set(sealed) - set(live))

    DETECTORS = [(sequtil, "find_orfs"), (domains, "scan_domains"),
                 (structural, "find_ltr"), (structural, "find_tir"), (structural, "find_polya"),
                 (structural, "find_pbs"), (structural, "find_ppt")]
    # Which seal key each detector parameter must appear under. Written out so coverage is tested by
    # KEY IDENTITY. Matching on value membership - "does this default equal some sealed number?" - is
    # what this check used to do, and it cannot fail usefully: the sealed set is full of small integers,
    # two keys already collide on 12 and two on 2, so 21 keys carry only 17 distinct values, and a
    # genuinely unsealed threshold defaulting to any of them was absorbed silently. A parameter absent
    # from this map is reported unaccounted for, which is the right treatment for a threshold nobody has
    # decided where to seal.
    SEAL_KEY_FOR = {
        ("find_orfs", "min_aa"): "orf_min_aa",
        ("scan_domains", "evalue"): "domain_evalue_max",
        ("scan_domains", "max_orfs"): "domain_orfs_scanned_max",
        ("find_ltr", "k"): "ltr_seed_k",
        ("find_ltr", "min_ltr"): "ltr_min_len",
        ("find_ltr", "min_anchors"): "ltr_min_anchors",
        ("find_tir", "k"): "tir_seed_k",
        ("find_tir", "min_tir"): "tir_min_len",
        ("find_tir", "max_tir"): "tir_max_len",
        ("find_tir", "min_anchors"): "tir_min_anchors",
        ("find_polya", "min_run"): "polya_min_run",
        ("find_pbs", "search"): "pbs_search_window",
        ("find_pbs", "min_ident"): "pbs_min_identity",
        ("find_ppt", "window"): "ppt_window",
        ("find_ppt", "min_len"): "ppt_min_len",
        ("find_ppt", "min_purine"): "ppt_min_purine",
        ("find_ppt", "max_defects"): "ppt_max_defects",
    }
    thresholds, unsealed = [], []
    for mod, fname in DETECTORS:
        fn = getattr(mod, fname, None)
        if fn is None:
            continue
        for pname, param in inspect.signature(fn).parameters.items():
            if param.default is inspect.Parameter.empty or isinstance(param.default, bool) \
                    or not isinstance(param.default, (int, float)):
                continue
            key = SEAL_KEY_FOR.get((fname, pname))
            entry = {"detector": f"{mod.__name__.split('.')[-1]}.{fname}", "parameter": pname,
                     "value": param.default, "expected_seal_key": key}
            thresholds.append(entry)
            if key is None or key not in sealed or sealed[key] != param.default:
                unsealed.append(entry)

    result["fidelity"] = {
        "sealed_parameters": len(sealed),
        "derivation": "engine._detector_parameters() via inspect.signature - never restated as a literal",
        "serialised_matches_live_derivation": not mismatches,
        "mismatches": mismatches,
        "keys_derived_but_not_sealed": missing_keys,
        "keys_sealed_but_not_derived": extra_keys,
        "detector_thresholds_found": len(thresholds),
        "detector_thresholds_unaccounted_for": unsealed,
        "coverage_complete": not unsealed,
        "detectors_examined": [f"{m.__name__.split('.')[-1]}.{f}" for m, f in DETECTORS],
    }

    # ---- 3b. the same coverage question, asked of the assay path -----------------------------------
    # The classification manifest is not the only one a reader has to trust. In-silico PCR carries its own
    # seal, and it had the weaker form of exactly the defect above: it recorded the parameters the caller
    # typed, so a run left at defaults sealed none of the numbers that chose its products. The test is
    # therefore run at defaults on purpose - that is the case where an under-specified seal is invisible.
    pcr_seal = engine.run_pcr({"sequence": INPUTS["A"], "fwd": _a[100:122], "rev": _pcr_rev})["provenance"]
    pcr_sealed = pcr_seal["parameters"]
    pcr_unsealed = []
    for pname, param in inspect.signature(primers.in_silico_pcr).parameters.items():
        if param.default is inspect.Parameter.empty or isinstance(param.default, bool) \
                or not isinstance(param.default, (int, float)):
            continue
        if pname not in pcr_sealed:
            pcr_unsealed.append({"detector": "primers.in_silico_pcr", "parameter": pname,
                                 "value": param.default})
    result["assay_fidelity"] = {
        "sealed_parameters": len(pcr_sealed),
        "run_at_defaults": True,
        "matcher_thresholds_found": sum(
            1 for _n, _p in inspect.signature(primers.in_silico_pcr).parameters.items()
            if _p.default is not inspect.Parameter.empty and not isinstance(_p.default, bool)
            and isinstance(_p.default, (int, float))),
        "matcher_thresholds_unaccounted_for": pcr_unsealed,
        "coverage_complete": not pcr_unsealed,
        "min_anneal_sealed_value": pcr_sealed.get("min_anneal"),
        "min_anneal_live_default": inspect.signature(primers.in_silico_pcr).parameters["min_anneal"].default,
    }

    # ---- seal contents, for the manuscript's description of what a seal carries ---------------------
    p = seals["A"]
    result["seal_contents"] = {
        "software_recorded": [s["name"] for s in p.get("software", [])],
        "databases_recorded": [d["name"] for d in p.get("databases", [])],
        "database_checksums": [d.get("sha256", "")[:12] for d in p.get("databases", [])],
        "parameters_recorded": len(p.get("parameters", {})),
        "references_recorded": len(p.get("references", [])),
        "reference_dois": [r.get("doi") for r in p.get("references", []) if r.get("doi")],
        "top_level_keys": sorted(p),
    }
    result["elapsed_s"] = round(time.time() - t0, 1)
    result["generated"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    tmp = OUT + ".part"
    json.dump(result, open(tmp, "w", encoding="utf-8"), indent=1, sort_keys=True)
    os.replace(tmp, OUT)

    s, sn, f = result["stability"], result["sensitivity"], result["fidelity"]
    print(f"TEagle {result['teagle_version']} - {REPLICATES} replicates per input\n")
    print("1. STABILITY")
    for n, v in s.items():
        print(f"   {n:22s} {v['distinct_seal_hashes']} distinct seal(s) over {REPLICATES} runs -> "
              f"{'STABLE' if v['stable'] else 'UNSTABLE'}")
    print("\n2. SENSITIVITY")
    print(f"   {sn['distinct_seals']} distinct seals over {sn['inputs']} inputs -> "
          f"{'all distinct' if sn['all_distinct'] else 'COLLISION'}")
    print(f"   single-base change alters the seal: {sn['one_base_change_alters_seal']}")
    print("\n3. FIDELITY")
    print(f"   {f['sealed_parameters']} thresholds sealed, each derived by inspect.signature")
    print(f"   serialised seal matches live derivation: {f['serialised_matches_live_derivation']}")
    for m in f["mismatches"]:
        print(f"     ! {m['key']}: seal says {m['sealed']}, live derivation gives {m['live']}")
    for k in f["keys_derived_but_not_sealed"]:
        print(f"     ! derived but absent from the seal: {k}")
    print(f"   coverage: {f['detector_thresholds_found']} numeric thresholds across "
          f"{len(f['detectors_examined'])} detectors, "
          f"{len(f['detector_thresholds_unaccounted_for'])} unaccounted for")
    for u in f["detector_thresholds_unaccounted_for"]:
        print(f"     ! {u['detector']}({u['parameter']}={u['value']}) is applied but not sealed")
    a = result["assay_fidelity"]
    print("\n4. FIDELITY, ASSAY PATH (run at defaults)")
    print(f"   {a['sealed_parameters']} parameters sealed; {a['matcher_thresholds_found']} numeric "
          f"thresholds in the matcher, {len(a['matcher_thresholds_unaccounted_for'])} unaccounted for")
    for u in a["matcher_thresholds_unaccounted_for"]:
        print(f"     ! primers.in_silico_pcr({u['parameter']}={u['value']}) is applied but not sealed")
    print(f"\n-> {OUT}")
    ok = (all(v["stable"] for v in s.values()) and sn["all_distinct"]
          and f["serialised_matches_live_derivation"] and f["coverage_complete"]
          and a["coverage_complete"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
