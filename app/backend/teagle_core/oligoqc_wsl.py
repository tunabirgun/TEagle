"""Out-of-process ViennaRNA engine for primer secondary-structure QC (optional WSL component).

ViennaRNA is deliberately not bundled: its licence forbids redistribution for a fee, which a copyleft
licence cannot carry, so it ships as a one-click backend component in its OWN micromamba env
(`teagle-vrna`) rather than inside the installer.

The remote worker mirrors `oligoqc`'s in-process computation exactly — same ViennaRNA Python API, same
DNA Mathews-2004 parameters, same temperature and monovalent salt — and the component pins the same
ViennaRNA minor version as the optional pip package. A primer therefore reports the SAME ΔG whichever
route supplied it; a version or parameter skew here would make a reported number depend on how the user
happened to install, which the provenance record could not express.

Solved against the shared `te` env the resolver picks viennarna 2.4.7 py36 — unusable bindings, wrong
version — and could re-solve hmmer/rmblast, which annotate runs seal. Hence the separate env.
"""
from __future__ import annotations
import json

from . import wsl

_ENV = "teagle-vrna"
_avail = None

# Mirrors oligoqc._vrna_md / _vrna_mfe / _vrna_binding. Kept as source text (not imported) because it
# executes under the WSL env's interpreter, not this one.
_WORKER = r'''
import json, sys
import RNA

job = json.load(sys.stdin)
cond = job["conditions"]
RNA.params_load_DNA_Mathews2004()          # process-global, matching the in-process path


def _md():
    m = RNA.md()
    m.temperature = cond["temp_c"]
    try:
        m.salt = cond["mv_conc"] / 1000.0   # mol/L monovalent (ViennaRNA >=2.6)
    except Exception:
        pass
    return m


def _mfe(s):
    return RNA.fold_compound(s, _md()).mfe()[1]


out = {"version": RNA.__version__, "mfe": {}, "binding": {}}
for s in job.get("mfe", []):
    out["mfe"][s] = _mfe(s)
for a, b in job.get("binding", []):
    out["binding"][a + "&" + b] = _mfe(a + "&" + b) - _mfe(a) - _mfe(b)
json.dump(out, sys.stdout)
'''


def available(refresh: bool = False) -> bool:
    """True when the optional ViennaRNA backend component is installed and importable."""
    global _avail
    if _avail is None or refresh:
        try:
            rc, out, _ = wsl._wsl_script(
                f'{wsl._MM} run -n {_ENV} python -c "import RNA;print(\'VRNA_OK\')" 2>/dev/null\n', timeout=60)
            _avail = rc == 0 and "VRNA_OK" in out
        except Exception:
            _avail = False
    return bool(_avail)


def version():
    try:
        rc, out, _ = wsl._wsl_script(
            f'{wsl._MM} run -n {_ENV} python -c "import RNA;print(RNA.__version__)" 2>/dev/null\n', timeout=60)
        return out.strip().splitlines()[-1].strip() if rc == 0 and out.strip() else None
    except Exception:
        return None


def compute(mfe_seqs, binding_pairs, conditions, timeout: int = 180):
    """One round trip for every fold a primer batch needs.

    Returns {("mfe", seq): dg, ("bind", a, b): dg} — the key shape oligoqc's remote cache expects —
    or {} when the component is absent or the worker fails. Batched deliberately: a per-metric call
    would be ~9 WSL round trips per primer pair, which is unusable for a full design.
    """
    mfe_seqs = sorted({s for s in mfe_seqs if s})
    binding_pairs = sorted({(a, b) for a, b in binding_pairs if a and b})
    if not mfe_seqs and not binding_pairs:
        return {}
    payload = json.dumps({"conditions": conditions, "mfe": mfe_seqs,
                          "binding": [list(p) for p in binding_pairs]})
    # both the worker and its payload travel as heredocs inside one STDIN-delivered script: wsl.exe
    # rebuilds an inline `bash -lc` command line and mangles quotes, $() and newlines.
    script = (f"cat > /tmp/teagle_vrna_worker.py <<'TEAGLE_PY'\n{_WORKER}\nTEAGLE_PY\n"
              f"cat > /tmp/teagle_vrna_job.json <<'TEAGLE_JSON'\n{payload}\nTEAGLE_JSON\n"
              f'{wsl._MM} run -n {_ENV} python /tmp/teagle_vrna_worker.py < /tmp/teagle_vrna_job.json\n')
    try:
        rc, out, _ = wsl._wsl_script(script, timeout=timeout)
        if rc != 0:
            return {}
        start = out.find("{")
        data = json.loads(out[start:out.rfind("}") + 1]) if start >= 0 else {}
    except Exception:
        return {}
    finally:
        # the job.json holds the primer sequences; clear it even when the round trip times out (a timeout
        # kills the script before an in-script rm would run), mirroring wsl.py's staged-dir cleanup
        try:
            wsl._wsl_script('rm -f /tmp/teagle_vrna_worker.py /tmp/teagle_vrna_job.json\n', timeout=30)
        except Exception:
            pass
    cache = {("mfe", s): dg for s, dg in (data.get("mfe") or {}).items()}
    for key, dg in (data.get("binding") or {}).items():
        a, _, b = key.partition("&")
        cache[("bind", a, b)] = dg
    return cache
