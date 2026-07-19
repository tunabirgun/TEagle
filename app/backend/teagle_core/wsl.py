"""Managed WSL2 backend — runs the Linux-only annotation stack (RepeatMasker + Dfam)
for family-level TE classification (Layer A homology).

Security: native Windows never runs a Linux shell built from user input. Commands go to
WSL via fixed argument vectors; the user's sequence is piped in as STDIN (data, never
part of a command); the only interpolated values are a self-generated run-id and a
strictly-validated species token. No shell-string concatenation of untrusted input.
"""
from __future__ import annotations
import os, subprocess, re, secrets, threading
from . import appdirs

# Suppress the console-window flash when the windowed GUI spawns wsl.exe (Windows only; 0 elsewhere).
_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0

_ENV = "$HOME/micromamba/envs/te"
_MM = "$HOME/bin/micromamba"
_SPECIES_RE = re.compile(r"^[A-Za-z0-9 _.-]{1,60}$")
_distro_cache = None
_DISTRO_FILE = os.path.join(appdirs.user_data_dir(), "wsl_distro.txt")   # the installer reads this to clean the RIGHT distro


def _persist_distro(name: str | None):
    try:                                            # best-effort record of the distro the backend actually uses
        if name:
            with open(_DISTRO_FILE, "w", encoding="utf-8") as f:
                f.write(name)
    except Exception:
        pass


def _decode(b: bytes) -> str:
    # wsl.exe list commands emit UTF-16LE; normal command output is UTF-8
    if b[:1] == b"\xff" or (len(b) > 1 and b[1] == 0):
        return b.decode("utf-16-le", "ignore")
    return b.decode("utf-8", "ignore")


def resolve_distro():
    """Return the DEFAULT WSL distro (marked '*' in `wsl -l -v`), else the first listed."""
    global _distro_cache
    if _distro_cache is not None:
        return _distro_cache
    try:
        out = subprocess.run(["wsl.exe", "-l", "-v"], capture_output=True, timeout=15, creationflags=_NO_WINDOW)
        txt = _decode(out.stdout).replace("\x00", "")
        default, names = None, []
        for line in txt.splitlines()[1:]:            # skip the header row
            s = line.strip()
            if not s:
                continue
            star = s.startswith("*")
            parts = s.lstrip("*").split()
            if parts:
                names.append(parts[0])
                if star:
                    default = parts[0]
        _distro_cache = default or (names[0] if names else None)
    except Exception:
        _distro_cache = None
    return _distro_cache


def _wsl(script: str, stdin: bytes | None = None, timeout: int = 600):
    """Run a bash -lc script inside the WSL distro. `script` must contain no untrusted input."""
    distro = resolve_distro()
    if not distro:
        raise RuntimeError("no WSL distribution found")
    cmd = ["wsl.exe", "-d", distro, "--", "bash", "-lc", script]
    p = subprocess.run(cmd, input=stdin, capture_output=True, timeout=timeout, creationflags=_NO_WINDOW)
    return p.returncode, _decode(p.stdout), _decode(p.stderr)


def _wsl_script(script: str, timeout: int = 90):
    """Run a multi-line bash script delivered via STDIN (not as a `-c` argument).

    wsl.exe rebuilds the Windows command line and mangles embedded double-quotes / $()/newlines
    in an inline `bash -lc <script>` argument, so a probe using "$VAR" or command substitution
    silently misbehaves. Feeding the script as STDIN bytes to a login shell avoids that entirely —
    the same reason the install script is delivered by `cat > file`."""
    distro = resolve_distro()
    if not distro:
        raise RuntimeError("no WSL distribution found")
    cmd = ["wsl.exe", "-d", distro, "--", "bash", "-l", "-s"]
    p = subprocess.run(cmd, input=script.encode(), capture_output=True, timeout=timeout, creationflags=_NO_WINDOW)
    return p.returncode, _decode(p.stdout), _decode(p.stderr)


def available() -> dict:
    global _distro_cache
    distro = resolve_distro()
    if not distro:
        return {"wsl2": False, "distro": None, "error": "WSL not installed / no distro"}
    try:
        rc, out, err = _wsl("echo ok", timeout=30)
        ok = rc == 0 and "ok" in out
        if not ok:
            _distro_cache = None                     # cached distro no longer usable -> re-resolve on the next probe
        else:
            _persist_distro(distro)                  # record the working distro so uninstall/clean targets the right one
        return {"wsl2": ok, "distro": distro, "error": None if rc == 0 else err.strip()[:120]}
    except Exception as e:
        _distro_cache = None                          # distro changed/removed mid-session -> drop the stale cache
        return {"wsl2": False, "distro": distro, "error": str(e)[:120]}


def env_status() -> dict:
    """Report the annotation stack state inside WSL (RepeatMasker version, Dfam library, minimap2)."""
    av = available()
    st = {**av, "repeatmasker": None, "engine": None, "dfam": False, "minimap2": None, "ready": False}
    if not av["wsl2"]:
        return st
    try:
        rc, out, err = _wsl(
            f'{_MM} run -n te RepeatMasker -v 2>/dev/null | head -1; '
            f'ls {_ENV}/share/RepeatMasker/Libraries/famdb/*.h5 2>/dev/null | head -1; '
            f'[ -x "{_ENV}/bin/minimap2" ] && {_ENV}/bin/minimap2 --version 2>/dev/null',
            timeout=60)
        m = re.search(r"RepeatMasker version ([\w.]+)", out)
        st["repeatmasker"] = m.group(1) if m else None
        st["dfam"] = ".h5" in out
        mm = re.search(r"^(\d+\.\d+[\w.-]*)$", out.strip().splitlines()[-1]) if out.strip() else None
        st["minimap2"] = mm.group(1) if mm else None
        st["ready"] = bool(st["repeatmasker"] and st["dfam"])
    except Exception as e:
        st["error"] = str(e)[:120]
    return st


def parse_out(text: str):
    """Parse a RepeatMasker .out table into structured hits."""
    hits = []
    for line in text.splitlines():
        f = line.split()
        if len(f) < 11 or not f[0].isdigit():
            continue                                  # skip headers / blank / non-data
        strand = "-" if f[8] in ("C", "-") else "+"
        try:
            q_start, q_end = int(f[5]), int(f[6])
            score, divergence = int(f[0]), float(f[1])       # one malformed line drops its own hit, not the whole table
        except ValueError:
            continue
        hits.append({
            "score": score, "divergence": divergence,
            "query": f[4], "q_start": q_start - 1, "q_end": q_end, "strand": strand,
            "family": f[9], "class_family": f[10],
        })
    return hits


# ============================ Managed install (component-wise) ============================
# Design goals (obs: "flawless for all users"): run from a SAFE cwd (WSL starts in the Windows
# DrvFs mount, where relative writes are permission-denied — every step cd's to $HOME first);
# pin the VERSIONED Dfam path (not the moving 'current' pointer) with EMBEDDED md5 trust anchors
# (R-SEC3: pinned checksums, not runtime-fetched); resumable downloads; each step idempotent so a
# repair re-run is safe. Steps share one prelude and are composed for either "install all" or a
# single-component "repair".
_DFAM_BASE = "https://www.dfam.org/releases/Dfam_4.0/families/FamDB"   # pinned & versioned
# (filename, md5-of-the-.gz) — trust anchors captured from Dfam_4.0/*.md5 (verified 2026-07-19)
_DFAM_FILES = {
    "dfam_root":    ("dfam40.0.h5", "234d177775f1bf3445b1fe146bc6e65e"),
    "dfam_curated": ("dfam40.curated.consensus.0.h5", "7892e18016fc820264e625cbb9ec607b"),
}

_PRELUDE = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || { echo "[teagle] FAILED: cannot cd to HOME"; exit 1; }
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"
FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
LOG="$HOME/teagle_wsl_install.log"; : > "$LOG"; exec > >(tee -a "$LOG") 2>&1
LOCK="$HOME/.teagle_install.lock"
# reap a lock orphaned by a crash / reboot / wsl --shutdown / timeout: its recorded PID is dead.
if [ -d "$LOCK" ] && { [ ! -f "$LOCK/pid" ] || ! kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; }; then
  rm -rf "$LOCK" 2>/dev/null
fi
if ! mkdir "$LOCK" 2>/dev/null; then echo "[teagle] FAILED: install already running"; exit 1; fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK" 2>/dev/null' EXIT
fail(){ echo "[teagle] FAILED: $1"; exit 1; }
echo "[teagle] START $(date -u +%FT%TZ)"
'''


def _dfam_step(key: str) -> str:
    fname, md5 = _DFAM_FILES[key]
    return f'''
echo "[teagle] STEP {key} START"
mkdir -p "$FAMDIR" || fail "mkdir famdb"
cd "$FAMDIR" || fail "cd famdb"
if [ -f "{fname}" ]; then echo "[teagle] {fname} already present"; else
  avail_gb=$(df -BG --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9'); avail_gb=${{avail_gb:-99}}
  [ "$avail_gb" -ge 2 ] || fail "insufficient disk space (${{avail_gb}}G free, need ~2G)"
  echo "[teagle] downloading {fname}.gz (resumable)"
  # if the .gz is already complete, a resumed request returns HTTP 416 and --fail errors — don't abort;
  # fall through to the md5 gate, which validates a good file or removes a bad one for a clean retry.
  curl -L --fail -C - --retry 5 --retry-delay 5 -o "{fname}.gz" "{_DFAM_BASE}/{fname}.gz" || {{ [ -f "{fname}.gz" ] || fail "download {fname}"; }}
  echo "{md5}  {fname}.gz" | md5sum -c - || {{ rm -f "{fname}.gz"; fail "md5 mismatch {fname} (removed; re-run to retry)"; }}
  echo "[teagle] md5 OK {fname}"
  gunzip -f "{fname}.gz" || fail "gunzip {fname}"
fi
cd "$HOME"
echo "[teagle] STEP {key} OK"
'''


# key -> idempotent bash body (prelude sets MM/ENV/FAMDIR/fail; every step is safe to re-run)
_STEP = {
    "micromamba": r'''
echo "[teagle] STEP micromamba START"
if [ -x "$MM" ]; then echo "[teagle] micromamba already present"; else
  mkdir -p "$HOME/bin" || fail "cannot create $HOME/bin"
  curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/latest" | tar -xj -C "$HOME" bin/micromamba || fail "micromamba download/extract"
  [ -x "$MM" ] || fail "micromamba missing after extract"
fi
echo "[teagle] STEP micromamba OK"
''',
    "repeatmasker": r'''
echo "[teagle] STEP repeatmasker START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if "$MM" run -n te RepeatMasker -v >/dev/null 2>&1; then echo "[teagle] RepeatMasker already present"; else
  [ -d "$ENV/conda-meta" ] || "$MM" create -y -n te -c conda-forge -c bioconda || fail "create te env"
  if ! "$MM" install -y -n te -c conda-forge -c bioconda repeatmasker; then
    # never `env remove` here: the Dfam .h5 libraries (multi-GB) live INSIDE this env prefix,
    # so recreating the env would silently wipe them. Force a clean package reinstall instead.
    echo "[teagle] install failed — forcing a clean reinstall of repeatmasker (env + Dfam preserved)"
    "$MM" install --force-reinstall -y -n te -c conda-forge -c bioconda repeatmasker || fail "install repeatmasker"
  fi
  "$MM" run -n te RepeatMasker -v >/dev/null 2>&1 || fail "RepeatMasker not runnable after install"
fi
echo "[teagle] STEP repeatmasker OK"
''',
    "minimap2": r'''
echo "[teagle] STEP minimap2 START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/minimap2" ]; then echo "[teagle] minimap2 already present"; else
  [ -d "$ENV/conda-meta" ] || "$MM" create -y -n te -c conda-forge -c bioconda || fail "create te env"
  "$MM" install -y -n te -c conda-forge -c bioconda minimap2 || fail "install minimap2"
  [ -x "$ENV/bin/minimap2" ] || fail "minimap2 missing after install"
fi
echo "[teagle] STEP minimap2 OK"
''',
    "dfam_root": _dfam_step("dfam_root"),
    "dfam_curated": _dfam_step("dfam_curated"),
    "famdb_conf": r'''
echo "[teagle] STEP famdb_conf START"
mkdir -p "$FAMDIR" || fail "mkdir famdb"
FDB=$(ls -d "$ENV"/share/famdb-* 2>/dev/null | head -1)
[ -n "$FDB" ] || FDB="$ENV/share/RepeatMasker"
[ -d "$FDB" ] || fail "famdb tool dir not found (repair RepeatMasker first)"
printf '[famdb]\nFAMDB_DATA_DIR = %s\n' "$FAMDIR" > "$FDB/famdb.conf" || fail "write famdb.conf"
if "$MM" run -n te famdb.py info >/tmp/teagle_famdb.txt 2>&1; then
  grep -iE "version|consensus|families" /tmp/teagle_famdb.txt | head -3
else
  echo "[teagle] note: famdb.py did not validate yet (Dfam libraries may be incomplete — repair them, then re-check integrity)"
fi
echo "[teagle] STEP famdb_conf OK"
''',
}

_ALL_STEPS = ["micromamba", "repeatmasker", "minimap2", "dfam_root", "dfam_curated", "famdb_conf"]

# component metadata surfaced to the install dialog (order = install order)
_COMP_META = [
    ("wsl2",         "WSL2 + Linux distro",       False, "Windows Subsystem for Linux — hosts the Dfam / RepeatMasker stack."),
    ("micromamba",   "micromamba (conda)",        True,  "Small conda package manager, installed under your Linux home."),
    ("repeatmasker", "RepeatMasker",              True,  "Homology-based TE annotator that names Dfam families."),
    ("minimap2",     "minimap2",                  True,  "Splice-aware aligner for de-novo exon / intron detection."),
    ("dfam_root",    "Dfam 4.0 root library",     True,  "Dfam root partition (dfam40.0.h5)."),
    ("dfam_curated", "Dfam 4.0 curated library",  True,  "Dfam curated consensus partition."),
    ("famdb_conf",   "FamDB configuration",       True,  "Points RepeatMasker at the downloaded Dfam library."),
]


def _build_script(keys) -> str:
    return _PRELUDE + "".join(_STEP[k] for k in keys) + '\necho "[teagle] DONE $(date -u +%FT%TZ)"\n'


# LIVE = a lock held by a still-running install; a lock whose recorded PID is dead reads FREE (reaped on next run)
_LOCK_LIVE = ('L="$HOME/.teagle_install.lock"; '
              'if [ -d "$L" ] && [ -f "$L/pid" ] && kill -0 "$(cat "$L/pid" 2>/dev/null)" 2>/dev/null; '
              'then echo LIVE; else echo FREE; fi\n')

_install_thread = None


def _run_attached():
    # Hold the WSL session open for the whole run: a detached (nohup) process is reaped by WSL2
    # when the launching session exits, so the script must run inside a long-lived attached call.
    try:
        _wsl('bash "$HOME/teagle_wsl_install.sh"', timeout=7200)
    except Exception:
        pass


def _start_thread(script: str) -> dict:
    global _install_thread
    if _install_thread is not None and _install_thread.is_alive():
        return {"started": False, "error": "install already running"}
    try:
        # only a LIVE install blocks a new one; a lock whose PID is dead (crash/reboot) is reaped by _PRELUDE
        _, out, _ = _wsl_script(_LOCK_LIVE, timeout=20)
        if "LIVE" in out:
            return {"started": False, "error": "install already running"}
        _wsl('cat > "$HOME/teagle_wsl_install.sh"', stdin=script.encode(), timeout=30)
    except Exception as e:
        return {"started": False, "error": str(e)[:160]}
    _install_thread = threading.Thread(target=_run_attached, daemon=True)
    _install_thread.start()
    return {"started": True}


def start_install() -> dict:
    """Install the whole annotation stack in a session-holding thread. Returns immediately;
    progress is read via install_log() / components_status()."""
    if not available()["wsl2"]:
        return {"started": False, "error": "WSL2 not available"}
    return _start_thread(_build_script(_ALL_STEPS))


def repair_component(key: str) -> dict:
    """Re-run a single component's idempotent install step (used by the dialog's per-package repair)."""
    if key not in _STEP:
        return {"started": False, "error": f"unknown component: {key}"}
    if not available()["wsl2"]:
        return {"started": False, "error": "WSL2 not available"}
    return _start_thread(_build_script([key]))


_STATUS_PROBE = r'''cd "$HOME" 2>/dev/null
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"; FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
echo "micromamba=$([ -x "$MM" ] && echo 1 || echo 0)"
rmv=$("$MM" run -n te RepeatMasker -v 2>/dev/null | grep -oiE 'version [0-9][0-9.]*' | head -1 | awk '{print $2}'); echo "repeatmasker=${rmv:-0}"
mmv=$([ -x "$ENV/bin/minimap2" ] && "$ENV/bin/minimap2" --version 2>/dev/null); echo "minimap2=${mmv:-0}"
echo "dfam_root=$([ -f "$FAMDIR/dfam40.0.h5" ] && echo 1 || echo 0)"
echo "dfam_curated=$([ -f "$FAMDIR/dfam40.curated.consensus.0.h5" ] && echo 1 || echo 0)"
echo "famdb_conf=$( { ls "$ENV"/share/famdb-*/famdb.conf >/dev/null 2>&1 || [ -f "$ENV/share/RepeatMasker/famdb.conf" ]; } && echo 1 || echo 0)"
if [ -d "$HOME/.teagle_install.lock" ] && [ -f "$HOME/.teagle_install.lock/pid" ] && kill -0 "$(cat "$HOME/.teagle_install.lock/pid" 2>/dev/null)" 2>/dev/null; then echo "installing=1"; else echo "installing=0"; fi
echo "disk_free_gb=$(df -BG --output=avail "$HOME" 2>/dev/null | tail -1 | tr -dc '0-9')"
'''


def components_status() -> dict:
    """Per-component state for the install dialog: WSL2, micromamba, RepeatMasker, minimap2,
    Dfam root/curated, FamDB config. One WSL round-trip; each component idempotently repairable."""
    av = available()
    comp = {c[0]: {"key": c[0], "name": c[1], "repairable": c[2], "desc": c[3], "ok": False, "detail": "—"}
            for c in _COMP_META}
    comp["wsl2"]["ok"] = bool(av["wsl2"])
    comp["wsl2"]["detail"] = av.get("distro") or av.get("error") or "not installed"
    if not av["wsl2"]:
        comp["wsl2"]["guide"] = "Open PowerShell as Administrator, run: wsl --install  (then restart Windows)"
        return {"wsl2": False, "installing": False, "ready": False,
                "components": [comp[c[0]] for c in _COMP_META]}
    try:
        rc, out, err = _wsl_script(_STATUS_PROBE, timeout=90)
    except Exception as e:
        comp["wsl2"]["detail"] = f"probe failed: {str(e)[:80]}"
        return {"wsl2": True, "installing": False, "ready": False,
                "components": [comp[c[0]] for c in _COMP_META], "error": str(e)[:120]}
    kv = {}
    for line in out.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()

    def present(key, ok, detail):
        comp[key]["ok"] = ok
        comp[key]["detail"] = detail
    present("micromamba", kv.get("micromamba") == "1", "installed" if kv.get("micromamba") == "1" else "missing")
    rm = kv.get("repeatmasker", "0")
    present("repeatmasker", rm not in ("0", ""), (f"v{rm}" if rm not in ("0", "") else "missing"))
    mm = kv.get("minimap2", "0")
    present("minimap2", mm not in ("0", ""), (mm if mm not in ("0", "") else "missing"))
    present("dfam_root", kv.get("dfam_root") == "1", "present" if kv.get("dfam_root") == "1" else "missing")
    present("dfam_curated", kv.get("dfam_curated") == "1", "present" if kv.get("dfam_curated") == "1" else "missing")
    present("famdb_conf", kv.get("famdb_conf") == "1", "configured" if kv.get("famdb_conf") == "1" else "missing")
    ready = all(comp[k]["ok"] for k in ("repeatmasker", "dfam_root", "dfam_curated"))
    return {"wsl2": True, "installing": kv.get("installing") == "1", "ready": ready,
            "disk_free_gb": kv.get("disk_free_gb"), "components": [comp[c[0]] for c in _COMP_META]}


_INTEGRITY_PROBE = r'''cd "$HOME" 2>/dev/null
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"; FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
echo "=RM="; "$MM" run -n te RepeatMasker -v 2>&1 | head -1
echo "=MM="; "$ENV/bin/minimap2" --version 2>&1 | head -1
echo "=FAMDB="; "$MM" run -n te famdb.py info 2>&1 | grep -iE "version|families|consensus" | head -3
echo "=FILES="; ls -l "$FAMDIR"/dfam40.0.h5 "$FAMDIR"/dfam40.curated.consensus.0.h5 2>&1 | awk '{print $5, $NF}'
'''


def integrity_check() -> dict:
    """Deep functional verification: does each installed tool actually run and does FamDB load?
    Complements components_status (existence) with a runs-clean test."""
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available", "checks": []}
    try:
        rc, out, err = _wsl_script(_INTEGRITY_PROBE, timeout=180)
    except Exception as e:
        return {"ok": False, "error": str(e)[:160], "checks": []}
    sec = {}
    cur = None
    for line in out.splitlines():
        m = re.match(r"=(\w+)=", line.strip())
        if m:
            cur = m.group(1); sec[cur] = []
        elif cur:
            sec[cur].append(line)
    rm_txt = " ".join(sec.get("RM", []))
    mm_txt = " ".join(sec.get("MM", [])).strip()
    fam_txt = " ".join(sec.get("FAMDB", []))
    files = [l for l in sec.get("FILES", []) if l.strip()]
    checks = [
        {"name": "RepeatMasker runs", "ok": "RepeatMasker version" in rm_txt, "detail": rm_txt.strip()[:80] or "no version reported"},
        {"name": "minimap2 runs", "ok": bool(re.match(r"^\d+\.\d+", mm_txt)), "detail": mm_txt[:60] or "no version reported"},
        {"name": "FamDB loads", "ok": bool(re.search(r"(?i)version|families|consensus", fam_txt)), "detail": fam_txt.strip()[:90] or "famdb.py info returned nothing"},
        {"name": "Dfam library files present", "ok": len(files) >= 2 and all("No such" not in f for f in files), "detail": "; ".join(files)[:90] or "missing"},
    ]
    return {"ok": all(c["ok"] for c in checks), "checks": checks, "raw": out[-1500:]}


def install_log(tail: int = 40) -> str:
    try:
        _, out, _ = _wsl(f'tail -n {int(tail)} "$HOME/teagle_wsl_install.log" 2>/dev/null || echo "(no install log yet)"', timeout=30)
        return out
    except Exception as e:
        return f"(log unavailable: {e})"


def annotate(fasta_text: str, species: str | None = None, threads: int = 4, timeout: int = 600) -> dict:
    """Run RepeatMasker (WSL) on the sequence and return family-level hits (Layer A)."""
    sp = ""
    if species:
        if not _SPECIES_RE.match(species):           # validate untrusted input first (hermetic, fast)
            return {"ok": False, "error": "invalid species token"}
        sp = f'-species "{species}"'
    st = env_status()
    if not st["ready"]:
        return {"ok": False, "error": "WSL annotation backend not ready "
                f"(RepeatMasker={st['repeatmasker']}, Dfam={st['dfam']})", "status": st}
    rid = "teagle_" + secrets.token_hex(6)
    try:
        rc, out, err = _wsl(f'mkdir -p /tmp/{rid} && cat > /tmp/{rid}/q.fa',
                            stdin=fasta_text.encode(), timeout=60)
        if rc != 0:
            return {"ok": False, "error": "failed to stage sequence: " + err.strip()[:200]}
        script = (f'cd /tmp/{rid} && {_MM} run -n te RepeatMasker -pa {int(threads)} {sp} -qq q.fa '
                  f'>rm.log 2>&1; echo "EXIT $?"; cat q.fa.out 2>/dev/null')
        rc, out, err = _wsl(script, timeout=timeout)
        hits = parse_out(out)
        # fetch tool versions for provenance
        rcv, ver, _ = _wsl(f'{_MM} run -n te RepeatMasker -v 2>/dev/null | head -1', timeout=30)
        _wsl(f'rm -rf /tmp/{rid}', timeout=30)
        return {"ok": True, "hits": hits, "n_hits": len(hits),
                "repeatmasker_version": st["repeatmasker"], "raw_out": out[-4000:],
                "species": species or "(all installed families)"}
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"RepeatMasker timed out after {timeout}s"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


_SPLICE_CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}    # U2 / minor U12 canonical splice sites
_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _parse_sam_splice(sam: str):
    """Primary alignment from SAM -> exon blocks and intron (CIGAR 'N') skips, in 0-based ref coords."""
    for line in sam.splitlines():
        if line.startswith("@") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        flag = int(f[1])
        if flag & 4 or flag & 0x100 or flag & 0x800:          # unmapped / secondary / supplementary
            continue
        pos = int(f[3]) - 1
        strand = "-" if (flag & 16) else "+"
        exons, introns, ref, cur = [], [], pos, pos
        for n, op in _CIGAR_RE.findall(f[5]):
            n = int(n)
            if op in "M=XD":
                ref += n
            elif op == "N":                                   # intron gap in the reference
                exons.append({"start": cur, "end": ref})
                introns.append({"start": ref, "end": ref + n})
                ref += n
                cur = ref
        exons.append({"start": cur, "end": ref})
        return {"exons": exons, "introns": introns, "strand": strand, "ref_start": pos, "ref_end": ref}
    return None


def splice_align(genomic_fasta: str, transcript_fasta: str, timeout: int = 180) -> dict:
    """Splice-aware alignment of a transcript/cDNA to the genomic sequence (minimap2 -x splice).
    Exons = aligned blocks; introns = CIGAR 'N' skips; splice sites compared to canonical motifs.
    Sequences are staged as files/STDIN (data), never interpolated into the command."""
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rid = "teagle_" + secrets.token_hex(6)
    try:
        rc, out, err = _wsl(f'[ -x "{_ENV}/bin/minimap2" ] && echo yes || echo no', timeout=20)
        if "yes" not in out:
            return {"ok": False, "error": "minimap2 not installed in the WSL backend (re-run Install backend)"}
        rc, out, err = _wsl(f'mkdir -p /tmp/{rid} && cat > /tmp/{rid}/ref.fa', stdin=genomic_fasta.encode(), timeout=60)
        if rc != 0:
            return {"ok": False, "error": "failed to stage genomic sequence: " + err.strip()[:200]}
        script = f'cd /tmp/{rid} && {_MM} run -n te minimap2 -a -x splice --secondary=no ref.fa - 2>/dev/null'
        rc, sam, err = _wsl(script, stdin=transcript_fasta.encode(), timeout=timeout)
        _rcv, ver, _ = _wsl(f'{_ENV}/bin/minimap2 --version 2>/dev/null', timeout=30)
        _wsl(f'rm -rf /tmp/{rid}', timeout=30)
        res = _parse_sam_splice(sam)
        if not res:
            return {"ok": False, "error": "transcript did not align to the genomic sequence (check they correspond)"}
        g = "".join(l.strip() for l in genomic_fasta.splitlines() if not l.startswith(">")).upper()
        minus = res["strand"] == "-"
        for it in res["introns"]:                             # transcribed-strand splice-site motifs (donor..acceptor)
            if minus:                                         # read the motif from the reverse-complement strand
                donor, acceptor = _revcomp(g[it["end"] - 2:it["end"]]), _revcomp(g[it["start"]:it["start"] + 2])
            else:
                donor, acceptor = g[it["start"]:it["start"] + 2], g[it["end"] - 2:it["end"]]
            it["donor"], it["acceptor"] = donor, acceptor
            it["canonical"] = (donor, acceptor) in _SPLICE_CANON
        res.update(ok=True, minimap2_version=ver.strip(),
                   counts={"exons": len(res["exons"]), "introns": len(res["introns"])},
                   canonical_introns=sum(1 for i in res["introns"] if i.get("canonical")))
        return res
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"minimap2 timed out after {timeout}s"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
