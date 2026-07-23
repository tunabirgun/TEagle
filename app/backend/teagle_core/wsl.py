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
_GENOMES = "$HOME/teagle_genomes"                       # per-assembly downloaded genome cache (kept after first prepare)
_SPECIES_RE = re.compile(r"^[A-Za-z0-9 _.-]{1,60}$")
_ACC_RE = re.compile(r"^GC[AF]_\d+\.\d+$")              # RefSeq/GenBank assembly accession (the reproducibility pin)
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
# A truncated/corrupt repodata shard (interrupted or concurrent fetch) makes every solve die with
# "Could not load repodata.json ... after retry" and stays stuck across re-runs. On the first failure
# purge the index cache + pkgs cache (mamba's own advice: `clean -a`) and retry once — never touch the
# env prefix, the multi-GB Dfam .h5 libraries live inside it.
mm_reset_cache(){
  echo "[teagle] purging ALL conda caches (corrupted/incompatible repodata recovery), then retrying"
  "$MM" clean --all --yes >/dev/null 2>&1 || "$MM" clean --index-cache --yes >/dev/null 2>&1 || true
  # --index-cache does NOT clear the newer SHARDED repodata cache; remove every known cache dir by hand.
  rm -rf "$MAMBA_ROOT_PREFIX/pkgs/cache" "$HOME/.cache/mamba" "$HOME/.cache/rattler" \
         "$HOME/.cache/conda" 2>/dev/null || true
}
mm_create(){   # create the shared 'te' env if absent; clean-index + retry once on a solve failure
  [ -d "$ENV/conda-meta" ] && return 0
  "$MM" create -y -n te -c conda-forge -c bioconda && return 0
  mm_reset_cache; "$MM" create -y -n te -c conda-forge -c bioconda
}
mm_install(){  # install package(s) into 'te'; clean-index + force-reinstall retry once on a solve failure
  "$MM" install -y -n te -c conda-forge -c bioconda "$@" && return 0
  mm_reset_cache; "$MM" install --force-reinstall -y -n te -c conda-forge -c bioconda "$@"
}
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
    # A fresh minimal Ubuntu WSL ships neither curl nor bzip2, so the old single `curl | tar -xj`
    # failed on other users' PCs (every downstream step then died "micromamba required first").
    # Robust order: (1) reuse a micromamba already on the box — another app (e.g. BulkSeq Studio)
    # may have installed one at ~/.local/bin; (2) python3 stdlib download (ships on default Ubuntu,
    # bz2 is built into CPython — no curl/bzip2/apt/sudo); (3) curl/wget+bzip2; (4) passwordless apt.
    # The reused/installed binary is copied to $MM ($HOME/bin/micromamba) so every hardcoded ref works.
    "micromamba": r'''
echo "[teagle] STEP micromamba START"
mkdir -p "$HOME/bin" || fail "cannot create $HOME/bin"
MM_URL="https://micro.mamba.pm/api/micromamba/linux-64/latest"
mm_py3(){ python3 - "$MM" "$MM_URL" <<'PY'
import io, os, stat, sys, tarfile, urllib.request
dest, url = sys.argv[1], sys.argv[2]
data = urllib.request.urlopen(url, timeout=180).read()
with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as tf:
    m = tf.extractfile(tf.getmember("bin/micromamba"))
    if m is None: raise SystemExit("bin/micromamba not in archive")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    open(dest, "wb").write(m.read())
os.chmod(dest, os.stat(dest).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
PY
}
if [ -x "$MM" ]; then echo "[teagle] micromamba already present"; else
  FOUND=""
  for c in "$HOME/.local/bin/micromamba" "$HOME/micromamba/bin/micromamba" "$(command -v micromamba 2>/dev/null)"; do
    if [ -n "$c" ] && [ "$c" != "$MM" ] && [ -x "$c" ] && "$c" --version >/dev/null 2>&1; then FOUND="$c"; break; fi
  done
  if [ -n "$FOUND" ] && cp -f "$FOUND" "$MM" && chmod +x "$MM"; then
    echo "[teagle] reused existing micromamba from $FOUND"
  elif command -v python3 >/dev/null 2>&1 && mm_py3; then
    echo "[teagle] micromamba installed via python3"
  elif command -v bzip2 >/dev/null 2>&1 && command -v curl >/dev/null 2>&1 && curl -fL "$MM_URL" | tar -xj -C "$HOME" bin/micromamba; then
    echo "[teagle] micromamba installed via curl"
  elif command -v bzip2 >/dev/null 2>&1 && command -v wget >/dev/null 2>&1 && wget -qO- "$MM_URL" | tar -xj -C "$HOME" bin/micromamba; then
    echo "[teagle] micromamba installed via wget"
  elif sudo -n true 2>/dev/null && sudo apt-get update && sudo apt-get install -y python3 ca-certificates && mm_py3; then
    echo "[teagle] micromamba installed via apt+python3"
  else
    fail "micromamba: no python3, no curl/wget+bzip2, and sudo needs a password. Open a WSL terminal, run:  sudo apt-get update && sudo apt-get install -y python3 ca-certificates  then click Repair again"
  fi
  [ -x "$MM" ] || fail "micromamba missing after install"
fi
echo "[teagle] STEP micromamba OK"
''',
    "repeatmasker": r'''
echo "[teagle] STEP repeatmasker START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if "$MM" run -n te RepeatMasker -v >/dev/null 2>&1; then echo "[teagle] RepeatMasker already present"; else
  # never `env remove` here: the Dfam .h5 libraries (multi-GB) live INSIDE this env prefix.
  # mm_create / mm_install self-recover from a corrupted repodata shard (clean-index + retry once).
  mm_create || fail "create te env"
  mm_install repeatmasker || fail "install repeatmasker"
  "$MM" run -n te RepeatMasker -v >/dev/null 2>&1 || fail "RepeatMasker not runnable after install"
fi
echo "[teagle] STEP repeatmasker OK"
''',
    "minimap2": r'''
echo "[teagle] STEP minimap2 START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/minimap2" ]; then echo "[teagle] minimap2 already present"; else
  mm_create || fail "create te env"
  mm_install minimap2 || fail "install minimap2"
  [ -x "$ENV/bin/minimap2" ] || fail "minimap2 missing after install"
fi
echo "[teagle] STEP minimap2 OK"
''',
    "miniprot": r'''
echo "[teagle] STEP miniprot START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/miniprot" ]; then echo "[teagle] miniprot already present"; else
  mm_create || fail "create te env"
  mm_install miniprot || fail "install miniprot"
  [ -x "$ENV/bin/miniprot" ] || fail "miniprot missing after install"
fi
echo "[teagle] STEP miniprot OK"
''',
    "genomescan": r'''
echo "[teagle] STEP genomescan START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/isPcr" ] && [ -x "$ENV/bin/datasets" ]; then echo "[teagle] isPcr + datasets already present"; else
  mm_create || fail "create te env"
  mm_install ispcr ncbi-datasets-cli || fail "install ispcr + ncbi-datasets-cli"
  [ -x "$ENV/bin/isPcr" ] || fail "isPcr missing after install"
  [ -x "$ENV/bin/datasets" ] || fail "datasets missing after install"
fi
# faToTwoBit is best-effort: it makes cached genomes compact + fast to load, but scans still work on plain
# FASTA if it is unavailable, so a missing package must NOT fail the step.
[ -x "$ENV/bin/faToTwoBit" ] || "$MM" install -y -n te -c bioconda -c conda-forge ucsc-fatotwobit >/dev/null 2>&1 \
  || echo "[teagle] note: faToTwoBit unavailable (genomes cached as FASTA — larger, still functional)"
echo "[teagle] STEP genomescan OK"
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

# miniprot (homology tier) is intentionally NOT in the default install list while that tier is on hold;
# its step + parser stay in the code, dormant, ready to re-enable when the homology UI ships.
_ALL_STEPS = ["micromamba", "repeatmasker", "minimap2", "genomescan", "dfam_root", "dfam_curated", "famdb_conf"]

# component metadata surfaced to the install dialog (order = install order)
_COMP_META = [
    ("wsl2",         "WSL2 + Linux distro",       False, "Windows Subsystem for Linux — hosts the Dfam / RepeatMasker stack."),
    ("micromamba",   "micromamba (conda)",        True,  "Small conda package manager, installed under your Linux home."),
    ("repeatmasker", "RepeatMasker",              True,  "Homology-based TE annotator that names Dfam families."),
    ("minimap2",     "minimap2",                  True,  "Splice-aware aligner for de-novo exon / intron detection."),
    ("genomescan",   "isPcr + NCBI Datasets",     True,  "Local whole-genome in-silico PCR engine + genome downloader."),
    ("dfam_root",    "Dfam 4.0 root library",     True,  "Dfam root partition (dfam40.0.h5)."),
    ("dfam_curated", "Dfam 4.0 curated library",  True,  "Dfam curated consensus partition."),
    ("famdb_conf",   "FamDB configuration",       True,  "Points RepeatMasker at the downloaded Dfam library."),
]


def _build_script(keys) -> str:
    script = _PRELUDE + "".join(_STEP[k] for k in keys) + '\necho "[teagle] DONE $(date -u +%FT%TZ)"\n'
    # force LF: bash chokes on CRLF (a heredoc terminator "PY\r" never matches "PY"), and a
    # Windows checkout with core.autocrlf=true could otherwise deliver a \r-poisoned script.
    return script.replace("\r\n", "\n").replace("\r", "\n")


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
echo "genomescan=$([ -x "$ENV/bin/isPcr" ] && [ -x "$ENV/bin/datasets" ] && echo 1 || echo 0)"
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
        # distinguish absent (no distro) from registered-but-won't-start (broken ext4.vhdx) — different actions
        broken = bool(av.get("distro"))
        win = os.name == "nt"
        if broken:
            comp["wsl2"]["detail"] = f"'{av['distro']}' registered but won't start"
            comp["wsl2"]["guide"] = (f"In an Administrator PowerShell run:  wsl --unregister {av['distro']}  "
                                     "then click Install WSL to reinstall (a restart may be required).")
        else:
            comp["wsl2"]["guide"] = ("Click Install WSL to install WSL2 + Ubuntu (needs Administrator; "
                                     "a Windows restart may be required before first use).")
        comp["wsl2"]["installable"] = win           # the dialog shows an in-app Install WSL button on Windows
        return {"wsl2": False, "installing": _wsl2_installing(),
                "ready": False, "components": [comp[c[0]] for c in _COMP_META]}
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
    present("genomescan", kv.get("genomescan") == "1", "installed" if kv.get("genomescan") == "1" else "missing")
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
echo "=SCAN="; { [ -x "$ENV/bin/isPcr" ] && echo "isPcr present"; } || echo "isPcr MISSING"; { [ -x "$ENV/bin/datasets" ] && echo "datasets present"; } || echo "datasets MISSING"
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
    scan_txt = " ".join(sec.get("SCAN", []))
    checks = [
        {"name": "RepeatMasker runs", "ok": "RepeatMasker version" in rm_txt, "detail": rm_txt.strip()[:80] or "no version reported"},
        {"name": "minimap2 runs", "ok": bool(re.match(r"^\d+\.\d+", mm_txt)), "detail": mm_txt[:60] or "no version reported"},
        {"name": "FamDB loads", "ok": bool(re.search(r"(?i)version|families|consensus", fam_txt)), "detail": fam_txt.strip()[:90] or "famdb.py info returned nothing"},
        {"name": "Dfam library files present", "ok": len(files) >= 2 and all("No such" not in f for f in files), "detail": "; ".join(files)[:90] or "missing"},
        # the whole-genome off-target scan needs isPcr + NCBI datasets; verify them here so a "healthy" deep-check
        # never precedes a scan that fails at first use (same binding-truthfulness class as the fixed unzip gap)
        {"name": "isPcr + NCBI datasets (whole-genome scan)", "ok": "isPcr present" in scan_txt and "datasets present" in scan_txt,
         "detail": scan_txt.strip()[:80] or "not probed"},
    ]
    return {"ok": all(c["ok"] for c in checks), "checks": checks, "raw": out[-1500:]}


def install_log(tail: int = 40) -> str:
    try:
        _, out, _ = _wsl(f'tail -n {int(tail)} "$HOME/teagle_wsl_install.log" 2>/dev/null || echo "(no install log yet)"', timeout=30)
        return out
    except Exception as e:
        return f"(log unavailable: {e})"


# ---------- WSL2 itself (Windows-side, elevated) — install the distro when WSL is absent ----------
_WIN_WSL_LOG = "wsl_win_install.log"
_wsl2_thread = None


def wsl2_install_log(tail: int = 200) -> str:
    """Tail the Windows-side WSL-install log. The in-WSL log can't exist until WSL is up, so the
    elevated installer logs here; wsl.exe emits UTF-16LE, so NULs are stripped for readability."""
    p = os.path.join(appdirs.user_data_dir(), _WIN_WSL_LOG)
    try:
        with open(p, "rb") as f:
            txt = f.read().decode("utf-8", "ignore").replace("\x00", "")
        return "\n".join(txt.splitlines()[-int(tail):])
    except Exception:
        return ""


def _wsl2_installing() -> bool:
    """A WSL2 install is in progress = the Windows-side log exists but has no terminal marker yet."""
    log = wsl2_install_log(500)
    return bool(log) and "DONE-WSL" not in log and "[teagle] FAILED" not in log


def _wsl2_bat_script(log: str) -> str:
    """The elevated batch: install WSL2 + Ubuntu (no interactive launch), tolerant of pre-reboot state,
    logging everything (incl. a terminal DONE-WSL marker) to the Windows-side `log` the dialog polls."""
    return (
        "@echo off\r\n"
        f'> "{log}" echo [teagle] installing WSL2 + Ubuntu (elevated)\r\n'
        f'wsl.exe --install -d Ubuntu --no-launch >> "{log}" 2>&1\r\n'
        "if %ERRORLEVEL% NEQ 0 (\r\n"
        f'  >> "{log}" echo [teagle] first attempt returned %ERRORLEVEL%; retrying: wsl --install\r\n'
        f'  wsl.exe --install >> "{log}" 2>&1\r\n'
        ")\r\n"
        f'wsl.exe --set-default-version 2 >> "{log}" 2>&1 || rem VirtualMachinePlatform is inactive until reboot\r\n'
        f'>> "{log}" echo [teagle] wsl --status:\r\n'
        f'wsl.exe --status >> "{log}" 2>&1\r\n'
        f'>> "{log}" echo [teagle] DONE-WSL %ERRORLEVEL% (restart Windows if the distro is not usable yet, then reopen this installer)\r\n'
    )


def install_wsl2() -> dict:
    """Install WSL2 + Ubuntu from an ELEVATED helper (`wsl --install` requires Administrator).
    Fire-and-forget: an elevated .bat runs `wsl --install -d Ubuntu --no-launch` (falling back to a
    plain `wsl --install` on older wsl.exe) and logs to a Windows-side file the dialog polls. If the
    UAC prompt is declined or the user is not an admin, a terminal FAILED marker is written so the
    poller never hangs. A Windows restart may be required before the new distro is usable."""
    global _wsl2_thread
    if os.name != "nt":
        return {"started": False, "error": "WSL installation is only available on Windows"}
    if _wsl2_thread is not None and _wsl2_thread.is_alive():
        return {"started": False, "error": "WSL install already running"}
    d = appdirs.user_data_dir()
    bat = os.path.join(d, "install_wsl.bat")
    log = os.path.join(d, _WIN_WSL_LOG)
    try:
        with open(bat, "w", encoding="ascii", errors="ignore", newline="") as f:
            f.write(_wsl2_bat_script(log))
        with open(log, "w", encoding="utf-8", newline="") as f:
            f.write("[teagle] launching the elevated WSL installer - accept the Windows (UAC) prompt...\n")
    except Exception as e:
        return {"started": False, "error": f"could not stage the WSL installer: {type(e).__name__}: {e}"}

    def _run():
        try:
            import ctypes
            # ShellExecuteW verb 'runas' -> UAC. Return > 32 = launched; <= 32 = failed (e.g. 1223 = declined).
            rc = int(ctypes.windll.shell32.ShellExecuteW(None, "runas", bat, None, None, 0))
            if rc <= 32:
                with open(log, "a", encoding="utf-8") as f:
                    f.write(f"\n[teagle] FAILED - could not elevate (code {rc}); the UAC prompt was declined or you "
                            "are not an administrator.\n[teagle] Manual: open PowerShell as Administrator and run:  wsl --install\n")
        except Exception as e:
            try:
                with open(log, "a", encoding="utf-8") as f:
                    f.write(f"\n[teagle] FAILED - {type(e).__name__}: {e}\n")
            except Exception:
                pass

    _wsl2_thread = threading.Thread(target=_run, daemon=True)
    _wsl2_thread.start()
    return {"started": True, "windows_log": True}


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


# ---------- homology-based coding/intron recovery (WSL / miniprot) ----------
def _mp_attrs(field9: str) -> dict:
    d = {}
    for kv in field9.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def _parse_miniprot_gff(gff_text: str, genomic_seq: str, max_hits: int = 25):
    """miniprot --gff -> ranked gene models. Pure function of the GFF3 text + genomic sequence.
    Each mRNA is one hit; its CDS features are exons; introns are the gaps between genomically
    adjacent CDS (miniprot keeps frameshifts in a single CDS, so a CDS break is a real intron).
    Splice motifs are read strand-aware and flagged canonical, as in splice_align. Frameshift
    and in-frame-stop counts come from the mRNA attributes / the preceding ##PAF tag line."""
    g = (genomic_seq or "").upper()
    hits, order, paf = {}, [], None
    for line in gff_text.splitlines():
        if line.startswith("##PAF"):
            paf = line.split("\t")
            continue
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        _seqid, _src, ftype, start, end, score, strand, phase, attrs = f[:9]
        a = _mp_attrs(attrs)
        if ftype == "mRNA":
            hid = a.get("ID") or f"MP{len(order) + 1:06d}"
            tgt = a.get("Target", "").split()
            qlen = fs = st = None
            if paf and len(paf) >= 3:
                try:
                    qlen = int(paf[2])                       # ##PAF col2 = query (protein) length
                except ValueError:
                    qlen = None
                tags = "\t".join(paf)
                mfs = re.search(r"fs:i:(\d+)", tags); fs = int(mfs.group(1)) if mfs else None
                mst = re.search(r"st:i:(\d+)", tags); st = int(mst.group(1)) if mst else None
            try:
                sc = float(score)
            except ValueError:
                sc = None
            hits[hid] = {
                "protein": tgt[0] if tgt else a.get("Target", ""),
                "protein_start": int(tgt[1]) if len(tgt) >= 3 and tgt[1].isdigit() else None,
                "protein_end": int(tgt[2]) if len(tgt) >= 3 and tgt[2].isdigit() else None,
                "protein_len": qlen,
                "strand": strand, "score": sc,
                "identity": float(a["Identity"]) if a.get("Identity") else None,
                "positive": float(a["Positive"]) if a.get("Positive") else None,
                "rank": int(a["Rank"]) if a.get("Rank", "").isdigit() else None,
                "frameshifts": int(a["Frameshift"]) if a.get("Frameshift", "").isdigit() else (fs or 0),
                "inframe_stops": int(a["StopCodon"]) if a.get("StopCodon", "").isdigit() else (st if st is not None else 0),
                "ref_start": int(start) - 1, "ref_end": int(end),
                "cds": [],
            }
            order.append(hid)
            paf = None
        elif ftype == "CDS":
            hid = a.get("Parent")
            if hid in hits:
                hits[hid]["cds"].append({"start": int(start) - 1, "end": int(end)})
    out = []
    for hid in order:
        h = hits[hid]
        cds = sorted(h["cds"], key=lambda c: c["start"])
        minus = h["strand"] == "-"
        exons = [{"start": c["start"], "end": c["end"]} for c in cds]
        introns = []
        for prev, nxt in zip(cds, cds[1:]):
            istart, iend = prev["end"], nxt["start"]         # 0-based half-open gap between exons
            if iend <= istart:
                continue
            if minus:                                        # motif read on the transcribed (reverse) strand
                donor, acceptor = _revcomp(g[iend - 2:iend]), _revcomp(g[istart:istart + 2])
            else:
                donor, acceptor = g[istart:istart + 2], g[iend - 2:iend]
            introns.append({"start": istart, "end": iend, "length": iend - istart,
                            "donor": donor, "acceptor": acceptor,
                            "canonical": (donor, acceptor) in _SPLICE_CANON})
        cov = None
        if h["protein_len"] and h["protein_start"] is not None and h["protein_end"] is not None:
            cov = round((h["protein_end"] - h["protein_start"] + 1) / h["protein_len"], 4)
        out.append({
            "protein": h["protein"], "strand": h["strand"], "score": h["score"],
            "identity": h["identity"], "positive": h["positive"], "rank": h["rank"],
            "frameshifts": h["frameshifts"], "inframe_stops": h["inframe_stops"],
            "ref_start": h["ref_start"], "ref_end": h["ref_end"],
            "protein_start": h["protein_start"], "protein_end": h["protein_end"],
            "protein_len": h["protein_len"], "protein_coverage": cov,
            "exons": list(reversed(exons)) if minus else exons,   # display in protein/transcription order
            "introns": introns,
            "counts": {"exons": len(exons), "introns": len(introns),
                       "canonical_introns": sum(1 for i in introns if i["canonical"])},
        })
    out.sort(key=lambda x: (-(x["score"] or 0), -(x["identity"] or 0)))
    return out[:max_hits]


def protein_align(genomic_fasta: str, protein_fasta: str, timeout: int = 180, max_hits: int = 25) -> dict:
    """Spliced protein-to-genome alignment (miniprot --gff): recovers CDS/exon boundaries, introns,
    and frameshift/stop lesions from a bare genomic sequence WITHOUT a transcript. The reference
    protein(s) are external evidence. Sequences are staged as files (data), never interpolated."""
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rid = "teagle_" + secrets.token_hex(6)
    try:
        rc, out, err = _wsl(f'[ -x "{_ENV}/bin/miniprot" ] && echo yes || echo no', timeout=20)
        if "yes" not in out:
            return {"ok": False, "error": "miniprot not installed in the WSL backend (re-run Install backend)"}
        rc, out, err = _wsl(f'mkdir -p /tmp/{rid} && cat > /tmp/{rid}/ref.fa', stdin=genomic_fasta.encode(), timeout=60)
        if rc != 0:
            return {"ok": False, "error": "failed to stage genomic sequence: " + err.strip()[:200]}
        rc, out, err = _wsl(f'cat > /tmp/{rid}/prot.faa', stdin=protein_fasta.encode(), timeout=60)
        if rc != 0:
            return {"ok": False, "error": "failed to stage reference protein(s): " + err.strip()[:200]}
        script = f'cd /tmp/{rid} && {_MM} run -n te miniprot --gff ref.fa prot.faa 2>/dev/null'
        rc, gff, err = _wsl(script, timeout=timeout)
        _rcv, ver, _ = _wsl(f'{_ENV}/bin/miniprot --version 2>/dev/null', timeout=30)
        _wsl(f'rm -rf /tmp/{rid}', timeout=30)
        g = "".join(l.strip() for l in genomic_fasta.splitlines() if not l.startswith(">")).upper()
        hits = _parse_miniprot_gff(gff, g, max_hits=max_hits)
        if not hits:
            return {"ok": False, "error": "no reference protein aligned to the sequence "
                    "(too diverged, or the query is non-coding for these proteins)"}
        return {"ok": True, "miniprot_version": ver.strip(), "hits": hits,
                "counts": {"hits": len(hits),
                           "with_introns": sum(1 for h in hits if h["counts"]["introns"]),
                           "with_frameshift": sum(1 for h in hits if h["frameshifts"]),
                           "with_inframe_stop": sum(1 for h in hits if h["inframe_stops"])}}
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"miniprot timed out after {timeout}s"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


# ---------- local whole-genome in-silico PCR (WSL / isPcr against a downloaded RefSeq assembly) ----------
def _parse_meta(text: str) -> dict:
    """Parse the key=value meta.txt (accession/target/sha256/n_seqs/bytes) written by genome_prepare."""
    d = {}
    for line in text.splitlines():
        line = line.strip()
        if "=" in line and not line.startswith(("[", "FAIL", "PREPARED", "NOTPREP")):
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def _ispcr_ver(banner: str) -> str:
    m = re.search(r"v\s*([\w.]+)", banner or "")
    return m.group(1) if m else "unknown"


# download the assembly, convert to compact 2bit (best-effort), seal the SOURCE FASTA sha256 (machine-
# independent, unlike the version-dependent 2bit), mark .done for idempotent resume. __ACC__ is a
# validated accession (safe to interpolate); no untrusted input enters the script.
_GENOME_PREP_LOG = "$HOME/teagle_genome_prepare.log"    # milestones the UI tails for a liveness indicator

_PREP_SCRIPT = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || exit 1
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
ENV="$HOME/micromamba/envs/te"
ACC="__ACC__"
GDIR="$HOME/teagle_genomes/$ACC"
LOG="$HOME/teagle_genome_prepare.log"
plog(){ echo "[prep] $1"; echo "$1" >> "$LOG"; }
mkdir -p "$GDIR" || { echo "FAIL mkdir"; exit 1; }
if [ -f "$GDIR/.done" ] && [ -f "$GDIR/meta.txt" ]; then echo "PREPARED"; cat "$GDIR/meta.txt"; exit 0; fi
# atomic per-genome lock so two concurrent prepares of the same accession can't race on dl.zip / genome.2bit
# / meta.txt and leave a half-written cache; a lock orphaned by a crash (dead PID) is reaped.
LOCK="$GDIR/.lock"
if [ -d "$LOCK" ] && { [ ! -f "$LOCK/pid" ] || ! kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; }; then rm -rf "$LOCK" 2>/dev/null; fi
if ! mkdir "$LOCK" 2>/dev/null; then echo "FAIL a download for this genome is already running"; exit 1; fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK" 2>/dev/null' EXIT
[ -f "$GDIR/.done" ] && { echo "PREPARED"; cat "$GDIR/meta.txt"; exit 0; }   # another prepare finished while we waited
cd "$GDIR"
avail=$(df -BG --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9'); avail=${avail:-0}
# blanket floor sized for a mammalian genome's peak (extracted FASTA ~3G + 2bit ~0.8G); the zip is freed
# before the conversion to keep the peak below this.
[ "$avail" -ge 8 ] || { echo "FAIL insufficient disk (${avail}G free, need >=8G for a genome)"; exit 1; }
plog "downloading $ACC (NCBI Datasets) — can take several minutes for a large genome"
# a ~1 GB mammalian download can drop mid-transfer; datasets writes a fresh zip each time, so retry a few
# times (the partial is discarded) rather than failing the whole prepare on one transient network blip.
ok=0
for attempt in 1 2 3; do
  if "$ENV/bin/datasets" download genome accession "$ACC" --include genome --filename dl.zip >/dev/null 2>dl.err; then ok=1; break; fi
  plog "download attempt $attempt failed — retrying"; rm -f dl.zip; sleep 5
done
# surface the last attempt's stderr (invalid/withdrawn accession, DNS/API/rate-limit) instead of an opaque failure
[ "$ok" = 1 ] || { echo "FAIL download (after 3 attempts): $(tail -c 200 dl.err 2>/dev/null | tr '\n' ' ')"; exit 1; }
plog "extracting genome FASTA"
rm -rf ex; mkdir -p ex
# a fresh minimal Ubuntu WSL ships python3 (zipfile is built into CPython) but NOT unzip — extract with the
# python3 stdlib first (guaranteed present), and fall back to a system unzip only if one happens to exist. The
# te conda env is never activated in this login shell, so a conda unzip at $ENV/bin is not on PATH; the fallback
# deliberately probes the system PATH. python3 -m zipfile -e preserves the internal ncbi_dataset/.../*.fna path.
if python3 -m zipfile -e dl.zip ex/ 2>/dev/null; then :
elif command -v unzip >/dev/null 2>&1 && unzip -oq dl.zip -d ex; then :
else echo "FAIL extract (need python3 with zipfile, or unzip)"; exit 1; fi
rm -f dl.zip                                            # free the zip before the (large) 2bit conversion -> lower peak disk
FNA=$(find ex -name '*_genomic.fna' | head -1)
[ -n "$FNA" ] || { echo "FAIL no genomic fna in package"; exit 1; }
plog "checksumming source FASTA"
FSHA=$(sha256sum "$FNA" | cut -d' ' -f1)
[ -n "$FSHA" ] || { echo "FAIL checksum failed"; exit 1; }   # never write .done/meta with an empty seal hash
NSEQ=$(grep -c '^>' "$FNA")
[ "${NSEQ:-0}" -ge 1 ] || { echo "FAIL extracted FASTA has no sequences"; exit 1; }   # never seal a 0-contig genome
plog "building compact isPcr target (2bit)"
TARGET=""
if [ -x "$ENV/bin/faToTwoBit" ] && "$ENV/bin/faToTwoBit" "$FNA" genome.2bit 2>/dev/null; then
  TARGET="genome.2bit"
else
  mv "$FNA" genome.fna 2>/dev/null || cp "$FNA" genome.fna; TARGET="genome.fna"
fi
TBYTES=$(stat -c %s "$TARGET" 2>/dev/null || echo 0)
rm -rf ex
printf 'accession=%s\ntarget=%s\nsha256=%s\nn_seqs=%s\nbytes=%s\n' "$ACC" "$TARGET" "$FSHA" "$NSEQ" "$TBYTES" > meta.txt
touch .done
plog "done"
echo "PREPARED"
cat meta.txt
'''


def genome_prepare(accession: str, assembly_name: str = "", timeout: int = 3600) -> dict:
    """One-time: download the RefSeq assembly (NCBI Datasets), build a compact isPcr target, and record
    the source-FASTA sha256 + contig count. Idempotent (a completed prepare returns instantly). Slow for
    large genomes — run it off the UI thread. The cache is kept for later scans."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rc, out, _ = _wsl(f'[ -x "{_ENV}/bin/datasets" ] && [ -x "{_ENV}/bin/isPcr" ] && echo yes || echo no', timeout=20)
    if "yes" not in out:
        return {"ok": False, "error": "genome-scan tools not installed in the WSL backend (run Install backend)"}
    try:
        rc, out, err = _wsl_script(_PREP_SCRIPT.replace("__ACC__", accession), timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"genome preparation timed out after {timeout}s (large genome — try again to resume)"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    meta = _parse_meta(out)
    if "PREPARED" not in out or not meta.get("sha256"):
        fl = next((l for l in out.splitlines() if l.startswith("FAIL")), "") or err.strip()[:200]
        return {"ok": False, "error": "genome preparation failed: " + (fl.replace("FAIL", "").strip() or "unknown error")}
    return {"ok": True, "accession": accession, "assembly_name": assembly_name,
            "target": meta.get("target"), "sha256": meta.get("sha256"),
            "n_seqs": int(meta.get("n_seqs", 0) or 0), "bytes": int(meta.get("bytes", 0) or 0)}


def genome_prepare_log(tail: int = 1) -> str:
    """Tail the genome-prepare milestone log for a UI liveness indicator during a long download."""
    try:
        _, out, _ = _wsl(f'tail -n {int(tail)} "{_GENOME_PREP_LOG}" 2>/dev/null || true', timeout=15)
        return out.strip()
    except Exception:
        return ""


def genome_scan(accession: str, query_text: str, max_size: int = 4000, min_size: int = 0,
                min_perfect: int = 15, min_good: int = 15, timeout: int = 600) -> dict:
    """Run isPcr for a prepared assembly against the query file (name<TAB>fwd<TAB>rev rows, staged as
    STDIN data). Returns the raw isPcr FASTA + the sealed genome sha256 + isPcr version. If the genome
    is not prepared yet, returns need_prepare so the UI can offer to download it first."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rc, out, _ = _wsl(f'[ -f "{_GENOMES}/{accession}/.done" ] && cat "{_GENOMES}/{accession}/meta.txt" || echo NOTPREP',
                      timeout=20)
    if "NOTPREP" in out or "accession=" not in out:
        return {"ok": False, "error": "genome not prepared — download it first", "need_prepare": True}
    meta = _parse_meta(out)
    target = meta.get("target") or "genome.fna"
    rid = "teagle_" + secrets.token_hex(6)
    try:
        rc, _o, err = _wsl(f'mkdir -p /tmp/{rid} && cat > /tmp/{rid}/q.txt', stdin=query_text.encode(), timeout=60)
        if rc != 0:
            return {"ok": False, "error": "failed to stage the primer query: " + err.strip()[:200]}
        # NB: isPcr v33 lists -minSize in its help but the binary REJECTS it (exit 255) — a lower size bound
        # is applied downstream in the parser, not here. -maxSize / -minPerfect / -minGood are honoured.
        opts = f"-maxSize={int(max_size)} -minPerfect={int(min_perfect)} -minGood={int(min_good)}"
        # Verify the cached target exists AND is non-empty AND isPcr exits cleanly. A missing/corrupt genome
        # (or any isPcr failure) MUST surface as an error, never a silent empty result — otherwise "0 off-target
        # sites" would falsely certify a primer pair as specific when the scan never actually ran. An isPcr
        # exit 0 with no products is the one legitimate empty case and stays ok. Delivered via STDIN (not an
        # inline `bash -lc` arg) so the shell variable survives wsl.exe's command-line rebuild.
        run = (f'cd /tmp/{rid} || exit 1\n'
               f'T="{_GENOMES}/{accession}/{target}"\n'
               f'[ -x "{_ENV}/bin/isPcr" ] || exit 8\n'      # isPcr missing/broken -> repair the backend, NOT re-download the genome
               f'[ -s "$T" ] || exit 9\n'
               f'"{_ENV}/bin/isPcr" {opts} "$T" q.txt stdout\n')
        rc, raw, err = _wsl_script(run, timeout=timeout)
        if rc == 8:                                           # distinguish a broken tool from a missing genome (don't misdirect the fix)
            return {"ok": False, "error": "the isPcr tool is missing from the WSL backend — open ⚙ Backend installer and "
                    "repair “isPcr + NCBI Datasets”. The cached genome is fine and does not need re-downloading."}
        if rc != 0:
            return {"ok": False, "error": "genome scan failed — the cached genome may be missing or incomplete; "
                    "re-download it from ⚙ Manage genomes. " + err.strip()[:160]}
        ver = _wsl(f'"{_ENV}/bin/isPcr" 2>&1 | head -1', timeout=20)[1]
        return {"ok": True, "raw": raw, "isPcr_version": _ispcr_ver(ver), "target": target,
                "sha256": meta.get("sha256"), "n_seqs": int(meta.get("n_seqs", 0) or 0)}
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"isPcr timed out after {timeout}s"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        try:
            _wsl(f'rm -rf /tmp/{rid}', timeout=30)             # always clean the staged query dir, even on timeout
        except Exception:
            pass


def genome_list() -> dict:
    """List prepared (cached) genomes: accession, on-disk target, sealed sha256, contig count, bytes."""
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available", "genomes": []}
    # deliver the loop via STDIN, never an inline `bash -lc` arg: wsl.exe rebuilds the Windows command line
    # and mangles the loop variable $d to empty, so an inline for-loop reports EVERY cached genome as missing
    # (reproduced live: inline `for d ... echo [$d]` prints []). STDIN bytes to `bash -l -s` round-trip intact.
    rc, out, _ = _wsl_script('for d in "$HOME"/teagle_genomes/*/; do [ -f "$d/.done" ] && cat "$d/meta.txt" && echo "==="; done 2>/dev/null || true',
                             timeout=30)
    genomes = []
    for block in out.split("==="):
        m = _parse_meta(block)
        if m.get("accession"):
            genomes.append({"accession": m["accession"], "target": m.get("target"), "sha256": m.get("sha256"),
                            "n_seqs": int(m.get("n_seqs", 0) or 0), "bytes": int(m.get("bytes", 0) or 0)})
    return {"ok": True, "genomes": genomes}


def genome_remove(accession: str) -> dict:
    """Delete a cached genome to reclaim disk."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rc, out, _ = _wsl(f'rm -rf "{_GENOMES}/{accession}" && echo REMOVED', timeout=60)
    return {"ok": "REMOVED" in out, "accession": accession}
