"""Managed WSL2 backend — runs the Linux-only annotation stack (RepeatMasker + Dfam)
for family-level TE classification (Layer A homology).

Security: native Windows never runs a Linux shell built from user input. Commands go to
WSL via fixed argument vectors; the user's sequence is piped in as STDIN (data, never
part of a command); the only interpolated values are a self-generated run-id and a
strictly-validated species token. No shell-string concatenation of untrusted input.
"""
from __future__ import annotations
import os, subprocess, re, secrets, threading, hashlib
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


_FAMDB_B, _FAMDB_E = "===TEAGLE_FAMDB_BEGIN===", "===TEAGLE_FAMDB_END==="
_PARTITION_RE = re.compile(r"partition\s+(\d+)\s+\[([^\]]+)\]\s*:\s*(.*)")


def parse_famdb_info(text: str):
    """Parse `famdb.py info` into the identity of the library that is ACTUALLY installed.

    The version alone does not identify the library: famdb ships Dfam in partitions (curated vs
    uncurated, consensus vs HMM, root vs clade), and only the installed subset is searchable. Two
    machines on the same Dfam version with different partitions present return different families
    for the same query, so the partition list is part of the identity and belongs in the seal.
    Returns None when the text carries no Database/Version header (famdb absent or failed to open)."""
    if not text or not text.strip():
        return None
    name = re.search(r"^\s*Database\s*:\s*(.+?)\s*$", text, re.M)
    ver = re.search(r"^\s*Version\s*:\s*(\S+)", text, re.M)      # 'FamDB Creation Format Version' starts with FamDB, so it cannot match
    if not (name and ver):
        return None
    date = re.search(r"^\s*Date\s*:\s*(\S+)", text, re.M)        # likewise distinct from 'FamDB Creation Date'
    fmt = re.search(r"^\s*FamDB Creation Format Version\s*:\s*(\S+)", text, re.M)
    ncon = re.search(r"^\s*Total consensus sequences present\s*:\s*([\d,]+)", text, re.M)
    parts = [m.group(2) for m in _PARTITION_RE.finditer(text) if "not present" not in m.group(3).lower()]
    return {"name": name.group(1).strip(), "version": ver.group(1),
            "date": date.group(1) if date else None,
            "famdbFormat": fmt.group(1) if fmt else None,
            "consensusSequences": int(ncon.group(1).replace(",", "")) if ncon else None,
            "partitions": sorted(parts)}


def env_status() -> dict:
    """Report the annotation stack state inside WSL (RepeatMasker version, Dfam library, minimap2)."""
    av = available()
    st = {**av, "repeatmasker": None, "engine": None, "dfam": False, "minimap2": None,
          "dfam_library": None, "dfam_version": None, "ready": False}
    if not av["wsl2"]:
        return st
    try:
        _famdb = f"{_ENV}/share/RepeatMasker/Libraries/famdb"
        # delivered via STDIN (_wsl_script), NOT inline _wsl: the $()/nested-quote probe below is exactly the
        # class wsl.exe mangles on a `bash -lc <arg>` command line (it would collapse every [ -f ] to 0).
        # famdb.py info rides along in the SAME round trip (it is what resolves the sealed library identity)
        # and is fenced by markers so its output cannot be mistaken for the minimap2 version line.
        rc, out, err = _wsl_script(
            f'{_MM} run -n te RepeatMasker -v 2>/dev/null | head -1\n'
            f'echo "dfam_root=$([ -f "{_famdb}/dfam40.0.h5" ] && echo 1 || echo 0)"\n'
            f'echo "dfam_curated=$([ -f "{_famdb}/dfam40.curated.consensus.0.h5" ] && echo 1 || echo 0)"\n'
            f'echo "{_FAMDB_B}"\n'
            f'{_MM} run -n te famdb.py info 2>/dev/null\n'
            f'echo "{_FAMDB_E}"\n'
            f'[ -x "{_ENV}/bin/minimap2" ] && {_ENV}/bin/minimap2 --version 2>/dev/null\n',
            timeout=90)
        if _FAMDB_B in out and _FAMDB_E in out:            # lift the famdb block out before the other probes parse
            head, rest = out.split(_FAMDB_B, 1)
            fam_txt, tail = rest.split(_FAMDB_E, 1)
            out = head + tail
            st["dfam_library"] = parse_famdb_info(fam_txt)
            st["dfam_version"] = (st["dfam_library"] or {}).get("version")
        m = re.search(r"RepeatMasker version ([\w.]+)", out)
        st["repeatmasker"] = m.group(1) if m else None
        # require BOTH pinned Dfam partitions — a root-only library is incomplete and must not gate annotation
        st["dfam"] = ("dfam_root=1" in out) and ("dfam_curated=1" in out)
        mm = re.search(r"^(\d+\.\d+[\w.-]*)$", out.strip().splitlines()[-1]) if out.strip() else None
        st["minimap2"] = mm.group(1) if mm else None
        st["ready"] = bool(st["repeatmasker"] and st["dfam"])
    except Exception as e:
        st["error"] = str(e)[:120]
    return st


def _int_or_none(tok):
    """RepeatMasker parenthesises the remaining-length column, e.g. '(1234)'."""
    try:
        return int(str(tok).strip().lstrip("(").rstrip(")"))
    except (ValueError, AttributeError):
        return None


def parse_out(text: str):
    """Parse a RepeatMasker .out table into structured hits.

    Columns 11-13 are the CONSENSUS-side coordinates, and RepeatMasker reverses their order by strand:
    a '+' hit reads `begin end (left)` while a 'C' hit reads `(left) begin end`. Reading them positionally
    without that switch yields negative consensus lengths on roughly half of all real hits, so the order
    is resolved from the strand rather than assumed. Column 15 is the fragment ID: RepeatMasker assigns
    one ID to the several lines of a single interrupted alignment, which is what lets a fragmented element
    be reported as one hit instead of N unrelated ones.

    `divergence` is RepeatMasker's RAW substitution percentage — NOT Kimura-corrected and not CpG-adjusted."""
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
        # %del and %ins sit beside %div and were previously discarded; together they separate a diverged
        # copy from a deleted/inserted one, which raw divergence alone cannot.
        pct_del = float(f[2]) if len(f) > 2 and _is_float(f[2]) else None
        pct_ins = float(f[3]) if len(f) > 3 and _is_float(f[3]) else None
        c1, c2, c3 = (f[11], f[12], f[13]) if len(f) > 13 else (None, None, None)
        if strand == "-":                              # 'C' rows are (left) END begin — note the reversal
            c_left, c_end, c_begin = _int_or_none(c1), _int_or_none(c2), _int_or_none(c3)
        else:                                          # '+' rows are begin end (left)
            c_begin, c_end, c_left = _int_or_none(c1), _int_or_none(c2), _int_or_none(c3)
        cons_len = None
        if c_end is not None and c_left is not None:
            cons_len = c_end + c_left                  # consensus total = aligned end + remaining
        coverage = None
        if cons_len and c_begin is not None and c_end is not None and cons_len > 0:
            coverage = round(100.0 * (c_end - c_begin + 1) / cons_len, 1)
        hits.append({
            "score": score, "divergence": divergence, "pct_del": pct_del, "pct_ins": pct_ins,
            "query": f[4], "q_start": q_start - 1, "q_end": q_end, "strand": strand,
            "family": f[9], "class_family": f[10],
            "cons_start": c_begin, "cons_end": c_end, "cons_left": c_left,
            "cons_length": cons_len, "cons_coverage_pct": coverage,
            "fragment_id": f[14] if len(f) > 14 else None,
        })
    return merge_fragments(hits)


def _is_float(tok):
    try:
        float(tok)
        return True
    except (TypeError, ValueError):
        return False


def merge_fragments(hits):
    """Collapse the lines RepeatMasker gave the same fragment ID into one hit.

    An interrupted alignment — an element split by a younger insertion, or by an internal deletion — is
    emitted as several lines sharing one ID. Reported separately they read as N independent copies of the
    family, which is how a single fragmented L1 came to look like a pile of unrelated hits. The merged hit
    spans the outermost query coordinates, keeps the best score, reports summed consensus coverage, and
    records how many pieces it came from so the join is visible rather than silent."""
    out, groups = [], {}
    for h in hits:
        key = (h.get("fragment_id"), h.get("family"), h.get("query"), h.get("strand"))
        if h.get("fragment_id") is None:
            out.append(h)
            continue
        groups.setdefault(key, []).append(h)
    for key, g in groups.items():
        if len(g) == 1:
            out.append(g[0])
            continue
        best = max(g, key=lambda x: x["score"])
        covs = [x["cons_coverage_pct"] for x in g if x["cons_coverage_pct"] is not None]
        merged = dict(best)
        merged.update({
            "q_start": min(x["q_start"] for x in g),
            "q_end": max(x["q_end"] for x in g),
            "cons_start": min((x["cons_start"] for x in g if x["cons_start"] is not None), default=None),
            "cons_end": max((x["cons_end"] for x in g if x["cons_end"] is not None), default=None),
            "cons_coverage_pct": round(min(100.0, sum(covs)), 1) if covs else None,
            "n_fragments": len(g),
            "fragment_note": f"{len(g)} alignment blocks sharing RepeatMasker fragment ID {key[0]} — one "
                             f"interrupted element, not {len(g)} separate copies",
        })
        out.append(merged)
    return sorted(out, key=lambda x: (x["query"], x["q_start"]))


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
    # OPTIONAL uncurated consensus partitions. What they add is entirely lineage-dependent — measured
    # curated-only vs all: Arabidopsis 9 -> 512 models, yeast 9 -> 398, Drosophila 399 -> 998, human
    # 1439 -> 1439 (see _DFAM_CURATED_COVERAGE). For a plant or a fungus a curated-only search finds no
    # transposable element at all. Not installed by default: partition 1 alone is 3.85 GiB compressed
    # and 22.6 GiB unpacked. md5s from Dfam's own .md5 sidecars.
    "dfam_unc_root": ("dfam40.uncurated.consensus.0.h5", "e7092e3ba01d887d4e7f84c86fa2d2ba"),
    "dfam_unc_euk":  ("dfam40.uncurated.consensus.1.h5", "2803414c3420cd8b9ebbc077d78491c4"),
}

# Unpacked size of each partition, in bytes, MEASURED with `stat -c %s` after gunzip (Dfam 4.0,
# 2026-07-31). It cannot be derived at run time: a gzip trailer stores the original size modulo 4 GiB,
# and the eukaryote partition is 22.6 GiB. This is the single source for both the disk-space gate in
# _dfam_step and the "free space needed" figure the installer panel shows, so the two cannot drift.
# The panel previously said ~14 GB for a file that needs 22.6 GiB unpacked: a machine with 15 GB free
# passed the gate, downloaded for an hour, then failed at decompression.
_DFAM_UNPACKED_B = {
    "dfam_root":       447531894,
    "dfam_curated":    155463080,
    "dfam_unc_root":     1865736,
    "dfam_unc_euk":  24242422200,
}

# Compressed size of each archive, from the server's own Content-Length (measured 2026-07-31). The run
# takes the live HEAD when it can; this is the fallback, so a HEAD that fails cannot silently drop the
# download half of the disk-space gate and let a transfer start that has nowhere to land.
_DFAM_ARCHIVE_B = {
    "dfam_root":       60708386,
    "dfam_curated":    28260067,
    "dfam_unc_root":     264500,
    "dfam_unc_euk":  4130166278,
}


def _gib(n: int) -> str:
    """Bytes as GiB for user-facing text, matching how `df` and the install log report free space."""
    return f"{n / 1073741824:.1f} GiB"


# What the CURATED partitions alone cover, per lineage: (curated-only, curated + uncurated) family models.
# MEASURED 2026-07-31 on Dfam 4.0 with the same query dfam_lineage_families() makes — `famdb.py families
# -a -f summary "<species>"` counted over its "len=" lines — first with only dfam40.0.h5 +
# dfam40.curated.consensus.0.h5 present, then with both uncurated partitions added.
#
# This exists because the app used to tell a user "the curated library holds just 9 families for
# Drosophila melanogaster, and copia, gypsy, hobo and mdg1 are not among them". Both halves are false:
# curated-only Drosophila holds 399 models and does contain Copia_I, Copia_LTR, Gypsy_I, Gypsy_LTR, hobo
# and MDG1_I/LTR. The 9 was the yeast figure attached to the wrong organism. Coverage is not a property
# of Dfam as a whole — it is a property of the lineage — so the panels quote this table rather than one
# number, and dfam_lineage_families() reports the live count for the organism actually selected.
_DFAM_CURATED_COVERAGE = {
    "Homo sapiens":             (1439, 1439),
    "Drosophila melanogaster":  (399, 998),
    "Arabidopsis thaliana":     (9, 512),
    "Saccharomyces cerevisiae": (9, 398),
}


def curated_coverage_sentence() -> str:
    """One sentence naming the best- and worst-covered lineages in the measured table, built from that
    table so a figure shown on screen can never drift from the measurement it came from."""
    best = max(_DFAM_CURATED_COVERAGE.items(), key=lambda kv: kv[1][0] / max(kv[1][1], 1))
    worst = min(_DFAM_CURATED_COVERAGE.items(), key=lambda kv: kv[1][0] / max(kv[1][1], 1))
    (bs, (bc, bu)), (ws, (wc, wu)) = best, worst
    same = " — the same set either way" if bc == bu else ""
    return (f"{bs} has {bc} of {bu} family models in the curated partitions alone{same}, "
            f"while {ws} has only {wc} of {wu}.")

_PRELUDE = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || { echo "[teagle] FAILED: cannot cd to HOME"; exit 1; }
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"
# ViennaRNA lives in its OWN env, never in 'te'. Solved against te (python 3.14) the resolver picks
# viennarna 2.4.7 py36 — an old build whose Python bindings cannot load, and a DIFFERENT version from
# the in-process package, which would make the same primer report a different dG depending on install
# route. Letting it re-solve 'te' could also move hmmer/rmblast, which annotate runs SEAL.
VRNAENV="$HOME/micromamba/envs/teagle-vrna"
FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
LOG="$HOME/teagle_wsl_install.log"; : > "$LOG"; exec > >(tee -a "$LOG") 2>&1
LOCK="$HOME/.teagle_install.lock"
# reap a lock orphaned by a crash / reboot / wsl --shutdown / timeout: its recorded PID is dead.
if [ -d "$LOCK" ] && { [ ! -f "$LOCK/pid" ] || ! kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; }; then
  rm -rf "$LOCK" 2>/dev/null
fi
if ! mkdir "$LOCK" 2>/dev/null; then echo "[teagle] FAILED: install already running"; exit 1; fi
echo $$ > "$LOCK/pid"
# BG holds any background helper (a download and its progress watcher). Closing the app kills the WSL
# session and with it this script; without the kill here the watcher is orphaned and keeps appending to
# the log after the lock is gone, so the panel would report progress for a download that is no longer running.
BG=""
cleanup(){ [ -n "$BG" ] && kill $BG 2>/dev/null; rm -rf "$LOCK" 2>/dev/null; return 0; }
trap cleanup EXIT
trap 'cleanup; exit 143' TERM HUP INT
fail(){ echo "[teagle] FAILED: $1"; exit 1; }
# byte count as a size a reader can act on. Fixed MiB renders the 0.3 MiB Dfam root partition as "0 MiB";
# numfmt picks the unit from the number itself, and falls back to MiB where coreutils is too old for it.
hr(){ [ "${1:-0}" -eq 0 ] && { echo 0; return 0; }; numfmt --to=iec --format='%.1f' "$1" 2>/dev/null || echo "$(($1/1048576))M"; }
# Did a download finish? Answered ONLY from positive evidence: curl exited clean, or the file reached a
# known total. Everything else is "partial". This decides whether a failed checksum may delete the file,
# and deleting a multi-gigabyte partial costs the user the whole transfer again — so an unknown total
# (the size HEAD can fail, or a server can answer without Content-Length) must never license a delete.
#   $1 curl exit status   $2 total bytes, 0 = unknown   $3 bytes on disk
transfer_state(){
  [ "$1" -eq 0 ] && { echo complete; return 0; }
  [ "$2" -gt 0 ] && [ "$3" -ge "$2" ] && { echo complete; return 0; }
  echo partial
}
# Free bytes a partition still needs: what is LEFT to fetch (bytes already on disk are not demanded a
# second time on a resume), plus the library it unpacks to, plus a gigabyte of headroom — gunzip holds
# the archive and its output at once. An unknown live total falls back to the measured archive size, so
# a failed size request cannot silently drop the download half of the estimate.
#   $1 live total, 0 = unknown   $2 measured archive size   $3 bytes on disk   $4 unpacked size
space_need(){
  a=$1; [ "$a" -gt 0 ] || a=$2
  l=$(( a > $3 ? a - $3 : 0 ))
  echo $(( l + $4 + 1073741824 ))
}
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
  # A prefix that exists but holds no conda-meta is not an environment, and micromamba refuses it
  # outright ("Non-conda folder exists at prefix - aborting"), permanently. A Dfam partition installed
  # before the environment existed leaves exactly that: its famdb directory creates the prefix. Set the
  # stray content aside, create the environment, then move the libraries back — the download can be
  # several gigabytes, so it is moved rather than copied and never re-fetched.
  STASH=""
  if [ -d "$ENV" ]; then
    echo "[teagle] $ENV exists but is not a conda environment — setting it aside to create one"
    STASH="$ENV.stash.$$"
    mv "$ENV" "$STASH" || fail "cannot set aside the non-environment prefix at $ENV"
  fi
  if ! "$MM" create -y -n te -c conda-forge -c bioconda; then
    mm_reset_cache
    if ! "$MM" create -y -n te -c conda-forge -c bioconda; then
      [ -n "$STASH" ] && mv "$STASH" "$ENV"      # restore what was there; the caller reports the failure
      return 1
    fi
  fi
  if [ -n "$STASH" ]; then
    OLDFAM="$STASH/share/RepeatMasker/Libraries/famdb"
    if [ -d "$OLDFAM" ]; then
      mkdir -p "$FAMDIR"
      for f in "$OLDFAM"/*; do [ -e "$f" ] && mv -f "$f" "$FAMDIR/"; done
      echo "[teagle] kept the Dfam files that were downloaded before the environment existed"
    fi
    rm -rf "$STASH"
  fi
  return 0
}
mm_install(){  # install package(s) into 'te'; clean-index + force-reinstall retry once on a solve failure
  "$MM" install -y -n te -c conda-forge -c bioconda "$@" && return 0
  mm_reset_cache; "$MM" install --force-reinstall -y -n te -c conda-forge -c bioconda "$@"
}
echo "[teagle] START $(date -u +%FT%TZ)"
'''


def _dfam_step(key: str) -> str:
    fname, md5 = _DFAM_FILES[key]
    unpacked = _DFAM_UNPACKED_B[key]
    unpacked_h = _gib(unpacked)
    archive = _DFAM_ARCHIVE_B[key]
    return f'''
echo "[teagle] STEP {key} START"
mkdir -p "$FAMDIR" || fail "mkdir famdb"
cd "$FAMDIR" || fail "cd famdb"
if [ -f "{fname}" ]; then echo "[teagle] {fname} already present — delete it inside WSL to force a fresh download"; else
  # Total size comes from a separate HEAD. A resumed (-C -) GET answers 206 with Content-Length set to
  # the REMAINING bytes, so reading the total off the transfer itself would make every percentage wrong.
  # Retried like the transfer itself: an unknown total costs the percentage AND the completeness evidence
  # that decides whether a failed checksum may delete the file.
  total_b=$(curl -sIL --max-time 60 --retry 3 --retry-delay 2 "{_DFAM_BASE}/{fname}.gz" | tr -d '\\r' | awk 'tolower($1)=="content-length:"{{n=$2}} END{{print n+0}}')
  total_b=${{total_b:-0}}
  have_b=$(stat -c %s "{fname}.gz" 2>/dev/null || echo 0)
  # Gate on what this file still needs. The unpacked size is measured ({unpacked_h}) because a gzip
  # trailer records the original size only modulo 4 GiB and this library is far larger than that.
  need_b=$(space_need "$total_b" {archive} "$have_b" {unpacked})
  avail_b=$(df -B1 --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9')
  # an EMPTY reading means df could not answer and the gate cannot judge; a reading of 0 means a full
  # disk, which must fail rather than be waved through as "unknown"
  if [ -n "$avail_b" ] && [ "$avail_b" -lt "$need_b" ]; then
    fail "not enough disk space in WSL for {fname}: $(hr $avail_b) free, $(hr $need_b) needed (it unpacks to $(hr {unpacked}), and the archive sits beside it until it does)"
  fi
  if [ "$total_b" -gt 0 ]; then
    echo "[teagle] downloading {fname}.gz — $(hr $total_b) total, $(hr $have_b) already here (resumes where it stopped if interrupted)"
  else
    echo "[teagle] downloading {fname}.gz (resumes where it stopped if interrupted)"
  fi
  # curl's own meter is one endless carriage-return line, which the app's log panel cannot tail — a
  # multi-GB download looked frozen there. Silence it (-sS keeps real errors) and report whole lines
  # from a watcher instead, so progress is visible while the transfer runs.
  curl -sS -L --fail -C - --retry 5 --retry-delay 5 -o "{fname}.gz" "{_DFAM_BASE}/{fname}.gz" &
  dl_pid=$!
  (
    while kill -0 "$dl_pid" 2>/dev/null; do
      sleep 15
      kill -0 "$dl_pid" 2>/dev/null || break
      now_b=$(stat -c %s "{fname}.gz" 2>/dev/null || echo 0)
      # a percentage is only printed when the total is known AND consistent with what is on disk;
      # a redirect hop's Content-Length would otherwise produce a figure above 100%
      if [ "$total_b" -gt 0 ] && [ "$now_b" -le "$total_b" ]; then
        echo "[teagle] {fname}.gz  $(hr $now_b) / $(hr $total_b)  ($((now_b*100/total_b))%)"
      else
        echo "[teagle] {fname}.gz  $(hr $now_b) downloaded"
      fi
    done
  ) &
  mon_pid=$!
  BG="$dl_pid $mon_pid"
  wait "$dl_pid"; dl_rc=$?
  kill "$mon_pid" 2>/dev/null; wait "$mon_pid" 2>/dev/null; BG=""
  now_b=$(stat -c %s "{fname}.gz" 2>/dev/null || echo 0)
  [ "$now_b" -gt 0 ] || fail "download {fname} produced no file"
  # the two silent stages warn about their cost only when the file is actually large enough to have one
  slow=""; [ "$now_b" -gt 1073741824 ] && slow=" — a few minutes at this size"
  # An INTERRUPTED transfer must not reach the delete-on-mismatch gate below: that would throw away a
  # multi-GB partial the next run could have resumed, and the panel promises it resumes. Without positive
  # evidence the transfer completed, the archive is still CHECKED — a complete file whose size was never
  # advertised verifies fine, and refusing it outright would loop forever — but it is never deleted.
  verified=0
  if [ "$(transfer_state "$dl_rc" "$total_b" "$now_b")" = partial ]; then
    echo "[teagle] the transfer ended early at $(hr $now_b) — checking whether what is here is nonetheless complete$slow"
    if echo "{md5}  {fname}.gz" | md5sum -c - >/dev/null 2>&1; then
      echo "[teagle] md5 OK {fname} — the file was complete despite the transfer error"
      verified=1
    else
      fail "download interrupted at $(hr $now_b) — the partial file is KEPT; start it again and it resumes from there"
    fi
  fi
  if [ "$verified" -eq 0 ]; then
    echo "[teagle] downloaded $(hr $now_b) — verifying checksum$slow"
    echo "{md5}  {fname}.gz" | md5sum -c - || {{ rm -f "{fname}.gz"; fail "md5 mismatch {fname} — the file transferred completely but its contents are wrong, so it was removed; start it again for a clean download"; }}
    echo "[teagle] md5 OK {fname}"
  fi
  echo "[teagle] decompressing {fname}.gz$slow"
  gunzip -f "{fname}.gz" || fail "gunzip {fname}"
  echo "[teagle] decompressed {fname} — $(hr $(stat -c %s "{fname}" 2>/dev/null || echo 0)) on disk"
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
    # Pinned to the SAME major/minor as the optional in-process package so both routes compute the
    # cross-check with one implementation and one parameter set. A version skew here would silently
    # make the reported dG depend on how ViennaRNA was installed.
    "viennarna": r'''
echo "[teagle] STEP viennarna START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
# ViennaRNA is an OPTIONAL cross-check engine in its own env. A failure here must NOT abort the run:
# the core Dfam / RepeatMasker stack downloads in the steps that follow, and "Install all" promises each
# component installs independently. So it logs a note and continues (the faToTwoBit policy above), and
# the component reports "not installed (optional)" via its own live import probe rather than failing.
if "$MM" run -n teagle-vrna python -c "import RNA" >/dev/null 2>&1; then
  echo "[teagle] ViennaRNA already present"
elif "$MM" create -y -n teagle-vrna -c conda-forge -c bioconda "viennarna>=2.7,<2.8" python=3.12 \
     || { mm_reset_cache; "$MM" create -y -n teagle-vrna -c conda-forge -c bioconda "viennarna>=2.7,<2.8" python=3.12; }; then
  "$MM" run -n teagle-vrna python -c "import RNA" >/dev/null 2>&1 \
    || echo "[teagle] note: ViennaRNA not importable after install (optional primer-QC engine skipped — Primer3 still reports every structure)"
else
  echo "[teagle] note: could not create the teagle-vrna env (optional primer-QC engine skipped — Primer3 still reports every structure)"
fi
"$MM" run -n teagle-vrna python -c "import RNA; print('[teagle] ViennaRNA ' + RNA.__version__)" 2>/dev/null || true
echo "[teagle] STEP viennarna OK"
''',
    "dfam_root": _dfam_step("dfam_root"),
    "dfam_unc_root": _dfam_step("dfam_unc_root"),
    "dfam_unc_euk": _dfam_step("dfam_unc_euk"),
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
_ALL_STEPS = ["micromamba", "repeatmasker", "minimap2", "genomescan", "viennarna", "dfam_root", "dfam_curated", "famdb_conf"]

# component metadata surfaced to the install dialog (order = install order)
_COMP_META = [
    ("wsl2",         "WSL2 + Linux distro",       False, "Windows Subsystem for Linux — hosts the Dfam / RepeatMasker stack."),
    ("micromamba",   "micromamba (conda)",        True,  "Small conda package manager, installed under your Linux home."),
    ("repeatmasker", "RepeatMasker",              True,  "Homology-based TE annotator that names Dfam families."),
    ("minimap2",     "minimap2",                  True,  "Splice-aware aligner for de-novo exon / intron detection."),
    ("genomescan",   "isPcr + NCBI Datasets",     True,  "Local whole-genome in-silico PCR engine + genome downloader."),
    ("viennarna",    "ViennaRNA (primer QC)",     True,  "Second, independent primer secondary-structure engine. Optional: Primer3 alone still reports every structure, this adds the cross-check. Installed in its own environment; not bundled, because its licence forbids redistribution inside AGPL software."),
    ("dfam_root",    "Dfam 4.0 root library",     True,  "Dfam root partition (dfam40.0.h5)."),
    ("dfam_curated", "Dfam 4.0 curated library",  True,  "Dfam curated consensus partition."),
    ("famdb_conf",   "FamDB configuration",       True,  "Points RepeatMasker at the downloaded Dfam library."),
    ("dfam_unc_root", "Dfam uncurated · root (optional)", True,
     "Optional. Uncurated root partition (0.3 MiB). Curated-only search names very few families outside "
     "the root set."),
    ("dfam_unc_euk",  "Dfam uncurated · Eukaryota (optional)", True,
     f"Optional, 3.9 GiB download that unpacks to {_gib(_DFAM_UNPACKED_B['dfam_unc_euk'])} — leave about "
     f"{_gib(_DFAM_UNPACKED_B['dfam_unc_euk'] + 5 * 1073741824)} free in WSL, since the archive and the "
     "library it unpacks to are both on disk for a while. Adds the uncurated eukaryote "
     "families — copia, gypsy, hobo, Ac, Tnt1, Tc1, the yeast Ty elements and most plant/invertebrate TEs "
     "are uncurated in Dfam 4.0 and cannot be named without it. Installing it is only half of what is "
     "needed: RepeatMasker searches curated families ONLY unless it is asked for both, so choose the "
     "'include uncurated families' library option to actually use what this downloads. Measured on baker's "
     "yeast, curated-only searches 9 families and finds no transposable element at all, while including "
     "uncurated searches 421 more and finds Ty elements over 4.5% of the genome. Expect roughly 40-60 minutes: Dfam serves "
     "this file at about 2 MB/s per client, and measured here, splitting it across parallel connections "
     "made it SLOWER, not faster. The download resumes where it stopped if it is interrupted, so a failed "
     "attempt does not start over. There is no smaller alternative: RepeatMasker's RMBLAST engine reads the "
     "consensus component, which Dfam 4.0 ships as only a root partition and this one."),
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
vrv=$("$MM" run -n teagle-vrna python -c "import RNA;print(RNA.__version__)" 2>/dev/null); echo "viennarna=${vrv:-0}"
echo "dfam_root=$([ -f "$FAMDIR/dfam40.0.h5" ] && echo 1 || echo 0)"
echo "dfam_curated=$([ -f "$FAMDIR/dfam40.curated.consensus.0.h5" ] && echo 1 || echo 0)"
echo "dfam_unc_root=$([ -f "$FAMDIR/dfam40.uncurated.consensus.0.h5" ] && echo 1 || echo 0)"
echo "dfam_unc_euk=$([ -f "$FAMDIR/dfam40.uncurated.consensus.1.h5" ] && echo 1 || echo 0)"
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
    vr = kv.get("viennarna", "0")
    # optional: its absence never blocks anything, so it is reported as such rather than as "missing"
    present("viennarna", vr not in ("0", ""), (f"v{vr}" if vr not in ("0", "") else "not installed (optional)"))
    present("dfam_root", kv.get("dfam_root") == "1", "present" if kv.get("dfam_root") == "1" else "missing")
    present("dfam_curated", kv.get("dfam_curated") == "1", "present" if kv.get("dfam_curated") == "1" else "missing")
    for _k in ("dfam_unc_root", "dfam_unc_euk"):      # optional: absent is a normal state, not "missing"
        present(_k, kv.get(_k) == "1", "present" if kv.get(_k) == "1" else "not installed (optional)")
    present("famdb_conf", kv.get("famdb_conf") == "1", "configured" if kv.get("famdb_conf") == "1" else "missing")
    ready = all(comp[k]["ok"] for k in ("repeatmasker", "dfam_root", "dfam_curated"))
    return {"wsl2": True, "installing": kv.get("installing") == "1", "ready": ready,
            "disk_free_gb": kv.get("disk_free_gb"), "components": [comp[c[0]] for c in _COMP_META]}


# which Dfam partitions the default install ships, and which are installed only on request — read off
# the step list rather than restated, so adding a partition cannot leave the integrity probe behind
_DFAM_REQUIRED = [k for k in _DFAM_FILES if k in _ALL_STEPS]
_DFAM_OPTIONAL = [k for k in _DFAM_FILES if k not in _ALL_STEPS]

_INTEGRITY_PROBE = r'''cd "$HOME" 2>/dev/null
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"; FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
echo "=RM="; "$MM" run -n te RepeatMasker -v 2>&1 | head -1
echo "=MM="; "$ENV/bin/minimap2" --version 2>&1 | head -1
echo "=FAMDB="; "$MM" run -n te famdb.py info 2>&1 | grep -iE "version|families|consensus" | head -3
echo "=FILES="; for f in ''' + " ".join(_DFAM_FILES[k][0] for k in _DFAM_REQUIRED) + r'''; do if [ -f "$FAMDIR/$f" ]; then echo "present $f $(stat -c %s "$FAMDIR/$f" 2>/dev/null)"; else echo "MISSING $f"; fi; done
echo "=OPTFILES="; for f in ''' + " ".join(_DFAM_FILES[k][0] for k in _DFAM_OPTIONAL) + r'''; do if [ -f "$FAMDIR/$f" ]; then echo "present $f"; else echo "absent $f"; fi; done
echo "=SCAN="; { [ -x "$ENV/bin/isPcr" ] && echo "isPcr present"; } || echo "isPcr MISSING"; { [ -x "$ENV/bin/datasets" ] && echo "datasets present"; } || echo "datasets MISSING"
echo "=VRNA="; "$MM" run -n teagle-vrna python -c "import RNA;print('ViennaRNA '+RNA.__version__)" 2>&1 | head -1
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
    # optional partitions are REPORTED, never graded: absent is a normal state, and a user who has just
    # spent 40 minutes downloading one should be able to see it named here rather than infer it
    opt = [l.split() for l in sec.get("OPTFILES", []) if l.strip()]
    opt_txt = ", ".join(f"{f} {'present' if s == 'present' else 'not installed'}" for s, f in opt if f)
    scan_txt = " ".join(sec.get("SCAN", []))
    checks = [
        {"name": "RepeatMasker runs", "ok": "RepeatMasker version" in rm_txt, "detail": rm_txt.strip()[:80] or "no version reported"},
        {"name": "minimap2 runs", "ok": bool(re.match(r"^\d+\.\d+", mm_txt)), "detail": mm_txt[:60] or "no version reported"},
        {"name": "FamDB loads", "ok": bool(re.search(r"(?i)version|families|consensus", fam_txt)), "detail": fam_txt.strip()[:90] or "famdb.py info returned nothing"},
        {"name": "Dfam library files present",
         "ok": len(files) >= len(_DFAM_REQUIRED) and all(f.startswith("present") for f in files),
         "detail": ("; ".join(files)[:90] or "missing")
                   + (f" | optional: {opt_txt}" if opt_txt else "")},
        # the whole-genome off-target scan needs isPcr + NCBI datasets; verify them here so a "healthy" deep-check
        # never precedes a scan that fails at first use (same binding-truthfulness class as the fixed unzip gap)
        {"name": "isPcr + NCBI datasets (whole-genome scan)", "ok": "isPcr present" in scan_txt and "datasets present" in scan_txt,
         "detail": scan_txt.strip()[:80] or "not probed"},
    ]
    return {"ok": all(c["ok"] for c in checks), "checks": checks, "raw": out[-1500:]}


def install_log(tail: int = 200) -> str:
    """Tail the in-WSL install log, normalised for a plain-text panel.

    A tool that draws a carriage-return progress meter (micromamba's solver, and curl before the Dfam
    step was changed to report whole lines) writes its entire meter as ONE line. `tail -n` then returns
    a single line of hundreds of kilobytes that no log panel can render usefully, so the read is bounded
    by bytes and each line is collapsed to the text after its last CR — the state the meter ended on."""
    try:
        _, out, _ = _wsl('tail -c 262144 "$HOME/teagle_wsl_install.log" 2>/dev/null || echo "(no install log yet)"', timeout=30)
    except Exception as e:
        return f"(log unavailable: {e})"
    lines = [ln.rpartition("\r")[2] for ln in out.replace("\r\n", "\n").split("\n")]
    return "\n".join(lines[-int(tail):])


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


def resolve_species(species: str, timeout: int = 90) -> dict:
    """Check a species name against famdb BEFORE handing it to RepeatMasker.

    RepeatMasker delegates the lineage lookup to famdb, and famdb rejects an ambiguous name — 'drosophila'
    matches 141 taxa — by printing its own usage text and exiting non-zero. RepeatMasker then exits 255,
    so the user saw a wall of famdb help as the error for what is really 'that name is not specific
    enough'. Resolving first costs one short probe and turns that into an answerable message.
    Returns {ok: True} or {ok: False, error, suggestions}."""
    if not _SPECIES_RE.match(species or ""):
        return {"ok": False, "error": "invalid species token"}
    try:
        rc, out, _ = _wsl_script(f'{_MM} run -n te famdb.py lineage -ad "{species}" 2>&1 | head -20\n',
                                 timeout=timeout)
    except Exception as e:
        return {"ok": True, "unchecked": f"{type(e).__name__}"}   # probe failed: let RepeatMasker try anyway
    low = out.lower()
    if "ambiguous search term" in low:
        n = re.search(r"found (\d+) results", out)
        return {"ok": False, "ambiguous": True,
                "error": (f"“{species}” matches {n.group(1) if n else 'many'} taxa in the Dfam library, so "
                          f"RepeatMasker cannot pick a lineage. Use the full scientific name — for example "
                          f"“Drosophila melanogaster” rather than “drosophila”, or “Homo sapiens”. A common "
                          f"name that maps to exactly one taxon (“human”) also works.")}
    if "no species" in low or "not found" in low:
        return {"ok": False, "error": f"“{species}” was not found in the installed Dfam library. Check the "
                                      f"spelling, or leave the organism blank to search all installed families."}
    return {"ok": True}


def annotate(fasta_text: str, species: str | None = None, threads: int = 4, timeout: int = 600,
             include_uncurated: bool = False) -> dict:
    """Run RepeatMasker (WSL) on the sequence and return family-level hits (Layer A).

    `include_uncurated` decides WHICH families are reachable, so it is part of the result's identity, not
    a convenience. RepeatMasker searches curated families only unless asked for both, and for most
    lineages nearly every family is uncurated — so without this flag the optional uncurated partitions a
    user downloaded are on disk but unreachable, and a blank result means "not searched", not "not there"."""
    sp = ""
    if species:
        if not _SPECIES_RE.match(species):           # validate untrusted input first (hermetic, fast)
            return {"ok": False, "error": "invalid species token"}
        chk = resolve_species(species)
        if not chk.get("ok"):
            return {"ok": False, "error": chk["error"], "ambiguous_species": chk.get("ambiguous", False)}
        sp = f'-species "{species}"' + (" -uncurated" if include_uncurated else "")
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
        # delivered via STDIN (_wsl_script): a multi-word -species value ("Homo sapiens") interpolated into an
        # inline `bash -lc <arg>` is re-split by wsl.exe into `-species Homo` + a stray `sapiens`; STDIN delivery
        # preserves the quotes. The species token is already validated (_SPECIES_RE), so no injection risk.
        # capture RepeatMasker's own exit code (not the trailing cat's) inside the SAME STDIN script — a nonzero
        # exit (unknown lineage, corrupt Dfam shard, OOM) must fail loudly, not be sealed as a clean "0 families".
        script = (f'cd /tmp/{rid} && {_MM} run -n te RepeatMasker -pa {int(threads)} {sp} -qq q.fa '
                  f'>rm.log 2>&1; ec=$?; echo "EXIT $ec"; [ "$ec" = "0" ] && cat q.fa.out 2>/dev/null || tail -8 rm.log')
        rc, out, err = _wsl_script(script, timeout=timeout)
        m = re.search(r"^EXIT (\d+)", out, re.M)
        if m and m.group(1) != "0":
            tail = out.split(f"EXIT {m.group(1)}", 1)[-1].strip()[:250]
            return {"ok": False, "error": f"RepeatMasker exited {m.group(1)}: {tail or 'see backend log'}", "status": st}
        hits = parse_out(out)                                 # repeatmasker version already resolved by env_status() above
        lib = st.get("dfam_library")                          # resolved from famdb in the SAME env_status probe
        return {"ok": True, "hits": hits, "n_hits": len(hits),
                "repeatmasker_version": st["repeatmasker"], "raw_out": out[-4000:],
                "dfam_version": (lib or {}).get("version"), "dfam_library": lib,
                # which family universe was searched — the same library with and without this flag
                # answers a different question, so it travels with the result and into the seal
                "include_uncurated": bool(include_uncurated),
                "library_kind": ("installed Dfam partitions, curated + uncurated" if include_uncurated
                                 else "installed Dfam partitions, curated families only"),
                "species": species or "(all installed families)"}
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"RepeatMasker timed out after {timeout}s"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    finally:
        try:
            _wsl(f'rm -rf /tmp/{rid}', timeout=30)             # always clear the staged query FASTA, even on timeout
        except Exception:
            pass


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
        if rc != 0:                                           # a minimap2 tool failure is not a genuine no-alignment result
            return {"ok": False, "error": f"minimap2 alignment failed (exit {rc}) — not a no-alignment result"}
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
    finally:
        try:
            _wsl(f'rm -rf /tmp/{rid}', timeout=30)             # always clear the staged genomic FASTA, even on timeout
        except Exception:
            pass


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
    finally:
        try:
            _wsl(f'rm -rf /tmp/{rid}', timeout=30)             # always clear the staged genomic + protein FASTA, even on timeout
        except Exception:
            pass


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
# NCBI Datasets transfers of large genomes routinely drop mid-stream, and the API rate-limits (harder when the app
# has also been fetching sequences from NCBI). datasets can't resume — it writes a fresh zip each time — so a blip
# discards the partial and the attempt is retried. Retry 5x with EXPONENTIAL backoff (5,10,20,40s) so a transient
# rate-limit window has time to clear, instead of hammering NCBI on a fixed short interval and exhausting the retries.
ok=0; delay=5
for attempt in 1 2 3 4 5; do
  if "$ENV/bin/datasets" download genome accession "$ACC" --include genome --filename dl.zip >/dev/null 2>dl.err; then ok=1; break; fi
  [ "$attempt" = 5 ] && break
  plog "download attempt $attempt did not complete (transient NCBI/network) — retrying in ${delay}s"; rm -f dl.zip; sleep "$delay"; delay=$((delay*2))
done
# surface the last attempt's stderr (invalid/withdrawn accession, DNS/API/rate-limit) instead of an opaque failure
[ "$ok" = 1 ] || { echo "FAIL download (after 5 attempts): $(tail -c 200 dl.err 2>/dev/null | tr '\n' ' ')"; exit 1; }
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


# ---------------- whole-genome TE annotation (RepeatMasker over a cached assembly) ----------------
# Design constraints this script exists to satisfy:
#   * a mammalian genome yields millions of hits and ~1 GB of .out text, so the per-hit rows NEVER cross
#     into the app — the summary is aggregated here, in awk, and the detail stays on disk in WSL until the
#     user exports it;
#   * work is done in contig CHUNKS so the UI can show N-of-M progress, a cancel can keep what finished,
#     and a re-run resumes instead of restarting;
#   * every chunk's RepeatMasker exit code is checked, so a failed chunk fails the run loudly instead of
#     being summarised as "no repeats found";
#   * the run is only ever sealed as COMPLETE when every chunk finished (see genome_annotate).
_GENOME_ANNOT_LOG = "$HOME/teagle_genome_annotate.log"

_ANNOT_SCRIPT = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || exit 1
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
ENV="$HOME/micromamba/envs/te"
ACC="__ACC__"
GDIR="$HOME/teagle_genomes/$ACC"
ADIR="$GDIR/annot"
LOG="$HOME/teagle_genome_annotate.log"
CHUNKBP=__CHUNKBP__
PA=__PA__
SENS="__SENS__"
plog(){ echo "[annot] $1"; echo "$1" >> "$LOG"; }
[ -f "$GDIR/.done" ] || { echo "FAIL genome not prepared"; exit 1; }
mkdir -p "$ADIR" || { echo "FAIL mkdir"; exit 1; }
LOCK="$ADIR/.lock"
if [ -d "$LOCK" ] && { [ ! -f "$LOCK/pid" ] || ! kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; }; then rm -rf "$LOCK" 2>/dev/null; fi
if ! mkdir "$LOCK" 2>/dev/null; then echo "FAIL an annotation for this genome is already running"; exit 1; fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK" 2>/dev/null' EXIT

# 1) a FASTA is required (the cache keeps a compact 2bit for isPcr; RepeatMasker needs the sequence).
FNA="$GDIR/genome.fna"
if [ ! -s "$FNA" ]; then
  plog "fetching genome FASTA (one time; kept for later runs)"
  cd "$GDIR" || exit 1
  ok=0; delay=5
  for attempt in 1 2 3 4 5; do
    if "$ENV/bin/datasets" download genome accession "$ACC" --include genome --filename adl.zip >/dev/null 2>adl.err; then ok=1; break; fi
    [ "$attempt" = 5 ] && break
    plog "download attempt $attempt did not complete — retrying in ${delay}s"; rm -f adl.zip; sleep "$delay"; delay=$((delay*2))
  done
  [ "$ok" = 1 ] || { echo "FAIL genome download (after 5 attempts): $(tail -c 200 adl.err 2>/dev/null | tr '\n' ' ')"; exit 1; }
  rm -rf aex; mkdir -p aex
  if python3 -m zipfile -e adl.zip aex/ 2>/dev/null; then :
  elif command -v unzip >/dev/null 2>&1 && unzip -oq adl.zip -d aex; then :
  else echo "FAIL extract (need python3 with zipfile, or unzip)"; exit 1; fi
  rm -f adl.zip
  SRC=$(find aex -name '*_genomic.fna' | head -1)
  [ -n "$SRC" ] || { echo "FAIL no genomic fna in package"; exit 1; }
  # the sealed identity is the SOURCE FASTA sha256 recorded at prepare time; refuse to annotate a file
  # that is not the sequence this genome was sealed from, rather than silently annotating something else.
  WANT=$(sed -n 's/^sha256=//p' "$GDIR/meta.txt")
  GOT=$(sha256sum "$SRC" | cut -d' ' -f1)
  if [ -n "$WANT" ] && [ "$WANT" != "$GOT" ]; then rm -rf aex; echo "FAIL fetched FASTA does not match the sealed genome checksum"; exit 1; fi
  mv "$SRC" "$FNA" 2>/dev/null || cp "$SRC" "$FNA"; rm -rf aex
fi

# 2) disk: the chunk set is written up front and each chunk is deleted as it is consumed, so the peak
# extra requirement is one FASTA. Derived from the actual file, never a fixed number.
FBYTES=$(stat -c %s "$FNA" 2>/dev/null || echo 0)
NEEDG=$(( FBYTES / 1073741824 + 2 ))
availg=$(df -BG --output=avail "$GDIR" 2>/dev/null | tail -1 | tr -dc '0-9'); availg=${availg:-0}
[ "$availg" -ge "$NEEDG" ] || { echo "FAIL insufficient disk (${availg}G free, need >=${NEEDG}G to annotate this genome)"; exit 1; }

# 3) split into contig chunks of ~CHUNKBP bases (one streaming pass; a contig is never split in half, so
# every RepeatMasker coordinate stays a real coordinate on a real contig).
# A resumed run MUST use the same library, species, sensitivity and chunk size as the chunks already on
# disk. Otherwise half the assembly is searched one way and half another, and the seal — which records a
# single value for each — would describe a run that never happened. The parameters are written beside the
# chunks at first split and compared on every resume; a mismatch stops the run and says what differs.
RUNSIG="species=__SPSIG__|sens=$SENS|chunkbp=$CHUNKBP"
if [ -f "$ADIR/.runsig" ]; then
  OLDSIG=$(cat "$ADIR/.runsig")
  if [ "$OLDSIG" != "$RUNSIG" ]; then
    echo "FAIL this genome has a part-finished annotation run with different settings ($OLDSIG) than the ones requested ($RUNSIG). Finish or discard that run before starting a different one."
    exit 1
  fi
fi
printf '%s' "$RUNSIG" > "$ADIR/.runsig"
if [ ! -f "$ADIR/.split_done" ]; then
  plog "splitting genome into work chunks"
  rm -f "$ADIR"/chunk_*.fa
  awk -v dir="$ADIR" -v lim="$CHUNKBP" '
    /^>/ { if (n > 0 && bp >= lim) { close(f); i++; bp = 0 }
           if (f == "") { i = 1 } ; f = sprintf("%s/chunk_%04d.fa", dir, i); n++ }
    f == "" { next }                      # sequence before any header: malformed, skip rather than write to ""
    { print >> f }
    !/^>/ { bp += length($0) }
  ' "$FNA" || { echo "FAIL split"; exit 1; }
  touch "$ADIR/.split_done"
fi
# a finished chunk has its .fa deleted and its .out kept, so the work total is pending + finished — counting
# only the .fa files made a fully-completed genome look like "no chunks produced" on a re-run.
PEND_N=$(ls "$ADIR"/chunk_*.fa 2>/dev/null | wc -l)
DONE_N=$(ls "$ADIR"/chunk_*.out 2>/dev/null | wc -l)
TOTAL=$((PEND_N + DONE_N))
[ "$TOTAL" -ge 1 ] || { echo "FAIL no chunks produced"; exit 1; }
plog "chunks: $TOTAL total, $DONE_N already finished"
echo "TOTAL_CHUNKS $TOTAL"

# 4) RepeatMasker per chunk. A chunk whose .out already exists is skipped (resume). Each exit code is
# checked: a nonzero exit aborts the run instead of contributing an empty .out to the summary.
cd "$ADIR" || exit 1
idx=0
for c in chunk_*.fa; do
  [ -e "$c" ] || continue
  idx=$((idx+1))
  base="${c%.fa}"
  [ -s "$base.out" ] && continue
  plog "RepeatMasker chunk $idx/$TOTAL"
  # invoked through micromamba (as annotate() does), not the bare binary: RepeatMasker needs the env's
  # Perl/library variables set, which `micromamba run` provides and a direct exec does not.
  "$HOME/bin/micromamba" run -n te RepeatMasker -pa "$PA" $SENS __SP__ "$c" >"$base.rmlog" 2>&1
  ec=$?
  if [ "$ec" != "0" ]; then echo "FAIL RepeatMasker exited $ec on chunk $idx: $(tail -3 "$base.rmlog" | tr '\n' ' ')"; exit 1; fi
  # RepeatMasker writes <chunk>.fa.out; a chunk with no repeats legitimately produces a header-only file.
  if [ -f "$c.out" ]; then mv "$c.out" "$base.out"; else printf '' > "$base.out"; fi
  rm -f "$c" "$c.masked" "$c.cat" "$c.cat.gz" "$c.tbl" "$c.ori.out" 2>/dev/null
  echo "CHUNK_DONE $idx/$TOTAL"
done

# 5) combine + summarise. The summary is computed HERE so the app never reads the per-hit rows: for each
# repeat class/family, the hit count, the bases covered, and the divergence averaged weighted by length.
plog "summarising"
cat chunk_*.out > all.out 2>/dev/null || true
GENOME_BP=$(awk '!/^>/ { t += length($0) } END { print t+0 }' "$FNA")
echo "GENOME_SHA256 $(sed -n 's/^sha256=//p' "$GDIR/meta.txt")"
# Coverage is computed on MERGED intervals, per family and overall. RepeatMasker legitimately emits
# overlapping rows for the same bases (a younger element inserted into an older one, and its own
# fragmented re-alignments), so summing row lengths double-counts sequence and can push the reported
# coverage above what is actually masked. Rows are sorted by family and start, then merged before the
# bases are counted; divergence stays length-weighted over the raw rows, where each alignment's own
# divergence belongs. The row COUNT is reported as alignment rows, which is what it is.
# sort by family, THEN contig, THEN start — merging needs rows grouped per (family, contig) and ordered
# by position; sorting on family+start alone interleaves contigs and merges intervals that never touch.
sort -k11,11 -k5,5 -k6,6n all.out 2>/dev/null | awk -v gbp="$GENOME_BP" '
  $1 ~ /^[0-9]+$/ {
    div=$2+0; qs=$6+0; qe=$7+0; seqid=$5; fam=$11
    if (qe < qs) { t=qs; qs=qe; qe=t }
    n[fam]++; dv[fam]+=div*(qe-qs+1); wt[fam]+=(qe-qs+1); tot_n++
    key=fam SUBSEP seqid
    if (key != prev) { if (prev != "") { bp[pf]+=(pe-ps+1); gtot+=(pe-ps+1) } ; prev=key; pf=fam; ps=qs; pe=qe }
    else if (qs > pe+1) { bp[pf]+=(pe-ps+1); gtot+=(pe-ps+1); ps=qs; pe=qe }
    else if (qe > pe) { pe=qe }
  }
  END {
    if (prev != "") { bp[pf]+=(pe-ps+1); gtot+=(pe-ps+1) }
    printf "GENOME_BP %d\n", gbp
    printf "TOTAL_HITS %d\nTOTAL_BP %d\n", tot_n+0, gtot+0
    for (f in n) printf "FAM\t%s\t%d\t%d\t%.2f\n", f, n[f], bp[f], (wt[f]>0 ? dv[f]/wt[f] : 0)
  }
'
echo "ANNOT_OK"
'''


_LIBDIR = "$HOME/micromamba/envs/te/share/RepeatMasker/Libraries/famdb"

# RepeatMasker's class/family column mixes transposable elements with tandem and non-TE repeats. Reporting
# one merged "% repeats" as a TE landscape would overstate TE content — a yeast run masks ~1.3% of the
# genome while finding ZERO transposable elements (measured), all of it simple/low-complexity sequence.
_TE_PREFIXES = ("LTR", "LINE", "SINE", "DNA", "RC", "Retroposon")
_TANDEM = ("Simple_repeat", "Low_complexity", "Satellite")


def repeat_kind(family_class: str) -> str:
    """TE / tandem / other for a RepeatMasker class-family token. 'other' covers rRNA, tRNA, snRNA and
    unclassified entries — real repeats, but not evidence of transposable-element content."""
    fam = (family_class or "").strip()
    head = fam.split("/")[0]
    if head in _TANDEM:
        return "tandem"
    if head in _TE_PREFIXES:
        return "TE"
    return "other"


def dfam_lineage_families(species: str) -> dict:
    """How many family models the INSTALLED Dfam partitions can actually supply for this lineage.

    This is the preflight that stops a multi-hour run from quietly producing a scientifically empty result.
    Dfam is partitioned; a machine that installed only the curated partitions holds deep coverage for a few
    heavily curated lineages and almost nothing elsewhere. Curated-only, measured 2026-07-31 (see
    _DFAM_CURATED_COVERAGE): Homo sapiens 1439 models, Drosophila melanogaster 399, Saccharomyces
    cerevisiae 9, Arabidopsis thaliana 9 — so a yeast or Arabidopsis genome annotation cannot find their
    transposable elements no matter how long it runs, and the user has to be told BEFORE starting."""
    if not species or not _SPECIES_RE.match(species):
        return {"ok": False, "error": "invalid species token"}
    # the family lines end with "len=<n>"; everything above them is the CC0 banner, which must not be counted
    script = (f'"{_MM}" run -n te famdb.py -i "{_LIBDIR}" families -a -f summary "{species}" 2>/dev/null '
              f'| grep -c "len="\n')
    try:
        rc, out, _ = _wsl_script(script, timeout=300)
        n = int((out.strip().splitlines() or ["0"])[-1] or 0)
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    return {"ok": True, "species": species, "families": n}


def annotate_budget(genome_bytes: int = 0) -> dict:
    """CPU / RAM / disk the annotation may use, measured INSIDE WSL.

    WSL2 caps cores and memory independently of the Windows host (.wslconfig), so the host's numbers would
    over-promise. RepeatMasker's -pa starts N parallel RMBlast jobs, each needing its own memory, so the
    thread count is the LOWER of what the cores allow and what the RAM allows — the same rule the project
    applies to every other parallel step."""
    try:
        rc, out, _ = _wsl_script(
            'echo "cores=$(nproc 2>/dev/null || echo 0)"\n'
            'echo "memkb=$(awk \'/MemTotal/ {print $2}\' /proc/meminfo 2>/dev/null || echo 0)"\n'
            'echo "availg=$(df -BG --output=avail "$HOME" 2>/dev/null | tail -1 | tr -dc \'0-9\')"\n', timeout=45)
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    d = dict(l.split("=", 1) for l in out.splitlines() if "=" in l)
    cores = int(d.get("cores") or 0)
    mem_gb = round(int(d.get("memkb") or 0) / 1048576, 1)
    avail_gb = int(d.get("availg") or 0)
    # ~2 GB per RMBlast job is the working figure for a Dfam-sized library; keep one core for the UI.
    by_cpu = max(1, cores - 1)
    by_ram = max(1, int(mem_gb // 2))
    threads = max(1, min(by_cpu, by_ram))
    need_gb = int(genome_bytes / 1073741824) + 2 if genome_bytes else None
    return {"ok": True, "cores": cores, "mem_gb": mem_gb, "avail_gb": avail_gb,
            "recommended_threads": threads, "limited_by": ("memory" if by_ram < by_cpu else "cores"),
            "disk_needed_gb": need_gb,
            "disk_ok": (avail_gb >= need_gb) if need_gb else None}


def stage_library(path: str) -> dict:
    """Copy a user-supplied repeat library (FASTA) into WSL and checksum it.

    A custom library is the standard RepeatMasker escape hatch (`-lib`) for a lineage Dfam does not cover
    well. TEagle cannot vouch for its contents, so the run records the file's sha256 and is labelled
    user-supplied: the checksum makes the run reproducible without implying TEagle validated the library."""
    if not path or not os.path.isfile(path):
        return {"ok": False, "error": "library file not found"}
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except Exception as e:
        return {"ok": False, "error": f"could not read the library file: {type(e).__name__}"}
    if not data.lstrip()[:1] == b">":
        return {"ok": False, "error": "the library must be a FASTA file (it does not start with '>')"}
    sha = hashlib.sha256(data).hexdigest()
    dest = f"$HOME/teagle_libs/{sha[:16]}.fa"
    rc, out, err = _wsl(f'mkdir -p "$HOME/teagle_libs" && cat > "{dest}"', stdin=data, timeout=300)
    if rc != 0:
        return {"ok": False, "error": "could not stage the library into WSL: " + err.strip()[:160]}
    n = data.count(b"\n>") + (1 if data.lstrip()[:1] == b">" else 0)
    return {"ok": True, "path": dest, "sha256": sha, "sequences": n, "bytes": len(data),
            "name": os.path.basename(path)}


def genome_annotate(accession: str, species: str | None = None, threads: int = 4,
                    sensitivity: str = "default", chunk_mb: int = 40, timeout: int = 86400,
                    custom_lib: str | None = None, include_uncurated: bool = False) -> dict:
    """Annotate every transposable element RepeatMasker can place in a cached assembly, against the
    installed Dfam library. Chunked, resumable, and summarised inside WSL — the caller receives per-family
    aggregates, never the millions of individual hits (those stay on disk for export).

    `sensitivity` is one of default / quick / slow, mapping to RepeatMasker's own speed-vs-recall flags.
    The run is homology-bound: only families present in the installed library can be found, so the caller
    must state which library was searched next to any coverage claim."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    sens = {"default": "", "quick": "-q", "slow": "-s"}.get(sensitivity)
    if sens is None:
        return {"ok": False, "error": "invalid sensitivity (use default, quick or slow)"}
    # LIBRARY CHOICE. Either the installed Dfam partitions filtered to a lineage (-species), or a
    # user-supplied FASTA library (-lib). RepeatMasker treats these as alternatives, not additions: -lib
    # replaces the database search entirely, so the two are never combined here — combining them would
    # make the result's provenance unstateable.
    lib_info = None
    sp = ""
    if custom_lib:
        staged = stage_library(custom_lib)
        if not staged.get("ok"):
            return staged
        lib_info = staged
        sp = f'-lib "{staged["path"]}"'
    elif species:
        if not _SPECIES_RE.match(species):
            return {"ok": False, "error": "invalid species token"}
        chk = resolve_species(species)
        if not chk.get("ok"):
            return {"ok": False, "error": chk["error"], "ambiguous_species": chk.get("ambiguous", False)}
        # RepeatMasker searches CURATED families only unless told otherwise, so the optional uncurated
        # Dfam partitions are inert without this flag. Measured on S. cerevisiae: curated-only sees 9
        # ancestor families and no lineage-specific ones and finds ZERO transposable elements, while
        # -uncurated sees the same 9 plus 421 lineage-specific families and recovers the Ty1/Copia and
        # Gypsy elements the genome actually carries. Uncurated families are auto-generated and lower
        # confidence, so this stays the caller's decision rather than a silent default.
        sp = f'-species "{species}"' + (" -uncurated" if include_uncurated else "")
    st = env_status()
    if not st["ready"]:
        return {"ok": False, "error": "WSL annotation backend not ready "
                f"(RepeatMasker={st['repeatmasker']}, Dfam={st['dfam']})", "status": st}
    # the signature identifies WHICH library/lineage was searched, so a resume with a different one is
    # caught; the thread count is deliberately NOT part of it (it changes speed, never which families hit)
    # The INSTALLED PARTITION SET belongs in the signature too. Adding a Dfam partition changes which
    # families exist to be found, so chunks finished before an install and chunks finished after it were
    # searched against different libraries — the same "one run, two searches" defect as a changed species,
    # and just as invisible in the finished summary.
    _parts = ",".join((st.get("dfam_library") or {}).get("partitions") or [])
    spsig = ((f"lib:{lib_info['sha256'][:16]}" if lib_info else f"species:{species or 'all'}")
             + f"|parts:{hashlib.sha256(_parts.encode()).hexdigest()[:12]}"
             + f"|unc:{int(bool(include_uncurated))}")     # curated-only vs +uncurated is a different search
    script = (_ANNOT_SCRIPT.replace("__ACC__", accession).replace("__PA__", str(max(1, int(threads))))
              .replace("__CHUNKBP__", str(max(1, int(chunk_mb)) * 1_000_000))
              .replace("__SENS__", sens).replace("__SPSIG__", spsig).replace("__SP__", sp))
    try:
        rc, out, err = _wsl_script(script, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"annotation timed out after {timeout}s — re-run to resume from the "
                                      "chunks that already finished"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    if "ANNOT_OK" not in out:
        fl = next((l for l in out.splitlines() if l.startswith("FAIL")), "") or err.strip()[:200]
        return {"ok": False, "error": "genome annotation failed: " + (fl.replace("FAIL", "").strip() or "unknown error"),
                "partial": True}
    fams, totals = [], {}
    for line in out.splitlines():
        if line.startswith("FAM\t"):
            _, fam, n, bp, div = line.split("\t")
            fams.append({"family": fam, "n": int(n), "bp": int(bp), "divergence": float(div)})
        elif line.startswith(("TOTAL_HITS ", "TOTAL_BP ", "GENOME_BP ", "TOTAL_CHUNKS ")):
            k, v = line.split(); totals[k.lower()] = int(v)
        elif line.startswith("GENOME_SHA256 "):
            totals["genome_sha256"] = line.split(None, 1)[1].strip()
    gbp = totals.get("genome_bp", 0)
    genome_sha = totals.get("genome_sha256")
    by_kind = {"TE": 0, "tandem": 0, "other": 0}
    for f in fams:
        f["percent"] = round(100.0 * f["bp"] / gbp, 4) if gbp else None
        f["kind"] = repeat_kind(f["family"])
        by_kind[f["kind"]] += f["bp"]
    fams.sort(key=lambda f: -f["bp"])
    lib = st.get("dfam_library")
    te_bp = by_kind["TE"]
    n_te_fams = sum(1 for f in fams if f["kind"] == "TE")
    # with a custom library the lineage count is meaningless — what was searched is the user's file
    avail = dfam_lineage_families(species) if (species and not lib_info) else {"ok": False}
    return {"ok": True, "accession": accession, "families": fams, "genome_sha256": genome_sha,
            "total_hits": totals.get("total_hits", 0), "total_bp": totals.get("total_bp", 0),
            "genome_bp": gbp, "chunks": totals.get("total_chunks", 0),
            # masked_percent is EVERY repeat RepeatMasker placed; te_percent is the transposable-element
            # subset. They are reported separately because they are different claims — a genome can mask
            # several percent while containing no detected TE at all.
            "masked_percent": round(100.0 * totals.get("total_bp", 0) / gbp, 3) if gbp else None,
            "te_bp": te_bp, "te_percent": round(100.0 * te_bp / gbp, 3) if gbp else None,
            "tandem_bp": by_kind["tandem"], "other_bp": by_kind["other"],
            "te_family_count": n_te_fams,
            "library_families_for_species": avail.get("families") if avail.get("ok") else None,
            # the honest headline when a run completes but the installed library had nothing to find with
            "library_kind": ("user-supplied FASTA library" if lib_info else
                             ("installed Dfam partitions, curated + uncurated" if include_uncurated
                              else "installed Dfam partitions, curated only")),
            "custom_library": lib_info, "include_uncurated": bool(include_uncurated),
            "coverage_warning": (
                None if n_te_fams else
                (("No transposable-element family was found. The supplied library "
                  f"({lib_info['name']}, {lib_info['sequences']} sequences) placed no TE in this assembly — "
                  "check that the library matches this organism before drawing any conclusion.")
                 if lib_info else
                 ("No transposable-element family was found. This is a limit of what was SEARCHED, not a "
                  "finding about the genome"
                  + (f" — the installed Dfam partitions hold {avail['families']} family model(s) for "
                     f"{species or 'this lineage'}" if avail.get("ok") else "")
                  + (", and this run used only the CURATED subset of them. Most families outside a few "
                     "heavily curated species are uncurated in Dfam, so re-running with uncurated families "
                     "included is the first thing to try: measured on S. cerevisiae, curated-only sees 9 "
                     "families and finds no transposable element, while including uncurated sees 421 more "
                     "and recovers the Ty elements the genome carries."
                     if not include_uncurated else
                     ". Uncurated families were already included, so consider whether Dfam covers this "
                     "lineage at all, or supply a species-specific library.")))),
            "repeatmasker_version": st["repeatmasker"], "dfam_version": (lib or {}).get("version"),
            "dfam_library": lib, "species": species or "(all installed families)",
            "sensitivity": sensitivity, "threads": int(threads), "chunk_mb": int(chunk_mb),
            "complete": True}


def genome_annotate_reset(accession: str) -> dict:
    """Discard a part-finished annotation so a run with different settings can start clean.

    The refusal that sends a user here exists because finished chunks and new chunks would otherwise have
    been searched differently; the only safe way to change the settings is to drop the finished work.
    Removes the annotation working directory only — the cached genome and its FASTA are untouched."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    rc, out, err = _wsl_script(f'A="{_GENOMES}/{accession}/annot"\n'
                              f'if [ -d "$A/.lock" ] && [ -f "$A/.lock/pid" ] && kill -0 "$(cat "$A/.lock/pid")" 2>/dev/null; then\n'
                              f'  echo RUNNING; exit 3\n'
                              f'fi\n'
                              f'rm -rf "$A" && echo RESET\n', timeout=120)
    if rc == 3 or "RUNNING" in out:
        return {"ok": False, "error": "an annotation is still running for this genome — let it finish or "
                                      "close it before discarding the results"}
    if "RESET" not in out:
        return {"ok": False, "error": "could not clear the previous annotation: " + err.strip()[:160]}
    return {"ok": True, "accession": accession}


def genome_annotate_log(tail: int = 1) -> str:
    """Tail the annotation milestone log — the UI's N-of-M progress line during a multi-hour run."""
    try:
        _, out, _ = _wsl(f'tail -n {int(tail)} "{_GENOME_ANNOT_LOG}" 2>/dev/null || true', timeout=15)
        return out.strip()
    except Exception:
        return ""


def genome_annotate_status(accession: str) -> dict:
    """How much of an annotation is already on disk (drives resume + the cost dialog's honesty)."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    rc, out, _ = _wsl_script(
        f'A="{_GENOMES}/{accession}/annot"\n'
        f'echo "chunks=$(ls "$A"/chunk_*.fa 2>/dev/null | wc -l)"\n'
        f'echo "done=$(ls "$A"/chunk_*.out 2>/dev/null | wc -l)"\n'
        f'echo "fasta=$([ -s "{_GENOMES}/{accession}/genome.fna" ] && echo 1 || echo 0)"\n', timeout=30)
    d = dict(l.split("=", 1) for l in out.splitlines() if "=" in l)
    return {"ok": True, "pending_chunks": int(d.get("chunks", 0) or 0),
            "finished_chunks": int(d.get("done", 0) or 0), "fasta_cached": d.get("fasta") == "1"}


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
