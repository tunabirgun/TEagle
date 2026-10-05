"""Managed WSL2 backend — runs the Linux-only annotation stack (RepeatMasker + Dfam)
for optional family-level homology naming and genome annotation.

Security: native Windows never runs a Linux shell built from user input. Commands go to
WSL via fixed argument vectors; the user's sequence is piped in as STDIN (data, never
part of a command); the only interpolated values are a self-generated run-id and a
strictly-validated species token. No shell-string concatenation of untrusted input.
"""
from __future__ import annotations
import os, subprocess, re, secrets, threading, hashlib, base64, shutil, math, time
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


def _remembered_distro() -> str | None:
    try:
        with open(_DISTRO_FILE, encoding="utf-8") as f:
            return f.read().strip() or None
    except OSError:
        return None


def set_distro(name: str):
    """Pin the distro the backend uses (for a user choice) and drop the cached resolution."""
    global _distro_cache
    _persist_distro(name)
    _distro_cache = None


def _decode(b: bytes) -> str:
    # wsl.exe list commands emit UTF-16LE; normal command output is UTF-8
    if b[:1] == b"\xff" or (len(b) > 1 and b[1] == 0):
        return b.decode("utf-16-le", "ignore")
    return b.decode("utf-8", "ignore")


# Container runtimes register their own minimal distros: they have no bash, are not a place to install
# conda packages, and Docker Desktop makes one the Windows default.
_RUNTIME_DISTROS = ("docker-desktop", "docker-desktop-data", "rancher-desktop", "rancher-desktop-data",
                    "podman-machine-default")


def is_runtime_distro(name: str) -> bool:
    return (name or "").strip().lower() in _RUNTIME_DISTROS


def parse_distro_table(text: str):
    """(names, default) from the text of `wsl -l -v`; only rows of the NAME/STATE/VERSION table count."""
    default, names = None, []
    for line in (text or "").replace("\x00", "").splitlines():
        s = line.strip()
        if not s:
            continue
        row = re.fullmatch(r"(.+?)\s+\S+\s+[12]", s.lstrip("*").strip())
        if row:
            name = row.group(1).strip()
            names.append(name)
            if s.startswith("*"):
                default = name
    return names, default


def _distro_usable(name: str):
    """True when the distro has bash and python3 (what the backend scripts need), False when it lacks
    either, None when it could not be asked in time (a cold start): unknown is not a verdict."""
    try:
        p = subprocess.run(["wsl.exe", "-d", name, "--", "sh", "-c",
                            "command -v bash >/dev/null && command -v python3 >/dev/null"],
                           capture_output=True, timeout=60, creationflags=_NO_WINDOW)
    except subprocess.TimeoutExpired:
        return None
    except Exception:
        return False
    return p.returncode == 0


_last_found: list[str] = []


def resolve_distro():
    """The WSL distro the backend runs in: the remembered one, else the Windows default, else the first
    listed, skipping container-runtime distros and any distro that lacks bash or python3."""
    global _distro_cache, _last_found
    if _distro_cache is not None:
        return _distro_cache
    try:
        out = subprocess.run(["wsl.exe", "-l", "-v"], capture_output=True, timeout=15, creationflags=_NO_WINDOW)
        if out.returncode != 0:
            _distro_cache = None
            return None
        names, default = parse_distro_table(_decode(out.stdout))
        _last_found = names
        eligible = [n for n in names if not is_runtime_distro(n)]
        order = []
        for n in (_remembered_distro(), default, *eligible):
            if n in eligible and n not in order:
                order.append(n)
        _distro_cache = next((n for n in order if _distro_usable(n) is not False), None)
    except Exception:
        _distro_cache = None
    return _distro_cache


def _windows_backing_free_bytes(distro: str | None = None) -> int | None:
    """Free bytes on the Windows volume that must grow the selected WSL2 virtual disk.

    WSL's ext4 filesystem can advertise nearly a terabyte while the sparse VHDX sits on a Windows
    volume with only a few gigabytes left. The registry BasePath is the authoritative location. If an
    older/custom registration does not expose it, the user-profile volume is a conservative fallback.
    Paths never leave this helper or enter logs/provenance.
    """
    if os.name != "nt":
        return None
    distro = distro or resolve_distro()
    try:
        import winreg
        root_path = r"Software\Microsoft\Windows\CurrentVersion\Lxss"
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, root_path) as root:
            index = 0
            while True:
                try:
                    child_name = winreg.EnumKey(root, index)
                except OSError:
                    break
                index += 1
                try:
                    with winreg.OpenKey(root, child_name) as child:
                        name = str(winreg.QueryValueEx(child, "DistributionName")[0])
                        if distro and name.casefold() != distro.casefold():
                            continue
                        base_path = os.path.expandvars(str(winreg.QueryValueEx(child, "BasePath")[0]))
                        return int(shutil.disk_usage(base_path).free)
                except (OSError, ValueError):
                    continue
    except (ImportError, OSError):
        pass
    try:
        return int(shutil.disk_usage(os.path.expanduser("~")).free)
    except OSError:
        return None


def _effective_free_bytes(wsl_free: int, host_free: int | None) -> int:
    """The space an expanding WSL2 VHDX can really use: the smaller positive limit."""
    wsl_free = max(0, int(wsl_free or 0))
    if host_free is None or int(host_free) < 0:
        return wsl_free
    return min(wsl_free, int(host_free))


def _disk_capacity(wsl_free: int) -> dict:
    host_free = _windows_backing_free_bytes()
    effective = _effective_free_bytes(wsl_free, host_free)
    gib = 1073741824
    wsl_gb = round(wsl_free / gib, 1)
    host_gb = round(host_free / gib, 1) if host_free is not None else None
    effective_gb = round(effective / gib, 1)
    limited_by = ("Windows backing volume" if host_free is not None and host_free < wsl_free
                  else "WSL filesystem")
    note = f"Usable backend storage: {effective_gb:.1f} GiB ({limited_by}). WSL reports {wsl_gb:.1f} GiB"
    if host_gb is not None:
        note += f"; its Windows backing volume has {host_gb:.1f} GiB free."
    else:
        note += ". Windows backing-volume capacity could not be measured."
    return {"wsl_free_bytes": wsl_free, "windows_free_bytes": host_free, "effective_free_bytes": effective,
            "wsl_disk_free_gb": wsl_gb, "windows_disk_free_gb": host_gb,
            "disk_free_gb": effective_gb, "disk_limited_by": limited_by, "disk_note": note}


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


_WSL_ERROR_CODE = re.compile(r"0x8[0-9a-f]{7}|Wsl/\w+|WSL_E_\w+|ERROR_\w+", re.I)


def _classify_probe_failure(rc: int, out: str, err: str) -> str:
    """Why `bash -lc 'echo ok'` failed. Only a WSL registration/start error code confirms a broken distro;
    a slow cold start, a restarting VM and a distro without bash are different problems with different fixes."""
    text = f"{out} {err}".lower()
    if rc == 127 or ("bash" in text and ("not found" in text or "no such file" in text)):
        return "no_bash"
    if "terminated" in text or "shut down" in text:
        return "starting"
    if _WSL_ERROR_CODE.search(f"{out} {err}"):
        return "start_failed"
    return "unknown"


def available() -> dict:
    global _distro_cache
    distro = resolve_distro()
    if not distro:
        if _last_found:
            return {"wsl2": False, "distro": None, "reason": "no_suitable_distro", "found": list(_last_found),
                    "error": "no usable Linux distribution (found: " + ", ".join(_last_found) + ")"}
        return {"wsl2": False, "distro": None, "error": "WSL not installed / no distro"}
    try:
        rc, out, err = _wsl("echo ok", timeout=30)
        ok = rc == 0 and "ok" in out
        if not ok:
            _distro_cache = None                     # cached distro no longer usable -> re-resolve on the next probe
            return {"wsl2": False, "distro": distro, "reason": _classify_probe_failure(rc, out, err),
                    "error": (err.strip() or out.strip() or f"exit {rc}")[:120]}
        _persist_distro(distro)                      # record the working distro so uninstall/clean targets the right one
        return {"wsl2": True, "distro": distro, "error": None}
    except subprocess.TimeoutExpired:
        _distro_cache = None
        return {"wsl2": False, "distro": distro, "reason": "timeout", "error": "the distro did not answer in 30 s"}
    except Exception as e:
        _distro_cache = None                          # distro changed/removed mid-session -> drop the stale cache
        return {"wsl2": False, "distro": distro, "reason": "unknown", "error": str(e)[:120]}


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
    st = {**av, "repeatmasker": None, "engine": None, "dfam": False, "dfam_uncurated": False,
          "minimap2": None, "dfam_library": None, "dfam_version": None, "ready": False}
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
            f'echo "dfam_unc_root=$([ -f "{_famdb}/dfam40.uncurated.consensus.0.h5" ] && echo 1 || echo 0)"\n'
            f'echo "dfam_unc_euk=$([ -f "{_famdb}/dfam40.uncurated.consensus.1.h5" ] && echo 1 || echo 0)"\n'
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
        # The optional search is truthful only when both uncurated consensus partitions are present.
        st["dfam_uncurated"] = ("dfam_unc_root=1" in out) and ("dfam_unc_euk=1" in out)
        mm = re.search(r"^(\d+\.\d+[\w.-]*)$", out.strip().splitlines()[-1]) if out.strip() else None
        st["minimap2"] = mm.group(1) if mm else None
        st["ready"] = bool(st["repeatmasker"] and st["dfam"] and st["dfam_library"])
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
_CONDA_PINS = {
    "repeatmasker": "4.2.4",
    "minimap2": "2.31",
    "miniprot": "0.18",
    "ispcr": "33",
    "ncbi-datasets-cli": "18.34.0",
    "ucsc-fatotwobit": "482",
    "viennarna": "2.7.2",
}
_MICROMAMBA_VERSION = "2.8.1"
_REPEATMASKER_VERSION = "4.2.5"
_REPEATMASKER_SOURCE_URL = "https://github.com/Dfam-consortium/RepeatMasker/archive/refs/tags/v4.2.5.tar.gz"
_REPEATMASKER_SOURCE_SHA256 = "2a863a62c01d631f733ee0e9d01278241a96dfc168dc8eb974c43a44bcf1d570"
_MICROMAMBA_URL = ("https://github.com/mamba-org/micromamba-releases/releases/download/"
                   "2.8.1-0/micromamba-linux-64.tar.bz2")
_MICROMAMBA_SHA256 = "a934c3709c997feae403a27fd1e321c106d26ffa4f294800ffb11cbc9a3e8515"

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


def dfam_coverage_counts(species: str) -> tuple[int, int, int]:
    """Curated, total, and additional-uncurated family counts for a measured lineage."""
    curated, total = _DFAM_CURATED_COVERAGE[species]
    return curated, total, total - curated


def curated_coverage_sentence() -> str:
    """One sentence naming the best- and worst-covered lineages in the measured table, built from that
    table so a figure shown on screen can never drift from the measurement it came from."""
    best = max(_DFAM_CURATED_COVERAGE.items(), key=lambda kv: kv[1][0] / max(kv[1][1], 1))
    worst = min(_DFAM_CURATED_COVERAGE.items(), key=lambda kv: kv[1][0] / max(kv[1][1], 1))
    (bs, (bc, bu)), (ws, (wc, wu)) = best, worst
    same = " — the same set either way" if bc == bu else ""
    return (f"{bs} has {bc} of {bu} family models in the curated partitions alone{same}, "
            f"while {ws} has only {wc} of {wu}.")


_DFAM_YEAST_CURATED, _DFAM_YEAST_TOTAL, _DFAM_YEAST_ADDITIONAL = \
    dfam_coverage_counts("Saccharomyces cerevisiae")

_PRELUDE = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || { echo "[teagle] FAILED: cannot cd to HOME"; exit 1; }
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"
HOST_FREE_B=__HOST_FREE_B__
effective_free_b(){
  local wsl_b=${1:-0}
  if [ "$HOST_FREE_B" -ge 0 ] && [ "$HOST_FREE_B" -lt "$wsl_b" ]; then echo "$HOST_FREE_B"; else echo "$wsl_b"; fi
}
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
# 'te' lives in micromamba's shared default root, so a same-named environment may belong to someone else.
# It is adopted only when TEagle made it: the marker, or the pinned RepeatMasker an earlier TEagle install left.
OWNED="$ENV/.teagle-owned"
env_owned(){
  [ -f "$OWNED" ] && return 0
  compgen -G "$ENV/conda-meta/repeatmasker-__RM_PIN__-*.json" >/dev/null || return 1
  : > "$OWNED" 2>/dev/null; echo "[teagle] adopted the existing 'te' environment from an earlier TEagle install"
}
mm_create(){   # create the shared 'te' env if absent; clean-index + retry once on a solve failure
  if [ -d "$ENV/conda-meta" ]; then
    env_owned && return 0
    echo "[teagle] an environment named 'te' already exists at $ENV and was not created by TEagle, so it was left untouched. Rename or remove it (micromamba env remove -n te), then run this step again."
    return 1
  fi
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
  : > "$OWNED"
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
mm_install(){  # called only after a functional/version probe fails, so restore files even if metadata survives
  "$MM" install --force-reinstall -y -n te -c conda-forge -c bioconda "$@" && return 0
  mm_reset_cache; "$MM" install --force-reinstall -y -n te -c conda-forge -c bioconda "$@"
}
echo "[teagle] START $(date -u +%FT%TZ)"
[ -d "$ENV/conda-meta" ] && { env_owned || true; }
'''.replace("__RM_PIN__", _CONDA_PINS["repeatmasker"])


def _dfam_step(key: str) -> str:
    fname, md5 = _DFAM_FILES[key]
    unpacked = _DFAM_UNPACKED_B[key]
    unpacked_h = _gib(unpacked)
    archive = _DFAM_ARCHIVE_B[key]
    return f'''
echo "[teagle] STEP {key} START"
mkdir -p "$FAMDIR" || fail "mkdir famdb"
cd "$FAMDIR" || fail "cd famdb"
dfam_file_valid() {{
  [ -f "$1" ] || return 1
  [ "$(stat -c %s "$1" 2>/dev/null || echo 0)" = "$2" ] || return 1
  # Return 0 only for a verified file, 1 for a definitively invalid file, and 2 when the verifier itself
  # is unavailable/broken. A verifier failure must preserve the multi-gigabyte destination for repair.
  local marker rc
  marker=$(mktemp "${{TMPDIR:-/tmp}}/teagle-dfam-verify.XXXXXX") || return 2
  if [ -x "$MM" ] && [ -d "$ENV/conda-meta" ]; then
    "$MM" run -n te python - "$1" "$marker" <<'PY'
import sys
from pathlib import Path

def invalid(message):
    print(message, file=sys.stderr)
    Path(sys.argv[2]).write_text("invalid", encoding="ascii")
    raise SystemExit(1)

try:
    import h5py
except Exception as exc:
    print(f"HDF5 verifier unavailable: {{exc}}", file=sys.stderr)
    raise SystemExit(2)
try:
    with h5py.File(sys.argv[1], "r") as db:
        if not list(db.keys()):
            invalid("invalid HDF5 file: no top-level groups")
except (OSError, ValueError) as exc:
    invalid(f"invalid HDF5 file: {{exc}}")
PY
  else
    python3 - "$1" "$marker" <<'PY'
import sys
from pathlib import Path
try:
    with open(sys.argv[1], "rb") as stream:
        header = stream.read(8)
except OSError as exc:
    print(f"HDF5 header verifier unavailable: {{exc}}", file=sys.stderr)
    raise SystemExit(2)
if header != bytes.fromhex("894844460d0a1a0a"):
    Path(sys.argv[2]).write_text("invalid", encoding="ascii")
    raise SystemExit(1)
PY
  fi
  rc=$?
  if [ "$rc" -eq 0 ]; then rm -f "$marker"; return 0; fi
  if [ -s "$marker" ]; then rm -f "$marker"; return 1; fi
  rm -f "$marker"
  return 2
}}
need_download=1
if [ -f "{fname}" ]; then
  dfam_file_valid "{fname}" {unpacked}; dfam_state=$?
  if [ "$dfam_state" -eq 0 ]; then
    echo "[teagle] {fname} verified — pinned name, exact size, and readable HDF5 structure"
    need_download=0
  elif [ "$dfam_state" -eq 1 ]; then
    echo "[teagle] {fname} is corrupt or does not match the pinned Dfam 4.0 partition — replacing it"
    rm -f "{fname}" || fail "remove invalid {fname}"
  else
    fail "could not verify {fname} because the HDF5 verifier is unavailable; the existing file was preserved — repair micromamba/FamDB, then run Repair again"
  fi
fi
if [ "$need_download" -eq 1 ]; then
  # A genuinely fresh Ubuntu WSL has python3 but often has neither curl nor wget. Keep curl as the
  # preferred path when present, but use only the Python standard library otherwise. The fallback sends
  # a Range request when a partial archive exists; if a server ignores it, a completed full response
  # replaces the partial atomically instead of corrupting it by append.
  dfam_head_bytes(){{
    if command -v curl >/dev/null 2>&1; then
      curl -sIL --max-time 60 --retry 3 --retry-delay 2 "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{{n=$2}} END{{print n+0}}'
    elif command -v python3 >/dev/null 2>&1; then
      python3 - "$1" <<'PY'
import sys
from urllib.request import Request, urlopen
try:
    with urlopen(Request(sys.argv[1], method="HEAD"), timeout=60) as response:
        print(int(response.headers.get("Content-Length") or 0))
except Exception:
    print(0)
PY
    else
      echo 0
    fi
  }}
  dfam_download(){{
    dest=$1; url=$2
    if command -v curl >/dev/null 2>&1; then
      curl -sS -L --fail -C - --retry 5 --retry-delay 5 -o "$dest" "$url"
    elif command -v python3 >/dev/null 2>&1; then
      python3 - "$dest" "$url" <<'PY'
import os, sys, time
from urllib.request import Request, urlopen

dest, url = sys.argv[1:3]
attempts = 3
last_error = None
for attempt in range(1, attempts + 1):
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {{"Range": f"bytes={{have}}-"}} if have else {{}}
    restart = None
    try:
        with urlopen(Request(url, headers=headers), timeout=180) as response:
            status = response.getcode()
            content_range = response.headers.get("Content-Range", "")
            if have and status == 206:
                if not content_range.startswith(f"bytes {{have}}-"):
                    raise ValueError(f"malformed Content-Range for resume: {{content_range!r}}")
                target, mode = dest, "ab"
            elif have and status == 200:
                # A server that ignores Range can still supply a clean archive. Download it beside the
                # retained partial and replace atomically only after the full response has been read.
                restart, target, mode = dest + ".restart", dest + ".restart", "wb"
            elif not have and status == 200:
                target, mode = dest, "wb"
            else:
                raise ValueError(f"unexpected HTTP status {{status}} for {{'resume' if have else 'download'}}")
            advertised = response.headers.get("Content-Length")
            expected = int(advertised) if advertised is not None else None
            written = 0
            with open(target, mode) as out:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    out.write(block)
                    written += len(block)
            if expected is not None and written != expected:
                raise EOFError(f"short HTTP response: received {{written}} of {{expected}} bytes")
        if restart:
            os.replace(restart, dest)
        break
    except Exception as exc:
        last_error = exc
        if restart:
            try:
                os.remove(restart)
            except FileNotFoundError:
                pass
        if attempt < attempts:
            time.sleep(attempt)
else:
    raise SystemExit(f"download failed after {{attempts}} attempts; partial archive was kept: {{last_error}}")
PY
    elif command -v wget >/dev/null 2>&1; then
      wget -q --continue -O "$dest" "$url"
    else
      echo "python3, curl, and wget are unavailable" >&2
      return 127
    fi
  }}
  # Total size comes from a separate HEAD. A resumed (-C -) GET answers 206 with Content-Length set to
  # the REMAINING bytes, so reading the total off the transfer itself would make every percentage wrong.
  # Retried like the transfer itself: an unknown total costs the percentage AND the completeness evidence
  # that decides whether a failed checksum may delete the file.
  total_b=$(dfam_head_bytes "{_DFAM_BASE}/{fname}.gz")
  total_b=${{total_b:-0}}
  have_b=$(stat -c %s "{fname}.gz" 2>/dev/null || echo 0)
  # Gate on what this file still needs. The unpacked size is measured ({unpacked_h}) because a gzip
  # trailer records the original size only modulo 4 GiB and this library is far larger than that.
  need_b=$(space_need "$total_b" {archive} "$have_b" {unpacked})
  wsl_avail_b=$(df -B1 --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9')
  avail_b=$(effective_free_b "${{wsl_avail_b:-0}}")
  # an EMPTY reading means df could not answer and the gate cannot judge; a reading of 0 means a full
  # disk, which must fail rather than be waved through as "unknown"
  if [ -n "$avail_b" ] && [ "$avail_b" -lt "$need_b" ]; then
    fail "not enough usable disk for {fname}: $(hr $avail_b) free after applying the smaller WSL/Windows backing-volume limit, $(hr $need_b) needed (it unpacks to $(hr {unpacked}), and the archive sits beside it until it does)"
  fi
  if [ "$total_b" -gt 0 ]; then
    echo "[teagle] downloading {fname}.gz — $(hr $total_b) total, $(hr $have_b) already here (resumes where it stopped if interrupted)"
  else
    echo "[teagle] downloading {fname}.gz (resumes where it stopped if interrupted)"
  fi
  # curl's own meter is one endless carriage-return line, which the app's log panel cannot tail — a
  # multi-GB download looked frozen there. Silence it (-sS keeps real errors) and report whole lines
  # from a watcher instead, so progress is visible while the transfer runs.
  dfam_download "{fname}.gz" "{_DFAM_BASE}/{fname}.gz" &
  dl_pid=$!
  (
    while kill -0 "$dl_pid" 2>/dev/null; do
      sleep 15
      kill -0 "$dl_pid" 2>/dev/null || break
      active="{fname}.gz"; [ -f "{fname}.gz.restart" ] && active="{fname}.gz.restart"
      now_b=$(stat -c %s "$active" 2>/dev/null || echo 0)
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
  dfam_file_valid "{fname}" {unpacked}; dfam_state=$?
  if [ "$dfam_state" -eq 1 ]; then
    rm -f "{fname}"
    fail "decompressed {fname} failed the pinned size/HDF5 identity check and was removed"
  elif [ "$dfam_state" -ne 0 ]; then
    fail "decompressed {fname} could not be verified because the HDF5 verifier is unavailable; the file was preserved — repair micromamba/FamDB, then run Repair again"
  fi
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
MM_VERSION="__MM_VERSION__"
MM_URL="__MM_URL__"
MM_SHA256="__MM_SHA256__"
mm_py3(){ python3 - "$MM" "$MM_URL" "$MM_SHA256" <<'PY'
import hashlib, io, os, stat, sys, tarfile, urllib.request
dest, url, expected = sys.argv[1:4]
data = urllib.request.urlopen(url, timeout=180).read()
got = hashlib.sha256(data).hexdigest()
if got != expected: raise SystemExit(f"micromamba sha256 mismatch: {got}")
with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as tf:
    m = tf.extractfile(tf.getmember("bin/micromamba"))
    if m is None: raise SystemExit("bin/micromamba not in archive")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    open(dest, "wb").write(m.read())
os.chmod(dest, os.stat(dest).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
PY
}
mm_fetch(){
  arc="$HOME/micromamba-$MM_VERSION.tar.bz2"
  rm -f "$arc"
  if command -v curl >/dev/null 2>&1; then curl -fL --retry 3 -o "$arc" "$MM_URL"
  elif command -v wget >/dev/null 2>&1; then wget -qO "$arc" "$MM_URL"
  else return 1; fi
  echo "$MM_SHA256  $arc" | sha256sum -c - || { rm -f "$arc"; return 1; }
  tar -xj -C "$HOME" -f "$arc" bin/micromamba || { rm -f "$arc"; return 1; }
  rm -f "$arc"; chmod +x "$MM"
}
if [ -x "$MM" ] && [ "$("$MM" --version 2>/dev/null)" = "$MM_VERSION" ]; then
  echo "[teagle] micromamba $MM_VERSION already present"
else
  # ~/bin/micromamba is where the official manual install puts it, so an existing one may be the user's own
  if [ -e "$MM" ] && [ ! -e "$MM.pre-teagle" ] && cp -p "$MM" "$MM.pre-teagle"; then
    echo "[teagle] kept the micromamba already at $MM as $MM.pre-teagle (it is not version $MM_VERSION)"
  fi
  FOUND=""
  for c in "$HOME/.local/bin/micromamba" "$HOME/micromamba/bin/micromamba" "$(command -v micromamba 2>/dev/null)"; do
    if [ -n "$c" ] && [ "$c" != "$MM" ] && [ -x "$c" ] && [ "$("$c" --version 2>/dev/null)" = "$MM_VERSION" ]; then FOUND="$c"; break; fi
  done
  if [ -n "$FOUND" ] && cp -f "$FOUND" "$MM" && chmod +x "$MM"; then
    echo "[teagle] reused existing micromamba from $FOUND"
  elif command -v python3 >/dev/null 2>&1 && mm_py3; then
    echo "[teagle] micromamba installed via python3"
  elif command -v bzip2 >/dev/null 2>&1 && mm_fetch; then
    echo "[teagle] micromamba installed via verified archive"
  elif sudo -n true 2>/dev/null && sudo apt-get update && sudo apt-get install -y python3 ca-certificates && mm_py3; then
    echo "[teagle] micromamba installed via apt+python3"
  else
    fail "micromamba: no python3, no curl/wget+bzip2, and sudo needs a password. Open a WSL terminal, run:  sudo apt-get update && sudo apt-get install -y python3 ca-certificates  then click Repair again"
  fi
  [ -x "$MM" ] && [ "$("$MM" --version 2>/dev/null)" = "$MM_VERSION" ] || fail "micromamba $MM_VERSION missing after install"
fi
echo "[teagle] STEP micromamba OK"
'''.replace("__MM_VERSION__", _MICROMAMBA_VERSION).replace("__MM_URL__", _MICROMAMBA_URL).replace("__MM_SHA256__", _MICROMAMBA_SHA256),
    "repeatmasker": r'''
echo "[teagle] STEP repeatmasker START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if "$MM" run -n te RepeatMasker -v 2>&1 | grep -q "RepeatMasker version __RM_VERSION__$"; then echo "[teagle] RepeatMasker __RM_VERSION__ already present"; else
  # never `env remove` here: the Dfam .h5 libraries (multi-GB) live INSIDE this env prefix.
  # mm_create / mm_install self-recover from a corrupted repodata shard (clean-index + retry once).
  mm_create || fail "create te env"
  # Bioconda supplies the dependency stack; the verified source supplies the corrected runtime.
  mm_install repeatmasker=__RM_BOOTSTRAP__ || fail "install RepeatMasker dependency stack"
  "$MM" run -n te python - "$ENV" '__RM_VERSION__' '__RM_URL__' '__RM_SHA__' <<'PY' || fail "RepeatMasker source setup failed; preserved files can be repaired by retrying this component"
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

env = Path(sys.argv[1])
version, url, expected = sys.argv[2:]
destination = env / "share" / ("teagle-repeatmasker-" + version)
marker = destination / ".teagle-source-owned.json"
if destination.exists() and not marker.is_file():
    raise SystemExit("RepeatMasker source directory is not TEagle-owned; it was left untouched")
ownership = {"version": version, "source_url": url, "source_sha256": expected}
if marker.is_file() and json.loads(marker.read_text(encoding="utf-8")) != ownership:
    raise SystemExit("RepeatMasker ownership marker is unrecognized; the directory was left untouched")
cache = env / ".teagle-source-downloads"
cache.mkdir(exist_ok=True)
archive = cache / ("RepeatMasker-" + version + ".tar.gz")
if not archive.exists():
    with tempfile.NamedTemporaryFile(dir=cache, delete=False) as stream:
        partial = Path(stream.name)
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                shutil.copyfileobj(response, stream)
        except Exception as exc:
            raise SystemExit("RepeatMasker source download failed; retry this component: " + type(exc).__name__)
    if hashlib.sha256(partial.read_bytes()).hexdigest() != expected:
        raise SystemExit("RepeatMasker source checksum mismatch; the untrusted download was retained for inspection")
    partial.replace(archive)
if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
    raise SystemExit("Cached RepeatMasker source checksum mismatch; inspect the cache before retrying")
stage = Path(tempfile.mkdtemp(prefix=".teagle-repeatmasker-", dir=env / "share"))
with tarfile.open(archive) as source:
    source.extractall(stage, filter="data")
candidate = stage / ("RepeatMasker-" + version)
if not (candidate / "RepeatMasker").is_file() or not (candidate / "configure").is_file():
    raise SystemExit("RepeatMasker source archive layout is invalid; staging was preserved")
famdb = sorted(env.glob("share/famdb-*/famdb.py"))
if not famdb:
    raise SystemExit("FamDB tool is missing; repair the dependency stack")
previous = None
if destination.exists():
    backup = Path(tempfile.mkdtemp(prefix=".teagle-repeatmasker-backup-", dir=env / "share"))
    previous = backup / destination.name
    destination.rename(previous)
candidate.rename(destination)
stage.rmdir()
(destination / ".teagle-source-owned.json").write_text(json.dumps(ownership) + "\n", encoding="utf-8")
command = [str(env / "bin/perl"), "configure", "-perlbin", str(env / "bin/perl"), "-trf_prgm", str(env / "bin/trf"), "-rmblast_dir", str(env / "bin"), "-hmmer_dir", str(env / "bin"), "-famdb_dir", str(famdb[0].parent), "-default_search_engine", "rmblast"]
configured = subprocess.run(command, cwd=destination, stdin=subprocess.DEVNULL)
if configured.returncode:
    if previous is not None:
        failed = Path(tempfile.mkdtemp(prefix=".teagle-repeatmasker-failed-", dir=env / "share"))
        destination.rename(failed / destination.name)
        previous.rename(destination)
    raise SystemExit("RepeatMasker configuration failed; source evidence was preserved")
binary = env / "bin/RepeatMasker"
if binary.exists() or binary.is_symlink():
    backup = Path(tempfile.mkdtemp(prefix="RepeatMasker-launcher-", dir=cache)) / "RepeatMasker"
    shutil.copy2(binary, backup, follow_symlinks=False)
temporary = binary.with_name(".RepeatMasker-source-link")
if temporary.exists() or temporary.is_symlink():
    raise SystemExit("RepeatMasker staging link already exists; inspect it before retrying")
temporary.symlink_to(destination / "RepeatMasker")
os.replace(temporary, binary)
PY
  "$MM" run -n te RepeatMasker -v 2>&1 | grep -q "RepeatMasker version __RM_VERSION__$" || fail "RepeatMasker __RM_VERSION__ not runnable after install"
fi
echo "[teagle] STEP repeatmasker OK"
'''.replace("__RM_VERSION__", _REPEATMASKER_VERSION).replace("__RM_BOOTSTRAP__", _CONDA_PINS["repeatmasker"]).replace("__RM_URL__", _REPEATMASKER_SOURCE_URL).replace("__RM_SHA__", _REPEATMASKER_SOURCE_SHA256),
    "minimap2": r'''
echo "[teagle] STEP minimap2 START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/minimap2" ] && [ "$("$ENV/bin/minimap2" --version 2>/dev/null | cut -d- -f1)" = "__MM2_VERSION__" ]; then echo "[teagle] minimap2 __MM2_VERSION__ already present"; else
  mm_create || fail "create te env"
  mm_install minimap2=__MM2_VERSION__ || fail "install minimap2 __MM2_VERSION__"
  [ -x "$ENV/bin/minimap2" ] && [ "$("$ENV/bin/minimap2" --version 2>/dev/null | cut -d- -f1)" = "__MM2_VERSION__" ] || fail "minimap2 __MM2_VERSION__ missing after install"
fi
echo "[teagle] STEP minimap2 OK"
'''.replace("__MM2_VERSION__", _CONDA_PINS["minimap2"]),
    "miniprot": r'''
echo "[teagle] STEP miniprot START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if [ -x "$ENV/bin/miniprot" ] && [ "$("$ENV/bin/miniprot" --version 2>&1 | head -1 | cut -d- -f1)" = "__MP_VERSION__" ]; then echo "[teagle] miniprot __MP_VERSION__ already present"; else
  mm_create || fail "create te env"
  mm_install miniprot=__MP_VERSION__ || fail "install miniprot __MP_VERSION__"
  [ -x "$ENV/bin/miniprot" ] && [ "$("$ENV/bin/miniprot" --version 2>&1 | head -1 | cut -d- -f1)" = "__MP_VERSION__" ] || fail "miniprot __MP_VERSION__ missing after install"
fi
echo "[teagle] STEP miniprot OK"
'''.replace("__MP_VERSION__", _CONDA_PINS["miniprot"]),
    "genomescan": r'''
echo "[teagle] STEP genomescan START"
[ -x "$MM" ] || fail "micromamba required first (repair micromamba)"
if compgen -G "$ENV/conda-meta/ispcr-__ISPCR_VERSION__-*.json" >/dev/null && compgen -G "$ENV/conda-meta/ncbi-datasets-cli-__DATASETS_VERSION__-*.json" >/dev/null && [ -x "$ENV/bin/isPcr" ] && [ -x "$ENV/bin/datasets" ]; then echo "[teagle] isPcr __ISPCR_VERSION__ + NCBI Datasets __DATASETS_VERSION__ already present"; else
  mm_create || fail "create te env"
  mm_install ispcr=__ISPCR_VERSION__ ncbi-datasets-cli=__DATASETS_VERSION__ || fail "install pinned ispcr + ncbi-datasets-cli"
  [ -x "$ENV/bin/isPcr" ] && compgen -G "$ENV/conda-meta/ispcr-__ISPCR_VERSION__-*.json" >/dev/null || fail "isPcr __ISPCR_VERSION__ missing after install"
  [ -x "$ENV/bin/datasets" ] && compgen -G "$ENV/conda-meta/ncbi-datasets-cli-__DATASETS_VERSION__-*.json" >/dev/null || fail "NCBI Datasets __DATASETS_VERSION__ missing after install"
fi
# faToTwoBit is best-effort: it makes cached genomes compact + fast to load, but scans still work on plain
# FASTA if it is unavailable, so a missing package must NOT fail the step.
[ -x "$ENV/bin/faToTwoBit" ] && compgen -G "$ENV/conda-meta/ucsc-fatotwobit-__FATOTWOBIT_VERSION__-*.json" >/dev/null || "$MM" install --force-reinstall -y -n te -c bioconda -c conda-forge ucsc-fatotwobit=__FATOTWOBIT_VERSION__ >/dev/null 2>&1 \
  || echo "[teagle] note: faToTwoBit unavailable (genomes cached as FASTA — larger, still functional)"
echo "[teagle] STEP genomescan OK"
'''.replace("__ISPCR_VERSION__", _CONDA_PINS["ispcr"]).replace("__DATASETS_VERSION__", _CONDA_PINS["ncbi-datasets-cli"]).replace("__FATOTWOBIT_VERSION__", _CONDA_PINS["ucsc-fatotwobit"]),
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
vrna_install(){
  if [ -d "$VRNAENV/conda-meta" ]; then
    "$MM" install --force-reinstall -y -n teagle-vrna -c conda-forge -c bioconda viennarna=__VRNA_VERSION__ python=3.12
  else
    "$MM" create -y -n teagle-vrna -c conda-forge -c bioconda viennarna=__VRNA_VERSION__ python=3.12
  fi
}
if "$MM" run -n teagle-vrna python -c "import RNA,sys;sys.exit(0 if RNA.__version__ == '__VRNA_VERSION__' else 1)" >/dev/null 2>&1; then
  echo "[teagle] ViennaRNA __VRNA_VERSION__ already present"
elif vrna_install || { mm_reset_cache; vrna_install; }; then
  "$MM" run -n teagle-vrna python -c "import RNA,sys;sys.exit(0 if RNA.__version__ == '__VRNA_VERSION__' else 1)" >/dev/null 2>&1 \
    || echo "[teagle] note: ViennaRNA not importable after install (optional primer-QC engine skipped — Primer3 still reports every structure)"
else
  echo "[teagle] note: could not create the teagle-vrna env (optional primer-QC engine skipped — Primer3 still reports every structure)"
fi
"$MM" run -n teagle-vrna python -c "import RNA; print('[teagle] ViennaRNA ' + RNA.__version__)" 2>/dev/null || true
echo "[teagle] STEP viennarna OK"
'''.replace("__VRNA_VERSION__", _CONDA_PINS["viennarna"]),
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

_DEFAULT_DFAM_DOWNLOAD_B = sum(_DFAM_ARCHIVE_B[k] for k in _DFAM_FILES if k in _ALL_STEPS)
_DEFAULT_DFAM_INSTALLED_B = sum(_DFAM_UNPACKED_B[k] for k in _DFAM_FILES if k in _ALL_STEPS)
_DEFAULT_INSTALL_COST = (
    f"The default Dfam data download is {_DEFAULT_DFAM_DOWNLOAD_B / 1048576:.0f} MiB and occupies "
    f"{_gib(_DEFAULT_DFAM_INSTALLED_B)} after extraction. Pinned package downloads, solver time, and "
    "their additional disk use vary by connection and the existing WSL environment; the installer checks "
    "Dfam free space before each transfer."
)
_WSLCONFIG_GUIDE = (
    "WSL2 shares Windows resources. To set global limits, edit %UserProfile%\\.wslconfig under [wsl2] "
    "(for example memory=8GB and processors=4), then run wsl --shutdown only after installs and analyses "
    "have stopped."
)

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
     "Optional. Uncurated root partition (0.3 MiB). How much a curated-only search finds depends on the "
     "lineage: " + curated_coverage_sentence()),
    ("dfam_unc_euk",  "Dfam uncurated · Eukaryota (optional)", True,
     f"Optional, 3.9 GiB download that unpacks to {_gib(_DFAM_UNPACKED_B['dfam_unc_euk'])} — leave about "
     f"{_gib(_DFAM_UNPACKED_B['dfam_unc_euk'] + 5 * 1073741824)} free in WSL, since the archive and the "
     "library it unpacks to are both on disk for a while. Adds the uncurated eukaryote families. Many "
     "lineage-specific families may be unavailable in curated-only mode; consult the live lineage-specific "
     "count before deciding whether to add this partition. Installing it is only half of what is "
     "needed: RepeatMasker searches curated families ONLY unless it is asked for both, so choose the "
     "'include uncurated families' library option to actually use what this downloads. Measured on baker's "
     f"yeast, curated-only searches {_DFAM_YEAST_CURATED} of {_DFAM_YEAST_TOTAL} family models and finds "
     f"no transposable element at all, while including uncurated adds {_DFAM_YEAST_ADDITIONAL} models and "
      "recovers the yeast Ty families. Expect roughly 40-60 minutes: Dfam serves "
     "this file at about 2 MB/s per client, and measured here, splitting it across parallel connections "
     "made it SLOWER, not faster. The download resumes where it stopped if it is interrupted, so a failed "
     "attempt does not start over. There is no smaller alternative: RepeatMasker's RMBLAST engine reads the "
     "consensus component, which Dfam 4.0 ships as only a root partition and this one."),
]


# Components installed independently: each step runs in its own subshell, so one failure ends only that
# step. A step whose prerequisite did not complete is skipped, and the run ends with a single terminal line.
_STEP_NEEDS = {"repeatmasker": ("micromamba",), "minimap2": ("micromamba",), "miniprot": ("micromamba",),
               "genomescan": ("micromamba",), "viennarna": ("micromamba",), "famdb_conf": ("repeatmasker",)}

_RUN_STEP = r'''
STEPS_OK=""; STEPS_FAILED=""; STEPS_SKIPPED=""
bg_cleanup(){ [ -n "$BG" ] && kill $BG 2>/dev/null; return 0; }
run_step(){
  local key=$1 needs=$2 d rc
  for d in $needs; do
    case " $STEPS_FAILED $STEPS_SKIPPED " in
      *" $d "*) echo "[teagle] STEP $key SKIPPED (needs $d, which did not complete)"; STEPS_SKIPPED="$STEPS_SKIPPED $key"; return 0;;
    esac
  done
  (
    fail(){ echo "[teagle] STEP $key FAILED: $1"; exit 1; }
    trap bg_cleanup EXIT
    trap 'bg_cleanup; exit 143' TERM HUP INT
    "step_$key"
  )
  rc=$?
  if [ "$rc" -eq 0 ]; then STEPS_OK="$STEPS_OK $key"
  else echo "[teagle] STEP $key ended with exit code $rc"; STEPS_FAILED="$STEPS_FAILED $key"; fi
  return 0
}
'''

_RUN_SUMMARY = r'''
if [ -n "$STEPS_FAILED$STEPS_SKIPPED" ]; then
  echo "[teagle] SUMMARY installed:${STEPS_OK:- none}; failed:${STEPS_FAILED:- none}; skipped:${STEPS_SKIPPED:- none}"
  echo "[teagle] FAILED: not installed:$STEPS_FAILED$STEPS_SKIPPED (the other components were installed; see the lines above)"
  exit 1
fi
echo "[teagle] DONE $(date -u +%FT%TZ)"
'''


def _build_script(keys) -> str:
    keys = list(keys)
    body = "".join(
        f"\nstep_{k}(){{\n{_STEP[k]}\n}}\nrun_step {k} '{' '.join(d for d in _STEP_NEEDS.get(k, ()) if d in keys)}'\n"
        for k in keys)
    script = _PRELUDE + _RUN_STEP + body + _RUN_SUMMARY
    host_free = _windows_backing_free_bytes()
    script = script.replace("__HOST_FREE_B__", str(host_free if host_free is not None else -1))
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
        _wsl_script('bash "$HOME/teagle_wsl_install.sh"\n', timeout=7200)
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
        payload = base64.b64encode(script.encode("utf-8")).decode("ascii")
        rc, _, err = _wsl_script(
            f"printf '%s' '{payload}' | base64 -d > \"$HOME/teagle_wsl_install.sh\"\n", timeout=30)
        if rc != 0:
            return {"started": False, "error": f"could not stage install script: {err.strip()[:120]}"}
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
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"; VRNAENV="$HOME/micromamba/envs/teagle-vrna"; FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
echo "micromamba=$([ -x "$MM" ] && "$MM" --version 2>/dev/null || echo 0)"
rmv=$("$MM" run -n te RepeatMasker -v 2>/dev/null | grep -oiE 'version [0-9][0-9.]*' | head -1 | awk '{print $2}'); echo "repeatmasker=${rmv:-0}"
mmv=$([ -x "$ENV/bin/minimap2" ] && "$ENV/bin/minimap2" --version 2>/dev/null); echo "minimap2=${mmv:-0}"
echo "genomescan=$([ -x "$ENV/bin/isPcr" ] && [ -x "$ENV/bin/datasets" ] && echo 1 || echo 0)"
echo "ispcr_pin=$(compgen -G "$ENV/conda-meta/ispcr-__ISPCR_VERSION__-*.json" >/dev/null && echo 1 || echo 0)"
echo "datasets_pin=$(compgen -G "$ENV/conda-meta/ncbi-datasets-cli-__DATASETS_VERSION__-*.json" >/dev/null && echo 1 || echo 0)"
if [ -x "$ENV/bin/faToTwoBit" ]; then
  compgen -G "$ENV/conda-meta/ucsc-fatotwobit-__FATOTWOBIT_VERSION__-*.json" >/dev/null && echo "fatotwobit=__FATOTWOBIT_VERSION__" || echo "fatotwobit=version drift"
else
  compgen -G "$ENV/conda-meta/ucsc-fatotwobit-__FATOTWOBIT_VERSION__-*.json" >/dev/null && echo "fatotwobit=broken (pinned metadata, executable missing)" || echo "fatotwobit=absent"
fi
if [ ! -d "$VRNAENV" ]; then
  echo "viennarna=absent"
elif vrv=$("$MM" run -n teagle-vrna python -c "import RNA;print(RNA.__version__)" 2>/dev/null) && [ -n "$vrv" ]; then
  echo "viennarna=$vrv"
else
  echo "viennarna=broken/unusable"
fi
echo "dfam_root=$([ "$(stat -c %s "$FAMDIR/dfam40.0.h5" 2>/dev/null || echo 0)" = "__DFAM_ROOT_SIZE__" ] && echo 1 || echo 0)"
echo "dfam_curated=$([ "$(stat -c %s "$FAMDIR/dfam40.curated.consensus.0.h5" 2>/dev/null || echo 0)" = "__DFAM_CURATED_SIZE__" ] && echo 1 || echo 0)"
echo "dfam_unc_root=$([ -f "$FAMDIR/dfam40.uncurated.consensus.0.h5" ] && echo 1 || echo 0)"
echo "dfam_unc_euk=$([ -f "$FAMDIR/dfam40.uncurated.consensus.1.h5" ] && echo 1 || echo 0)"
echo "famdb_conf=$( { ls "$ENV"/share/famdb-*/famdb.conf >/dev/null 2>&1 || [ -f "$ENV/share/RepeatMasker/famdb.conf" ]; } && echo 1 || echo 0)"
echo "famdb_loads=$("$MM" run -n te famdb.py info >/dev/null 2>&1 && echo 1 || echo 0)"
if [ -d "$HOME/.teagle_install.lock" ] && [ -f "$HOME/.teagle_install.lock/pid" ] && kill -0 "$(cat "$HOME/.teagle_install.lock/pid" 2>/dev/null)" 2>/dev/null; then echo "installing=1"; else echo "installing=0"; fi
echo "disk_free_bytes=$(df -B1 --output=avail "$HOME" 2>/dev/null | tail -1 | tr -dc '0-9')"
'''.replace("__ISPCR_VERSION__", _CONDA_PINS["ispcr"]).replace("__DATASETS_VERSION__", _CONDA_PINS["ncbi-datasets-cli"]).replace("__FATOTWOBIT_VERSION__", _CONDA_PINS["ucsc-fatotwobit"]).replace("__DFAM_ROOT_SIZE__", str(_DFAM_UNPACKED_B["dfam_root"])).replace("__DFAM_CURATED_SIZE__", str(_DFAM_UNPACKED_B["dfam_curated"]))


def components_status() -> dict:
    """Per-component state for the install dialog: WSL2, micromamba, RepeatMasker, minimap2,
    Dfam root/curated, FamDB config. One WSL round-trip; each component idempotently repairable."""
    av = available()
    comp = {c[0]: {"key": c[0], "name": c[1], "repairable": c[2], "desc": c[3], "ok": False, "detail": "—"}
            for c in _COMP_META}
    comp["wsl2"]["ok"] = bool(av["wsl2"])
    comp["wsl2"]["detail"] = av.get("distro") or av.get("error") or "not installed"
    if not av["wsl2"]:
        # Diagnose before advising: a cold start, a restarting VM, a distro without bash and a registration
        # failure need different actions, and only the last is ever answered by `wsl --unregister`.
        d, reason = av.get("distro"), av.get("reason")
        log_state = _wsl2_log_state()
        reboot = log_state == "done"
        win = os.name == "nt"
        installable = win
        if reboot:
            comp["wsl2"]["detail"] = "installed — restart Windows to finish setup, then reopen this installer"
            comp["wsl2"]["guide"] = ("Nothing is lost by restarting. If Ubuntu then asks for a Linux user name, "
                                     "create one.")
        elif reason in ("timeout", "starting"):
            comp["wsl2"]["detail"] = f"WSL is still starting ('{d}' did not answer yet)"
            comp["wsl2"]["guide"] = ("Wait a minute, then click Refresh. The first start after boot or after "
                                     "`wsl --shutdown` can take longer than this check allows.")
            installable = False
        elif reason == "no_bash":
            comp["wsl2"]["detail"] = f"'{d}' is not a usable Linux distribution (it has no bash)"
            comp["wsl2"]["guide"] = ("Click Install WSL to add Ubuntu. If Ubuntu is already installed, make it "
                                     "the default with  wsl --set-default Ubuntu  and click Refresh.")
        elif reason == "no_suitable_distro":
            comp["wsl2"]["detail"] = av.get("error") or "no usable Linux distribution"
            comp["wsl2"]["guide"] = ("Container-runtime distributions cannot host the backend. Click Install WSL "
                                     "to add Ubuntu, then click Refresh.")
        elif reason == "start_failed":
            comp["wsl2"]["detail"] = f"'{d}' is registered but could not be started"
            comp["wsl2"]["guide"] = (f"Click Refresh first. If it still fails, in an Administrator PowerShell run:  "
                                     f"wsl --unregister {d}  then click Install WSL to reinstall. This permanently "
                                     "deletes everything in that distribution: its home directory, the Dfam libraries "
                                     "and the cached genomes (a restart may be required).")
        elif d:
            comp["wsl2"]["detail"] = f"'{d}' did not respond ({av.get('error') or 'no error text'})"
            comp["wsl2"]["guide"] = (f"Click Refresh. If it persists, run  wsl -d {d} -- echo ok  in a Windows "
                                     "terminal to see the error.")
            installable = False
        else:
            comp["wsl2"]["guide"] = ("Click Install WSL to install WSL2 + Ubuntu (needs Administrator; "
                                     "a Windows restart may be required before first use).")
        comp["wsl2"]["installable"] = installable   # the dialog shows an in-app Install WSL button on Windows
        return {"wsl2": False, "installing": False, "wsl2_installing": log_state == "running",
                "wsl2_install_stale": log_state == "stale", "reboot_pending": reboot, "ready": False,
                "install_cost": _DEFAULT_INSTALL_COST, "resource_guide": _WSLCONFIG_GUIDE,
                "components": [comp[c[0]] for c in _COMP_META]}
    try:
        rc, out, err = _wsl_script(_STATUS_PROBE, timeout=90)
    except Exception as e:
        comp["wsl2"]["detail"] = f"probe failed: {str(e)[:80]}"
        return {"wsl2": True, "installing": False, "ready": False,
                "install_cost": _DEFAULT_INSTALL_COST, "resource_guide": _WSLCONFIG_GUIDE,
                "components": [comp[c[0]] for c in _COMP_META], "error": str(e)[:120]}
    kv = {}
    for line in out.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()

    def present(key, ok, detail):
        comp[key]["ok"] = ok
        comp[key]["detail"] = detail
    mmgr = kv.get("micromamba", "0")
    present("micromamba", mmgr == _MICROMAMBA_VERSION,
            (f"v{mmgr}" if mmgr not in ("0", "") else "missing"))
    rm = kv.get("repeatmasker", "0")
    present("repeatmasker", rm == _REPEATMASKER_VERSION, (f"v{rm}" if rm not in ("0", "") else "missing"))
    mm = kv.get("minimap2", "0")
    mm_base = mm.split("-", 1)[0]
    present("minimap2", mm_base == _CONDA_PINS["minimap2"], (mm if mm not in ("0", "") else "missing"))
    scan_ok = kv.get("genomescan") == "1" and kv.get("ispcr_pin") == "1" and kv.get("datasets_pin") == "1"
    fa = kv.get("fatotwobit", "absent")
    fa_detail = ("faToTwoBit not installed (optional; genomes remain FASTA)" if fa == "absent" else
                 (f"faToTwoBit v{fa} installed" if fa == _CONDA_PINS["ucsc-fatotwobit"] else
                  f"faToTwoBit {fa} (optional; compact-cache conversion unavailable)"))
    scan_detail = (("installed at pinned isPcr/Datasets versions" if scan_ok else "isPcr or NCBI Datasets missing/version drift")
                   + "; " + fa_detail)
    present("genomescan", scan_ok, scan_detail)
    vr = kv.get("viennarna", "absent")
    # Optional absence is normal, but a present environment that cannot import RNA is a repairable fault.
    vr_detail = (f"v{vr}" if vr == _CONDA_PINS["viennarna"] else
                 ("not installed (optional)" if vr == "absent" else
                  ("broken/unusable (repair optional component)" if vr == "broken/unusable" else
                   f"v{vr} (version drift; repair optional component)")))
    present("viennarna", vr == _CONDA_PINS["viennarna"], vr_detail)
    present("dfam_root", kv.get("dfam_root") == "1", "present" if kv.get("dfam_root") == "1" else "missing")
    present("dfam_curated", kv.get("dfam_curated") == "1", "present" if kv.get("dfam_curated") == "1" else "missing")
    for _k in ("dfam_unc_root", "dfam_unc_euk"):      # optional: absent is a normal state, not "missing"
        present(_k, kv.get(_k) == "1", "present" if kv.get(_k) == "1" else "not installed (optional)")
    famdb_ok = kv.get("famdb_conf") == "1" and kv.get("famdb_loads") == "1"
    famdb_detail = ("configured; library opens successfully" if famdb_ok else
                    ("configuration present, but FamDB cannot open the library" if kv.get("famdb_conf") == "1"
                     else "missing configuration"))
    present("famdb_conf", famdb_ok, famdb_detail)
    ready = all(comp[k]["ok"] for k in ("repeatmasker", "dfam_root", "dfam_curated", "famdb_conf"))
    wsl_free_bytes = int(kv.get("disk_free_bytes") or 0)
    if not wsl_free_bytes and kv.get("disk_free_gb"):
        wsl_free_bytes = int(float(kv["disk_free_gb"]) * 1073741824)
    capacity = _disk_capacity(wsl_free_bytes)
    return {"wsl2": True, "installing": kv.get("installing") == "1", "ready": ready,
            "install_cost": _DEFAULT_INSTALL_COST, "resource_guide": _WSLCONFIG_GUIDE,
            "components": [comp[c[0]] for c in _COMP_META], **capacity}


# which Dfam partitions the default install ships, and which are installed only on request — read off
# the step list rather than restated, so adding a partition cannot leave the integrity probe behind
_DFAM_REQUIRED = [k for k in _DFAM_FILES if k in _ALL_STEPS]
_DFAM_OPTIONAL = [k for k in _DFAM_FILES if k not in _ALL_STEPS]

_INTEGRITY_PROBE = r'''cd "$HOME" 2>/dev/null
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"; ENV="$HOME/micromamba/envs/te"; VRNAENV="$HOME/micromamba/envs/teagle-vrna"; FAMDIR="$ENV/share/RepeatMasker/Libraries/famdb"
echo "=RM="; "$MM" run -n te RepeatMasker -v 2>&1 | head -1
echo "=MM="; "$ENV/bin/minimap2" --version 2>&1 | head -1
echo "=FAMDB="; famdb_out=$("$MM" run -n te famdb.py info 2>&1); echo "FAMDB_RC=$?"; printf '%s\n' "$famdb_out" | head -60
echo "=FILES="; for spec in ''' + " ".join(f"{_DFAM_FILES[k][0]}:{_DFAM_UNPACKED_B[k]}" for k in _DFAM_REQUIRED) + r'''; do f=${spec%%:*}; expected=${spec##*:}; got=$(stat -c %s "$FAMDIR/$f" 2>/dev/null || echo 0); if [ "$got" = "$expected" ]; then echo "present $f $got"; elif [ -f "$FAMDIR/$f" ]; then echo "INVALID $f size=$got expected=$expected"; else echo "MISSING $f"; fi; done
echo "=OPTFILES="; for f in ''' + " ".join(_DFAM_FILES[k][0] for k in _DFAM_OPTIONAL) + r'''; do if [ -f "$FAMDIR/$f" ]; then echo "present $f"; else echo "absent $f"; fi; done
echo "=SCAN="
# The pinned UCSC isPcr/faToTwoBit builds deliberately return 255 for a no-argument usage request.
# Their exact first-line banners plus exact conda metadata distinguish that healthy behavior from a
# missing ELF interpreter or a wrong binary without coupling this backend check to assay fixtures.
if [ -x "$ENV/bin/isPcr" ]; then
  ispcr_out=$("$ENV/bin/isPcr" 2>&1); ispcr_rc=$?; ispcr_first=$(printf '%s\n' "$ispcr_out" | head -1); printf '%s\n' "$ispcr_first"
  case "$ispcr_first" in "isPcr - Standalone v 33x2"*) [ "$ispcr_rc" -eq 255 ] && echo "isPcr signature 33x2 OK" || echo "isPcr BROKEN rc=$ispcr_rc";; *) echo "isPcr BROKEN rc=$ispcr_rc wrong signature";; esac
else echo "isPcr MISSING"; fi
if [ -x "$ENV/bin/datasets" ]; then
  datasets_out=$("$ENV/bin/datasets" version 2>&1); datasets_rc=$?; datasets_first=$(printf '%s\n' "$datasets_out" | head -1); printf '%s\n' "$datasets_first"
  if [ "$datasets_rc" -eq 0 ] && printf '%s\n' "$datasets_out" | grep -Eq '(^|[^0-9.])__DATASETS_VERSION__([^0-9.]|$)'; then echo "datasets version __DATASETS_VERSION__ OK"; else echo "datasets BROKEN rc=$datasets_rc wrong version/signature"; fi
else echo "datasets MISSING"; fi
compgen -G "$ENV/conda-meta/ispcr-__ISPCR_VERSION__-*.json" >/dev/null && echo "ispcr package __ISPCR_VERSION__" || echo "ispcr package MISMATCH"
compgen -G "$ENV/conda-meta/ncbi-datasets-cli-__DATASETS_VERSION__-*.json" >/dev/null && echo "datasets package __DATASETS_VERSION__" || echo "datasets package MISMATCH"
echo "=FATOTWOBIT="
if [ -x "$ENV/bin/faToTwoBit" ]; then
  fa_out=$("$ENV/bin/faToTwoBit" 2>&1); fa_rc=$?; fa_first=$(printf '%s\n' "$fa_out" | head -1); printf '%s\n' "$fa_first"
  case "$fa_first" in "faToTwoBit - Convert DNA"*) [ "$fa_rc" -eq 255 ] && echo "faToTwoBit signature OK" || echo "faToTwoBit BROKEN rc=$fa_rc";; *) echo "faToTwoBit BROKEN rc=$fa_rc wrong signature";; esac
  compgen -G "$ENV/conda-meta/ucsc-fatotwobit-__FATOTWOBIT_VERSION__-*.json" >/dev/null && echo "present __FATOTWOBIT_VERSION__" || echo "present MISMATCH"
else
  compgen -G "$ENV/conda-meta/ucsc-fatotwobit-__FATOTWOBIT_VERSION__-*.json" >/dev/null && echo "MISSING executable package __FATOTWOBIT_VERSION__" || echo "absent"
fi
echo "=VRNA="
if [ -d "$VRNAENV" ]; then
  v=$("$MM" run -n teagle-vrna python -c "import RNA;print(RNA.__version__)" 2>/dev/null) || v="UNUSABLE"
  echo "present ${v:-UNUSABLE}"
else
  echo "absent"
fi
'''.replace("__ISPCR_VERSION__", _CONDA_PINS["ispcr"]).replace("__DATASETS_VERSION__", _CONDA_PINS["ncbi-datasets-cli"]).replace("__FATOTWOBIT_VERSION__", _CONDA_PINS["ucsc-fatotwobit"])


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
    fam_lines = sec.get("FAMDB", [])
    fam_rc = next((l.strip().partition("=")[2] for l in fam_lines if l.startswith("FAMDB_RC=")), None)
    fam_body = "\n".join(l for l in fam_lines if not l.startswith("FAMDB_RC="))
    fam_ok = fam_rc == "0" and parse_famdb_info(fam_body) is not None     # exit status AND the header a real library prints
    fam_txt = " ".join(l.strip() for l in fam_body.splitlines() if l.strip())
    files = [l for l in sec.get("FILES", []) if l.strip()]
    # optional partitions are REPORTED, never graded: absent is a normal state, and a user who has just
    # spent 40 minutes downloading one should be able to see it named here rather than infer it
    opt = [l.split() for l in sec.get("OPTFILES", []) if l.strip()]
    opt_txt = ", ".join(f"{f} {'present' if s == 'present' else 'not installed'}" for s, f in opt if f)
    scan_lines = [line.strip() for line in sec.get("SCAN", []) if line.strip()]
    scan_txt = " ".join(scan_lines)
    scan_fault = next((line for line in scan_lines if any(flag in line for flag in ("BROKEN", "MISSING", "MISMATCH"))), None)
    scan_detail = ((scan_fault + (f" | {scan_lines[0]}" if scan_lines and scan_lines[0] != scan_fault else ""))
                   if scan_fault else scan_txt)
    fatotwobit_lines = [line.strip() for line in sec.get("FATOTWOBIT", []) if line.strip()]
    fatotwobit_txt = " ".join(fatotwobit_lines)
    fatotwobit_fault = next((line for line in fatotwobit_lines if any(flag in line for flag in ("BROKEN", "MISSING", "MISMATCH"))), None)
    fatotwobit_detail = ((fatotwobit_fault + (f" | {fatotwobit_lines[0]}" if fatotwobit_lines and fatotwobit_lines[0] != fatotwobit_fault else ""))
                         if fatotwobit_fault else fatotwobit_txt)
    vrna_txt = " ".join(sec.get("VRNA", [])).strip()
    fa_absent = fatotwobit_txt == "absent"
    fa_exact = ("faToTwoBit signature OK" in fatotwobit_txt
                and f"present {_CONDA_PINS['ucsc-fatotwobit']}" in fatotwobit_txt
                and "BROKEN" not in fatotwobit_txt and "MISSING" not in fatotwobit_txt)
    vrna_absent = vrna_txt == "absent"
    vrna_exact = vrna_txt == f"present {_CONDA_PINS['viennarna']}"
    checks = [
        {"name": "RepeatMasker runs", "ok": f"RepeatMasker version {_REPEATMASKER_VERSION}" in rm_txt, "detail": rm_txt.strip()[:80] or "no version reported"},
        {"name": "minimap2 runs", "ok": mm_txt.split("-", 1)[0] == _CONDA_PINS["minimap2"], "detail": mm_txt[:60] or "no version reported"},
        {"name": "FamDB loads", "ok": fam_ok,
         "detail": (f"exit {fam_rc}: " if fam_rc not in (None, "0") else "") + (fam_txt[:90] or "famdb.py info returned nothing")},
        {"name": "Dfam library files present",
         "ok": len(files) >= len(_DFAM_REQUIRED) and all(f.startswith("present") for f in files),
         "detail": ("; ".join(files)[:90] or "missing")
                   + (f" | optional: {opt_txt}" if opt_txt else "")},
        # the whole-genome off-target scan needs isPcr + NCBI datasets; verify them here so a "healthy" deep-check
        # never precedes a scan that fails at first use (same binding-truthfulness class as the fixed unzip gap)
        {"name": "isPcr + NCBI datasets (whole-genome scan)",
         "ok": ("isPcr signature 33x2 OK" in scan_txt
                and f"datasets version {_CONDA_PINS['ncbi-datasets-cli']} OK" in scan_txt
                and "MISSING" not in scan_txt and "BROKEN" not in scan_txt
                and f"ispcr package {_CONDA_PINS['ispcr']}" in scan_txt
                and f"datasets package {_CONDA_PINS['ncbi-datasets-cli']}" in scan_txt),
         "detail": scan_detail[:120] or "not probed"},
        {"name": "faToTwoBit (compact genome cache; optional)",
         "ok": fa_absent or fa_exact,
         "detail": ("not installed (optional)" if fa_absent else fatotwobit_detail[:120] or "not probed")},
        {"name": "ViennaRNA (primer-QC cross-check; optional)",
         "ok": vrna_absent or vrna_exact,
         "detail": ("not installed (optional)" if vrna_absent else vrna_txt[:90] or "not probed")},
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


_WSL2_STALE_AFTER_S = 2 * 3600


def _boot_time():
    try:
        import psutil
        return psutil.boot_time()
    except Exception:
        return None


def _wsl2_log_state() -> str:
    """none / running / stale / done / failed, from the Windows-side install log. A log with no terminal
    marker is only 'running' while it is recent and written since the last boot: a batch killed by a restart
    or abandoned for hours never writes one. A 'done' log from before the last boot is a restart already
    taken, so it no longer means one is pending."""
    try:
        mtime = os.stat(os.path.join(appdirs.user_data_dir(), _WIN_WSL_LOG)).st_mtime
    except OSError:
        return "none"
    log = wsl2_install_log(500)
    if not log:
        return "none"
    boot = _boot_time()
    before_boot = boot is not None and mtime < boot
    if "[teagle] FAILED" in log:
        return "failed"
    if "DONE-WSL" in log:
        return "none" if before_boot else "done"
    return "stale" if (before_boot or time.time() - mtime > _WSL2_STALE_AFTER_S) else "running"


def _wsl2_installing() -> bool:
    """A WSL2 install is in progress = the Windows-side log is recent and has no terminal marker yet."""
    return _wsl2_log_state() == "running"


def _wsl2_elevated_command(log: str, wsl_exe: str = r"%SystemRoot%\System32\wsl.exe"):
    """The elevated command line: install WSL2 + Ubuntu (no interactive launch), tolerant of the pre-reboot
    state, logging everything (incl. a terminal DONE-WSL marker) to the file the dialog polls.

    The command travels as ShellExecute parameters, not as a script staged in a user-writable folder, so
    nothing that could be rewritten between staging and the UAC consent is ever executed with admin rights.
    `wsl --install` is only ever called with --no-launch: a plain one can open the distro's interactive first
    run in a window nobody can see, which stalls the install until the dialog gives up."""
    q = f'"{log}"'
    steps = [
        f"echo [teagle] installing WSL2 + Ubuntu (elevated) >> {q}",
        f'"{wsl_exe}" --install -d Ubuntu --no-launch >> {q} 2>&1',
        "set TEAGLE_RC=!errorlevel!",
        # VirtualMachinePlatform is inactive until the reboot, so this may legitimately fail
        f'"{wsl_exe}" --set-default-version 2 >> {q} 2>&1',
        f"echo [teagle] wsl --status: >> {q}",
        f'"{wsl_exe}" --status >> {q} 2>&1',
        f"echo [teagle] DONE-WSL !TEAGLE_RC! (restart Windows if the distro is not usable yet, then reopen "
        f"this installer) >> {q}",
    ]
    cmd = os.path.join(os.environ.get("SystemRoot", r"C:\Windows"), "System32", "cmd.exe")
    return cmd, '/v:on /d /s /c "' + " & ".join(steps) + '"'


def _shell_execute_runas(exe: str, params: str) -> int:
    import ctypes
    # verb 'runas' -> UAC. Return > 32 = launched; <= 32 = failed (e.g. 1223 = declined).
    return int(ctypes.windll.shell32.ShellExecuteW(None, "runas", exe, params, None, 0))


def _account_is_administrator() -> bool:
    """True when this account is in Administrators (a filtered token still lists the group). A standard user
    reaching the UAC prompt must type another account's credentials, and WSL distributions are registered
    per user, so the distro would land in that other account's profile. Unknown reads as True: no warning."""
    try:
        p = subprocess.run(["whoami", "/groups"], capture_output=True, timeout=15, creationflags=_NO_WINDOW)
        return "S-1-5-32-544" in p.stdout.decode("utf-8", "ignore")
    except Exception:
        return True


_OVER_THE_SHOULDER = ("This Windows account is not an administrator, so the prompt asks for a different account's "
                      "credentials, and the Linux distribution is registered to the account that approves it, "
                      "not to yours. If TEagle still reports no distribution afterwards, ask an administrator to "
                      "enable WSL, then run  wsl --install -d Ubuntu  from your own account.")


def install_wsl2() -> dict:
    """Install WSL2 + Ubuntu from an ELEVATED helper (`wsl --install` requires Administrator).
    Fire-and-forget: an elevated cmd.exe runs `wsl --install -d Ubuntu --no-launch` and logs to a Windows-side
    file the dialog polls. If the UAC prompt is declined or the user is not an admin, a terminal FAILED marker
    is written so the poller never hangs. A Windows restart may be required before the distro is usable."""
    global _wsl2_thread
    if os.name != "nt":
        return {"started": False, "error": "WSL installation is only available on Windows"}
    if _wsl2_thread is not None and _wsl2_thread.is_alive():
        return {"started": False, "error": "WSL install already running"}
    log = os.path.join(appdirs.user_data_dir(), _WIN_WSL_LOG)
    warning = None if _account_is_administrator() else _OVER_THE_SHOULDER
    try:
        with open(log, "w", encoding="utf-8", newline="") as f:
            f.write("[teagle] launching the elevated WSL installer - accept the Windows (UAC) prompt...\n")
            if warning:
                f.write("[teagle] NOTE: " + warning + "\n")
    except Exception as e:
        return {"started": False, "error": f"could not stage the WSL installer: {type(e).__name__}: {e}"}
    exe, params = _wsl2_elevated_command(log)

    def _run():
        try:
            rc = _shell_execute_runas(exe, params)
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
    return {"started": True, "windows_log": True, **({"warning": warning} if warning else {})}


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
                                      f"spelling, or leave the organism blank to use RepeatMasker's default "
                                      f"lineage (Homo sapiens)."}
    return {"ok": True}


_RM_DEFAULT_LINEAGE = "Homo sapiens"        # RepeatMasker 4.2.5 retains this default when no lineage is supplied.


def _backend_not_ready(st: dict, need_dfam: bool = True) -> str:
    """Name the condition that failed, in words, and where to fix it."""
    if st.get("wsl2") is False:
        why = "WSL is not installed or not running — install or start it from the Backend installer"
    elif not st.get("repeatmasker"):
        why = "RepeatMasker is not installed in the WSL backend — install it from the Backend installer"
    elif need_dfam and not st.get("dfam"):
        why = "the Dfam root and curated libraries are not both installed — install them from the Backend installer"
    elif need_dfam and not st.get("dfam_library"):
        why = "FamDB could not open the Dfam library — run Check integrity in the Backend installer"
    else:
        why = "the backend status could not be read — try again"
    return f"WSL annotation backend not ready: {why}."


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
    elif include_uncurated:
        return {"ok": False, "error": "select an organism/lineage before including uncurated Dfam families",
                "uncurated_available": True}
    st = env_status()
    if include_uncurated and not st.get("dfam_uncurated"):
        return {"ok": False, "error": "the optional uncurated Dfam partitions are not installed; install "
                "both uncurated partitions in the backend manager or use curated families only",
                "status": st, "uncurated_available": False}
    if not st["ready"]:
        return {"ok": False, "error": _backend_not_ready(st), "status": st}
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
        # The script's own exit status is RepeatMasker's, and RepeatMasker prints the taxon it actually
        # searched (the line after 'Species/Taxa Search:'), which is read back so the result can say so.
        script = (f'cd /tmp/{rid} || exit 97\n'
                  f'{_MM} run -n te RepeatMasker -pa {int(threads)} {sp} -qq q.fa >rm.log 2>&1; ec=$?\n'
                  'echo "EXIT $ec"\n'
                  "sed -n '/^Species\\/Taxa Search:/{n;p;q}' rm.log 2>/dev/null | sed 's/^/TEAGLE_SEARCHED /'\n"
                  'if [ "$ec" = "0" ]; then cat q.fa.out 2>/dev/null; else tail -8 rm.log; fi\n'
                  'exit "$ec"\n')
        rc, out, err = _wsl_script(script, timeout=timeout)
        m = re.search(r"^EXIT (\d+)", out, re.M)
        if not m:                                             # wsl.exe died (VM shut down, crash): nothing ran to completion
            return {"ok": False, "error": "the WSL session ended before RepeatMasker reported, so nothing was "
                    "searched. Check that WSL is running, then run the annotation again.", "status": st}
        if m.group(1) != "0":
            tail = out.split(f"EXIT {m.group(1)}", 1)[-1].strip()[:250]
            return {"ok": False, "error": f"RepeatMasker exited {m.group(1)}: {tail or 'see backend log'}", "status": st}
        if rc != 0:
            return {"ok": False, "error": f"RepeatMasker reported success but the WSL session returned exit {rc}, "
                    "so the result was discarded; run the annotation again.", "status": st}
        sm = re.search(r"^TEAGLE_SEARCHED\s+(.+?)\s*(?:\[NCBI Taxonomy ID:\s*(\d+)\])?\s*$", out, re.M)
        lineage, taxid = (sm.group(1), sm.group(2)) if sm else (None, None)
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
                # With no organism RepeatMasker searches ITS default lineage (Homo sapiens), not every
                # installed family, so that is what is named here and sealed.
                "species": species or f"{lineage or _RM_DEFAULT_LINEAGE} (RepeatMasker default lineage)",
                "species_requested": species, "species_defaulted": not species,
                "lineage_searched": lineage, "lineage_taxid": taxid}
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
        # FLAG 16 says which way the READ lies on the reference; minimap2's ts:A tag (-ub, splice preset) gives
        # the transcript strand relative to the read. The gene's genomic strand combines the two, so an
        # antisense read (3' ESTs, first-strand cDNA) is not mistaken for a minus-strand gene.
        aln = "-" if (flag & 16) else "+"
        ts = next((t[5:] for t in f[11:] if t.startswith("ts:A:")), None)
        ts = ts if ts in ("+", "-") else None
        strand = aln if ts in (None, "+") else ("+" if aln == "-" else "-")
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
        return {"exons": exons, "introns": introns, "strand": strand, "alignment_strand": aln,
                "transcript_strand": ts, "strand_determined": ts is not None, "ref_start": pos, "ref_end": ref}
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


_CACHE_CORRUPT_MARKER = "CACHE_CORRUPT"
_CACHE_RECOVERY = "Open Manage genomes, remove this cached genome, then download it again."


def _operation_failure(prefix: str, output: str, stderr: str = "", *, partial: bool = False) -> dict:
    """Translate a shell failure while preserving a recoverable cache-integrity diagnosis."""
    failure = next((line for line in output.splitlines() if line.startswith("FAIL")), "")
    detail = failure[4:].strip() if failure else (stderr or "").strip()[:200]
    cache_corrupt = _CACHE_CORRUPT_MARKER in detail
    detail = detail.replace(_CACHE_CORRUPT_MARKER, "").strip(" :-") or "unknown error"
    result = {"ok": False, "error": f"{prefix}: {detail}"}
    if cache_corrupt:
        result["cache_corrupt"] = True
        result["error"] = result["error"].rstrip(".") + ". " + _CACHE_RECOVERY
    if partial:
        result["partial"] = True
    return result


# download the assembly, convert to compact 2bit (best-effort), seal the SOURCE FASTA sha256 (machine-
# independent, unlike the version-dependent 2bit), mark .done for idempotent resume. __ACC__ is a
# validated accession (safe to interpolate); no untrusted input enters the script.
_GENOME_PREP_LOG = "$HOME/teagle_genome_prepare.log"    # milestones the UI tails for a liveness indicator

_PREP_SCRIPT = r'''#!/usr/bin/env bash
set -uo pipefail
cd "$HOME" || exit 1
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
ENV="$HOME/micromamba/envs/te"
HOST_FREE_B=__HOST_FREE_B__
effective_free_b(){
  local wsl_b=${1:-0}
  if [ "$HOST_FREE_B" -ge 0 ] && [ "$HOST_FREE_B" -lt "$wsl_b" ]; then echo "$HOST_FREE_B"; else echo "$wsl_b"; fi
}
ACC="__ACC__"
GDIR="$HOME/teagle_genomes/$ACC"
LOG="$HOME/teagle_genome_prepare.log"
plog(){ echo "[prep] $1"; echo "$1" >> "$LOG"; }
mkdir -p "$GDIR" || { echo "FAIL mkdir"; exit 1; }
# A completion marker is only a hint. Recompute the cached target checksum before every reuse so a
# truncated or modified cache can never be returned as PREPARED. Version-3.8.0 FASTA caches remain
# verifiable because their target is the sealed source FASTA; legacy 2bit caches have no target checksum
# and must be downloaded once more rather than being trusted without evidence.
verify_prepared(){
  TARGET=$(sed -n 's/^target=//p' "$GDIR/meta.txt" | head -1)
  case "$TARGET" in genome.fna|genome.2bit) ;; *) echo "FAIL CACHE_CORRUPT cached genome metadata has an invalid target"; return 1 ;; esac
  WANT=$(sed -n 's/^target_sha256=//p' "$GDIR/meta.txt" | head -1)
  [ -n "$WANT" ] || { [ "$TARGET" = "genome.fna" ] && WANT=$(sed -n 's/^sha256=//p' "$GDIR/meta.txt" | head -1); }
  [ -n "$WANT" ] || { echo "FAIL CACHE_CORRUPT cached genome has no sealed target checksum"; return 1; }
  [ -s "$GDIR/$TARGET" ] || { echo "FAIL CACHE_CORRUPT cached genome target is missing or empty"; return 1; }
  GOT=$(sha256sum "$GDIR/$TARGET" | cut -d' ' -f1)
  [ "$GOT" = "$WANT" ] || { echo "FAIL CACHE_CORRUPT cached genome target failed its sealed checksum"; return 1; }
  echo "PREPARED"
  cat "$GDIR/meta.txt"
}
if [ -f "$GDIR/.done" ] && [ -f "$GDIR/meta.txt" ]; then verify_prepared || exit 1; exit 0; fi
# atomic per-genome lock so two concurrent prepares of the same accession can't race on dl.zip / genome.2bit
# / meta.txt and leave a half-written cache; a lock orphaned by a crash (dead PID) is reaped.
LOCK="$GDIR/.lock"
if [ -d "$LOCK" ] && { [ ! -f "$LOCK/pid" ] || ! kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; }; then rm -rf "$LOCK" 2>/dev/null; fi
if ! mkdir "$LOCK" 2>/dev/null; then echo "FAIL a download for this genome is already running"; exit 1; fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK" 2>/dev/null' EXIT
if [ -f "$GDIR/.done" ] && [ -f "$GDIR/meta.txt" ]; then verify_prepared || exit 1; exit 0; fi   # another prepare finished while we waited
cd "$GDIR"
wsl_avail_b=$(df -B1 --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9')
avail_b=$(effective_free_b "${wsl_avail_b:-0}")
# blanket floor sized for a mammalian genome's peak (extracted FASTA ~3G + 2bit ~0.8G); the zip is freed
# before the conversion to keep the peak below this.
[ "$avail_b" -ge 8589934592 ] || { echo "FAIL insufficient usable disk ($((avail_b/1073741824))G free after applying the smaller WSL/Windows backing-volume limit, need >=8G for a genome)"; exit 1; }
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
FABYTES=$(stat -c %s "$FNA" 2>/dev/null || echo 0)       # the source FASTA size is what annotation must re-fetch and split
plog "building compact isPcr target (2bit)"
TARGET=""
if [ -x "$ENV/bin/faToTwoBit" ] && "$ENV/bin/faToTwoBit" "$FNA" genome.2bit 2>/dev/null; then
  TARGET="genome.2bit"
else
  mv "$FNA" genome.fna 2>/dev/null || cp "$FNA" genome.fna; TARGET="genome.fna"
fi
TBYTES=$(stat -c %s "$TARGET" 2>/dev/null || echo 0)
TSHA=$(sha256sum "$TARGET" | cut -d' ' -f1)
[ -n "$TSHA" ] || { echo "FAIL target checksum failed"; exit 1; }
rm -rf ex
printf 'accession=%s\ntarget=%s\nsha256=%s\ntarget_sha256=%s\nn_seqs=%s\nbytes=%s\nfasta_bytes=%s\n' "$ACC" "$TARGET" "$FSHA" "$TSHA" "$NSEQ" "$TBYTES" "$FABYTES" > meta.txt
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
    try:
        rc, out, _ = _wsl(f'[ -x "{_ENV}/bin/datasets" ] && [ -x "{_ENV}/bin/isPcr" ] && echo yes || echo no', timeout=20)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "WSL did not answer the check for the genome-scan tools in time (it may "
                "still be starting); try again in a minute"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    if "yes" not in out:
        return {"ok": False, "error": "genome-scan tools not installed in the WSL backend (run Install backend)"}
    try:
        host_free = _windows_backing_free_bytes()
        script = (_PREP_SCRIPT.replace("__ACC__", accession)
                  .replace("__HOST_FREE_B__", str(host_free if host_free is not None else -1)))
        rc, out, err = _wsl_script(script, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": (f"genome preparation timed out after {timeout}s; try again, but the "
                                       "NCBI genome ZIP restarts because Datasets does not resume partial ZIPs")}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    meta = _parse_meta(out)
    if "PREPARED" not in out or not meta.get("sha256"):
        return _operation_failure("genome preparation failed", out, err)
    return {"ok": True, "accession": accession, "assembly_name": assembly_name,
            "target": meta.get("target"), "sha256": meta.get("sha256"),
            "n_seqs": int(meta.get("n_seqs", 0) or 0), "bytes": int(meta.get("bytes", 0) or 0),
            "fasta_bytes": int(meta.get("fasta_bytes", 0) or 0)}


def genome_prepare_log(tail: int = 1) -> str:
    """Tail the genome-prepare milestone log for a UI liveness indicator during a long download."""
    try:
        _, out, _ = _wsl(f'tail -n {int(tail)} "{_GENOME_PREP_LOG}" 2>/dev/null || true', timeout=15)
        return out.strip()
    except Exception:
        return ""


# The checksum reads the whole cached target. 30 s timed out on a multi-GB genome read cold from the VHDX;
# 1800 s covers a 25 GB target at 15 MB/s, and the probe before it (available(), 30 s) already filters a hung WSL.
_HASH_PROBE_TIMEOUT_S = 1800


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
    # Verify the cached search target against its sealed target checksum before isPcr reads it. The source
    # FASTA checksum identifies the assembly; target_sha256 separately protects a derived 2bit file. Legacy
    # FASTA caches can use the source checksum, while legacy 2bit caches must be downloaded once more because
    # no trustworthy target baseline exists for them.
    probe = f'''GDIR="{_GENOMES}/{accession}"
if [ ! -f "$GDIR/.done" ] || [ ! -f "$GDIR/meta.txt" ]; then echo NOTPREP; exit 0; fi
TARGET=$(sed -n 's/^target=//p' "$GDIR/meta.txt" | head -1)
case "$TARGET" in genome.fna|genome.2bit) ;; *) echo "CACHE_BAD invalid target metadata"; exit 0 ;; esac
SOURCE_SHA=$(sed -n 's/^sha256=//p' "$GDIR/meta.txt" | head -1)
WANT=$(sed -n 's/^target_sha256=//p' "$GDIR/meta.txt" | head -1)
[ -n "$WANT" ] || {{ [ "$TARGET" = "genome.fna" ] && WANT="$SOURCE_SHA"; }}
[ -n "$SOURCE_SHA" ] && [ -n "$WANT" ] || {{ echo "CACHE_BAD missing sealed checksum"; exit 0; }}
[ -s "$GDIR/$TARGET" ] || {{ echo "CACHE_BAD target missing or empty"; exit 0; }}
GOT=$(sha256sum "$GDIR/$TARGET" | cut -d' ' -f1)
[ "$GOT" = "$WANT" ] || {{ echo "CACHE_BAD checksum mismatch"; exit 0; }}
echo CACHE_OK
cat "$GDIR/meta.txt"
'''
    try:
        rc, out, err = _wsl_script(probe, timeout=_HASH_PROBE_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"checking the cached genome's checksum timed out after {_HASH_PROBE_TIMEOUT_S} s "
                "(WSL may be busy or still starting, or the disk is slow); try again"}
    except Exception as e:
        return {"ok": False, "error": f"could not verify the cached genome: {type(e).__name__}: {e}"}
    if rc != 0:
        return {"ok": False, "error": "could not verify the cached genome: " + err.strip()[:160]}
    if "CACHE_BAD" in out:
        detail = next((line.partition("CACHE_BAD")[2].strip() for line in out.splitlines()
                       if line.startswith("CACHE_BAD")), "integrity check failed")
        return {"ok": False, "error": "cached genome failed its integrity check (" + detail
                + "); remove it in Manage genomes, then download it again", "cache_corrupt": True}
    if "NOTPREP" in out or "accession=" not in out:
        return {"ok": False, "error": "genome not prepared — download it first", "need_prepare": True}
    meta = _parse_meta(out)
    if not re.fullmatch(r"[0-9a-fA-F]{64}", meta.get("sha256", "")):
        return {"ok": False, "error": "cached genome metadata has no valid sealed source checksum; "
                "remove it in Manage genomes, then download it again", "cache_corrupt": True}
    target = meta.get("target")
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
HOST_FREE_B=__HOST_FREE_B__
effective_free_b(){
  local wsl_b=${1:-0}
  if [ "$HOST_FREE_B" -ge 0 ] && [ "$HOST_FREE_B" -lt "$wsl_b" ]; then echo "$HOST_FREE_B"; else echo "$wsl_b"; fi
}
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
WANT=$(sed -n 's/^sha256=//p' "$GDIR/meta.txt" | head -1)
[ -n "$WANT" ] || { echo "FAIL CACHE_CORRUPT prepared genome has no sealed source FASTA checksum"; exit 1; }
if [ ! -s "$FNA" ]; then
  # Check free space BEFORE the download, not after it: the zip and the extracted FASTA coexist (at most twice
  # the FASTA), and the stage-2 check below only runs once that space has already been spent.
  FB=$(sed -n 's/^fasta_bytes=//p' "$GDIR/meta.txt" | head -1)
  # caches prepared before the FASTA size was recorded: a 2bit packs four bases per byte, so the FASTA is at least 4x it
  [ -n "$FB" ] || { BY=$(sed -n 's/^bytes=//p' "$GDIR/meta.txt" | head -1 | tr -dc '0-9'); FB=$(( ${BY:-0} * 4 )); }
  wsl_avail_b=$(df -B1 --output=avail "$GDIR" 2>/dev/null | tail -1 | tr -dc '0-9')
  avail_b=$(effective_free_b "${wsl_avail_b:-0}")
  need_b=$(( 2 * FB + 1073741824 ))
  [ "$avail_b" -ge "$need_b" ] || { echo "FAIL insufficient usable disk to re-fetch the genome FASTA ($((avail_b/1073741824))G free after applying the smaller WSL/Windows backing-volume limit, need >=$((need_b/1073741824))G for the download and extraction)"; exit 1; }
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
  GOT=$(sha256sum "$SRC" | cut -d' ' -f1)
  if [ "$WANT" != "$GOT" ]; then rm -rf aex; echo "FAIL CACHE_CORRUPT fetched FASTA does not match the sealed genome checksum"; exit 1; fi
  mv "$SRC" "$FNA" 2>/dev/null || cp "$SRC" "$FNA"; rm -rf aex
fi
# Do not trust an existing genome.fna merely because it was checked when first downloaded. Rehash it on
# every annotation and fail before reading or reusing any chunk if its bytes no longer match the seal.
GENOME_SHA=$(sha256sum "$FNA" | cut -d' ' -f1)
[ "$GENOME_SHA" = "$WANT" ] || { echo "FAIL CACHE_CORRUPT cached genome FASTA failed its sealed checksum"; exit 1; }

# 2) disk: the chunk set is written up front and each chunk is deleted as it is consumed, so the peak
# extra requirement is one FASTA. Derived from the actual file, never a fixed number.
FBYTES=$(stat -c %s "$FNA" 2>/dev/null || echo 0)
NEEDG=$(( FBYTES / 1073741824 + 2 ))
wsl_avail_b=$(df -B1 --output=avail "$GDIR" 2>/dev/null | tail -1 | tr -dc '0-9')
avail_b=$(effective_free_b "${wsl_avail_b:-0}")
need_b=$(( NEEDG * 1073741824 ))
[ "$avail_b" -ge "$need_b" ] || { echo "FAIL insufficient usable disk ($((avail_b/1073741824))G free after applying the smaller WSL/Windows backing-volume limit, need >=${NEEDG}G to annotate this genome)"; exit 1; }

# 3) split into contig chunks of ~CHUNKBP bases (one streaming pass; a contig is never split in half, so
# every RepeatMasker coordinate stays a real coordinate on a real contig).
# A resumed run MUST use the same library, species, sensitivity and chunk size as the chunks already on
# disk. Otherwise half the assembly is searched one way and half another, and the seal — which records a
# single value for each — would describe a run that never happened. The parameters are written beside the
# chunks at first split and compared on every resume; a mismatch stops the run and says what differs.
RUNSIG="genome_sha256=$GENOME_SHA|species=__SPSIG__|sens=$SENS|chunkbp=$CHUNKBP"
if [ -f "$ADIR/.runsig" ]; then
  OLDSIG=$(cat "$ADIR/.runsig")
  if [ "$OLDSIG" != "$RUNSIG" ]; then
    if [ -f "$ADIR/.complete" ]; then
      echo "FAIL a completed annotation of this genome exists with different settings ($OLDSIG) than the ones requested ($RUNSIG). Discard it to run again with the new settings."
    else
      echo "FAIL this genome has a part-finished annotation run with different settings ($OLDSIG) than the ones requested ($RUNSIG). Finish or discard that run before starting a different one."
    fi
    exit 1
  fi
fi
printf '%s' "$RUNSIG" > "$ADIR/.runsig"
if [ ! -f "$ADIR/.split_done" ]; then
  plog "splitting genome into work chunks"
  # No completion marker means the split was interrupted. Rebuild that scratch set from the verified
  # genome instead of mixing it with unsealed partial inputs or outputs from the aborted attempt.
  rm -f "$ADIR"/chunk_*.fa "$ADIR"/chunk_*.input.sha256 "$ADIR"/chunk_*.out \
        "$ADIR"/chunk_*.out.sha256 "$ADIR"/all.out
  awk -v dir="$ADIR" -v lim="$CHUNKBP" '
    /^>/ { if (n > 0 && bp >= lim) { close(f); i++; bp = 0 }
           if (f == "") { i = 1 } ; f = sprintf("%s/chunk_%04d.fa", dir, i); n++ }
    f == "" { next }                      # sequence before any header: malformed, skip rather than write to ""
    { print >> f }
    !/^>/ { bp += length($0) }
  ' "$FNA" || { echo "FAIL split"; exit 1; }
  SPLIT_N=0
  for c in "$ADIR"/chunk_*.fa; do
    [ -e "$c" ] || continue
    SPLIT_N=$((SPLIT_N+1))
    sha256sum "$c" | cut -d' ' -f1 > "${c%.fa}.input.sha256" \
      || { echo "FAIL could not seal work chunk $c"; exit 1; }
  done
  [ "$SPLIT_N" -ge 1 ] || { echo "FAIL no chunks produced"; exit 1; }
  touch "$ADIR/.split_done"
fi

# A resume may reuse work written hours or days earlier. Verify every pending FASTA before RepeatMasker
# reads it and every completed output before it is credited to the final summary. Runs created before
# these seals existed stop explicitly: there is no defensible baseline from which to infer integrity.
verify_chunk_inputs(){
  for c in "$ADIR"/chunk_*.fa; do
    [ -e "$c" ] || continue
    seal="${c%.fa}.input.sha256"
    [ -f "$seal" ] || { echo "FAIL pending work chunk lacks its integrity seal; discard this partial run"; return 1; }
    want_chunk=$(head -1 "$seal")
    got_chunk=$(sha256sum "$c" | cut -d' ' -f1)
    [ "$got_chunk" = "$want_chunk" ] \
      || { echo "FAIL pending work chunk failed its integrity seal; discard this partial run"; return 1; }
  done
}
verify_finished_outputs(){
  for o in "$ADIR"/chunk_*.out; do
    [ -e "$o" ] || continue
    seal="$o.sha256"
    [ -f "$seal" ] || { echo "FAIL completed annotation output lacks its integrity seal; discard this partial run"; return 1; }
    want_out=$(head -1 "$seal")
    got_out=$(sha256sum "$o" | cut -d' ' -f1)
    [ "$got_out" = "$want_out" ] \
      || { echo "FAIL completed annotation output failed its integrity seal; discard this partial run"; return 1; }
  done
}
verify_chunk_inputs || exit 1
verify_finished_outputs || exit 1
# A process can die after committing a valid output and its seal but before deleting the corresponding
# input. Both files then describe one completed chunk, not two chunks. Because both sides were verified
# above, finish that interrupted commit now; otherwise progress and sealed chunk counts double-count it.
for c in "$ADIR"/chunk_*.fa; do
  [ -e "$c" ] || continue
  base="${c%.fa}"
  [ -f "$base.out" ] || continue
  rm -f "$c" "$base.input.sha256" "$c.masked" "$c.cat" "$c.cat.gz" "$c.tbl" "$c.ori.out" \
    || { echo "FAIL could not reconcile a committed annotation chunk"; exit 1; }
done
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
idx=$DONE_N        # progress counts the chunks already finished, so a resumed run still ends at N of N
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
  # Commit the seal before moving the output into its resumable name. A crash between those moves leaves
  # the input intact and can safely rerun; it never leaves an apparently completed, unsealed output.
  [ -f "$c.out" ] || printf '' > "$c.out"
  sha256sum "$c.out" | cut -d' ' -f1 > "$base.out.sha256.tmp" \
    || { echo "FAIL could not seal completed annotation output for chunk $idx"; exit 1; }
  mv "$base.out.sha256.tmp" "$base.out.sha256" \
    || { echo "FAIL could not commit annotation output seal for chunk $idx"; exit 1; }
  mv "$c.out" "$base.out" \
    || { echo "FAIL could not commit completed annotation output for chunk $idx"; exit 1; }
  rm -f "$c" "$base.input.sha256" "$c.masked" "$c.cat" "$c.cat.gz" "$c.tbl" "$c.ori.out" 2>/dev/null
  echo "CHUNK_DONE $idx/$TOTAL"
done

# 5) combine + summarise. The summary is computed HERE so the app never reads the per-hit rows: for each
# RepeatMasker repeat name/model, its class/family label, the alignment-row count, the bases covered, and divergence
# averaged weighted by length.
verify_finished_outputs || exit 1
FINAL_GENOME_SHA=$(sha256sum "$FNA" | cut -d' ' -f1)
[ "$FINAL_GENOME_SHA" = "$GENOME_SHA" ] \
  || { echo "FAIL CACHE_CORRUPT genome FASTA changed during annotation; no result was sealed"; exit 1; }
plog "summarising"
cat chunk_*.out > all.out 2>/dev/null || true
GENOME_BP=$(awk '!/^>/ { t += length($0) } END { print t+0 }' "$FNA")
echo "GENOME_SHA256 $GENOME_SHA"
# the taxon RepeatMasker actually searched (its default lineage when no -species was given), from a chunk log
RML=$(ls chunk_*.rmlog 2>/dev/null | head -1)
if [ -n "$RML" ]; then SPL=$(sed -n '/^Species\/Taxa Search:/{n;p;q}' "$RML" 2>/dev/null); [ -n "$SPL" ] && echo "SPECIES_SEARCHED $SPL"; fi
# BEGIN_COVERAGE_SUMMARY -- tests execute this exact production block on synthetic RepeatMasker rows.
# Coverage is computed on interval unions. RepeatMasker legitimately emits
# overlapping rows for the same bases (a younger element inserted into an older one, and its own
# fragmented re-alignments), so summing row lengths double-counts sequence and can push the reported
# coverage above what is actually masked. Per-name/model, all-repeat, and biological-kind unions are kept
# separate so overlaps across names/models count once in each applicable headline. Divergence stays
# length-weighted over raw rows. The row COUNT is reported as alignment rows, which is what it is.
# RepeatMasker column 10 is the matched repeat name/model; column 11 is its broader class/family. Sort by
# both, then contig and start. Grouping on column 11 alone collapses distinct names/models in the same class.
sort -k10,10 -k11,11 -k5,5 -k6,6n all.out 2>/dev/null | awk -v gbp="$GENOME_BP" '
  function flush_name() { if (have) bp[pn]+=(pe-ps+1) }
  $1 ~ /^[0-9]+$/ {
    div=$2+0; qs=$6+0; qe=$7+0; seqid=$5; repeat_name=$10; classfam=$11
    if (qe < qs) { t=qs; qs=qe; qe=t }
    nkey=repeat_name SUBSEP classfam
    names[nkey]=repeat_name; classes[nkey]=classfam
    n[nkey]++; dv[nkey]+=div*(qe-qs+1); wt[nkey]+=(qe-qs+1); tot_n++
    key=nkey SUBSEP seqid
    if (!have || key != prev) { flush_name(); prev=key; pn=nkey; ps=qs; pe=qe; have=1 }
    else if (qs > pe+1) { flush_name(); ps=qs; pe=qe }
    else if (qe > pe) { pe=qe }
  }
  END {
    flush_name()
    printf "GENOME_BP %d\nTOTAL_HITS %d\n", gbp, tot_n+0
    for (k in n) printf "REP\t%s\t%s\t%d\t%d\t%.2f\n", names[k], classes[k], n[k], bp[k], (wt[k]>0 ? dv[k]/wt[k] : 0)
  }
' || { echo "FAIL per-repeat-name/model coverage summary"; exit 1; }

TOTAL_BP=$(sort -k5,5 -k6,6n all.out 2>/dev/null | awk '
  function flush_interval() { if (have) total+=(pe-ps+1) }
  $1 ~ /^[0-9]+$/ {
    seqid=$5; qs=$6+0; qe=$7+0; if (qe < qs) { t=qs; qs=qe; qe=t }
    if (!have || seqid != prev) { flush_interval(); prev=seqid; ps=qs; pe=qe; have=1 }
    else if (qs > pe+1) { flush_interval(); ps=qs; pe=qe }
    else if (qe > pe) { pe=qe }
  }
  END { flush_interval(); print total+0 }
') || { echo "FAIL overall coverage summary"; exit 1; }
printf "TOTAL_BP %d\n" "$TOTAL_BP"

awk '
  $1 ~ /^[0-9]+$/ {
    head=$11; sub(/\/.*/, "", head)
    if (head=="LTR" || head=="LINE" || head=="SINE" || head=="DNA" || head=="RC" || head=="Retroposon" || head=="PLE") kind="TE"
    else if (head=="Simple_repeat" || head=="Low_complexity" || head=="Satellite") kind="tandem"
    else kind="other"
    qs=$6+0; qe=$7+0; if (qe < qs) { t=qs; qs=qe; qe=t }
    printf "%s\t%s\t%d\t%d\n", kind, $5, qs, qe
  }
' all.out | sort -k1,1 -k2,2 -k3,3n | awk '
  function flush_interval() { if (have) bp[pk]+=(pe-ps+1) }
  {
    kind=$1; seqid=$2; qs=$3+0; qe=$4+0; key=kind SUBSEP seqid
    if (!have || key != prev) { flush_interval(); prev=key; pk=kind; ps=qs; pe=qe; have=1 }
    else if (qs > pe+1) { flush_interval(); ps=qs; pe=qe }
    else if (qe > pe) { pe=qe }
  }
  END { flush_interval(); for (k in bp) printf "KIND\t%s\t%d\n", k, bp[k] }
' || { echo "FAIL repeat-kind coverage summary"; exit 1; }
# END_COVERAGE_SUMMARY
FINAL_GENOME_SHA=$(sha256sum "$FNA" | cut -d' ' -f1)
[ "$FINAL_GENOME_SHA" = "$GENOME_SHA" ] \
  || { echo "FAIL CACHE_CORRUPT genome FASTA changed during annotation; no result was sealed"; exit 1; }
printf '%s' "$RUNSIG" > "$ADIR/.complete"      # marks the run finished, so a later settings clash is worded for a completed run
echo "ANNOT_OK"
'''


_LIBDIR = "$HOME/micromamba/envs/te/share/RepeatMasker/Libraries/famdb"

# RepeatMasker's class/family column mixes transposable elements with tandem and non-TE repeats. Reporting
# one merged "% repeats" as a TE landscape would overstate TE content: curated-only yeast annotation
# can mask non-TE simple/low-complexity sequence while finding no TE-classified repeat name/model.
_TE_PREFIXES = ("LTR", "LINE", "SINE", "DNA", "RC", "Retroposon", "PLE")
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


def dfam_lineage_families(species: str, include_uncurated: bool = False) -> dict:
    """How many family models the INSTALLED Dfam partitions can actually supply for this lineage.

    This is the preflight that stops a multi-hour run from quietly producing a scientifically empty result.
    Dfam is partitioned; a machine that installed only the curated partitions holds deep coverage for a few
    heavily curated lineages and almost nothing elsewhere. Curated-only, measured 2026-07-31 (see
    _DFAM_CURATED_COVERAGE): Homo sapiens 1439 models, Drosophila melanogaster 399, Saccharomyces
    cerevisiae 9, Arabidopsis thaliana 9 — so a yeast or Arabidopsis genome annotation cannot find their
    transposable elements no matter how long it runs, and the user has to be told BEFORE starting."""
    if not species or not _SPECIES_RE.match(species):
        return {"ok": False, "error": "invalid species token"}
    # `famdb families` includes both curated and uncurated records unless -c is present. Match the
    # preflight count to the library RepeatMasker will actually search: curated-only uses -c, while the
    # selected curated+uncurated search deliberately omits it. The family lines end with "len=<n>";
    # everything above them is the CC0 banner, which must not be counted.
    curated_flag = "" if include_uncurated else " -c"
    availability_gate = (
        f'[ -f "{_LIBDIR}/dfam40.uncurated.consensus.0.h5" ] && '
        f'[ -f "{_LIBDIR}/dfam40.uncurated.consensus.1.h5" ] || '
        '{ echo "TEAGLE_UNCURATED_MISSING"; exit 3; }\n'
        if include_uncurated else "")
    # famdb's own exit status decides: an ambiguous or unknown name makes it print usage and fail, and a dead
    # WSL session prints nothing, and neither is a lineage with zero family models.
    script = (availability_gate
              + f'o=$("{_MM}" run -n te famdb.py -i "{_LIBDIR}" families -a{curated_flag} -f summary "{species}" 2>/dev/null)\n'
              'frc=$?\n'
              '[ "$frc" = "0" ] || { echo "TEAGLE_FAMDB_FAILED $frc"; exit 4; }\n'
              'n=$(printf \'%s\\n\' "$o" | grep -c "len=")\n'
              'echo "$n"\n')
    try:
        rc, out, _ = _wsl_script(script, timeout=300)
        if "TEAGLE_UNCURATED_MISSING" in out:
            return {"ok": False, "species": species, "include_uncurated": True,
                    "uncurated_available": False,
                    "error": "the optional uncurated Dfam partitions are not installed"}
        if rc != 0 or "TEAGLE_FAMDB_FAILED" in out:
            return {"ok": False, "species": species,
                    "error": f"the Dfam library could not list the families for {species} "
                             "(the name may be ambiguous, or WSL was not running); no count is available"}
        n = int((out.strip().splitlines() or [""])[-1])
        if n == 0:
            # famdb exits 0 with no families for an ambiguous or unknown name ('drosophila' matches 141 taxa),
            # so a zero is only a count once the name is known to resolve to one taxon
            chk = resolve_species(species)
            if not chk.get("ok"):
                return {"ok": False, "species": species, "error": chk["error"],
                        "ambiguous_species": chk.get("ambiguous", False)}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    return {"ok": True, "species": species, "families": n,
            "include_uncurated": bool(include_uncurated),
            "uncurated_available": True if include_uncurated else None}


def annotate_budget(genome_bytes: int = 0, fasta_bytes: int | None = None, fasta_cached: bool | None = None) -> dict:
    """CPU / RAM / disk the annotation may use, bounded by WSL and its Windows backing volume.

    Disk is sized from the SOURCE FASTA (`fasta_bytes`, recorded at prepare time). `genome_bytes` is the cached
    search target, which is about a quarter of that when it is a 2bit file, so it under-states the need and is
    used only when the FASTA size is unknown. A FASTA that is not cached must be downloaded and extracted first
    (the zip and the FASTA coexist, at most twice the FASTA); the chunk copy of it comes after.

    WSL2 caps cores and memory independently of the Windows host (.wslconfig), so the host's numbers would
    over-promise. RepeatMasker's -pa starts N parallel RMBlast jobs, each needing its own memory, so the
    thread count is the LOWER of what the cores allow and what the RAM allows — the same rule the project
    applies to every other parallel step."""
    try:
        rc, out, _ = _wsl_script(
            'echo "cores=$(nproc 2>/dev/null || echo 0)"\n'
            'echo "memkb=$(awk \'/MemTotal/ {print $2}\' /proc/meminfo 2>/dev/null || echo 0)"\n'
            'echo "availb=$(df -B1 --output=avail "$HOME" 2>/dev/null | tail -1 | tr -dc \'0-9\')"\n', timeout=45)
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    d = dict(l.split("=", 1) for l in out.splitlines() if "=" in l)
    cores = int(d.get("cores") or 0)
    mem_gb = round(int(d.get("memkb") or 0) / 1048576, 1)
    wsl_free = int(d.get("availb") or 0)
    capacity = _disk_capacity(wsl_free)
    avail_gb = capacity["disk_free_gb"]
    # ~2 GB per RMBlast job is the working figure for a Dfam-sized library; keep one core for the UI.
    by_cpu = max(1, cores - 1)
    by_ram = max(1, int(mem_gb // 2))
    threads = max(1, min(by_cpu, by_ram))
    if fasta_bytes:
        need_gb = math.ceil((1 if fasta_cached else 2) * fasta_bytes / 1073741824) + 2
    else:
        need_gb = int(genome_bytes / 1073741824) + 2 if genome_bytes else None
    return {"ok": True, "cores": cores, "mem_gb": mem_gb, "avail_gb": avail_gb,
            "recommended_threads": threads, "limited_by": ("memory" if by_ram < by_cpu else "cores"),
            "disk_needed_gb": need_gb,
            "disk_ok": (avail_gb >= need_gb) if need_gb else None,
            "wsl_avail_gb": capacity["wsl_disk_free_gb"],
            "windows_avail_gb": capacity["windows_disk_free_gb"],
            "disk_limited_by": capacity["disk_limited_by"], "disk_note": capacity["disk_note"]}


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
    """Run RepeatMasker against installed Dfam models or a replacing custom FASTA. The annotation is
    chunked, resumable, and summarised inside WSL — the caller receives per-repeat-name/model aggregates, never the
    millions of individual hits (those stay on disk for export).

    `sensitivity` is one of default / quick / slow, mapping to RepeatMasker's own speed-vs-recall flags.
    Interspersed-repeat calls are homology-bound by the selected library, so the caller must state which
    library was searched next to any coverage claim. RepeatMasker's separate simple/low-complexity screen
    can still emit non-TE rows during a custom-library run."""
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
        # Dfam partitions are inert without this flag. Measured on S. cerevisiae, curated-only sees the
        # curated subset and finds ZERO transposable elements, while -uncurated adds the lineage's
        # remaining family models and recovers the Ty1/Copia and Gypsy elements the genome actually
        # carries. Uncurated families are auto-generated and lower
        # confidence, so this stays the caller's decision rather than a silent default.
        sp = f'-species "{species}"' + (" -uncurated" if include_uncurated else "")
    elif not lib_info and include_uncurated:
        return {"ok": False, "error": "select an organism/lineage before including uncurated Dfam families",
                "uncurated_available": True}
    st = env_status()
    # RepeatMasker's -lib mode searches only the supplied FASTA and does not read FamDB. Requiring Dfam
    # here would make a valid custom-library run unavailable because of a database it never uses.
    backend_ready = bool(st.get("repeatmasker") and (lib_info or st.get("dfam")))
    if not backend_ready:
        return {"ok": False, "error": _backend_not_ready(st, need_dfam=not lib_info), "status": st}
    if include_uncurated and not lib_info and not st.get("dfam_uncurated"):
        return {"ok": False, "error": "the optional uncurated Dfam partitions are not installed; install "
                "both uncurated partitions in the backend manager or use curated families only",
                "status": st, "uncurated_available": False}
    # the signature identifies WHICH library/lineage was searched, so a resume with a different one is
    # caught; the thread count is deliberately NOT part of it (it changes speed, never which families hit)
    # The INSTALLED PARTITION SET belongs in the signature too. Adding a Dfam partition changes which
    # families exist to be found, so chunks finished before an install and chunks finished after it were
    # searched against different libraries — the same "one run, two searches" defect as a changed species,
    # and just as invisible in the finished summary.
    if lib_info:
        # Dfam partitions and -uncurated are not read in -lib mode. An unused Dfam install/change must not
        # invalidate resumable chunks produced from the same supplied FASTA.
        spsig = f"lib:{lib_info['sha256']}"
        include_uncurated = False
    else:
        _parts = ",".join((st.get("dfam_library") or {}).get("partitions") or [])
        spsig = (f"species:{species or 'all'}"
                 + f"|parts:{hashlib.sha256(_parts.encode()).hexdigest()[:12]}"
                 + f"|unc:{int(bool(include_uncurated))}")  # curated-only vs +uncurated is a different search
    host_free = _windows_backing_free_bytes()
    script = (_ANNOT_SCRIPT.replace("__ACC__", accession).replace("__PA__", str(max(1, int(threads))))
              .replace("__CHUNKBP__", str(max(1, int(chunk_mb)) * 1_000_000))
              .replace("__SENS__", sens).replace("__SPSIG__", spsig).replace("__SP__", sp)
              .replace("__HOST_FREE_B__", str(host_free if host_free is not None else -1)))
    try:
        rc, out, err = _wsl_script(script, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"annotation timed out after {timeout}s — re-run to resume from the "
                                      "chunks that already finished"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}
    if "ANNOT_OK" not in out:
        failure = _operation_failure("genome annotation failed", out, err, partial=True)
        if "a completed annotation" in failure["error"]:
            failure["annotation_state"] = "complete"          # a finished run, not part-finished work, is what would be discarded
        elif "part-finished" in failure["error"]:
            failure["annotation_state"] = "part_finished"
        return failure
    repeat_names, totals = [], {}
    lineage = taxid = None
    for line in out.splitlines():
        if line.startswith("REP\t"):
            _, repeat_name, class_family, n, bp, div = line.split("\t")
            repeat_names.append({"repeat_name": repeat_name, "class_family": class_family, "n": int(n),
                                 "bp": int(bp), "divergence": float(div)})
        elif line.startswith(("TOTAL_HITS ", "TOTAL_BP ", "GENOME_BP ", "TOTAL_CHUNKS ")):
            k, v = line.split(); totals[k.lower()] = int(v)
        elif line.startswith("KIND\t"):
            _, kind, bp = line.split("\t")
            totals[kind.lower() + "_bp"] = int(bp)
        elif line.startswith("GENOME_SHA256 "):
            totals["genome_sha256"] = line.split(None, 1)[1].strip()
        elif line.startswith("SPECIES_SEARCHED "):
            sm = re.match(r"SPECIES_SEARCHED\s+(.+?)\s*(?:\[NCBI Taxonomy ID:\s*(\d+)\])?\s*$", line)
            if sm:
                lineage, taxid = sm.group(1), sm.group(2)
    gbp = totals.get("genome_bp", 0)
    genome_sha = totals.get("genome_sha256")
    for row in repeat_names:
        row["percent"] = round(100.0 * row["bp"] / gbp, 4) if gbp else None
        row["kind"] = repeat_kind(row["class_family"])
    repeat_names.sort(key=lambda row: -row["bp"])
    # Do not attach an installed-but-unused Dfam identity to a custom-library result.
    lib = None if lib_info else st.get("dfam_library")
    # Headline kind coverage comes from contig-wise unions emitted by the shell summary. Summing repeat-name/model
    # totals would count a base twice whenever two TE repeat names/models overlap.
    te_bp = totals.get("te_bp", 0)
    n_te_names = len({row["repeat_name"] for row in repeat_names if row["kind"] == "TE"})
    # with a custom library the lineage count is meaningless — what was searched is the user's file
    avail = (dfam_lineage_families(species, include_uncurated=include_uncurated)
             if (species and not lib_info) else {"ok": False})
    return {"ok": True, "accession": accession, "repeat_names": repeat_names, "genome_sha256": genome_sha,
            "total_hits": totals.get("total_hits", 0), "total_bp": totals.get("total_bp", 0),
            "genome_bp": gbp, "chunks": totals.get("total_chunks", 0),
            # masked_percent is EVERY repeat RepeatMasker placed; te_percent is the transposable-element
            # subset. They are reported separately because they are different claims — a genome can mask
            # several percent while containing no detected TE at all.
            "masked_percent": round(100.0 * totals.get("total_bp", 0) / gbp, 3) if gbp else None,
            "te_bp": te_bp, "te_percent": round(100.0 * te_bp / gbp, 3) if gbp else None,
            "tandem_bp": totals.get("tandem_bp", 0), "other_bp": totals.get("other_bp", 0),
            "te_repeat_name_count": n_te_names,
            "library_families_for_species": avail.get("families") if avail.get("ok") else None,
            # the honest headline when a run completes but the installed library had nothing to find with
            "library_kind": ("user-supplied FASTA library" if lib_info else
                             ("installed Dfam partitions, curated + uncurated" if include_uncurated
                              else "installed Dfam partitions, curated only")),
            "custom_library": lib_info, "include_uncurated": bool(include_uncurated),
            "coverage_warning": (
                None if n_te_names else
                ((("The supplied library produced homology matches, but none carried a recognized "
                   "RepeatMasker TE class/family label. Those rows are retained as unclassified/other and "
                   "excluded from TE coverage because their classification is unavailable, not because they "
                   "were demonstrated to be non-TE. Use FASTA headers such as >name#LTR/Gypsy when a reference "
                   "has a supported TE class.")
                   if repeat_names else
                  ("No TE-classified repeat name/model was found. The supplied library "
                   f"({lib_info['name']}, {lib_info['sequences']} sequences) placed no repeat in this assembly — "
                   "check that the library matches this organism before drawing any conclusion."))
                 if lib_info else
                 ("No TE-classified repeat name/model was found. This is a limit of what was SEARCHED, not a "
                  "finding about the genome"
                  + (f" — the installed Dfam partitions hold {avail['families']} family model(s) for "
                     f"{species or 'this lineage'}" if avail.get("ok") else "")
                  + (", and this run used only the CURATED subset of them. Most families outside a few "
                     "heavily curated species are uncurated in Dfam, so re-running with uncurated families "
                     f"included is the first thing to try: measured on S. cerevisiae, curated-only sees "
                     f"{_DFAM_YEAST_CURATED} of {_DFAM_YEAST_TOTAL} family models and finds no transposable "
                     f"element, while including uncurated adds {_DFAM_YEAST_ADDITIONAL} models and recovers "
                     "the Ty elements the genome carries."
                     if not include_uncurated else
                     ". Uncurated families were already included, so consider whether Dfam covers this "
                     "lineage at all, or supply a species-specific library.")))),
            "repeatmasker_version": st["repeatmasker"], "dfam_version": (lib or {}).get("version"),
            "dfam_library": lib,
            # no organism means RepeatMasker's default lineage, not every installed family
            "species": (None if lib_info else
                        (species or f"{lineage or _RM_DEFAULT_LINEAGE} (RepeatMasker default lineage)")),
            "species_requested": species, "species_defaulted": bool(not lib_info and not species),
            "lineage_searched": None if lib_info else lineage, "lineage_taxid": None if lib_info else taxid,
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
                            "n_seqs": int(m.get("n_seqs", 0) or 0), "bytes": int(m.get("bytes", 0) or 0),
                            "fasta_bytes": int(m.get("fasta_bytes", 0) or 0)})
    return {"ok": True, "genomes": genomes}


def genome_remove(accession: str) -> dict:
    """Delete a cached genome to reclaim disk.

    Refuses while an annotation is running against it. The annotation's working directory lives INSIDE the
    genome directory ({_GENOMES}/<acc>/annot), so an unguarded rm -rf here deletes the working set out from
    under a live multi-hour run — the same hazard genome_annotate_reset already guards, checked the same way."""
    if not _ACC_RE.match(accession or ""):
        return {"ok": False, "error": "invalid assembly accession"}
    av = available()
    if not av["wsl2"]:
        return {"ok": False, "error": "WSL2 not available"}
    rc, out, err = _wsl_script(f'A="{_GENOMES}/{accession}/annot"\n'
                               f'if [ -d "$A/.lock" ] && [ -f "$A/.lock/pid" ] && kill -0 "$(cat "$A/.lock/pid")" 2>/dev/null; then\n'
                               f'  echo RUNNING; exit 3\n'
                               f'fi\n'
                               f'rm -rf "{_GENOMES}/{accession}" && echo REMOVED\n', timeout=60)
    if rc == 3 or "RUNNING" in out:
        return {"ok": False, "error": "an annotation is still running against this genome — let it finish, "
                                      "or discard the annotation first, before deleting the genome"}
    return {"ok": "REMOVED" in out, "accession": accession,
            **({} if "REMOVED" in out else {"error": "could not delete the genome: " + err.strip()[:160]})}
