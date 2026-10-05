# PyInstaller spec for the native PySide6 TEagle — one-click Windows build (onedir).
# Build:  pyinstaller installer/teagle_native.spec --noconfirm
# Everything the science needs (primer3, pyhmmer, the CC0 Pfam HMM profiles) is bundled;
# a non-coder installs one thing and runs it. The optional Dfam/splice features install from inside the app.
import os, re, sys
from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata

ROOT = os.path.abspath(os.path.join(SPECPATH, ".."))         # project root (SPECPATH = installer/)

# single source of truth: parse __version__ from teagle_core (no import -> no heavy deps at build time)
_vsrc = open(os.path.join(ROOT, "app", "backend", "teagle_core", "__init__.py"), encoding="utf-8").read()
VERSION = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', _vsrc).group(1)
_vt = tuple(int(x) for x in (re.findall(r"\d+", VERSION) + [0, 0, 0, 0])[:4])

os.makedirs(workpath, exist_ok=True)                         # PyInstaller's build workpath: the generated file stays out of the tree
_verinfo = os.path.join(workpath, "_version_info.txt")
open(_verinfo, "w", encoding="utf-8").write(f"""VSVersionInfo(
  ffi=FixedFileInfo(filevers={_vt}, prodvers={_vt}, mask=0x3f, flags=0x0,
                    OS=0x40004, fileType=0x1, subtype=0x0, date=(0, 0)),
  kids=[
    StringFileInfo([StringTable('040904B0', [
      StringStruct('CompanyName', 'Tuna Birgun'),
      StringStruct('FileDescription', 'TEagle - Transposable Elements Assay Terminal'),
      StringStruct('FileVersion', '{VERSION}'),
      StringStruct('LegalCopyright', 'Copyright (C) 2026 Tuna Birgun. Licensed under AGPL-3.0-or-later.'),
      StringStruct('ProductName', 'TEagle'),
      StringStruct('ProductVersion', '{VERSION}')])]),
    VarFileInfo([VarStruct('Translation', [1033, 1200])])])
""")

datas = [
    (os.path.join(ROOT, "app", "backend", "data", "te_domains.hmm"), "data"),   # CC0 Pfam TE-domain profiles
    (os.path.join(ROOT, "app", "backend", "requirements.txt"), "."),            # for the in-app environment panel
    (os.path.join(ROOT, "app", "native", "assets", "teagle-mark.svg"), "native/assets"),      # header brand mark
    (os.path.join(ROOT, "app", "native", "assets", "teagle-wordmark.svg"), "native/assets"),  # header wordmark (Cascadia Code)
    (os.path.join(SPECPATH, "teagle.ico"), "."),                                # runtime window/taskbar QIcon (same crisp ICO the exe embeds)
]
# bundled per-assembly chromosome maps for coordinate fetch -> teagle_core/data/assemblies/*.json
_asmdir = os.path.join(ROOT, "app", "backend", "teagle_core", "data", "assemblies")
datas += [(os.path.join(_asmdir, f), "teagle_core/data/assemblies") for f in os.listdir(_asmdir) if f.endswith(".json")]
# bundled example elements (real published specimens; ~40 KB) -> teagle_core/data/examples
_exdir = os.path.join(ROOT, "app", "backend", "teagle_core", "data", "examples")
datas += [(os.path.join(_exdir, f), "teagle_core/data/examples") for f in os.listdir(_exdir) if f.endswith(".fasta")]
# bundle the Cascadia fonts (SIL OFL) so the UI renders identically without them installed -> native/assets/fonts
_fontdir = os.path.join(ROOT, "app", "native", "assets", "fonts")
datas += [(os.path.join(_fontdir, fn), "native/assets/fonts") for fn in os.listdir(_fontdir)]
# metadata so envcheck's importlib.metadata.version(...) resolves for the bundled deps (env panel stays green)
for dist in ("PySide6", "shiboken6", "primer3-py", "pyhmmer", "psutil", "openpyxl", "et_xmlfile", "certifi", "pypdf"):
    datas += copy_metadata(dist)                              # a missing distribution must fail the build, not ship a red Environment panel

binaries = []
hiddenimports = [
    # native UI package (app/native) + engine adapter (app/backend)
    "main", "engine", "engine_worker", "figures", "widgets", "sample", "theme", "fonts", "install_dialog",
    "coords", "report", "envcheck",
    "ui_common", "panel_chrome", "panel_family", "panel_genome", "panel_genome_annotate",
    "panel_pcr", "panel_primers", "panel_provenance", "panel_records", "panel_report",
    "panel_specimen", "panel_splice", "panel_structure",
    "teagle_core", "teagle_core.appdirs", "teagle_core.sequtil", "teagle_core.structural",
    "teagle_core.classify", "teagle_core.domains", "teagle_core.primers", "teagle_core.oligoqc",
    "teagle_core.fetch", "teagle_core.refs", "teagle_core.provenance", "teagle_core.timing", "teagle_core.wsl",
    "teagle_core.examples", "teagle_core.oligoqc_wsl", "teagle_core.retroviral", "teagle_core.genomepcr",
    "teagle_core.artifact_metadata", "teagle_core.dotplot", "teagle_core.gff3",
    # Qt modules actually used (PySide6's hook collects the matching libs + plugins, incl. QtSvg)
    "PySide6.QtCore", "PySide6.QtGui", "PySide6.QtWidgets", "PySide6.QtSvg",
]

# openpyxl.utils.dataframe is a pandas-interop shim TEagle never imports; collecting it manufactured the
# edge openpyxl -> pandas -> scipy -> torch -> transformers (~528 MB of dead weight). Drop it at the source.
_SHIM = "openpyxl.utils.dataframe"

# collect_all also returns tests shipped inside some wheels. They are development payload, not runtime
# data: filter them by complete path/module component so names such as ``latest`` remain untouched.
_TEST_COMPONENTS = {"test", "tests", "testing", "testdata", "test_data", "fixture", "fixtures"}


def _is_test_payload(parts):
    return any(part.lower() in _TEST_COMPONENTS for part in parts)


def _has_test_filename(path):
    stem = os.path.splitext(os.path.basename(path))[0].lower()
    return "test" in re.split(r"[^a-z0-9]+", stem)


def _collect_runtime_package(package):
    package_datas, package_binaries, package_imports = collect_all(package)
    package_datas = [
        item for item in package_datas
        if not _is_test_payload(item[1].replace("\\", "/").split("/"))
        and not _has_test_filename(item[0])
    ]
    package_imports = [
        module for module in package_imports
        if not _is_test_payload(module.split("."))
    ]
    package_imports += [
        module for module in collect_submodules(package)
        if not _is_test_payload(module.split("."))
    ]
    return package_datas, package_binaries, package_imports


def _is_forbidden_bundle_path(path):
    parts = str(path).replace("\\", "/").lower().split("/")
    name = parts[-1]
    return (
        "qml" in parts
        # Qt uses Windows ICU C APIs; an unrelated ICU found on PATH can shadow them.
        or name in {"icu.dll", "icuuc.dll", "icuin.dll"}
        or (name.startswith("icudt") and name.endswith(".dll"))
        or any("virtualkeyboard" in part for part in parts)
        or name == "qpdf.dll"
        or name.startswith("qt6pdf")
        or name.startswith("qt6quick")
        or name.startswith("qt6qml")
        or name.startswith("qtquick")
        or name.startswith("qtqml")
    )


def _without_forbidden_bundle_entries(entries):
    return [entry for entry in entries if not _is_forbidden_bundle_path(entry[0])]


def _assert_bundle_inventory_clean(paths, label):
    forbidden = sorted({str(path) for path in paths if _is_forbidden_bundle_path(path)})
    if forbidden:
        raise SystemExit(
            f"teagle_native.spec BUNDLE GUARD: forbidden unused Qt payload entered {label}:\n  "
            + "\n  ".join(forbidden)
        )


def _bundle_inventory(root):
    return [
        os.path.relpath(os.path.join(directory, filename), root)
        for directory, _, filenames in os.walk(root)
        for filename in filenames
    ]


# C-extension packages + openpyxl (XLSX table export): pull binaries, runtime data and dynamic submodules
for pkg in ("primer3", "pyhmmer", "psutil", "openpyxl", "et_xmlfile", "certifi", "pypdf"):
    d, b, h = _collect_runtime_package(pkg)
    datas += d; binaries += b; hiddenimports += h
hiddenimports = [m for m in hiddenimports if m != _SHIM]
hiddenimports += ["pyhmmer.platform", "pyhmmer.platform.win32", "openpyxl", "et_xmlfile", "certifi", "pypdf"]

# trim the bundle: exclude heavy Qt modules the app never imports (keeps the installer small)
excludes = [
    _SHIM,                                                   # pandas-interop shim -> pandas/scipy/torch/transformers
    # openpyxl's three try/except optional accelerators. TEagle's XLSX export writes str/int/float cells and a
    # bold header only: no ws.add_image (PIL), no numpy scalars, and the stdlib ElementTree writer is equivalent.
    "numpy", "PIL", "lxml",
    # pypdf carries version-guarded fallback imports for Python <3.11. The controlled build is exactly
    # Python 3.12 and executes the stdlib typing branches, so do not bundle an undeclared backport.
    "typing_extensions",
    # pypdf can optionally inspect and rewrite embedded font programs through fontTools. TEagle only
    # clones PDF objects and removes document metadata; the optional font-rewrite path is never used.
    "fontTools",
    # ViennaRNA is NOT redistributable under AGPL-3.0: its licence forbids redistribution for a fee and
    # asks that commercial inclusion be cleared with the authors, which GPL forbids adding. oligoqc
    # already guards `import RNA`, so an absent ViennaRNA degrades to primer3-only with the cross-check
    # reported as unavailable. Users who install it themselves still get the dual-engine check.
    "RNA",
    "tkinter", "test", "pytest", "playwright", "webview", "pywebview",
    "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineQuick",
    "PySide6.QtQml", "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtQuickWidgets",
    "PySide6.Qt3DCore", "PySide6.Qt3DRender", "PySide6.Qt3DExtras",
    "PySide6.QtMultimedia", "PySide6.QtMultimediaWidgets", "PySide6.QtCharts",
    "PySide6.QtDataVisualization", "PySide6.QtPdf", "PySide6.QtPdfWidgets",
    "PySide6.QtDesigner", "PySide6.QtBluetooth", "PySide6.QtPositioning",
    "PySide6.QtWebSockets", "PySide6.QtWebChannel", "PySide6.QtSql", "PySide6.QtTest",
    "PySide6.QtSensors", "PySide6.QtSerialPort", "PySide6.QtNfc", "PySide6.QtSpatialAudio",
    "PySide6.QtHelp", "PySide6.QtOpenGL", "PySide6.QtOpenGLWidgets", "PySide6.QtQuickControls2",
]

_icon = os.path.join(SPECPATH, "teagle.ico")

a = Analysis(
    [os.path.join(SPECPATH, "teagle_native.py")],
    pathex=[os.path.join(ROOT, "app", "backend"), os.path.join(ROOT, "app", "native")],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[], hooksconfig={}, runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)

# Qt's generic image/input hooks can add native payload even when the corresponding Python modules
# are excluded. Filter the resolved collection tables so unused PDF, Quick/QML and virtual-keyboard
# dependency chains cannot reach COLLECT.
a.binaries = _without_forbidden_bundle_entries(a.binaries)
a.datas = _without_forbidden_bundle_entries(a.datas)

# build guard: only the pinned deps (requirements.txt) + their real transitive deps + TEagle's own modules
# may enter the module graph. One stray import edge once dragged in ~528 MB; fail the build, don't ship it.
_ALLOWED_TOP = {
    "PySide6", "shiboken6",                                  # PySide6==6.11.1 (+ its shiboken6 runtime)
    "primer3", "pyhmmer", "psutil",                          # primer3-py, pyhmmer (+ pinned psutil)
    "openpyxl", "et_xmlfile",                                # openpyxl (+ its only dependency)
    "certifi", "pypdf",
    "teagle_core", "teagle_native", "main", "engine", "engine_worker", "envcheck",
    "figures", "widgets", "sample", "theme", "fonts", "install_dialog", "coords", "report",   # app/backend + app/native
    "ui_common", "panel_chrome", "panel_family", "panel_genome", "panel_genome_annotate",
    "panel_pcr", "panel_primers", "panel_provenance", "panel_records", "panel_report",
    "panel_specimen", "panel_splice", "panel_structure",
    "_pyi_rth_utils",                                        # PyInstaller runtime shim
}
_counts = {}
for _m in a.pure:
    _t = _m[0].split(".")[0]
    if _t not in _ALLOWED_TOP and _t not in sys.stdlib_module_names:
        _counts[_t] = _counts.get(_t, 0) + 1
if _counts:
    raise SystemExit(
        "teagle_native.spec BUILD GUARD: unlisted top-level packages entered the bundle:\n  "
        + "\n  ".join(f"{t} ({n} modules)" for t, n in sorted(_counts.items(), key=lambda kv: -kv[1]))
        + "\nEither exclude the import edge that pulled it in, or add it to _ALLOWED_TOP if it is a real dependency."
    )

pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True,
    name="TEagle", debug=False, strip=False, upx=False,
    console=False,                                           # GUI app, no console window
    icon=_icon if os.path.exists(_icon) else None,
    version=_verinfo,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="TEagle")
_assert_bundle_inventory_clean((entry[0] for entry in coll.toc), "the COLLECT table")
_assert_bundle_inventory_clean(_bundle_inventory(coll.name), "the assembled bundle")
