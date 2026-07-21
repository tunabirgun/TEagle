# PyInstaller spec for the native PySide6 TEagle — one-click Windows build (onedir).
# Build:  pyinstaller installer/teagle_native.spec --noconfirm
# Everything the science needs (primer3, pyhmmer, the CC0 Pfam HMM profiles) is bundled;
# a non-coder installs one thing and runs it. The optional Dfam/splice features install from inside the app.
import os, re
from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata

ROOT = os.path.abspath(os.path.join(SPECPATH, ".."))         # project root (SPECPATH = installer/)

# single source of truth: parse __version__ from teagle_core (no import -> no heavy deps at build time)
_vsrc = open(os.path.join(ROOT, "app", "backend", "teagle_core", "__init__.py"), encoding="utf-8").read()
VERSION = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', _vsrc).group(1)
_vt = tuple(int(x) for x in (re.findall(r"\d+", VERSION) + [0, 0, 0, 0])[:4])

_verinfo = os.path.join(SPECPATH, "_version_info.txt")
open(_verinfo, "w", encoding="utf-8").write(f"""VSVersionInfo(
  ffi=FixedFileInfo(filevers={_vt}, prodvers={_vt}, mask=0x3f, flags=0x0,
                    OS=0x40004, fileType=0x1, subtype=0x0, date=(0, 0)),
  kids=[
    StringFileInfo([StringTable('040904B0', [
      StringStruct('CompanyName', 'Tuna Birgun'),
      StringStruct('FileDescription', 'TEagle - Transposable Elements Assay Terminal'),
      StringStruct('FileVersion', '{VERSION}'),
      StringStruct('InternalName', 'TEagle'),
      StringStruct('LegalCopyright', 'MIT License'),
      StringStruct('OriginalFilename', 'TEagle.exe'),
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
# bundle the Cascadia fonts (SIL OFL) so the UI renders identically without them installed -> native/assets/fonts
_fontdir = os.path.join(ROOT, "app", "native", "assets", "fonts")
datas += [(os.path.join(_fontdir, fn), "native/assets/fonts") for fn in os.listdir(_fontdir)]
# metadata so envcheck's importlib.metadata.version(...) resolves for the bundled deps (env panel stays green)
for dist in ("PySide6", "shiboken6", "primer3-py", "pyhmmer"):
    try:
        datas += copy_metadata(dist)
    except Exception:
        pass

binaries = []
hiddenimports = [
    # native UI package (app/native) + engine adapter (app/backend)
    "main", "engine", "engine_worker", "figures", "widgets", "sample", "theme", "fonts", "install_dialog",
    "envcheck", "server",
    "teagle_core", "teagle_core.appdirs", "teagle_core.sequtil", "teagle_core.structural",
    "teagle_core.classify", "teagle_core.domains", "teagle_core.primers", "teagle_core.fetch",
    "teagle_core.refs", "teagle_core.provenance", "teagle_core.timing", "teagle_core.wsl",
    # Qt modules actually used (PySide6's hook collects the matching libs + plugins, incl. QtSvg)
    "PySide6.QtCore", "PySide6.QtGui", "PySide6.QtWidgets", "PySide6.QtSvg",
]

# C-extension packages + openpyxl (XLSX table export): pull binaries, data and dynamic submodules
for pkg in ("primer3", "pyhmmer", "openpyxl"):
    d, b, h = collect_all(pkg)
    datas += d; binaries += b; hiddenimports += h
    hiddenimports += collect_submodules(pkg)
hiddenimports += ["pyhmmer.platform", "pyhmmer.platform.win32", "openpyxl", "et_xmlfile"]

# trim the bundle: exclude heavy Qt modules the app never imports (keeps the installer small)
excludes = [
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
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True,
    name="TEagle", debug=False, strip=False, upx=False,
    console=False,                                           # GUI app, no console window
    icon=_icon if os.path.exists(_icon) else None,
    version=_verinfo,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="TEagle")
