# Changelog

All notable changes to TEagle are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The version is defined once in `app/backend/teagle_core/__init__.py` (`__version__`)
and propagates to the backend health endpoint, the UI header badge, every run
provenance manifest, the packaged executable's Windows file-version metadata, and
the LaTeX report title page.

## [2.0.0] — 2026-07-18

Native desktop rewrite. The user interface is now a **native PySide6/Qt application**;
the browser + WebView2 stack is retired. All scientific behaviour and results are
unchanged — the same validated engine drives both — but the app is now a true native
window with no embedded browser and no local web server on the core path.

### Added — native parity & install hardening (2026-07-19)
- **Right-click context menus** on structural, ORF, and domain rows (copy FASTA/DNA/
  coordinates/protein, design a primer here), matching the web UI's feature menu.
- **Interactive figures**: hover any genome-viewer feature or gel band for its size/type,
  and right-click it to copy FASTA/coordinates or design a primer.
- **Source-citation links** ("source ↗") on the classification, structural, domain, family,
  splice, and specimen panels — opening the verified DOI/record for Wicker 2007, Pfam,
  Dfam, RepeatMasker, Primer3, minimap2, and NCBI.
- **Complete glossary tooltips** on every results-table header.
- **Dedicated backend installer** dialog: per-component status (WSL2, micromamba,
  RepeatMasker, minimap2, Dfam root/curated, FamDB config), one-click install-all,
  per-component **repair**, and a **check-integrity** pass, with a live log — off the GUI thread.

### Fixed
- **Backend install failed at the first step for real users**: the WSL session starts in the
  Windows `/mnt/c` mount, so `tar -xj bin/micromamba` hit `Permission denied`. Every install
  step now runs from `$HOME` and extracts to an absolute path. The status probe is delivered
  via stdin (not an inline `bash -lc` argument), because `wsl.exe` mangled embedded quotes and
  made a healthy backend read as "missing".
- The Dfam download now uses a **pinned, versioned** URL with **embedded md5 trust anchors**
  and resumable transfers, instead of the moving `current` pointer and a runtime-fetched checksum.
- **No console-window flashes**: `wsl.exe` and pip helpers spawn with `CREATE_NO_WINDOW`, so the
  windowed app never pops a terminal.

### Changed
- **Native PySide6/Qt UI** replacing the HTML/JS UI hosted in WebView2. Feature parity:
  specimen intake, classification & structure with an interactive genome viewer, primer
  design, staged multi-lane in-silico PCR with a to-scale gel, Dfam family annotation,
  de-novo splice detection, run provenance, dark/light themes, glossary tooltips,
  context-menu FASTA/DNA/coordinate copy, and CSV/TSV + SVG/PNG export.
- **In-process engine.** A new `engine.py` is the single source of truth for every
  operation (request validation + science); the native app calls it directly on a
  background thread (three-tier result taxonomy: success / user error / unexpected fault).
  `server.py` is now a thin HTTP adapter over the same engine, kept only for the legacy
  `--server` web mode. Request-validation behaviour (non-finite parameters, non-string
  species, non-nucleotide input) is now covered by tests at the engine level.
- **Figures** (genome viewer, agarose gel) are rendered through Qt's SVG engine from the
  same builders used for export, so on-screen and exported figures match.
- **Packaging** now freezes the PySide6 app (`installer/teagle_native.spec`); the frozen
  bundle ships PySide6, Primer3, HMMER, and the Pfam profiles, and the build is gated on a
  headless self-test that proves the scientific stack imports and Qt SVG renders. The
  installer build is one command (`installer/build_installer.ps1`). A single-instance mutex
  and a kill-on-close Job Object clean up any WSL subprocess tree when the app closes.

### Verified
- Full hermetic suite (backend + engine + native Qt), golden-fixture classification routed
  through the new in-process adapter (copia→LTR/Copia, gypsy→LTR/Gypsy, L1→LINE, Tc1→DNA/Tc1,
  Ac→DNA/hAT), and explicit broken-WSL / no-WSL / stack-missing degradation trials.

### Branding & typography (2026-07-20)
- **Cascadia Mono across the whole app.** The UI, tables, data fields, and figure labels
  (gel ladder, genome ruler) now render in Cascadia Mono. The fonts are **bundled in the
  app** (`app/native/assets/fonts`, SIL OFL 1.1) and loaded at startup, so the interface
  looks identical without them installed. The frozen self-test fails if any bundled font
  is missing from the build.
- **Clean Cascadia Code wordmark.** The header/README wordmark is Cascadia Code Bold rendered
  from the static face — the earlier square notches (an artifact of synthetic-bolding the
  variable font) are gone. TE is white, AGLE takes the eagle mark's teal. The wordmark sits
  a step smaller than the eagle mark.
- **Eagle logo in the installer wizard.** The Inno Setup wizard now shows the eagle mark as
  its top-right and Welcome/Finish images (DPI ladder of BMPs).

## [1.0.0] — 2026-07-18

First public release. A self-contained Windows desktop application for
transposable-element annotation, gene-structure inspection, and TE-aware PCR
primer design, usable without a command line.

### Added
- **Sequence intake** — paste, FASTA upload, or NCBI accession fetch (E-utilities,
  content-addressed cache) with IUPAC validation and composition summary.
- **Structural detection** — LTR and terminal-inverted-repeat (TIR) finders with
  terminal-arm boundary detection at the identity cliff; ORF finder; TSD reporting.
- **Protein-domain scan** — native HMMER (pyhmmer) against a bundled CC0 Pfam
  TE-domain profile set; hits mapped back to nucleotide coordinates. Fully offline.
- **Classification** — evidence-backed superfamily/order assignment with transparent
  confidence and per-call reasoning.
- **Primer design** — Primer3 (primer3-py) with exposed parameters (size, Tm, GC,
  product size, poly-X, GC clamp) and domain/target region constraints.
- **In-silico PCR** — pair-aware, both-strand search with a strict 3' matching rule,
  multi-background specificity check, and gel-lane amplicon visualization.
- **Genome viewer** — interactive browser-style track view with wheel-zoom, crosshair,
  keyboard navigation, an overview minimap, and transparent-background figure export
  (PNG/SVG).
- **Exports** — CSV/TSV for every result table.
- **Reproducibility** — every result carries a provenance manifest packing exact tool
  and database versions, input checksum, parameters, and environment; identical inputs
  yield an identical content-addressed seal.
- **Optional WSL annotation backend** — RepeatMasker + Dfam curated + minimap2,
  auto-installed by the app; core analysis runs without it.

### Packaging
- One-click Windows build (PyInstaller onedir): bundled Python and C extensions, so no
  system Python, pip, or manual dependency downloads are required.
- Graceful degradation: a broken or missing primer3/pyhmmer disables only its feature;
  the rest of the engine keeps running.
- Real native application window (pywebview + Microsoft Edge WebView2): its own window and
  taskbar entry, no browser chrome. Falls back to a chromeless Edge/Chrome `--app` window if
  the WebView2 runtime is absent. A kill-on-close Job Object ties the whole process tree to the
  launcher, so an in-place upgrade never orphans a window.

[2.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.0.0
[1.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v1.0.0
