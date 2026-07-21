# Changelog

All notable changes to TEagle are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The version is defined once in `app/backend/teagle_core/__init__.py` (`__version__`)
and propagates to the backend health endpoint, the UI header badge, every run
provenance manifest, the packaged executable's Windows file-version metadata, and
the LaTeX report title page.

## [2.4.0] — 2026-07-21

Fetch by genomic coordinate, an explicit table-export format menu, and a fixed Windows taskbar icon.

### Added
- **Fetch by coordinate (UCSC-style).** Alongside accession fetch, the Specimen panel gains a **Fetch by coordinate** section: pick an organism from 17 curated reference assemblies (or **Other organism / assembly…** to resolve any species or GCF/GCA accession through NCBI Datasets), choose the strand, and paste one or more loci in browser notation — `chr13:33,016,423-33,066,143`, one region per line for multi-region. Coordinates are 1-based inclusive, identical to the UCSC/NCBI browser display, so the numbers pass through with no conversion. Organism-specific chromosome names (`2L`, roman numerals, `X`, `MT`) resolve against the assembly's own map. Each region is fetched from NCBI E-utilities as an exact base range; multi-region fetches concatenate all regions into the sequence box, and analysis runs on the first region. The pinned assembly accession, taxon id, resolved chromosome RefSeq accessions, coordinates, and strand are recorded in the run provenance seal, so a coordinate run is as reproducible as an accession run.
- **Explicit table-export format menu.** The **Export table** button and the table right-click menu now offer **Excel (.xlsx)**, **CSV**, and **TSV** as named choices instead of hiding the format behind a save-dialog filter. The save dialog opens pre-set to the chosen format and appends the extension if you omit it.

### Fixed
- **Windows taskbar icon.** The running app and its top-level window now set the bundled icon explicitly, so TEagle shows its eagle mark in the taskbar and Alt-Tab instead of a generic placeholder.

## [2.3.0] — 2026-07-21

Install WSL2 directly from the app, plus a crisper desktop-shortcut icon.

### Added
- **Install WSL2 from the backend installer.** When the Windows Subsystem for Linux is absent, the installer now offers an **Install WSL** action that installs WSL2 + Ubuntu through an elevated helper (Windows UAC) instead of only printing manual steps. It distinguishes an absent WSL from a registered-but-won't-start distro (guiding you to unregister and reinstall), surfaces progress in the log, and reports when a Windows restart is required. The **Install / update all** button routes through the WSL install first when WSL is missing.

### Changed
- **Sharper app and desktop-shortcut icon** — the eagle icon is rendered at higher supersampling with edge sharpening on the small frames, so the shortcut and taskbar read crisply instead of soft.

## [2.2.0] — 2026-07-21

In-silico PCR and gel-imaging upgrades for repeat-rich elements, plus spreadsheet-native table export.

### Added
- **Single-primer (self-priming) in-silico PCR products.** The amplicon search now also reports products a single primer makes by priming in both orientations across an inverted repeat (F+F / R+R) — common at TE terminal inverted repeats and LTRs. These are flagged distinctly, listed in the amplicon table, exported with a `singleprimer` tag, and drawn in their own gel colour.
- **Excel (XLSX) table export with a visible button.** The protein-domain and Dfam / RepeatMasker family tables gain an **Export table** button; every table now exports to CSV, TSV, and native XLSX (numbers typed as numbers, header frozen and bold, spreadsheet-formula injection neutralised) in addition to the existing right-click export.

### Changed
- **Gel imaging.** Equal-size amplicons co-migrate into a single band (a gel cannot resolve them) and the on-target colour always wins, so an on-target band is never painted over by an off-target of the same size. Band intensity now tracks priming efficiency (fewer mismatches → brighter). A lane carrying bands but no intended product is labelled **no on-target**, and the legend gains a single-primer swatch when relevant.

## [2.1.1] — 2026-07-20

Bugfix release: the optional WSL annotation backend failed to install on some users' machines.

### Fixed
- **WSL backend install failed on a freshly installed Linux distribution.** The micromamba
  bootstrap assumed `curl` and `bzip2` were present, but a fresh Ubuntu WSL ships neither, so the
  download silently failed and every dependent step (RepeatMasker, minimap2, Dfam) then reported
  *"micromamba required first."* Installation now tries, in order: a micromamba already present on
  the machine (reused if another tool installed one), a download using only the Python 3 standard
  library (no `curl`/`bzip2`/`apt`/`sudo` needed), then `curl`/`wget`, then a passwordless `apt`,
  and finally a clear message naming the one command to run. The generated install script is forced
  to Unix line endings so it cannot break on a Windows checkout.

## [2.1.0] — 2026-07-20

UI and workflow refinements for wet-lab use, plus fixes found by a two-round debugging swarm.

### Added
- **Sortable, centered result tables.** Click any header to sort; sorting is numeric-aware for
  score, E-value, aa, divergence, and coordinates (composite `F/R` cells sort by the combined
  value), and the default view keeps the engine's order. All cells and headers are centered.
- **Organism dropdown for family annotation.** RepeatMasker's `-species` is now chosen from a list
  of common model organisms, with an **Other…** free-text field for any lineage; it auto-selects
  from a fetched accession's organism.
- **Send to splice detection.** Right-click any feature (structural, ORF, domain, family, genome
  viewer, amplicon) to send its sequence to the splice tool as the transcript. The splice card now
  states the **genomic reference** it aligns against (the loaded specimen).
- **Amplicon FASTA export** writes a file; figure export (gel / genome) now writes the **currently
  selected background mode** (dark / light / UV / mono) instead of always transparent.

### Fixed
- **Right-click actions used the wrong sequence.** Family-annotation hits (when run on a pasted
  sequence) and custom-background PCR amplicons copied/derived DNA by re-slicing the panel-01
  specimen at coordinates that indexed a *different* sequence. Copies, coordinates, primer design,
  and send-to-splice now use each row's actual source sequence.
- Result-card headers kept their ALL-CAPS style when expanded, and now show their one-line
  description (previously dropped).
- The in-silico-PCR amplicon table no longer lets a long *Source* value push the on/off-target
  *Call* column off-screen (the source elides, full value on hover).
- The header tagline reads **TRANSPOSABLE ELEMENTS ASSAY TERMINAL**; the "copied to clipboard"
  message no longer sticks in the status area (brief cursor tooltip instead).
- **App icon quality.** The `.ico` embedded only a 16 px frame (Windows upscaled it); it now
  carries every size up to 256 px, rendered cleanly.

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

[2.4.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.4.0
[2.3.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.3.0
[2.2.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.2.0
[2.1.1]: https://github.com/tunabirgun/TEagle/releases/tag/v2.1.1
[2.1.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.1.0
[2.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.0.0
[1.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v1.0.0
