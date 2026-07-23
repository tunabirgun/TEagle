# Changelog

All notable changes to TEagle are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The version is defined once in `app/backend/teagle_core/__init__.py` (`__version__`)
and propagates to the backend health endpoint, the UI header badge, every run
provenance manifest, the packaged executable's Windows file-version metadata, and
the LaTeX report title page.

## [2.8.0] — 2026-07-23

Dual-engine, cross-checked primer secondary-structure analysis (hairpin / self-dimer / cross-dimer / 3′-end ΔG), a feature sub-region picker that routes only a chosen interval to primer design or splice, a global UI-scale setting with a collapsible specimen panel for small screens, and an explicit statement of every detection method — plus a whole-genome scan that reads a real on-target from the design locus.

### Added
- **Primer secondary-structure QC, cross-checked by two engines.** Every designed pair now carries hairpin, self-dimer, cross-(hetero-)dimer, and 3′-end anneal free energies (ΔG, kcal/mol) the way IDT OligoAnalyzer reports them. Because a single engine can disagree with the numbers seen elsewhere, TEagle computes each with two independent nearest-neighbor implementations — Primer3 (`thal`, SantaLucia 1998 parameters) as primary, cross-checked against ViennaRNA (RNAfold/RNAcofold, DNA parameters) at matched conditions — and shows both side by side, flagging a pair only when the two agree (amber ≤ −5, red ≤ −9 kcal/mol; the 3′ end weighted more strictly), with a ‡ marker and a caution when they diverge. A per-pair detail dialog shows the full breakdown. Validated against 12 published primer pairs (PrimerBank, Funakoshi 2017, Misak 2025): 11 of 12 carry no flag and the engines agree on 11 of 12. The screen is advisory — it never removes a designed pair, and its method references are recorded but kept out of the run seal.
- **Send a sub-region of a feature to primer design or splice.** Right-click any structural feature, ORF, protein domain, or amplicon and pick a coordinate sub-interval (e.g. bases 150–400 of the feature); only that subset is routed to primer design or splice detection, so you can target a specific part of a domain rather than the whole feature.
- **Explicit *methods and databases* in the classification card.** A one-click panel states exactly what defines each evidence layer: protein domains via HMMER against a bundled Pfam-A TE-domain profile set (the 14 accessions, ORFs ≥ 40 aa, E-value ≤ 1e-3); structural features via heuristic terminal-repeat detectors with their thresholds; superfamily via the Wicker 2007 scheme and the integrase-vs-RT order; and family naming via RepeatMasker + Dfam 4.0 — so the annotation is never a black box.
- **On-target vs off-target, from the design locus.** When the specimen was fetched at a known position in the scanned assembly, the whole-genome scan marks the product overlapping that position as the ON-TARGET and every other product as an off-target paralog, then leads with a specificity verdict over the split — *copy-specific*, *low-copy / paralogous*, *family-generic*, or *off-target-only*. A bare consensus with no genomic position yields neutral *genomic priming sites* rather than off-targets, so the gel, table, and verdict read neutrally. The scan renders in its own result card (06), separate from the in-silico PCR gel, and lists all products together, on-target first.

### Changed
- **Global UI scale and a collapsible specimen panel.** A new *⤢ Scale* control sets an overall interface scale (75–150%, persisted, applied on restart) so the whole window fits a small screen; the left specimen panel now collapses (Ctrl+B, or the header toggle) to give the results area full width and remove horizontal scrolling of the wider tables.
- **A co-migrating off-target is drawn as off-target.** When an on-target product shares a band size with one or more off-targets, a real gel cannot resolve them, so the band is drawn in the off-target colour as a specificity warning; the full on/off-target list stays in the table below the gel.

### Fixed
- **Duplicate section labels across successive runs.** Re-rendering a result card cleared only its widgets and left orphaned nested layouts behind, so headings such as *Structural evidence* and *Protein domains* could accumulate and appear twice when a second sequence was analysed. Card bodies are now cleared recursively, so each run renders exactly one of each section.

## [2.6.0] — 2026-07-23

Interpretable off-target results, honest backend health, progress on every long operation, IUPAC-degenerate primers, and a tabular genome manager — plus a critical fix that had left the cached-genome list empty.

### Added
- **Interpretable whole-genome off-target scan.** A scan now leads with a plain-language specificity verdict — *locus-specific*, *low-copy / paralogous*, or *family-generic* (expected for a TE-consensus pair) — computed over the true forward+reverse pair products, with a per-chromosome breakdown and a product-size cluster. Single-primer (F+F / R+R) artefacts are separated from real pair products, and the count is framed as a conservative floor (isPcr's ≥15 bp 3′-perfect rule does not count more-diverged copies). The genome-scan result table shows genome-specific columns (source, coordinates, length, strand, pair-vs-single-primer kind).
- **IUPAC-degenerate primer support in local in-silico PCR.** Consensus / wobble primers (R, Y, N, …) — standard for transposable-element work — now bind correctly against a plain-ACGT template, matching the ambiguity-aware genome-scan (isPcr) path instead of silently reporting the pair as non-binding.
- **Progress indicators on long operations.** Genome download, whole-genome scan, RepeatMasker family annotation, and minimap2 splice alignment now show an animated indeterminate progress bar, so a multi-minute backend call reads as working rather than hung.
- **Tabular Manage genomes panel.** The genome manager is now a table with one row per organism — assembly, accession, download status, cached size, contig count — with a per-row download / delete action.
- **Nested / composite element flag.** When a transposase domain co-occurs with reverse transcriptase (a nested or composite locus), the classification surfaces it as evidence and caps the confidence, instead of presenting a confident single-element call.

### Changed
- **The off-target scan organism menu lists only downloaded genomes.** Download an organism once from Manage genomes and it appears in the scan menu; this makes the "which genomes can I scan" state explicit and prevents a scan against a genome that is not on the machine.
- **Status banners carry a level.** Success, informational, warning, and error messages are now styled distinctly (a tick for success, amber for advisory, red only for real failures) — a completed download or scan is no longer shown as a red error with a warning triangle.
- **A single-exon splice result carries a caution.** A gapless alignment (0 introns) now always notes that it is consistent with either a genuine single-exon transcript or a genomic slice pasted as the transcript, so a common novice mistake is not read as a real biological finding.

### Fixed
- **Critical — the cached-genome list always came back empty.** `genome_list` ran its shell loop through an inline WSL command whose loop variable wsl.exe silently mangled to empty, so the whole-genome-scan organism menu and the Manage genomes panel reported zero cached genomes even when genomes were downloaded and present on disk. The loop is now delivered to the shell over STDIN, so cached genomes list correctly.
- **Genome download failed on a fresh WSL with "genome preparation failed: unzip".** A freshly installed Ubuntu WSL ships no `unzip`; extraction now uses python3's built-in `zipfile` (guaranteed present) with a system-`unzip` fallback. Verified end-to-end with a fresh *Drosophila* download.
- **The deep integrity check certified the backend healthy while the genome-scan tools were missing.** It now verifies isPcr and NCBI Datasets, so a "healthy" report is not immediately followed by a scan that fails at first use.
- **A broken isPcr binary was misattributed to a missing genome**, sending users to re-download a multi-gigabyte assembly for the wrong cause; the scan now reports a missing tool distinctly from a missing genome.
- **An edit to the sequence during an in-flight analysis could defeat the stale-sequence guard** and let primer design / PCR run against features indexed on the pre-edit sequence; the guard now compares against the analyzed snapshot, not the live box.
- **Manage-genomes row buttons could stay disabled** after a failed download, or after a scan that was started while the dialog was open.
- Corrected the "no timeouts" overstatement in the genome-scan documentation (a local safety timeout applies), the WSL-not-installed guidance (it now points at the in-app installer), and several banner-lifecycle and window-teardown-timer issues.

## [2.5.0] — 2026-07-22

A whole-genome off-target scan that runs entirely locally, primer design on flanks and gaps, automatic transcript-based exon/intron detection, and theme-following genome viewers.

### Added
- **Whole-genome off-target scan (local, no remote timeouts).** Right-click a designed primer pair → **⊕ Scan whole genome for off-targets**, pick an organism, and TEagle runs UCSC isPcr against that organism's RefSeq genome to report every candidate off-target priming site as a to-scale gel + coordinate table. The genome is downloaded once through NCBI Datasets and kept locally (a **⚙ Manage genomes** dialog lists cached genomes with sizes and lets you pre-download or delete them), so the first scan of an organism triggers a one-time download and every later scan is a fast local search — no remote query and no server-side queue. Validated end-to-end on the yeast, *Drosophila*, and human genomes. The scan is **reproducible and sealed**: the run provenance records the assembly accession *with version*, the source-genome SHA-256, the isPcr version, and the priming parameters, so an identical scan seals identically on any machine. Results are advisory — candidate priming sites under isPcr's ≥15 bp 3′-perfect-match rule, not wet-lab-validated amplicons.
- **Primer design on flanks and gaps.** The gene-model viewer now exposes the 5′/3′ flanking regions and interior gaps (neither exon nor intron) as clickable features — copy their FASTA or design primers there, the same as on exons, introns, and domains.
- **Automatic exon/intron detection with an annotation cross-check.** Splice detection aligns a transcript / cDNA / mRNA to the loaded genomic sequence (minimap2) and, when the specimen is a fetched record, reports an independent advisory cross-check of the alignment against the record's own feature-table annotation (matched / alignment-only / annotation-only introns).
- **Benchmarks.** Ten family-level naming specimens (Dfam / RepeatMasker) and de-novo splice / exon–intron benchmarks, each with the expected result verified against NCBI.
- **WSL backend components.** The managed Linux backend adds isPcr and NCBI Datasets (with a compact 2bit genome cache) for the local genome scan.

### Changed
- **Genome viewers follow the app theme.** Switching the app between dark and light now propagates to every open genome viewer and gel by default; a manual per-viewer background pick (including the gel's UV/mono modes) is kept and no longer reset by the next app-theme toggle. Pan and zoom are preserved.
- **Gene-model completion, honestly marked.** A coding exon the record implies through its CDS but omits from its exon annotation (e.g. the middle exon of insulin gene J00265) is now shown so the model is complete — rendered in a distinct lighter green with an `exon*` label and a legend entry, so a tool-inferred exon is never mistaken for a GenBank-annotated one.

### Fixed
- **WSL conda-cache recovery.** A corrupted or incompatible package-index shard (which made every environment solve fail with a repodata parse error) is now cleared with a full cache purge and the solve retried once, so installing the backend recovers instead of getting stuck.
- **HTTPS certificate verification.** certifi's CA bundle is bundled and used for NCBI / EBI / UCSC requests, fixing certificate-verification failures on Windows Python builds that lack a usable system trust store.

## [2.4.1] — 2026-07-21

Correctness and robustness release from a comprehensive multi-agent debug pass (three adversarial review loops with independent verification and advisor review). No feature changes.

### Fixed
- **Reproducibility seal.** The manifest hash is now invariant to which database served a fetched sequence (NCBI vs the ENA fallback) — two byte-identical runs of the same accession seal identically again. The in-silico splice manifest now seals the transcript as well as the genomic sequence, so two splice runs on the same locus with different transcripts no longer collide. In-silico PCR rejects a non-finite (NaN/Inf) target span instead of sealing it.
- **Specimen identity and stale state.** Editing a fetched sequence now clears the accession identity so a pasted edit is never sealed under the previous record; loading a new specimen clears the previous accession's gene model; and feature copy/design/send-to-splice now slice the *analysed* sequence rather than a later unanalysed edit, so copied bases always match the displayed coordinates. The splice "genomic reference" label refreshes when the specimen changes.
- **Classification and domains.** The aspartic-protease evidence line is now recorded; a 5′ poly-T tract is no longer mislabelled as a 3′ poly-A tail; and overlapping same-domain hits on opposite strands are no longer discarded (strand-aware de-duplication).
- **In-silico PCR.** A primer pair is only called on-target against the template it was actually designed on (never a false on-target from a pair designed on a different sequence); a multi-pair run still renders the gel for the successful lanes if one pair fails; and concurrent primer-design or fetch requests can no longer race.
- **Fetch robustness.** A transient NCBI response that is not valid JSON now surfaces as a clean, retry-suggesting message instead of an internal error; the ENA fallback also fires when NCBI raises a request error (not only on a non-FASTA body); a corrupt bundled assembly map falls through to a live resolve; and the served-database label links correctly for ENA-served records.
- **Miscellaneous.** RNA detection reads the sequence body only (a header containing the letter "U" no longer stamps a false RNA note); the genome viewer's coordinate mapping is correct at any panel width; a Primer3 environment fault is reported as an internal fault rather than a bad-input error; the local web server returns 400 (not 500) for a malformed Content-Length and tightens its static-file path guard; and the app version shown in the UI is taken from the single source of truth.

### Known limitations
- **WSL2-install-from-app path (unverified).** Three robustness issues on the elevated WSL2 installer — the completion marker can report a `wsl --install` failure as success, a non-ASCII Windows user path can be dropped from the install-log path, and a corrupted conda index shard is not auto-repaired on a Repair re-run — are identified but **not fixed in this release**, because this path cannot be verified without a WSL-less test machine (as noted for 2.4.0). The core app (classification, domains, primer design, in-silico PCR, coordinate/accession fetch) is unaffected and works fully offline.
- **Long LTRs (> ~1800 bp).** The terminal-repeat search window is capped at ~1800 bp, so an element with LTRs longer than that is reported with a truncated LTR length and an inward-shifted element span. This is a pre-existing limitation; a correct fix requires validation against long-LTR reference elements and is deferred.

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

[2.8.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.8.0
[2.6.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.6.0
[2.5.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.5.0
[2.4.1]: https://github.com/tunabirgun/TEagle/releases/tag/v2.4.1
[2.4.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.4.0
[2.3.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.3.0
[2.2.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.2.0
[2.1.1]: https://github.com/tunabirgun/TEagle/releases/tag/v2.1.1
[2.1.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.1.0
[2.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.0.0
[1.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v1.0.0
