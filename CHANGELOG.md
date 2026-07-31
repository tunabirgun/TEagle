# Changelog

All notable changes to TEagle are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The version is defined once in `app/backend/teagle_core/__init__.py` (`__version__`)
and propagates to the backend health endpoint, the UI header badge, every run
provenance manifest, the packaged executable's Windows file-version metadata, and
the LaTeX report title page.

## [3.3.0] — 2026-07-31

Two additions: an element's LTR now reports its polyadenylation-signal motif and a wider set of terminal motifs, and a downloaded genome can be annotated for its whole transposable-element landscape. Both are reported with the limits of the evidence stated on screen, because both are easy to over-read. **Results may differ from 3.2.1** for LTR elements: see *Fixed*, where a mis-placed detection window is corrected.

### Added
- **Whole-genome transposable-element annotation.** A genome already downloaded for the off-target scan can be run through RepeatMasker against the installed Dfam library, giving the families present, how many copies of each, how much of the assembly they cover, and how diverged those copies are — with XLSX/CSV/TSV table export, an SVG/PNG composition figure, a Markdown report and a JSON provenance manifest (GFF3/BED remains the per-record annotation export). The assembly is processed in contig chunks so progress is shown, a cancelled run keeps what finished, and a re-run resumes rather than restarting; the per-hit rows (millions, for a mammalian genome) are summarised inside the backend and never loaded into the window.
  - **Before the run starts** the app states the measured cost: cores, memory and free disk in WSL, the parallel job count it will use and whether cores or memory is the limit, the disk the run needs (derived from the genome's own size), a realistic time range, and — the part that decides whether the run is worth starting at all — **how many family models the installed Dfam partitions actually hold for that organism**.
  - **Transposable elements are reported separately from every other repeat.** Simple repeats, low-complexity sequence and satellites are tandem repeats, and rRNA/tRNA entries are not TEs; merging them into one "% repeats" would overstate TE content.
  - **You choose whether to search Dfam's uncurated families, and it matters more than any other setting.** RepeatMasker searches curated families only unless asked for both, and outside a few intensively studied species almost all families are uncurated. Measured on the same yeast assembly, changing nothing else: curated-only searches 9 families and reports **0.00% transposable elements**; including uncurated searches 421 more and reports **4.53%**, dominated by LTR/Copia — the Ty elements yeast actually carries. A run that finds nothing therefore says so as a statement about what was *searched*, not about the genome, and the library setting is printed with every result and sealed into its provenance.
  - **Library choice is the user's.** Either the installed Dfam partitions, or your own repeat library in FASTA, which RepeatMasker searches instead of Dfam (`-lib`). A custom library's checksum is recorded so the run is reproducible, and it is labelled user-supplied because TEagle cannot vouch for its contents.
  - **Sensitivity** is selectable — default, quick (`-q`) or slow (`-s`) — quoting RepeatMasker's own figures: quick is 2–5× faster and 5–10% less sensitive, slow is 2–3× slower and 0–5% more sensitive, where sensitivity means recovering older, more diverged copies.
- **The polyadenylation-signal motif of an LTR**, reported as advisory context. A hexamer from the published panel (AATAAA and its documented variants) counts only when a GU/U-rich downstream element follows it, because the hexamer alone occurs about once per 4 kb by chance and roughly half of genomic occurrences are never used. The same offset is checked in the other LTR copy: agreement corroborates, disagreement is shown as drift. It is a **motif and never a located cleavage site** — the U3–R–U5 boundaries and the transcript end need RNA evidence TEagle does not use, and the panel itself is derived from mammalian data, which is stated wherever the call appears.
- **Non-canonical LTR terminal motifs.** When the canonical TG…CA is absent, the seven non-TGCA termini LTR_retriever searches (TCCA, TGCT, TACA, TACT, TGGA, TATA, TGTA) are tested, exact-match only. It stays a badge that never moves a boundary — *gypsy* legitimately reads AG…TT, and snapping boundaries onto a motif would corrupt a correct call to satisfy a rule that element does not follow.
- **Upstream / downstream flanking sequence** can be taken from any feature — choose the side and the length, then copy it, export it as FASTA, or design primers on it. Offered for terminal repeats too, because amplifying *across* an insertion from its flanks is how an insertion is genotyped, and it is exactly the case where designing inside the element would be wrong. Sides are named in the record's orientation so their meaning does not silently change with strand.
- **Hover explanations** on the new scientific terms and parameters (sensitivity, library, TE % versus all-repeat %, family-model count, chunks, divergence), and context menus matched to what was clicked — a genome-landscape row offers aggregate actions and a plain-language explanation of divergence, not sequence actions it has no sequence for.
- **A BACKEND button in the window header** opens the backend installer from anywhere. It was previously reachable only from a secondary button inside panel 03, so a user who had not opened that card never met the thing that turns on family naming, splice detection and whole-genome scans. The in-card button stays; both open the same window.

### Performance
- **The application starts about twice as fast.** Its import graph drops from 537 ms to 151 ms, and a window is up in roughly 250 ms. Two costs were being paid on every launch for work most sessions never do: `openpyxl` (245 ms, and it pulls in numpy) was loaded for spreadsheet export, and the TLS certificate bundle (137 ms) was parsed even by a session that never opens a socket. Both now load on first use. A present-but-broken `openpyxl` still degrades to CSV exactly as an absent one does, rather than offering an export that fails.
- **Closing is clean.** The window stopped one named polling timer on close; it now stops every one, so a background poll can no longer fire into a window being destroyed.
- Interactive latency was measured and left alone: window build 30 ms, theme switch 64 ms, rendering a finished analysis 18 ms, self-similarity matrix 59 ms at 50 kb.
- **Downloads were measured, and one plausible optimisation was rejected on the evidence.** Dfam's server limits each client to about 2 MB/s, and splitting the 3.9 GB library across eight parallel connections measured *slower* (1.23 vs 1.81 MB/s), so no parallel download was added; the installer states the real time instead. NCBI, by contrast, is about twice as fast in parallel — recorded in `verification/perf_measurements.md` for a future change, since the genome download is one-time and cached.

### Fixed
- **The poly(A) downstream-element window was measured from the wrong point.** The hexamer lies 10–30 nt upstream of the cleavage site and the downstream element begins about 15 nt *after* that site, so relative to the hexamer it sits at roughly +25 to +60. The first implementation scored +10 to +30 — the gap between the hexamer and the cleavage site — and therefore judged the wrong stretch of sequence. Corrected to +20–60. Hits that only passed under the old window (measured: Ty1 and Ta1) are no longer reported.
- **The in-app methods panel claimed a 21-model Pfam panel** after the panel had grown to 30 — the application under-reporting its own method in a string no test covered. It is now derived from the profile table the scan loads, with a test that fails if the two ever disagree.
- **The genome-viewer legend no longer repeats its colours**; it reads them from the palette, so a legend can no longer drift from the bands it labels.
- **The optional uncurated Dfam partitions were unreachable from single-sequence family naming.** RepeatMasker reads the curated families only unless it is asked for both, and panel 03 never asked — so a user could download the 3.9 GB uncurated library, watch it install, and get exactly the same blank result, while the panel told them to install what they already had. Panel 03 now carries the same **Library** choice the whole-genome scan has, offered once the partitions are present, defaulting to curated-only as RepeatMasker itself does, and sealed into the run's provenance because the two settings answer different questions. Measured on this backend with yeast Ty1 (`M18706`): curated-only returns no hit at all; including the uncurated families returns an LTR/Copia family, which is what Ty1 is. A blank result now says which family set was actually searched, and the result header no longer claims "curated + uncurated" merely because those partitions are installed.
- **A false statement about what the curated Dfam partitions contain.** When a locus got no family match on a curated-only backend, the panel said "the curated library holds just 9 families for *Drosophila melanogaster*, and copia, gypsy, hobo and mdg1 are not among them". Measured against `famdb` on Dfam 4.0: curated-only *D. melanogaster* holds **399** family models and does contain `Copia_I`, `Copia_LTR`, `Gypsy_I`, `Gypsy_LTR`, `hobo` and `MDG1_I`/`MDG1_LTR`. The 9 was the yeast figure attached to the wrong organism. Coverage is a property of the lineage, not of Dfam as a whole — curated-only holds 1439 of 1439 models for human, 399 of 998 for *Drosophila*, 9 of 512 for *Arabidopsis* and 9 of 398 for yeast — so the panel now states that range from a measured table rather than one organism's number, and the cost estimate continues to report the live count for the organism actually selected.
- **Each opening of the backend installer left the window behind.** A dialog with a parent is owned by Qt and closing only hides it, so every open retained a dialog with its own worker pool — and the new header button makes reopening cheap. It is now destroyed on close.
- **A status probe that timed out replaced the progress line with a Python exception**, and a failed log read wiped the accumulated log. Both are transient — the probe runs several tools inside WSL and can time out while a download saturates the disk — so both now leave the display alone and wait for the next poll.
- **Only the component actually being installed shows the working marker.** The in-progress flag is global, so a single-component install marked every pending row as working. Which component is running is read back from the install script's own step markers, so it is right in a window opened *after* the run started — an install outlives the window that began it, and a reopened one cannot remember what a previous one clicked.
- **An integrity check finishing during an install no longer re-enables the buttons** while that install is still running.
- **A Dfam download showed no progress in the installer's log panel.** `curl` draws its progress meter with carriage returns and no newline, so a 40–60 minute transfer was a single line hundreds of kilobytes long: tailing the log by lines returned that one line, and the download read as frozen. The step now silences curl's meter and reports whole lines — bytes fetched, total, and percentage — every 15 seconds, with the total taken from a separate HEAD request because a resumed transfer reports only its *remaining* bytes. The log reader also collapses any carriage-return meter it still meets to the state that meter ended on, so a solver drawing one cannot flood the panel again.
- **An interrupted multi-gigabyte download is no longer discarded.** The checksum gate deletes a file that fails it, and a transfer cut short reached that gate, so an interrupted 3.9 GB download was deleted despite the panel promising it resumes where it stopped. Deleting now requires positive evidence that the transfer finished — the download exited cleanly, or the file reached its known size. Without that evidence the archive is still checked, so a complete file whose size a server never advertised is not rejected forever, but it is never deleted.
- **The free space a Dfam partition needs is now the space it actually needs.** The uncurated eukaryote partition was described as needing "~14 GB free while extracting" and gated on that figure; it unpacks to 22.6 GiB. A machine with 15 GB free passed the check, downloaded for an hour, then failed at decompression. The gate is now the measured unpacked size plus the archive plus headroom, and the figure the panel shows is computed from the same measurement, so the two cannot disagree.
- **A completed install could announce itself as stopped.** The install lock and the log are read by two separate calls into WSL at different moments, so a run finishing between the two read as "no lock, no completion marker" for one poll — and a 40-minute download that had just succeeded reported "the install stopped before it finished". Two consecutive observations are now required.
- **The two optional Dfam partitions offered "Repair" for a download that had never been run.** A component that is absent now reads **Install**; one that is present reads **Repair**, which is what re-running its idempotent step does. It no longer reads "Reinstall", which promised a fresh download the step does not perform on a file already on disk.
- **Opening the installer while an install was already running left its log blank** for the whole run, because polling started only when the user clicked something in that session. The window now attaches to a run in progress, and its status line no longer reports "ready" over the top of one.
- **The checksum and decompression stages announce themselves.** Both take minutes on a 3.9 GB library and both were silent.
- **A closed window no longer leaves a download reporting progress.** The install script kills its background transfer and progress watcher on exit, so an orphaned watcher cannot keep writing to a log whose run has ended.

### Verification
- Eight specimens across six organisms, each checked against the published or curator-annotated architecture: 6 of 8 match. Both mismatches are documented — a terminal repeat extended 1 bp past the annotated boundary in *Arabidopsis* Ta1-3, and an element embedded in host flanks (maize Ac9) where the terminal scan locks onto an internal repeat instead of the element's ends.
- A reference audit of that table corrected three sources before it could be trusted: a superseded accession (copia M11240 → X04456), an LTR length no publication states, and a TIR length cited to a paper that does not state it.
- Whole-genome runs: *D. melanogaster* 18.5% TE in 5.1 minutes with the family set the literature describes, and *S. cerevisiae* as the library-coverage control described above.

## [3.2.1] — 2026-07-30

A fix release: no new feature, one scientific-correctness fix and a set of robustness and documentation corrections found by a full review of 3.2.0. **A reported target-site duplication can differ from 3.2.0** for a flanked element whose superfamily has a literature target-site length and whose insertion happens to sit in a coincidental longer exact repeat — see *Fixed*. A bare pasted element is unaffected.
### Fixed
- **The target-site duplication is now measured with the classified superfamily's expected length preferred.** Structural detection runs before the superfamily is known, so it selected the longest exact flanking repeat; a coincidental longer repeat in an AT-rich flank could outrank the diagnostic short duplication — a Tc1/*mariner* element's 2 bp TA in particular — which then read as *incongruent* (discrediting genuinely complete termini) and exported the wrong duplication coordinates. Once the superfamily is resolved, the duplication is re-detected with its literature length preferred; the correction can only shorten a coincidental repeat to the diagnostic length when that length genuinely flanks, never lengthen or fabricate one.
- **Adding a custom organism no longer races on its store.** Two concurrent resolves could read the assembly store, each add its entry, and the second write clobber the first — an organism the app had reported as added would silently vanish. The read-modify-write is now serialised.
- **The "Add" control is disabled for the duration of its network resolve.** Only the text field was disabled before, so a double-click or Enter-then-click queued duplicate NCBI lookups; a single-flight guard now survives a manager close/reopen.
- **The "Resolving… against NCBI" status shows as a Notice, not an Error.** It had inherited the default error styling — an *Error*-titled, focus-stealing dialog for a routine in-progress message.
- **Documentation corrections.** The report cited the wrong *mariner* Mos1 accession (M14653 → X78906, the record actually benchmarked) and an unsubstantiated Tam3 specimen (removed). `app/README.md` pointed the provenance cross-reference at the wrong panel (05 → 07), listed a superfamily set that predated the 3.2.0 panel growth (now includes ERV, DIRS, CACTA, MULE and IS4-like/piggyBac), and omitted XLSX from the table-export formats.

## [3.2.0] — 2026-07-28

TEagle is now free software. The annotation can leave the tool as a genome-browser file, a first-time user opens a real published element instead of a synthetic null, a self-similarity view catches repeats the targeted detectors miss, and the primer secondary-structure cross-check no longer requires a manual ViennaRNA install. **Results may differ from 3.1.0**: the wider domain panel adds a DIRS-group class, and terminal-repeat lengths are now measured by diagonal extension — re-run any record whose classification or terminal-repeat length you have on file (see *Changed*).
### Added
- **The Pfam TE-domain panel grows from 21 to 30 models, widening coverage across all four TE classes.** LINE ORF1p is now a genuine multi-domain module (RNA-binding, trimerisation and dsRBD-like domains), the ORF2p apurinic-like endonuclease is tested (credited only in the EN–RT reading-frame arrangement), a tyrosine recombinase is modelled — which lets a **DIRS-group** element be called directly instead of misfiling as a LINE — and the Helitron helicase plus the CACTA, MULE and IS4-like/piggyBac transposases are added. Each carries a completeness ledger scoped to its class (a DIRS call is scored against RT + tyrosine recombinase + Gag, never the LTR ledger's DDE integrase it structurally lacks). The target-site duplication length is checked for congruence with the called superfamily's expected duplication.
- **GFF3 / BED export of an annotation.** The element, its terminal repeats (LTR/TIR), its cis-elements (TSD, PBS, PPT) and its protein domains export with verified Sequence Ontology terms and accessions, coordinates locus-relative unless a genome-anchored seqid is supplied, with the input sequence embedded after `##FASTA` so the file stands alone in a browser. Column 3 is gated on the completeness tier: a call whose defining evidence was not established degrades to a generic term rather than asserting a coding subclass the run does not support.
- **A self-similarity dot plot,** in its own resizable window. Exact k-mer matches — forward on the diagonals (direct repeats), reverse-complement on the anti-diagonals (inverted repeats) — surface a repeat of any arrangement without the scan being told what to expect, so a block the targeted structural detectors capped or missed still shows. A present diagonal is positive evidence; an absent one is not (a diverged repeat leaves a faint diagonal or none), and the panel states the asymmetry. The window opens fitted to its contents and supports zoom (buttons or Ctrl + wheel) with pan; the direct- and inverted-repeat mark colours are user-selectable (defaulting to a colour-vision-safe Okabe–Ito pair, with a reset), and the figure exports to SVG, PNG and vector PDF with a one-line method caveat embedded in the file. Pure stdlib, so it adds nothing to the frozen bundle.
- **Bundled real, published TE examples.** Five canonical specimens (Copia M11240, Gypsy M12927, LINE-1 M80343, Tc1 X01005, Ac X05424) chosen to exercise every branch of the classifier, each carrying its published identity so the call can be checked against the literature. The previous synthetic demo had a random ORF that could not hit a Pfam profile, so a first-time user's first result was a structural-only null.
- **An optional ViennaRNA primer-QC engine installable from the app.** The backend installer offers ViennaRNA as a one-click WSL component in its own `teagle-vrna` environment (its licence forbids redistribution inside an AGPL work, so it is not bundled). The out-of-process worker mirrors the in-process computation exactly — same ViennaRNA minor version, same DNA Mathews-2004 parameters, same temperature and salt — so a primer reports the same ΔG whichever route supplied it.
### Changed
- **Relicensed from proprietary to the GNU Affero General Public License, version 3 or later.** `LICENSE`, `THIRD-PARTY-NOTICES.md`, the README badge and prose, and the packaged executable's Windows file-version metadata all carry the new terms.
- **A locus with reverse transcriptase and a tyrosine recombinase but no DDE integrase is now classified as a DIRS-group retroelement.** The tyrosine-recombinase model is new to the panel, so 3.1.0 — unable to see it — read the same locus as a LINE. Re-run any record whose class you have on file; a DIRS-group call now carries its own completeness ledger (RT + tyrosine recombinase + Gag) and is never treated as an ERV.
- **Terminal-repeat (LTR/TIR) lengths are measured by X-drop diagonal extension across the whole sequence**, replacing the seeding-window diagonal cluster. The seed only locates the repeat; its extent is then extended base by base, so a reported length is no longer capped by the seeding window. Reported lengths may differ from 3.1.0 for some records.

## [3.1.0] — 2026-07-26

**Displayed structural-completeness tiers move for DNA transposons.** A Class II element no longer earns the *intact / autonomous-consistent* tier from a transposase hit alone. Re-run any DNA-transposon record whose tier you have recorded from an earlier release: the value may differ, and the new value is the one backed by the evidence ledger shown beside it.
### Changed
- **DNA transposons are now scored against a real completeness ledger.** Until now the Class II branch declared the strongest phrase in the application — *intact / autonomous-consistent* — whenever a transposase co-occurred with an inverted-repeat pair, while its own bookkeeping was hardcoded to *expected = present = transposase, nothing missing*: the tier could not be contradicted by the evidence it claimed to rest on. The branch now reasons the way the retroelement branch already did. An autonomous cut-and-paste transposon requires its transposase **and** both terminal inverted repeats — the cis ends its own product binds and excises (Wicker 2007) — so the ledger carries all three, the two arms are credited only when they actually enclose the transposase, and the tier follows the ledger instead of preceding it. A record whose ends were never recovered (5′- or 3′-truncated, internally deleted, or with termini too diverged for the scan) now reads **partial (transposase present)** with both arms listed as not detected, where it previously read a bare *transposase present*. An inverted repeat sitting elsewhere in the record no longer counts as that element's ends. The tier stays a categorical architecture call scoped to the tested Pfam panel — not a claim of transposition competence — and the classification card states in words which of the three components was recovered and why.
- **A named GAG domain no longer renders as an unnamed one.** The gag capsid/matrix bands were drawn in `#7A7A7A`, only 5.63 CAM02-UCS from the `#888` fallback used for domains with no palette entry — and, both being achromatic, colour-vision deficiency moved that separation by 0.00. GAG now has its own hue, measured at 3.02:1 on the dark track and 5.61:1 on the light track, with a worst-case separation of 18.4 CAM02-UCS from every other domain hue (13.8 once the cis-element, ORF and gene-model bands that can share a render are included), under normal, deuteranomalous and protanomalous vision alike.

## [3.0.0] — 2026-07-26

The installer drops by four fifths, the interface is rebuilt on a single source of colour truth, and every place where colour alone, a silent omission, or an unconfirmed click was carrying scientific meaning now says what it means in words. No detection method, threshold, database, panel or hedge changed.
### Added
- **A build guard on the bundle contents.** The PyInstaller spec now walks the resolved module graph and fails the build if any top-level package outside the pinned dependency set (and its real transitive dependencies, and TEagle's own modules) has entered it. A single stray import edge is what put half a gigabyte in the last installer; it cannot pass silently again.
- **Keyboard focus is visible on every control.** The native port had lost the coverage the web UI got from `:focus-visible`; buttons, links, card headers and the primary action each carry a focus ring tuned so the box never shifts and the ring never lands on the glyphs. Tables show a current-cell cue, citation links are Tab-reachable (and only when they actually contain a link, so no dead tab stops), and the disclosure toggles, close buttons and unlabelled icon controls carry accessible names and descriptions.
- **A confirmation before a destructive genome delete**, naming the organism and the cost of getting it back (~1 GB re-download, several minutes), defaulting to No.
### Changed
- **The installer is 147.6 MiB, down from 776.8 MiB (−81%); `TEagle.exe` is 4.4 MB, down from 71.2 MB.** One `collect_submodules("openpyxl")` call was pulling in openpyxl's pandas-interop shim, and through it pandas, scipy, scikit-learn, torch, transformers, PIL, lxml and numpy — roughly 528 MB of packages the code never imports. The shim is excluded at the source and the optional accelerators with it. Nothing the science uses was removed: Primer3, pyhmmer, ViennaRNA, the Pfam TE-domain profiles, openpyxl's XLSX export, certifi and PySide6 all ship as before.
- **One source of colour truth.** Every chrome token and every figure hue is defined once in `theme.py`; `figures.py` derives its palettes from it instead of restating hexes, so a colour used in both the interface and a figure exists in exactly one place.
- **A type scale in whole pixels.** Qt rounds a fractional QSS size up to the next pixel, so the old 10.5 / 11.5 / 12.5 px steps rendered as 11 / 12 / 13 px — the same pixels described by two different numbers. The scale is now written as the integers it actually produces.
- **Spacing snapped onto its scale.** The spacing helper rounds its argument onto the declared 6 / 10 / 16 px ladder rather than passing arbitrary values through, so the rhythm is enforced instead of merely documented (no call site moves more than 2 px).
- **Loading, error and empty states read differently.** An in-flight operation, a failure, a no-result and a genuinely empty panel each have their own treatment, so "still working" is never mistaken for "nothing found".
- **The disk requirement is stated before the download starts.** Manage genomes and the download action both state up front that a mammalian genome needs ≥ 8 GB free disk in WSL, with extraction peaking near 4 GB (FASTA ~3 GB plus the 2bit) before the temporary files are cleared — instead of the shortfall surfacing as a failure part-way through.
- **A backend failure is announced in words**, not by colour alone.
- **Both export entry points propose the same filename** for a given table, so the button and the right-click no longer disagree about what the file is called.
### Fixed
- **The frozen-bundle self-test gate never actually ran.** `build_installer.ps1` invoked the packaged executable with the call operator, but the exe is GUI-subsystem, so the call returned immediately and `$LASTEXITCODE` was never set — the gate could not fail. It now waits on the process and reads its real exit code, so a broken bundle stops the build.
- **ENV domains are no longer drawn in the "unknown domain" grey.** The env glycoprotein bands had no palette entry and fell through to the grey fallback, next to GAG's grey — the two most distinguishable retroviral domains rendered nearly alike. ENV now has its own hue, picked by measured contrast on both the dark and light tracks and checked for separation under deuteranomaly and protanomaly.
- **The gel distinguishes on- from off-target by shape as well as hue.** Each band carries a redundant non-colour mark (on-target, off-target, single-primer artefact, neutral priming site), so the specificity call survives greyscale printing and colour-blind reading.
- **Feature labels too narrow to letter in place are surfaced, not dropped.** A feature only a few pixels wide previously lost its label silently — which is exactly what was hiding the PBS tRNA-identity hedge, the one label that must never disappear. Narrow features now get a leader line and a caption line beneath the track; if more captions exist than fit, the overflow is counted rather than discarded.
- **A poly-A/T tail is no longer drawn as a terminal repeat.** It gets its own genome-viewer track, so a LINE's tail is not read as an LTR/TIR-style repeat.
- **A self-priming product is its own call.** An F+F or R+R product falling inside the target window was being labelled on-target; it is an artefact, never the intended amplicon. On-target, off-target and single-primer are now disjoint and exhaustive, so a derived off-target count (total − on − single) can no longer go negative.
- **Card 06 states isPcr's real coordinate convention.** Product coordinates are reported exactly as isPcr gives them — 1-based inclusive — and the column header says so, so length reads as end − start + 1.

## [2.12.0] — 2026-07-25

Post-overhaul refinements from direct use: the tool now carries you to the next step, warnings are proper dialogs you can dismiss, the genome manager shows everything at a fixed comfortable size, and download failures read plainly instead of as a crash. Driven by a user-directed fix loop with per-round triple-reviewer verification.
### Added
- **The tool moves you to the next phase.** "Design primer here" now brings the primer panel up (with a busy cue) and lands you on it; "send to in-silico PCR" and "send to splice" scroll their panel into view. A routed design no longer leaves you looking at the previous region's stale table.
- **Progressive disclosure in the classification banner.** The full domain-panel methodology folds behind a **Scope and methods** toggle; the one-line completeness caveat stays visible, so the card is less crowded without hiding the honest limitation.
### Changed
- **Warnings and confirmations are dismissable dialogs**, centred over the window, instead of a banner wedged above the panels. Errors and warnings persist until closed; confirmations auto-dismiss without stealing focus; the message text is selectable so an accession or seal can be copied.
- **Manage genomes is a fixed-size panel** sized to show every organism, column and Download button with no sideways scroll — wide enough that names never truncate, at a comfortable row height that keeps every button fully legible.
- **The backend installer opens larger** so each component row is comfortably readable.
- **Genome downloads retry with exponential backoff** (five attempts, 5→10→20→40 s) so a transient NCBI rate-limit or dropped transfer recovers instead of exhausting a short fixed retry.
### Fixed
- A genome download or scan that fails now shows one clear **warning** with the real reason, not a red "unexpected error" crash dialog stacked on top of the status line.
- The light-theme dialog and confidence-chip text is darkened to meet **WCAG AA** (≥ 4.5:1) against its tinted background.
- The **Manage genomes Download/Delete button no longer renders clipped or squeezed** — the label is centred at a comfortable height in every row, at every UI scale.
- The genome-manager dialog no longer collapses to a sliver when it rebuilds while open (rescale, add-organism, or a download completing).

## [2.11.0] — 2026-07-24

A UI/UX overhaul that keeps the assay-terminal identity but strips the noise: a legible body font, a calmer light default, honest and traceable readouts, and a whole-genome scan you can actually find — driven by a multi-persona design swarm and a five-round verification loop.
### Added
- **In-card whole-genome off-target scan.** Card 06 now carries an organism picker, a designed-pair picker, and a primary **Scan whole genome** button, so the scan is reachable without a right-click (the right-click stays as a power shortcut). Every route calls one sealed handler, so the isPcr job is byte-identical.
- **Add an organism to the genome manager** by name *or* assembly accession. The tool resolves it once, pins the versioned RefSeq accession, and unions it into every organism dropdown — a user-added genome seals byte-identically to a curated one.
- **Record's-own transcript picker for splice detection.** A fetched, annotated record lists its own mRNAs; one click loads a transcript. Aligning a record's own transcript back to its locus is labelled a **consistency check (same annotation source)**, never independent confirmation.
- **Traceable specimen readouts.** RECORDS links to the source accession; STRUCTURAL EVIDENCE and ORFS are in-app links to their tables, captioned "detected de novo — not database-retrieved" (a heuristic call is never dressed up as a database record).
- **Export-table buttons** on the primer, in-silico PCR, and off-target-scan tables, alongside the existing figure/FASTA exports.
### Changed
- **Roboto** is the body/UI font for legibility; sequences, accessions, coordinates and numeric tables stay in Cascadia Mono.
- **Light mode is the first-run default**, and the theme choice is remembered.
- **UI scale applies live** — no restart. The old "Restart now" (which only closed the app) is removed.
- **Calmer panels.** The classification banner is a three-tier CALL / WHY / SCOPE read with the reliability caveat kept visible; card titles are sentence-case; spacing, text tiers and padding are unified.
- **Gene model vs splice** are labelled by provenance — "NCBI annotation" vs "de novo alignment" — with a reciprocal cross-link, so an annotation is never conflated with an independent measurement.
- **Manage Genomes** is a designed panel (coloured status, aligned columns, an add-organism row) instead of a bare table.
- **Right-click menus are contextual** — a structural motif (LTR/TIR/TSD/PBS/PPT) offers only copy actions, not primer design or splice routing.
### Fixed
- Decorative glyphs and emoji are swept from the UI to plain text, while scientific notation (5′→3′, ΔG, ≤/≥), the ■ colour-keys, the ‡ engine-disagreement marker and the 01–07 card badges are kept.
- The light-mode accent is darkened to meet WCAG AA contrast for button and link text; the "source" citation link no longer renders a missing-glyph box.
- Content clipping / horizontal overflow at high UI scale and narrow windows is fixed — the page body never scrolls sideways and wide tables scroll inside their own viewport.
- The export right-click no longer overlaps a submenu arrow (flattened to a single action).

## [2.10.0] — 2026-07-24

Endogenous retroviruses are now read as retroviruses, not host genes: an explicit spliced-env transcript architecture plus the LTR cis-elements — the answer to "why does my HERV-K show one exon and no introns."

### Added
- **Retroviral transcript architecture (ERV).** A HERV-K (or any ERV) no longer shows a single host-style "exon = CDS". The analysis card draws the retroviral transcript model: **env is expressed from a spliced subgenomic mRNA** (a short 5′ leader exon joined to the env exon), with the whole **gag–pro–pol span removed as one intron** and labelled a **frameshift-fused polyprotein** — not a set of host exons. Junctions are anchored on the LTR and protein-domain positions and are **explicitly approximate**; TEagle never fabricates a single-base donor/acceptor from motif guessing (the exact env intron cannot be recovered from proviral DNA alone). HML-2 *rec/np9* sub-splicing is noted. Model after Löwer 1995 / Magin 1999 / Schmitt 2015.
- **LTR cis-elements.** The structural layer now detects the **primer-binding site (PBS)** in the leader just 3′ of the 5′ LTR — matched to a primer-tRNA panel and reported with the priming tRNA (tRNA-Lys for HERV-K), **hedged when the match is diverged** as endogenised proviruses are — and the **polypurine tract (PPT)** abutting the 3′ LTR. Both render as their own genome-viewer track and rows in the structural table.

### Changed
- For an endogenous retrovirus the host-style **gene-model view is de-emphasised** (collapsed behind a toggle) so the retroviral transcript architecture is the primary picture; for a non-ERV TE embedded in a host gene, the gene model stays visible.

### Fixed
- **Sub-panels follow the theme.** Every dialog — the secondary-structure detail, the sub-region picker, manage-genomes, the gel/PCR views, and the backend installer — now renders in the current dark/light theme instead of falling back to the operating-system default background.
- **The structural-evidence table shows both LTR copies.** An LTR (or TIR) pair now lists its 5′ and 3′ copy coordinates, matching the two blocks drawn in the genome viewer, rather than a single element-span row that read as one oversized repeat.
- **IDT reference reframed.** The Owczarzy 2008 (IDT SciTools) citation is labelled a comparability reference: TEagle matches OligoAnalyzer's ΔG convention and −9 kcal/mol threshold but computes every value independently with Primer3 + ViennaRNA and does not run IDT SciTools.

## [2.9.0] — 2026-07-24

Full retroviral GAG–POL–ENV domain coverage so endogenous retroviruses are read as complete elements, a two-axis reliability report (per-domain confidence + a scoped structural-completeness tier), a nine-specimen HERV benchmark across seven families, and the v2.8.x primer-QC polish.

### Added
- **Retroviral GAG and ENV domains.** The bundled Pfam profile set grows from 14 to 21 models (all CC0): gag matrix (PF02337), capsid (PF00607, PF19317) and nucleocapsid (PF14787), and env glycoprotein (PF13804), transmembrane (PF00517) and surface (PF00429). An endogenous retrovirus such as HERV-K now reports the full **GAG–PR–RT–RNaseH–INT–ENV** architecture in genomic order, not just the *pol* enzymes; the previous set had no env model and only a retrotransposon-gag profile that HERV-K capsid diverges from. An env domain flanked by paired LTRs flags the element as an **ERV**.
- **Reliability on two independent axes (no fabricated score).** A **per-domain call confidence** from the HMMER i-Evalue (high ≤ 1e-10 / moderate), shown in a new column; and a categorical **structural-completeness** tier — *intact / autonomous-consistent*, *near-complete*, *partial*, *structural-only* — mapped to the Wicker 2007, TEsorter and LTR_retriever criteria and always scoped to the models actually queried, so a missing-model result is never mistaken for element decay.
- **HERV benchmark.** Nine verified HERV proviruses across seven families (HERV-K/H/L/T/E/I/W) run through the tool: intact HML-2 proviruses recover the full architecture, while env-less HERV-L, gag-degraded HERV-H, and pseudogenised members are correctly reported partial with no spurious domains. Included as a report table and a reproducible test.

### Changed
- The classification card states the completeness tier, the domain architecture, and a plain note that a transposable element's coding organisation is its domain architecture (for an ERV, gag–pol–env), not a host exon–intron structure.
- The "Methods & databases" panel and the self-test now describe and verify the 21-model profile set; the self-test fails the build if a gag/env model is missing.
- Primer-QC polish (from v2.8.1): a visible legend for the ΔG colour flags and the ‡ engine-disagreement marker, per-theme flag colours that re-style on a theme toggle, a UI-scale "Restart now" option, and a sub-region picker that echoes the 0-based span.

### Fixed
- **Backend readiness verifies every pinned component.** The deep integrity check and the annotation gate could report an incomplete Dfam library as ready when only one of the two pinned partitions was present; both now require each pinned file individually. The installer status line claims splice detection or the whole-genome scan only when their tools (minimap2, isPcr + NCBI Datasets) are actually installed.
- **Honest completeness wording.** The structural-completeness tier is an architecture-level call, not a claim of functional or transposition/infection competence; *gag* counts toward an *intact* call only with a capsid or matrix domain, not the promiscuous nucleocapsid zinc-finger alone; a detected *env* (or other domain) is never reported as "no coding domain detected"; and the tier can no longer contradict the adjacent architecture / not-detected lines.
- **Fits small and high-DPI screens.** The window opens no larger than the display and keeps a usable minimum even at 125–150% interface scale, so it never overflows a 1366×768 laptop.

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

[3.1.0]: https://github.com/tunabirgun/TEagle/releases/tag/v3.1.0
[3.0.0]: https://github.com/tunabirgun/TEagle/releases/tag/v3.0.0
[2.12.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.12.0
[2.11.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.11.0
[2.10.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.10.0
[2.9.0]: https://github.com/tunabirgun/TEagle/releases/tag/v2.9.0
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
