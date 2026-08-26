# Screenshots

The present files were captured from the native TEagle 4.0.0 app in the light theme at 1920x1080, with a 1500x1150 window and UI scale 1.0. Every accepted figure is a real run: GenBank specimens (M11240 copia, M12927 gypsy, M80343 LINE, X01005 Tc1, X05424 hAT, AY037928 HERV-K113, AJ289709 HERV-H), a real NCBI coordinate fetch, a real insulin splice alignment (J00265 x NM_000207), and live WSL-backed runs (RepeatMasker/Dfam and isPcr against the cached *Drosophila* assembly). No figure contains synthesised results. The capture sequence was `python verification/recapture_shots.py` followed by `python verification/screenshot_runs.py`; the recapture sidecar and `screenshots/manifest.txt` bind the accepted files to their dimensions and SHA-256 hashes and record what each figure shows.

| File | Shows |
|---|---|
| `overview.png` | Main window: specimen rail, LTR/Copia classification banner, interactive genome viewer (terminal-repeat / cis-element / domain / ORF tracks), and the structural-evidence table with source-citation links. |
| `coord_fetch.png` | Panel 01 — coordinate fetch (organism + `chr:start-end`, UCSC-style 1-based) resolving chr13:33,016,423-33,066,143 on GRCh38.p14 to 49,721 bp. |
| `primers.png` | Panel 04 — Primer3 design with presets and advanced parameters; the ranked pair table with per-pair secondary-structure QC and a source citation. |
| `pcr_gel.png` | Panel 05 — staged in-silico PCR rendered as a to-scale agarose gel (MW ladder + sample lane) with the amplicon table. |
| `genome_scan.png` | Panel 06 — live whole-genome off-target scan (isPcr) of the designed pair against the cached *Drosophila* assembly: gel, verdict, and the full match table. |
| `installer.png` | The backend installer: each WSL/conda/Dfam component with a live status tick, a per-component Install (absent) or Repair (present) button, Install-all, Check-integrity, and a live log. |
| `screenshots/` | The full figure set used by the manual and the report — one per analysis outcome. See `screenshots/manifest.txt`. |
