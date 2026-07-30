<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/img/teagle-banner-dark.png">
  <img alt="TEagle" src="docs/img/teagle-banner-light.png" width="460">
</picture>

![Version](https://img.shields.io/badge/version-3.2.1-0A7259) ![Platform](https://img.shields.io/badge/platform-Windows%20x64-1FB89C) ![License](https://img.shields.io/badge/license-AGPL--3.0--or--later-0A7259) ![Built with](https://img.shields.io/badge/built%20with-PySide6%20%C2%B7%20Primer3%20%C2%B7%20HMMER-2B3740)
</div>

**TEagle** is a native Windows desktop tool that annotates transposable elements, reads their gene structure, and designs TE-aware PCR primers — all in one window, no command line, with every result reproducible from the exact database and software versions that produced it. The scientific core (structural detection, HMMER protein-domain scanning, superfamily classification, Primer3 design, primer secondary-structure QC, in-silico PCR, provenance) runs in-process; two optional features (Dfam family naming, de-novo splice) use a managed WSL backend the app installs for you.

![TEagle analysis view — classification, interactive genome viewer, and structural evidence for a Drosophila copia element](docs/img/overview.png)

> **License:** TEagle is free software under the **GNU Affero General Public License, version 3 or later** — see [LICENSE](LICENSE). You may use, study, modify and redistribute it under those terms. Because it is AGPL, anyone who runs a modified TEagle as a network service must also offer its source to the users of that service. Third-party components and their licences are listed in [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

## Install

1. Download **`TEagle-Setup-<version>.exe`** (≈ 40 MiB) from the [Releases](../../releases) page.
2. Run it (per-user, no admin) and launch **TEagle**.
3. Paste a sequence, open a FASTA, type an NCBI accession, or click **Load example element** for a bundled real, published TE → **Run analysis**.

Everything the core needs — Python, PySide6, Primer3, HMMER (pyhmmer), and the CC0 Pfam TE-domain profiles — is bundled. ViennaRNA is optional (its licence forbids redistribution inside a copyleft work); the backend installer can add it as a one-click WSL component, or install it yourself with pip/conda, to enable the second, independent primer-QC engine. Nothing to `pip install`, no command line. The optional Dfam/splice backend installs from within the app (**03 → Backend installer**, one click per component, with repair and integrity checks).

## What it does

- **Classify** an element from its structure (LTR/TIR/TSD/poly-A, ORFs) and its protein-domain architecture, into a superfamily under the Wicker 2007 scheme — Copia vs Gypsy by strand-aware integrase-vs-RT order, LINE, DNA transposons.
- **Coverage across all four TE classes** — the 30-model Pfam panel spans retroviral gag (matrix/capsid/nucleocapsid), pol (PR/RT/RNaseH/INT) and env, so an endogenous retrovirus (HERV-K, -W, -L, …) is read as a complete element and flagged as an ERV; the LINE modules ORF1p and the ORF2p apurinic-like endonuclease; a tyrosine recombinase for DIRS-group elements; the Helitron helicase; and transposases of the Tc1/Mariner, hAT, CACTA, MULE and IS4 groups.
- **Reliability, honestly** — a per-domain confidence (from the HMMER E-value) plus a categorical structural-completeness tier (*intact / near-complete / partial / structural-only*, after Wicker 2007 / TEsorter / LTR_retriever), always scoped to the models tested.
- **Retroviral transcript architecture, not a host gene model** — an endogenous retrovirus is read the way a retrovirus is expressed: env from a spliced subgenomic mRNA, the gag–pro–pol span drawn as the single frameshift-fused intron (junctions labelled approximate, never guessed from motifs), alongside the LTR **cis-elements** (primer-binding site, polypurine tract). The misleading host exon–intron view is de-emphasised for an ERV.
- **Cross-check the structure with a self-similarity dot plot** — in its own resizable window, exact k-mer matches (forward on the diagonals, reverse-complement on the anti-diagonals) surface direct and inverted repeats of any arrangement without being told what to expect, so a repeat the targeted detectors capped or missed still shows. A present diagonal is positive evidence; an absent one is not, and the panel says so. Zoom and pan; pick the mark colours; export to SVG, PNG or vector PDF.
- **Name the family** (optional, WSL) against the Dfam 4.0 library, and **resolve exon–intron structure** from a transcript with minimap2. The curated partitions install by default; most TE families of most organisms are *uncurated* in Dfam 4.0 — the curated set holds only nine families for *Drosophila melanogaster* — so the uncurated partitions ship as an optional component and the panel says which were searched.
- **Design and screen primers** — Primer3 with presets and full parameters; every pair carries a **secondary-structure check** (hairpin / self-dimer / cross-dimer / 3′-end ΔG) from Primer3's thermodynamics, cross-checked against a second independent engine when ViennaRNA is installed; pair-aware **in-silico PCR** as a to-scale gel; and a local **whole-genome off-target scan** against a downloaded RefSeq genome, with on-target/off-target framing.
- **Export for a genome browser** — the annotation leaves as GFF3 or BED with verified Sequence Ontology terms; column 3 is gated on the completeness tier, so a structural-only or evidence-incomplete call degrades to a generic term rather than asserting a coding subclass the run did not establish. Tables export to XLSX/CSV/TSV, figures to SVG/PNG, sequence to FASTA.
- **Reproducible by construction** — every result carries a provenance manifest sealing the exact tool/database versions, parameters, and input checksums; fetched sequences and genomes are content-addressed.

![HERV-K endogenous retrovirus: the full GAG–PR–RT–RNaseH–INT–ENV architecture with per-domain confidence and an intact structural-completeness tier](docs/img/screenshots/herv_k_domains.png)

![Primer design with secondary-structure columns (hairpin / self-dimer / cross-dimer / 3′-end ΔG), colour-flagged, with the optional ViennaRNA cross-check installed](docs/img/primers.png)

![Whole-genome off-target scan: copia primers against the Drosophila genome give one on-target at the design locus and off-target paralogs, led by a specificity verdict](docs/img/genome_scan.png)

## Full guide

The complete, illustrated user guide — every panel, every option, and how to read every result — ships as **`TEagle-User-Manual.pdf`** with each release. Fetch-by-coordinate, the backend installer, splice detection, the whole-genome scan, and the reproducibility record are all documented there.

## Develop / build

```powershell
python app/teagle.py             # native window (first run auto-installs pinned deps)
python app/teagle.py --selftest  # headless bundle self-test (imports + QtSvg + a real analysis)
python -m pytest tests/ -q       # test suite (500+ hermetic tests; @wsl/@network gated separately)
powershell -File installer/build_installer.ps1   # freeze + bundle guard + self-test gate + Inno Setup
```

## Reproducibility

Every analysis packs the databases and package versions plus input checksums that produced it, so a run reproduces byte-for-byte on another machine. The seal excludes volatile fields (retrieval timestamps, unused-tool versions), and derived/advisory annotations (the primer QC, the on/off-target labelling) are recorded but kept out of the seal so they never change a result's identity.
