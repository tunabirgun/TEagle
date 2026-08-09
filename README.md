<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/img/teagle-banner-dark.png">
  <img alt="TEagle" src="docs/img/teagle-banner-light.png" width="460">
</picture>

![Version](https://img.shields.io/badge/version-3.8.0-0A7259) ![Platform](https://img.shields.io/badge/platform-Windows%20x64-1FB89C) ![License](https://img.shields.io/badge/license-AGPL--3.0--or--later-0A7259) ![Built with](https://img.shields.io/badge/built%20with-PySide6%20%C2%B7%20Primer3%20%C2%B7%20HMMER-2B3740)
</div>

**TEagle** takes one transposable element and tells you what it is, how it is built, and how to amplify it — in one window, without a command line, and with every result carrying the exact software and database versions that produced it.

![TEagle analysis view — classification, interactive genome viewer, and structural evidence for a Drosophila copia element](docs/img/overview.png)

## Why it exists

Naming an element, reading its internal structure, and designing a PCR assay for it are normally three separate jobs across three tools, each with its own input format and none of which records how it reached its answer. TEagle joins them, and is built on one rule: **a result must not claim more than its evidence supports.** Where the evidence is thin, the tool says so on screen rather than rounding up to a confident-sounding answer.

It is written for someone who runs benches — not a bioinformatician, not a programmer.

## What it is not

Being clear about this saves time:

- **Not a genome assembler or a variant caller.** It analyses sequence you give it.
- **Not a discovery tool for new TE families.** It recognises families that are already in the library it searches. Finding families nobody has described yet is what RepeatModeler2 and EDTA are for.
- **Not a replacement for wet-lab validation.** Primers are designed and screened *in silico*; the gel it draws is a prediction.
- **Not a web service.** Everything runs on your machine. It fetches sequences from NCBI when you ask it to, and never uploads yours.
- **Not a claim of completeness.** A family absent from the installed library cannot be reported, so absence is never evidence of absence — and the panel says which library was searched.

## Scope

| | |
|---|---|
| **Input** | One element or locus: pasted sequence, a FASTA file, an NCBI accession, or genomic coordinates. Up to roughly tens of kb. |
| **Also** | A whole downloaded genome, for a transposable-element landscape (families, copy numbers, coverage). |
| **Platform** | Windows 10/11, 64-bit. The optional Linux-only steps run in a WSL2 environment the app installs for you. |
| **Offline** | Everything except NCBI fetch and the one-time backend/genome downloads. |

## Getting started

1. Download **`TEagle-Setup-<version>.exe`** (≈ 40 MiB) from the [Releases](../../releases) page.
2. Run it — per-user, no administrator rights — and launch **TEagle**.
   *Windows will show a blue **“Windows protected your PC”** SmartScreen warning: the installer is not code-signed, so Windows has no publisher reputation for it. Click **More info → Run anyway**. On a managed or institutional laptop this prompt may be blocked outright, in which case your IT team has to approve the executable.*
3. Click **Load example element** for a bundled, published TE (or paste a sequence, open a FASTA, or type an NCBI accession) and press **Run analysis**.

That is the whole first run. Everything the core needs — Python, Qt, Primer3, HMMER and the CC0 Pfam TE-domain profiles — is inside the installer; there is nothing to `pip install`. The optional Linux backend (Dfam family naming, splice detection, whole-genome scans) installs from inside the app: click **BACKEND** in the window header (or the same button inside panel 03), one click per component, with repair and integrity checks. The backend installer uses exact package versions and a checksum-pinned micromamba bootstrap; a repair rejects version drift and restores missing package files. Components that are not there yet offer **Install**; each reports its download live in the window's log.

The illustrated guide to every panel ships with each release as **`TEagle-User-Manual.pdf`**. The version-locked benchmark evidence is archived separately at [Zenodo](https://doi.org/10.5281/zenodo.21848814).

> **License:** TEagle is free software under the **GNU Affero General Public License, version 3 or later** — see [LICENSE](LICENSE). You may use, study, modify and redistribute it under those terms. Because it is AGPL, anyone who runs a modified TEagle as a network service must also offer its source to the users of that service. Third-party components and their licences are listed in [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

## How you use it

The window is a numbered sequence of panels. You work down it, and each panel fills in from the one above.

| | Panel | You do | You get |
|---|---|---|---|
| 01 | **Specimen** | Paste, upload, fetch an accession, or load an example → **Run analysis** | Length, GC%, validity, feature counts |
| 02 | **Classification & structure** | (fills in) | Superfamily call with the evidence behind it, an interactive map, terminal repeats, cis-elements, ORFs, protein domains |
| 03 | **Dfam family** | Set the species, choose the family library → **Run family annotation** *(needs the Linux backend)* | The community family name, with divergence and score, and which library was searched |
| — | **Splice detection** | Paste a transcript *(needs the backend)* | Exons and introns resolved by alignment |
| 04 | **Primer design** | Pick a preset or open the parameters → **Design primers** | Ranked pairs with Tm, GC, product size and a secondary-structure check |
| 05 | **In-silico PCR** | Load pairs → **Run** | A to-scale gel and an amplicon table |
| 06 | **Whole-genome scan** | Choose a downloaded genome | On- and off-target priming across the genome |
| 07 | **Run provenance** | (fills in) | The versions, checksums and parameters behind the result |

Right-click almost anything — a table cell or row, a feature in the map, a band on the gel — for actions that fit *that* item: copy the selected cell or row, copy its sequence or coordinates, take its flanks, design a primer there, or send it onward.

Separately, **Manage genomes → Annotate TE landscape** runs a whole downloaded genome and reports which TE families live in it and how much of it they occupy.

## What it does

Short version — the manual explains each in full.

**Reads an element.** Terminal repeats (LTR/TIR), target-site duplications, tails, ORFs, and protein domains from a 30-model Pfam panel spanning major retrotransposon and DNA-transposon groups. For an LTR element it also reads the cis-elements: primer-binding site, polypurine tract, the terminal motif, and the polyadenylation-signal motif.

**Calls a superfamily, with its reasons — and declines when it cannot.** Copia versus Gypsy from strand-aware integrase-vs-RT order, LINE, DIRS, ERV, and the DNA-transposon superfamilies — under the Wicker 2007 scheme, with the domains that support the call listed beside it. Where the evidence cannot carry a call the element is reported as undetermined rather than guessed: an integrase and RT on opposite strands, or on overlapping spans, leave the order unreadable, so no Copia/Gypsy call is made.

**States how much to trust it.** A per-domain confidence from the HMMER E-value, and a categorical domain-completeness class (*intact / near-complete / partial / structural-only*) — always scoped to the models actually tested. An endogenous retrovirus is read as a retrovirus, with env from a spliced subgenomic mRNA rather than a misleading host gene model.

**Names the family** against Dfam, and **resolves exons and introns** from a transcript (both optional, via the Linux backend).

**Designs and screens primers.** Primer3 with presets or full parameters, a secondary-structure check on every pair (cross-checked against a second engine when ViennaRNA is installed), pair-aware in-silico PCR drawn as a to-scale gel, and a whole-genome off-target scan against a genome you have downloaded.

**Annotates a whole genome's TE landscape** — which families, how many copies, how much of the assembly, how diverged — counting transposable elements separately from tandem and non-TE repeats, and telling you before the run whether the installed library can find TEs in your organism at all.

**Exports everything.** Tables to XLSX/CSV/TSV, figures to SVG/PNG (the self-similarity dot plot also exports vector PDF), sequence to FASTA, annotation to GFF3/BED with verified Sequence Ontology terms, and a provenance manifest sealing the versions, parameters and checksums behind the result.

![HERV-K endogenous retrovirus: the full GAG–PR–RT–RNaseH–INT–ENV architecture with per-domain confidence and intact domain completeness](docs/img/screenshots/herv_k_domains.png)

![Primer design with secondary-structure columns (hairpin / self-dimer / cross-dimer / 3′-end ΔG), colour-flagged, with the optional ViennaRNA cross-check installed](docs/img/primers.png)

![Whole-genome off-target scan: copia primers against the Drosophila genome give one on-target at the design locus and off-target paralogs, led by a specificity verdict](docs/img/genome_scan.png)

## Where to read more

- **`TEagle-User-Manual.pdf`** (ships with each release) — the illustrated guide: every panel, every option, what each number means and how to read it. Start here.
- **[Benchmark evidence](https://doi.org/10.5281/zenodo.21848814)** — the sanitized v3.8.0 corpora, aggregate results, provenance audit, method-audit scripts, tables and figures.
- **[CHANGELOG.md](CHANGELOG.md)** — what changed in each version, including any change that can move a reported value.

## Develop / build

```powershell
python app/teagle.py             # native window (first run auto-installs pinned deps)
python app/teagle.py --selftest  # headless bundle self-test (imports + QtSvg + a real analysis)
python -m pytest tests/ -q       # hermetic test suite (@wsl / @network sets gated separately)
powershell -File installer/build_installer.ps1   # freeze + bundle guard + self-test gate + Inno Setup
```

## Reproducibility

Every analysis records the database and package versions, parameters, and input checksums that produced it. Identical retained or refetched inputs analysed with the same versions and parameters produce the same seal. The seal excludes volatile fields (retrieval timestamps, unused-tool versions), and derived/advisory annotations (the primer QC, the on/off-target labelling) are recorded but kept out of the seal so they never change a result's identity.

The separately archived v3.8.0 benchmark evidence includes the measured ten-run seal-stability result, input-sensitivity checks and a compact audit of the recorded thresholds. The archive is a sanitized evidence package rather than a source release; exact reruns require the tagged TEagle source and refetching the cited public records.

## Benchmarks

The version-locked v3.8.0 benchmark evidence is deposited separately at [Zenodo](https://doi.org/10.5281/zenodo.21848814). It contains sanitized corpus definitions, scored aggregates, compact sequence-free provenance evidence, audit scripts, tables and figures. It excludes fetched records and sequence-bearing per-case detector output. Exact aggregate regeneration requires the tagged TEagle source, refetching the cited public accessions and rerunning the tools.

Measured on a held-out validation corpus under the documented bundled-example exclusion rule: 43 records
carrying 36 distinct organism annotations. Each label is sourced from the publication or authority that
described it, independently of TEagle. TEagle named the Wicker order for 34 of the 35 transposable elements
and was correct on all 34, and declined to name an order on all 8 negative controls. The broader benchmark
is a curated, non-random development corpus and does not estimate population-level accuracy or broad
utility. Its primary distinct-input summary contains 93 gradable content hashes at class and order rank:
TEagle answered 66 and all 66 were correct (answered-case accuracy 1.000; Wilson 95% CI 0.945–1.000).
Among 82 gradable element hashes, it answered 66 and withheld the order on 16; all 11 negative-control
hashes received no order call. The row-level sensitivity summary is 67/67 answered among 94 gradable
curated rows. At superfamily rank, 31 of 33 answered calls were correct among 59 gradable hashes
(answered-case accuracy 0.939; Wilson 95% CI 0.804–0.983). Compared against TEsorter 1.5.1 on
byte-identical input, neither tool made an incorrect order call on any of the 82 paired accession-level
inputs; they differ only in how often they commit. These figures describe this fixed curated evidence set.

The assay half is measured the same way, against 17 published assays from 8 sources whose product sizes
were determined at the bench. Every one returned a product on its cited template; of the 16 targeting a
single locus, all 16 fell within tolerance and 14 reproduced the published size exactly, with a median
absolute difference of 0 bp and all 16 binding with no mismatch in either primer.

The deposit documents what each corpus can and cannot support, including the curated sampling, repeated-accession structure, abstention denominators, paired-comparator exclusions and absence of prospective wet-lab validation.
