# Try it yourself — five elements, five organisms

A guided set of representative workflows built around five published transposable-element accessions. Each primary accession is a public GenBank record, and each expected biological feature is stated in the cited paper or annotated in the record. A separate bundled synthetic construct is used only as a negative control.

The expected biological features below come from the cited publications or public records. They are interpretive checks, not a claim that the current release has reproduced every value. Clean-machine release validation is recorded separately, and this guide itself claims no benchmark or prospective wet-lab result. The TEagle 4.0.0 prospective benchmark is recorded separately ([results package](https://doi.org/10.5281/zenodo.22159977)) and applies only to that release. Protocol 2.3.6 results remain version-locked to TEagle 3.9.0. TEagle 4.1.0 has not been separately evaluated: the statements below about what this release reports are single-record checks, not an estimate of accuracy. Two boundary examples deliberately do not use the published length as an exact acceptance criterion; see §6.

Work down the list. Each specimen is chosen to drive a different branch of the tool, and together with the extras that follow they cover panels 01, 02 and 04–07, the 03b splice-structure card and the whole-genome mode. Panel 03's single-sequence Dfam family naming is the one step not walked through here; the whole-genome run in §5 searches the same Dfam library from the genome manager instead.

Button labels are written as they appear in the code; the application draws action buttons in capitals. The **Start here** strip above the results has one-click buttons for the bundled examples, **Help ▾** in the header opens the quick start and the glossary (F1), and hovering a table header gives its definition. The **Coords** button in the header switches the on-screen coordinates between 0-based half-open (the default) and 1-based inclusive; it changes the display only.

---

## 1. *Drosophila melanogaster* — **copia** — the canonical LTR retrotransposon

Accession: **`X04456`**  (5,146 bp. Note: `M11240` is the same element but NCBI has **replaced** it, so prefer X04456.)

**Run it:** type `X04456` into the accession box of panel 01 → **Fetch** → **Run analysis**.

| Expect | Value | Source |
|---|---|---|
| Class | LTR/Copia | Wicker 2007 scheme |
| LTR length | **276 bp**, two copies, ~100 % identity | Levis, Dunsmuir & Rubin 1980 *Cell* 21:581 |
| Terminal motif | **TG…CA**, canonical | integrase att convention |
| Superfamily logic | integrase **N-terminal to** RT | the diagnostic that separates Copia from Gypsy |
| Cis-elements | PPT present (~9 bp) | — |

Under the call, the Wicker rank line reads `Class I › LTR › Copia (RLC)` with the RepeatMasker class `LTR/Copia`, a display mapping of the call. This sequence is identical to the bundled `M11240` example, so a **Published as** note also sets the published class beside TEagle's call.

**What this tests:** accession fetch, structural detection, HMMER domain scan, superfamily classification by domain order, the genome viewer, and the provenance seal.

**Then also try:** right-click the LTR row, which copies either arm or the whole element and offers **Flanks of the whole element…** (the flanks are what you would amplify across to genotype an insertion; this record is the element alone, so there is no host flank to take, whereas `K01964` in §4 carries one). Press **Self-similarity plot…** — the two LTRs appear as short off-diagonal segments.

---

## 2. *Homo sapiens* — **HERV-K113** — an endogenous retrovirus, read as a retrovirus

Accession: **`AY037928`**  (9,472 bp)

**Run it:** fetch → **Run analysis**.

| Expect | Value |
|---|---|
| Flagged as | **ERV** (env-bearing LTR retroelement) |
| Domain architecture | **GAG–PR–RT–RNaseH–INT–ENV**, the complete retroviral set |
| Structural completeness | *intact / autonomous-consistent* |
| PBS | reported, but the priming tRNA **hedged as undetermined** — HERV-K113 matches the canonical tRNA-Lys3 template at only ~56 %, which is normal for an endogenised, diverged provirus |
| Transcript view | *env* from a spliced subgenomic mRNA, with gag–pro–pol drawn as one frameshift-fused intron — **not** a host exon–intron model |

TEagle 4.1.0 stops the Wicker rank line at the order for this element (`Class I › LTR › superfamily not resolved — env-bearing (RLX)`), because Gypsy, Bel-Pao and Retrovirus/ERV share the integrase-after-RT pol order and the bundled panel cannot separate them. A Gypsy or Copia label is therefore not expected here, and the call confidence is Moderate.

**What this tests:** the ERV branch, the full gag/pol/env panel, the two-axis reliability report, LTR cis-elements, and — importantly — that the tool **hedges rather than guesses**. A tool that named the priming tRNA confidently here would be overclaiming.

---

## 3. *Caenorhabditis elegans* — **Tc1** — a DNA transposon with terminal inverted repeats

Accession: **`X01005`**  (1,610 bp)

**Run it:** fetch → **Run analysis**.

| Expect | Value | Source |
|---|---|---|
| Class | DNA/Tc1-Mariner | — |
| TIR length | **54 bp**, perfect inverted repeat | Rosenzweig, Liao & Hirsh 1983 *Nucleic Acids Res* 11:4201 (and annotated in the record) |
| Evidence line | transposase enclosed by the inverted-repeat pair | — |

**What this tests:** the Class II branch, terminal-inverted-repeat detection (a different detector from LTR), and the completeness ledger for a cut-and-paste transposon.

---

## 4. *Zea mays* — **Activator (Ac)** — the classic maize element

Accession: **`X05424`**  (4,565 bp)

**Run it:** fetch → **Run analysis**.

| Expect | Value | Source |
|---|---|---|
| Class | DNA/hAT | — |
| TIR length | **11 bp**, imperfect terminal repetition | Pohlman, Fedoroff & Messing 1984 *Cell* 37:635 |

**What this tests:** the hAT branch and the fact that TIR detection is not tuned to one length.

> **Instructive variant.** Try **`K01964`** as well — the *same element* sequenced with its host *waxy* flanking DNA, so the terminal detectors must find the termini inside the record rather than at its ends. The record annotates the two 11 bp inverted-repeat copies at 122–132 and 4672–4682 (1-based) with 8 bp insertion-target copies beside them. TEagle 4.1.0 reports an 11 bp TIR at the same positions (shown as 121–132 and 4671–4682 in the default 0-based half-open numbering) and an 8 bp target-site duplication; TEagle 4.0.0 reported a 13 bp repeat at different coordinates for this record. This is one record, not an evaluation, and a different flanked record can still differ from its annotation: compare the coordinates your run reports with the record.

---

## 5. *Saccharomyces cerevisiae* — **Ty1** and the whole-genome mode

Accession: **`M18706`**  (5,918 bp, strain Ty1-H3)

**Run it (single element):** fetch → **Run analysis**.

| Expect | Value | Source |
|---|---|---|
| Class | LTR/Copia | — |
| LTR (delta) length | **334 bp** | Boeke, Eichinger, Castrillon & Fink 1988 *Mol Cell Biol* 8:1432; annotated in the record |

**Then run the genome — this is the most important test in the document.**

Panel 06 → **Manage genomes** → download *Saccharomyces cerevisiae* (small, ~12 Mb) → **Annotate TE landscape…**

Run it **twice**, changing only the **Library** setting:

| Library | Models available | Expected interpretation |
|---|---|---|
| Installed Dfam library — curated families only | Live lineage count shown by TEagle | The searched set can omit the yeast Ty families; a zero TE call is not evidence that the genome lacks TEs. |
| Installed Dfam library — include uncurated families | Live combined count; both shards required | Ty models become searchable; compare the matched repeat names/models, class/family labels, and interval-union coverage reported by your run. |

A curated-only zero means that the searched curated set produced no TE-classified match. Because *S. cerevisiae* carries Ty elements, that result must not be interpreted as biological absence. TEagle therefore shows the available model count before starting and records the searched library with the result. The known Ty1/Ty2 (LTR/Copia) and Ty3 (LTR/Gypsy) composition provides qualitative context rather than an exact numerical acceptance threshold.

The second run needs both optional uncurated Dfam partitions. Click **Backend** in the window header and install *Dfam uncurated · root (optional)* and *Dfam uncurated · Eukaryota (optional)*. The Eukaryota shard is the large download; the log reports its progress and an interrupted transfer can resume. TEagle blocks the curated-plus-uncurated run if either partition is absent; it never relabels a curated-only search. In the result table, **Repeat name/model** is RepeatMasker output column 10 and **Class/family** is the broader label in column 11. **Copies** counts raw alignment rows rather than insertion events. The TE count is the number of distinct TE repeat names/models. Mean divergence is weighted by the aligned spans of those raw rows, including overlaps and fragments; it is not a mutation rate or age. Every percentage divides by all submitted FASTA bases, including ambiguous bases. TE, tandem, and other unions may overlap and are not additive.

---

## Also worth trying

**Fetch by coordinate** (panel 01 → *Fetch by coordinate*): choose *Homo sapiens · GRCh38.p14*, paste `chr13:33,016,423-33,066,143` and press **Fetch region(s)**. Expect the panel to resolve chromosome 13 to `NC_000013.11` and report **49,721 bp** — the coordinates are 1-based inclusive, exactly as UCSC and NCBI display them, so numbers copied from a browser need no conversion.

**Primer design and in-silico PCR** (panels 04–05): with copia loaded, leave **Assay goal** on **Detect or quantify this element**, pick a domain, and press **Design primers**. Read the **Location** and **Notes** columns, select a pair, press **Send to in-silico PCR** (or **+ stage all designed** in panel 05), then **Run loaded pairs**. Primer3's target, included region, intended product, and primer coordinates are separate 0-based half-open quantities. The PCR table separately reports whether a physical product matches the intended endpoints and whether it covers the Primer3 target; without design-bound geometry those assessments are unavailable rather than negative. A non-templated 5′ tail contributes to the full product length shown for the predicted band but not to the template-coordinate span. To see multi-copy priming, paste the specimen sequence twice, back to back, into the **Optional custom background (FASTA)** box of the PCR panel and re-run the same pair: the pair now yields one product from the specimen and one from each pasted copy. (TEagle offers no *Design primer here* action on an LTR row, because a primer placed inside a repeat is non-specific by design.) The extra products are listed rather than hidden, and the table marks a product from a custom background as having no intended-product or target assessment; TEagle does not declare the pair specific. To check a pair of your own, type it under **Or check a primer pair you already have** (IUPAC letters, at least 18 bases each, an optional 5′ tail and name) and press **+ add pair**.

**Genotyping an insertion** (panel 04): **Assay goal → Genotype this insertion** designs across the 5′ junction, the 3′ junction, or flank to flank, and needs at least 18 bp of genomic flank on the side it uses. Copia from the accession above has no flank, so those options are disabled and the card says so; the flanked `K01964` record from §4 offers them. **Copy order rows** and **Export order sheet** list each unique oligo once, with empty Scale and Purification columns to fill in at the vendor.

**Several records and the assay record:** paste two FASTA records into panel 01 and run the analysis; click a row of the record table (or move with the arrow keys) to choose the record that panels 04 and 05 use. In panel 07, **Record this assay** offers **Copy notebook entry**, **Assay report ▾** (PDF or HTML) and **Export classification (TSV)**, and the **Runs** list shows each sealed run in sections.

**Export and invalidation checks:** in the native application, export a result table as XLSX, CSV, and TSV; export a figure as SVG, PNG, and PDF; export the record as FASTA (**Export full records (FASTA)**); and export the element annotation as GFF3 or BED (**Export annotation ▾**). A fetched or splice-derived gene model also exports as GFF3, BED, or FASTA. Then edit one base in the specimen: every dependent result and export must become stale until analysis is rerun. Manifest schema 0.2 and complete-analysis result schema 1.1 JSON are available from panel 07 (**Export analysis manifest (.json)** and **Export analysis result (.json)**) only after their bindings revalidate.

**Whole-genome off-target scan** (panel 06): scan the copia primers against a downloaded *Drosophila* genome with **Scan whole genome**. An accession-only specimen has no assembly coordinate, so its products are reported as neutral genomic priming sites. To exercise the on-target/off-target split, design from a specimen fetched by genomic coordinate in the same assembly.

**Splice detection** (needs the Linux backend; panel 03b, **Splice structure (transcript alignment)**): load `J00265` (human insulin), paste the `NM_000207` transcript sequence into the transcript box and press **Detect exons / introns** — expect **3 exons / 2 introns** with each junction checked against canonical GT–AG.

**Negative control:** *Load example element ▾* → **Synthetic construct (negative control)**, or the **negative control** button in the Start here strip. It has structural repeats but no real protein domains, so it should return a **Terminal direct-repeat candidate** with the order withheld and no superfamily call.

---

## 6. Boundary cases that do not use an exact-length acceptance rule

These two records are kept because an earlier release disagreed with the record's annotated repeat boundaries. They show what TEagle 4.1.0 reports in a single run for each record; the published length is still not an exact acceptance rule, and the coordinates in your own output are what to inspect.

**Arabidopsis Ta1-3 (`X13291`).** The record annotates the 5′ and 3′ LTR as 21–534 and 4725–5238 (1-based; **514 bp** each) with five-base target-site copies at 16–20 and 5239–5243. TEagle 4.1.0 reports a 514 bp LTR pair at the same positions (shown as 20–534 and 4724–5238 in the default 0-based half-open numbering), canonical TG…CA termini, and a 5 bp target-site duplication. TEagle 4.0.0 reported 515 bp, with both 3′ ends one base further and the canonical TG…CA badge lost.

**Maize Ac9 (`K01964`).** Described in §4 above.

A different record can still differ from its annotation by a base or more, so these two are not reduced to a pass/fail rule.

---

## What you should *not* expect

- **A Dfam family label without a matched backend search.** Superfamily classification works offline from protein domains; a label such as `Copia_I` or `L1HS` is reported only when RepeatMasker/Dfam in WSL returns a match.
- **A poly(A) signal to mean a located cleavage site.** When TEagle reports one it is a *motif* with its downstream element, labelled advisory. The U3–R–U5 boundaries and the transcript end need RNA evidence this tool does not use.
- **The priming tRNA to be named for a plant or fungal element.** The bundled panel is anchored on tRNA-Lys3; for *copia*, Ty1, Ta1 or Tnt1 the honest answer is "undetermined", and that is what you get.
- **A computational specificity result to establish wet-lab specificity.** TEagle reports the templates and criteria searched, the predicted products, and a scoped specificity verdict; experimental specificity still requires laboratory validation under the intended assay conditions.
