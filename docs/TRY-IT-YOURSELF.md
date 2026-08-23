# Try it yourself — five elements, five organisms

A guided set of representative workflows built around five published transposable-element accessions.
Each primary accession is a public GenBank record, and each expected biological feature is stated in the
cited paper or annotated in the record. A separate bundled synthetic construct is used only as a negative
control.

The expected biological features below come from the cited publications or public records. They are interpretive checks, not a claim that the 3.9.0 release has reproduced every value. Clean-machine release validation is recorded separately, and no new 3.9.0 classification, assay, or prospective benchmark result is claimed in this guide. Two boundary examples deliberately do *not* use the published length as an exact acceptance criterion; see §6.

Work down the list. Each specimen is chosen to drive a different branch of the tool, and together with
the extras that follow they cover panels 01, 02 and 04–07, the splice-detection card and the
whole-genome mode. Panel 03's single-sequence Dfam family naming is the one step not walked through
here; the whole-genome run in §5 searches the same Dfam library from the genome manager instead.

---

## 1. *Drosophila melanogaster* — **copia** — the canonical LTR retrotransposon

Accession: **`X04456`**  (5,146 bp. Note: `M11240` is the same element but NCBI has **replaced** it, so
prefer X04456.)

**Run it:** paste `X04456` into panel 01 → **FETCH** → **RUN ANALYSIS**.

| Expect | Value | Source |
|---|---|---|
| Class | LTR/Copia | Wicker 2007 scheme |
| LTR length | **276 bp**, two copies, ~100 % identity | Levis, Dunsmuir & Rubin 1980 *Cell* 21:581 |
| Terminal motif | **TG…CA**, canonical | integrase att convention |
| Superfamily logic | integrase **N-terminal to** RT | the diagnostic that separates Copia from Gypsy |
| Cis-elements | PPT present (~9 bp) | — |

**What this tests:** accession fetch, structural detection, HMMER domain scan, superfamily classification
by domain order, the genome viewer, and the provenance seal.

**Then also try:** right-click the LTR row → **Flanking sequence** (the flanks are what you would amplify
across to genotype an insertion), and **SELF-SIMILARITY PLOT** — the two LTRs appear as short off-diagonal
segments.

---

## 2. *Homo sapiens* — **HERV-K113** — an endogenous retrovirus, read as a retrovirus

Accession: **`AY037928`**  (9,472 bp)

**Run it:** fetch → **RUN ANALYSIS**.

| Expect | Value |
|---|---|
| Flagged as | **ERV** (env-bearing LTR retroelement) |
| Domain architecture | **GAG–PR–RT–RNaseH–INT–ENV**, the complete retroviral set |
| Domain completeness | *intact / autonomous-consistent* |
| PBS | reported, but the priming tRNA **hedged as undetermined** — HERV-K113 matches the canonical tRNA-Lys3 template at only ~56 %, which is normal for an endogenised, diverged provirus |
| Transcript view | *env* from a spliced subgenomic mRNA, with gag–pro–pol drawn as one frameshift-fused intron — **not** a host exon–intron model |

**What this tests:** the ERV branch, the full gag/pol/env panel, the two-axis reliability report, LTR
cis-elements, and — importantly — that the tool **hedges rather than guesses**. A tool that named the
priming tRNA confidently here would be overclaiming.

---

## 3. *Caenorhabditis elegans* — **Tc1** — a DNA transposon with terminal inverted repeats

Accession: **`X01005`**  (1,610 bp)

**Run it:** fetch → **RUN ANALYSIS**.

| Expect | Value | Source |
|---|---|---|
| Class | DNA/Tc1-Mariner | — |
| TIR length | **54 bp**, perfect inverted repeat | Rosenzweig, Liao & Hirsh 1983 *Nucleic Acids Res* 11:4201 (and annotated in the record) |
| Evidence line | transposase enclosed by the inverted-repeat pair | — |

**What this tests:** the Class II branch, terminal-inverted-repeat detection (a different detector from
LTR), and the completeness ledger for a cut-and-paste transposon.

---

## 4. *Zea mays* — **Activator (Ac)** — the classic maize element

Accession: **`X05424`**  (4,565 bp)

**Run it:** fetch → **RUN ANALYSIS**.

| Expect | Value | Source |
|---|---|---|
| Class | DNA/hAT | — |
| TIR length | **11 bp**, imperfect terminal repetition | Pohlman, Fedoroff & Messing 1984 *Cell* 37:635 |

**What this tests:** the hAT branch and the fact that TIR detection is not tuned to one length.

> **Instructive variant.** Try **`K01964`** as well — the *same element* sequenced with its host *waxy*
> flanking DNA. TEagle reports a **13 bp** repeat at different coordinates, not the annotated 11 bp
> termini. That is a real, documented limit: the terminal detectors are tuned for a record that *is* the
> element (how the tool is normally used); with an element embedded in host flanks the scan can lock onto
> an internal repeat instead. Seeing the tool's boundary honestly is more useful than seeing only its
> successes.

---

## 5. *Saccharomyces cerevisiae* — **Ty1** and the whole-genome mode

Accession: **`M18706`**  (5,918 bp, strain Ty1-H3)

**Run it (single element):** fetch → **RUN ANALYSIS**.

| Expect | Value | Source |
|---|---|---|
| Class | LTR/Copia | — |
| LTR (delta) length | **334 bp** | Boeke, Eichinger, Castrillon & Fink 1988 *Mol Cell Biol* 8:1432; annotated in the record |

**Then run the genome — this is the most important test in the document.**

Manage genomes → download *Saccharomyces cerevisiae* (small, ~12 Mb) → **Annotate TE landscape**.

Run it **twice**, changing only the **Library** setting:

| Library | Models available | Expected interpretation |
|---|---|---|
| Curated families only | Live lineage count shown by TEagle | The searched set can omit the yeast Ty families; a zero TE call is not evidence that the genome lacks TEs. |
| Include uncurated families | Live combined count; both shards required | Ty models become searchable; compare the matched repeat names/models, class/family labels, and interval-union coverage reported by your run. |

A curated-only zero means that the searched curated set produced no TE-classified match. Because
*S. cerevisiae* carries Ty elements, that result must not be interpreted as biological absence. TEagle
therefore shows the available model count before starting and records the searched library with the result.
The known Ty1/Ty2 (LTR/Copia) and Ty3 (LTR/Gypsy) composition provides qualitative context rather than an
exact numerical acceptance threshold.

The second run needs both optional uncurated Dfam partitions. Click **BACKEND** in the window header and
install *Dfam uncurated · root* and *Dfam uncurated · Eukaryota*. The Eukaryota shard is the large download;
the log reports its progress and an interrupted transfer can resume. TEagle blocks the
curated-plus-uncurated run if either partition is absent; it never relabels a curated-only search.
In the result table, **Repeat name/model** is RepeatMasker output column 10 and **Class/family** is the
broader label in column 11. **Copies** counts raw alignment rows rather than insertion events. The TE count
is the number of distinct TE repeat names/models. Mean divergence is weighted by the aligned spans of those
raw rows, including overlaps and fragments; it is not a mutation rate or age. Every percentage divides by
all submitted FASTA bases, including ambiguous bases. TE, tandem, and other unions may overlap and are not additive.

---

## Also worth trying

**Fetch by coordinate** (panel 01 → *Fetch by coordinate*): choose *Homo sapiens · GRCh38.p14* and paste
`chr13:33,016,423-33,066,143`. Expect the panel to resolve chromosome 13 to `NC_000013.11` and report
**49,721 bp** — the coordinates are 1-based inclusive, exactly as UCSC and NCBI display them, so numbers
copied from a browser need no conversion.

**Primer design and in-silico PCR** (panels 04–05): with copia loaded, design primers, load a pair, and run. Primer3's target, included region, intended product, and primer coordinates are separate 0-based half-open quantities. The PCR table separately reports whether a physical product matches the intended endpoints and whether it covers the Primer3 target; without design-bound geometry those assessments are unavailable rather than negative. A non-templated 5′ tail contributes to the full product length shown for the predicted band but not to the template-coordinate span. Then design a pair *inside* the LTR and re-run — because the LTR occurs twice, you should see it prime at **both** copies. TEagle reports the second physical product as an additional product rather than declaring the pair specific.

**Export and invalidation checks:** in the native application, export a result table as XLSX, CSV, and TSV; export a figure as SVG, PNG, and PDF; export the record as FASTA; and export the element annotation as GFF3 or BED. A fetched or splice-derived gene model also exports as GFF3, BED, or FASTA. Then edit one base in the specimen: every dependent result and export must become stale until analysis is rerun. Manifest schema 0.2 and complete-analysis result schema 1.0 JSON are available from the native Run provenance panel only after their bindings revalidate. The legacy browser deliberately exposes only CSV/TSV, SVG/PNG, FASTA, and gene/splice GFF3/BED/FASTA; provenance and complete-result JSON export are disabled there.

**Whole-genome off-target scan** (panel 06): scan the copia primers against a downloaded *Drosophila*
genome. An accession-only specimen has no assembly coordinate, so its products are reported as neutral
genomic priming sites. To exercise the on-target/off-target split, design from a specimen fetched by
genomic coordinate in the same assembly.

**Splice detection** (needs the Linux backend): load `J00265` (human insulin) and give it the transcript
`NM_000207` — expect **3 exons / 2 introns** with each junction checked against canonical GT–AG.

**Negative control:** *Load example element* → **Synthetic construct**. It has structural repeats but no
real protein domains, so it should return a **structural-only** result with no superfamily call.

---

## 6. Boundary cases that do not use an exact-length acceptance rule

**Arabidopsis Ta1-3 (`X13291`).** The record annotates the LTR as **514 bp**; TEagle reports **515**. The
start coordinates agree exactly; both 3′ ends run one base further, and that extra base also costs the
canonical TG…CA badge (the annotated boundary reads TG…**CA**, TEagle's reads TG…**AA**). The repeat
genuinely does extend one more matching base, so this is a boundary-convention difference, not a
misdetection. It is *not* fixed by snapping boundaries onto the motif, because *gypsy* legitimately reads
AG…TT — a rule that moved boundaries to satisfy TG…CA would corrupt a correct call to repair a rarer one.

**Maize Ac9 (`K01964`).** Described in §4 above.

These two examples expose boundary conventions that should be inspected in the current output rather than
reduced to an exact-length pass/fail rule.

---

## What you should *not* expect

- **A Dfam family label without a matched backend search.** Superfamily classification works offline from
  protein domains; a label such as `Copia_I` or `L1HS` is reported only when RepeatMasker/Dfam in WSL
  returns a match.
- **A poly(A) signal to mean a located cleavage site.** When TEagle reports one it is a *motif* with its
  downstream element, labelled advisory. The U3–R–U5 boundaries and the transcript end need RNA evidence
  this tool does not use.
- **The priming tRNA to be named for a plant or fungal element.** The bundled panel is anchored on
  tRNA-Lys3; for *copia*, Ty1, Ta1 or Tnt1 the honest answer is "undetermined", and that is what you get.
- **A computational specificity result to establish wet-lab specificity.** TEagle reports the templates
  and criteria searched, the predicted products, and a scoped specificity verdict; experimental specificity
  still requires laboratory validation under the intended assay conditions.
