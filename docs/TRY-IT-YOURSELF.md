# Try it yourself — five elements, five organisms

A guided test of every major function, using real, published transposable elements whose architecture is
described in the literature. Nothing here is synthetic: each accession is a GenBank record you can look up
independently, and each expected value is either stated in the cited paper or annotated in the record
itself.

**Every number in the "expect" columns was produced by TEagle on this machine and checked against the
source.** Two of them deliberately do *not* match the published figure, and those are the most informative
rows in the document — see §6.

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

| Library | Models available | Expect TE % | Families |
|---|---|---|---|
| Curated families only | 9 | **0.00 %** | 0 |
| Include uncurated families | 398 | **≈4.5 %** | 8, led by LTR/Copia ≈4.2 % and LTR/Gypsy ≈0.24 % |

*S. cerevisiae* genuinely carries Ty elements at a few percent of its genome, so **the first run is wrong
about the biology — and it completes normally while being wrong.** That is exactly why TEagle shows the
available model count before starting and states, when it finds nothing, that this reflects what was
searched rather than what the genome contains. LTR/Copia (Ty1/Ty2) dominating over LTR/Gypsy (Ty3) is the
known ordering for this genome.

The second run needs the optional uncurated Dfam partitions. Click **BACKEND** in the window header and
press **Install** on *Dfam uncurated · Eukaryota* — a 3.9 GB, ~40–60 minute download that reports its
progress in the log pane and resumes where it stopped if you close the app part-way.

---

## Also worth trying

**Fetch by coordinate** (panel 01 → *Fetch by coordinate*): choose *Homo sapiens · GRCh38.p14* and paste
`chr13:33,016,423-33,066,143`. Expect the panel to resolve chromosome 13 to `NC_000013.11` and report
**49,721 bp** — the coordinates are 1-based inclusive, exactly as UCSC and NCBI display them, so numbers
copied from a browser need no conversion.

**Primer design and in-silico PCR** (panels 04–05): with copia loaded, design primers, load a pair, and
run. Then design a pair *inside* the LTR and re-run — because the LTR occurs twice, you should see it
prime at **both** copies, which is the point: TEagle shows you the off-target product rather than
declaring the pair "specific".

**Whole-genome off-target scan** (panel 06): copia primers against the *Drosophila* genome give one
on-target product at the design locus plus off-target paralogues, led by a specificity verdict.

**Splice detection** (needs the Linux backend): load `J00265` (human insulin) and give it the transcript
`NM_000207` — expect **3 exons / 2 introns** with each junction checked against canonical GT–AG.

**Negative control:** *Load example element* → **Synthetic construct**. It has structural repeats but no
real protein domains, so it should return a **structural-only** result with no superfamily call. A tool
that invented a classification here would be telling you what you want to hear.

---

## 6. The two rows that do not match — read these

A benchmark where everything agrees has usually not been checked hard enough.

**Arabidopsis Ta1-3 (`X13291`).** The record annotates the LTR as **514 bp**; TEagle reports **515**. The
start coordinates agree exactly; both 3′ ends run one base further, and that extra base also costs the
canonical TG…CA badge (the annotated boundary reads TG…**CA**, TEagle's reads TG…**AA**). The repeat
genuinely does extend one more matching base, so this is a boundary-convention difference, not a
misdetection. It is *not* fixed by snapping boundaries onto the motif, because *gypsy* legitimately reads
AG…TT — a rule that moved boundaries to satisfy TG…CA would corrupt a correct call to repair a rarer one.

**Maize Ac9 (`K01964`).** Described in §4 above.

Both are documented in `verification/bench_cis_elements.md` with their full coordinates.

---

## What you should *not* expect

- **A family name without the Linux backend.** Superfamily classification works offline from protein
  domains; the community family name (`Copia_I`, `L1HS`) needs Dfam via WSL.
- **A poly(A) signal to mean a located cleavage site.** When TEagle reports one it is a *motif* with its
  downstream element, labelled advisory. The U3–R–U5 boundaries and the transcript end need RNA evidence
  this tool does not use.
- **The priming tRNA to be named for a plant or fungal element.** The bundled panel is anchored on
  tRNA-Lys3; for *copia*, Ty1, Ta1 or Tnt1 the honest answer is "undetermined", and that is what you get.
- **"Specific" primers.** TEagle never uses the word. It reports what was searched and what was found.
