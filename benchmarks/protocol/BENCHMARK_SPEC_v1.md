# TEagle benchmark specification v1

For a tool/evaluation paper. Target: BMC Bioinformatics or Bioinformatics Advances.
Written against TEagle **3.4.0** (commit `71f7f5b`, tag `v3.4.0`).

> **Amendment 2 — 2026-08-06.** The pin moved again, to **3.7.0**. The pre-registered design,
> thresholds, metrics and gate conditions below are **unchanged**; nothing was revised after seeing a
> result. Three things are recorded here.
>
> 1. A classifier defect was corrected. A lone reverse transcriptase, unaccompanied by any of the coding
>    modules that distinguish a retrotransposon from another reverse-transcribing gene, no longer
>    establishes a class or an order. Telomerase and group II intron maturases both present that way, and
>    a poly-A tail on an mRNA deposit was being credited as a LINE signature. Negative-control false
>    positives fell from two to zero. This changes classification output and is the reason for the pin.
>
> 2. A data-integrity defect in the harness was found and corrected. Per-case output is named by corpus
>    row index, so an inserted row shifts every index below it and a forced re-run left the files written
>    under the old names in place. The directory held 232 files for 135 corpus rows, and the scorer read
>    all of them, so cases were counted more than once and records produced under two different decision
>    rules were scored together. **Figures reported before this correction are superseded.** Class- and
>    order-level accuracy are identical under the old and the new decision rules once the duplicates are
>    removed, so the change in those figures is attributable to the de-duplication and not to item 1.
>
> 3. The monotonicity gate in section 6 is reported as follows, since the release makes it partly
>    untestable. At order rank no tier produced an error, so there is no error to distribute and the gate
>    cannot be evaluated. At superfamily rank — the only rank at which this release errs — the ordering
>    does **not** hold: both errors carry the *High* label. Neither outcome caused any threshold or tier
>    boundary to be refitted, which would be tuning on the evaluation corpus.

> **Amendment 1 — 2026-08-05.** The version this specification pins moved to **3.5.0** during execution.
> The pre-registered design, thresholds, metrics and gate conditions below are **unchanged**; nothing was
> revised after seeing a result. Two things changed in the software, both recorded here so the pin is
> honest rather than quietly updated:
>
> 1. Running the reproducibility panel revealed that one detector threshold (`find_tir(max_tir)`) was
>    applied but not recorded in the provenance manifest. It is now sealed, which changes every manifest
>    hash and therefore required a version bump. **No detection, classification or primer behaviour
>    changed**; the divergence panel was re-executed under 3.5.0 and its 600 records compared against the
>    3.4.0 run to demonstrate this rather than assert it.
> 2. The classification corpus was found not to isolate elements within their deposited records, so the
>    classification, calibration, specificity and head-to-head panels could not be scored as originally
>    assembled. See `benchmarks/CORPUS_DEFECT.md`. This is a shortfall against the pre-registered plan and
>    is reported as one; the affected panels are not presented as completed.

---

## 1. The claim the benchmark has to support

Most TE-annotation tools return a call for every input. TEagle's design rule is that a result must not
claim more than its evidence supports, so it **withholds** a call when the evidence cannot carry one.

> **Thesis.** On inputs where the diagnostic evidence is absent or unreadable, TEagle declines while
> comparable tools assert — and the declined cases are enriched for cases where those assertions are wrong.

This is the paper's contribution. Everything in §5 is table stakes a reviewer checks off; §4 is the part
that is genuinely new, and it is also the part most likely to fail. **If it fails, report that it failed.**
A tool that declines too often is a real finding and a publishable one.

### 1.1 The refusal classes, and whether each is benchmarkable

| Refusal | Trigger | Independent ground truth? | Use |
|---|---|---|---|
| **Superfamily undetermined** | integrase and RT on opposite strands, or overlapping spans | Only if a curated element exists whose deposit has that arrangement AND whose superfamily is fixed by Dfam/RepBase membership or phylogeny | §4.1 — **verify supply before committing**; if real cases < 5, drop the class from the headline claim |
| **Sub-threshold terminal repeat** | LTR pair at 73–80% identity, below the acceptance floor | **Yes, strongest case.** Dfam/RepBase family membership fixes the truth; old insertions have diverged LTRs by construction | §4.2 — primary evidence for the thesis |
| **Domain scan did not run** | `scan_domains` raised | No — an internal robustness property, not a biological one | §6, robustness. **Not** in the comparison |

If only §4.2 yields real cases, narrow the thesis to terminal-repeat age-censoring rather than padding with
synthetic constructs. Reviewers discount synthetic evidence precisely for claims of this shape.

---

## 2. Corpus rules

1. **Public accessions only.** Every case is an NCBI nucleotide accession (or a Dfam family + assembly
   coordinates). A reviewer must be able to re-fetch every input.
2. **Ground truth is independent of TEagle.** The label comes from the primary literature or from Dfam/
   RepBase family membership, recorded *before* TEagle is run on the case.
   *This is a change of direction from `verification/validation_matrix.txt`, whose entries are "the tool's
   actual output reconciled with the cited source". That is fine for regression testing and circular for
   benchmarking. Those cases may be reused only after their labels are re-derived source-first.*
3. **One record per case**, with: accession, organism, coordinates if a subrange, expected class/order/
   superfamily, the citation establishing it, and the evidence type (phylogeny / library membership /
   experimental).
4. **Held-out.** No case used to develop or tune a detector may appear in the benchmark. The five bundled
   examples (`M11240`, `M12927`, `M80343`, `X01005`, `X05424`) are **excluded** — they are development data.
5. **Target n ≥ 20** per major class: LTR/Copia, LTR/Gypsy, LINE, TIR DNA transposon, ERV. A class that
   cannot reach 20 is reported as a count, never as a percentage.

---

## 3. Comparator protocol

**Primary comparator: TEsorter.** Same input granularity (one element, classified from protein-domain
architecture via HMMER), so the comparison is like-for-like.

- Installed in its **own** environment, not the shipped `te` env that powers the RepeatMasker/Dfam backend.
- Default parameters, stated version, no per-case tuning. Any deviation is reported.
- Both tools receive **byte-identical input** (the same fetched FASTA), run on the same machine.

**Reference standard, not a competitor: RepeatMasker + Dfam.** Homology to a curated library is a
*different task* from de novo structural inference. It is used to establish ground truth where literature
does not, and is labelled as such. Presenting it as a competitor would be an unfair comparison.

**Explicitly out of scope: LTR_retriever, EDTA, RepeatModeler.** Whole-genome tools; comparing them on
single loci would misrepresent them. Named in the discussion as out-of-domain, not benchmarked.

---

## 4. Headline panels (the thesis)

### 4.1 Refusal precision — undetermined superfamily
- **Input:** elements whose deposits give an unreadable integrase/RT arrangement, superfamily known independently.
- **Measure:** for each tool, the call made. For TEagle: refusal rate. For TEsorter: call rate and accuracy on the same cases.
- **The claim holds if** TEsorter's calls on this subset are materially less accurate than on the readable subset, while TEagle abstains.
- **Report regardless of outcome**, including "TEsorter was correct on n of m, so abstention cost information."

### 4.2 Age censoring — sub-threshold terminal repeats
- **Input:** LTR elements spanning a divergence gradient, family membership known from Dfam/RepBase.
- **Measure:** detection as a function of LTR–LTR identity; the identity at which each tool stops detecting.
- **Purpose:** quantify TEagle's ~80% acceptance floor and ~72% seeding limit as an explicit sensitivity boundary, and show the advisory near-miss row recovers information that silent rejection destroyed.
- **This is the strongest panel** — the ground truth is unambiguous and the effect is measurable.

### 4.3 Calibration — does stated confidence track correctness?
- **Input:** the whole corpus.
- **Measure:** accuracy stratified by TEagle's own confidence tier (High / Moderate / Candidate) and by completeness tier.
- **Purpose:** a tool that reports confidence must show it means something. Report accuracy per tier with Wilson 95% intervals.
- **A monotone relationship is the claim.** If it is not monotone, that is a finding about the tier logic and must be reported.

---

## 5. Table-stakes panels

| Panel | Input | Metric |
|---|---|---|
| **5.1 Classification accuracy** | full corpus | per-class accuracy vs ground truth, Wilson 95% CI; confusion matrix; head-to-head with TEsorter on identical input |
| **5.2 Structural detection** | elements with literature-stated LTR/TIR/TSD coordinates | boundary agreement (bp offset distribution), detection rate |
| **5.3 Specificity / negative controls** | host genes, satellites, low-complexity tracts, shuffled TE sequence | false-positive rate. **Shuffled controls preserve base composition**, so a hit is structure, not composition |
| **5.4 Domain architecture** | ERVs and full-length elements with known architecture | per-domain precision/recall against the 30-model panel; explicitly scoped to what the panel can model |
| **5.5 Primer design** | published TE-specific primer pairs | do published pairs pass TEagle's QC? does in-silico PCR predict the published amplicon size? |
| **5.6 Off-target scan** | pairs with known genomic copy number | recovered copy count vs published, on a downloaded assembly |
| **5.7 Provenance reproducibility** | any case, run twice | identical `manifestSha256` for identical input; manifest sufficient to reconstruct every threshold |

---

## 6. Robustness (no comparator)

Truncated elements (5′ and 3′), N-runs, IUPAC-degenerate input, nested insertions, solo LTRs, sequences
below detector minimums, and a failed domain scan. **Measure: does the tool degrade honestly** — the
correct outcome is an explicit limitation, not a confident wrong answer. Also: runtime vs sequence length,
and the 12-ORF search budget's effect on domain recall.

---

## 7. Statistics

- Proportions with **Wilson score 95% intervals** (not normal-approximation — n is small).
- Paired tool comparison on identical inputs: **McNemar's test** on discordant pairs.
- Boundary offsets: median and IQR, not mean ± SD (offsets are not normal).
- **No percentage is reported for n < 20.** Counts only, stated as counts.
- Every panel reports its n alongside its estimate.

---

## 8. Reproducibility package

- `benchmarks/corpus.tsv` — one row per case: accession, organism, coords, expected label, citation, evidence type
- `benchmarks/run_teagle.py`, `benchmarks/run_tesorter.sh` — runners writing **raw** tool output per case
- `benchmarks/raw/` — every tool's unmodified output, committed
- `benchmarks/analyse.py` — raw output → tables and figures; no hand-entered numbers anywhere
- `benchmarks/RESULTS.md` — generated, with tool versions, database versions and checksums
- A single `make benchmark` path from empty directory to final tables

**Every number in the paper is produced by executing the real tools on the real corpus and is traceable to
a file in `benchmarks/raw/`.** No figure is transcribed by hand or reported from memory.

---

## 9. Known threats to validity — to state in the paper, not hide

1. **Panel scope.** The 30-model Pfam panel is retroviral-tuned. Plant Copia GAG/PR/RNaseH are not modelled,
   so completeness tiers for plant elements are systematically conservative. Already observed on `D83003`.
2. **Single-sequence inference.** No consensus, no multiple alignment, no copy number. Nested insertions,
   solo LTRs, MITEs and non-autonomous derivatives are the known weak cases and must be reported as such.
3. **Corpus bias.** Elements with literature-established labels are the well-studied ones — easier than a
   random element from an unstudied genome. Accuracy here is an upper bound.
4. **Ground-truth circularity.** Where Dfam supplies the label and TEagle's optional backend uses Dfam, the
   comparison is not independent. Those cases are marked and analysed separately.
5. **Version pinning.** Every result is tied to a stated TEagle version (3.4.0 as originally pinned;
   3.5.0 after Amendment 1, with the divergence panel re-executed and compared), a stated TEsorter
   version, a stated Dfam
   version with checksum, and the 30-profile panel hash `57f2b7881f35`.

---

## 10. Build order

1. Confirm §4.1 case supply. **If < 5 real cases exist, narrow the thesis to §4.2 now**, before collecting.
2. Assemble `corpus.tsv` to n ≥ 20 per major class, labels source-first.
3. Install TEsorter in an isolated environment; record its version.
4. Write runners; execute; commit raw outputs.
5. Analyse; report what the data shows, including where TEagle loses.
