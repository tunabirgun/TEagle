# Corpus scope and known limitations

Two corpora are used. This note records what each can and cannot support, so a reader of the raw output
reaches the same conclusions as a reader of the paper.

## The `record_scope` column

Both corpora carry a curated `record_scope` column with exactly two permitted values:

- `element` — the deposit IS the transposable element, or the negative-control gene, essentially in full.
  An element with short flanks is `element`.
- `containing_record` — the deposit is a clone, contig, chromosome, assembly, vector or genomic region
  that carries the element among other sequence.

The value is curated per row from that row's own evidence: the deposit's GenBank DEFINITION line, its
deposited length, the element coordinates the corpus or the record's feature table supplies, and the row's
`why_safe` note. It is not inferred at scoring time, because no mechanical test survives contact with the
records — four *Drosophila* P1 clones are titled "… DNA sequence (P1s …), complete sequence", and one Ty3
deposit leads with the tRNA gene it sits beside, so any keyword rule either keeps the clones or discards
the element. Length fails in the same way: element size and record size overlap across 0.4–20 kb.

`score.py` reads the column and aborts the run if a row is missing it or carries any other value. Nothing
defaults. Scope describes the DEPOSIT; what the engine actually analysed can be narrower, so `score.py`
also re-admits a `containing_record` case whose corpus coordinate was applied as a plain span before
analysis — the engine never saw the rest of the record in that case.

## Validation corpus (`corpus_holdout.tsv`) — 44 cases

44 distinct accessions, 36 species, no repeated records. All 44 rows are `element`, verified row by row:
the paper's claim that every entry is a deposit whose record IS the element holds, and no row contradicts
it. Three DEFINITION lines lead with host sequence and were checked individually rather than by keyword —
`M34549` (tRNA-Cys gene named first; the record's own feature table puts the two sigma LTRs at 121..460 and
5132..5471 of 5,510 bp, so the element is the record bar 159 bp of flank), `X17551` (white locus named
first; `mobile_element 1..4725` covers the whole deposit) and `AF435967` (a BAC isolate name in the title;
the deposit is the 8,911 bp provirus, not the clone). Labels are taken from the publication that described
each element and were verified against the deposited record before use. Eight entries are negative controls
chosen to carry a protein domain a TE classifier might match.

This is the corpus the paper's primary accuracy figures come from, and none of its numbers changed when
the column was introduced.

## Broader corpus (`corpus.tsv`) — 135 cases

135 cases across 122 distinct accessions and 50 species. Wider taxonomic and structural coverage, at the
cost of three limitations that the paper states and that any reuse should respect:

1. **Repeated records.** Sixteen cases share an accession with another case — eight elements within one
   maize bacterial artificial chromosome (`AF391808`), six within one maize contig (`AF123535`), and one
   pair (`AY037928`, once in the divergence panel and once in the ERV panel). Per-case counts are
   therefore not independent observations; the 95 cases that enter the primary figures come from 94
   distinct analysed inputs, and `scores.json` reports that count alongside the case count.
2. **The divergence-gradient panel cannot be scored as assembled.** 23 of its 25 rows are
   `containing_record` — BAC clones and contigs of 20–306 kb, and whole human chromosomes of 146–198 Mb,
   recorded without applicable element coordinates. Only `AY037928` and `AF164611`, both dedicated
   provirus deposits, are gradable. `extract_coords.py` recovers coordinates from the GenBank feature
   table where the record has one; the *Drosophila* `AC*` clones carry no `mobile_element` feature and
   would need the source publication's own table.
3. **Whole-record inputs.** 49 of the 135 rows are `containing_record`. Twelve of them supply a plain
   coordinate span that `run_teagle.py` applies before analysis, so the engine sees only the named
   sub-range and the case stays in the primary figures — six negative controls are deposited this way, as
   a chromosome or a whole bacterial genome carrying the control feature as a sub-range, and stratifying
   them out on deposit scope alone would have deleted the cases most able to expose a false positive.
   Four produced no output at all (fetch failures). The remaining 33 were analysed whole: the tool sees
   many elements at once and correctly declines to name a superfamily. Those 33 are stratified out of the
   primary figures and reported as their own stratum in `scores.json`
   (`stratum_containing_record`), not dropped.

   Stratified does not mean defective. The `refusal-supply` panel is a containing record BY DESIGN — its
   three rows exist to test what the tool does when handed the maize *adh1-F* nested-retrotransposon
   contig, and the answer, which is that it declines to name a superfamily, is the result the panel was
   built to obtain. All three are `containing_record`, so the panel now reports n = 0 in the accuracy
   table with its three cases in the stratum. That is the correct place for them: a whole-contig refusal
   is not comparable on the same scale as a call on an element deposit, and averaging the two would let a
   panel designed to elicit an abstention contribute to an accuracy figure.

## What changed when this was enforced

The stratification described in item 3 had never fired. `score.py` inferred scope from the DEFINITION line
using two keyword lists combined by OR, so the element list could only ever re-admit a row and never
exclude one, and the coordinate half of the test was vacuously true for every row carrying no coordinate —
101 of the 124 rows then being scored. 33 whole-record cases sat in the headline figures as a result.

Removing them moved the broader-corpus figures against the tool, which is the honest direction: class
accuracy 0.981 (95% CI 0.933–0.995) → 0.974 (0.910–0.993), order accuracy 0.990 (0.943–0.998) → 0.986
(0.922–0.997), with abstention rising from 0.148 to 0.181 (class) and 0.221 to 0.266 (order). None of the
33 carried an incorrect call — every one of them was answered correctly or abstained on — so the fall is
entirely the loss of easy correct answers from the numerator, and the intervals widen accordingly. The
`divergence-gradient` panel drops from 21 gradable cases at 1.00 to 2, which is the finding that matters:
its perfect score was built almost entirely on inputs with no discriminative power.

## A corpus defect corrected at the same time

`NC_000008.10` — human chromosome 8 — appeared twice, once labelled HERV-K115 and once HERV-K70. Neither
`chrN:` coordinate can be applied without an assembly, so both rows would have analysed the identical
146 Mb record under two different ground truths. The HERV-K115 row was removed: that element is already in
the corpus at element scope as `AY037929`, whose own DEFINITION line reads "Human endogenous retrovirus
K115 complete genome", so no ground truth was lost. The HERV-K70 row is kept. The removed row had produced
no output, so no reported number depended on it.

## Rule this establishes

A stratifier that removes cases from a benchmark must be justified against the code it claims to model,
and it must be able to fail. Length and open-reading-frame count both correlate with record size rather
than with whether a record is one element, so neither is a sound criterion — but the criterion that
replaced them was worse in a way that was harder to see, because it was written so that it could only ever
admit. A rule whose every branch says "keep" is not a rule, and it will pass review indefinitely while the
figures it is supposed to protect quietly include what it claims to exclude. The check that catches this
is not reading the criterion; it is counting how many rows it actually removed.
