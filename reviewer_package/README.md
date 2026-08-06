# TEagle — reviewer package

Everything needed to check a number in the manuscript, in one place. This is a curated copy of material
that also lives in `benchmarks/` and is assembled so that a reviewer does not have to navigate the working
tree. Nothing here is unique to this folder; it is a convenience, not a separate deposit.

TEagle 3.6.0. Corpora and results correspond to the manuscript as submitted.

## What is here

```
corpora/    the benchmark inputs, as tab-separated tables and one JSON
results/    the scored output every figure, table and number in the paper is computed from
scripts/    the programs that produced results/ from corpora/
figures/    the manuscript figures, as PNG and SVG
```

### corpora/

| file | contents |
|---|---|
| `corpus.tsv` | the development corpus: 135 rows, 122 distinct GenBank accessions, 50 species. One row per case, carrying the expected class, order and superfamily, the publication each label was taken from, and the curated `record_scope` field |
| `corpus_holdout.tsv` | the independent validation corpus: 44 records from 36 species, disjoint from `corpus.tsv` at accession level |
| `assay_corpus.json` | 17 published PCR assays from 8 sources, each with the primers as printed, the template accession, the published product size, and the passage in the source that supplies them |
| `CORPUS_SCOPE.md` | what the corpora do and do not support, including which panels are stratified out of the primary figures and why |

Every row records where its label came from. The corpora are built entirely from public GenBank
accessions, so they can be reconstructed independently of this repository.

### results/

Scored output, one file per experiment, plus the generated tables as CSV.

| file | experiment |
|---|---|
| `scores.json` | classification over the development corpus |
| `scores_holdout.json` | classification over the independent validation corpus |
| `head_to_head.json` | paired comparison against TEsorter 1.5.1, both database modes |
| `assay_scores.json` | amplicon prediction against published product sizes |
| `sim_divergence.json` | 600-run LTR divergence gradient |
| `reproducibility.json` | provenance seal: stability, sensitivity, threshold coverage |
| `diff_anchored_sites.json` | 6,000-trial matcher equivalence test |
| `chance_tir.json` | 4,500-trial chance inverted-repeat measurement |
| `testsuite.json` | test-suite counts |
| `tool_capability_matrix.json` | the 16-tool, 16-axis capability survey, every cell sourced |
| `table1..4`, `tableS1` `.csv` | the manuscript tables, generated from the files above |

### scripts/

The programs that produce `results/` from `corpora/`. Run from the repository root, not from this folder —
they resolve paths relative to the project.

```
python benchmarks/sim_divergence.py     # divergence gradient
python benchmarks/figures.py            # Figure 2
python benchmarks/run_teagle.py         # classification corpus
python benchmarks/score.py              # scoring: class, order and superfamily
python benchmarks/run_tesorter.py       # the comparator, both database passes
python benchmarks/head_to_head.py       # paired comparison, McNemar exact
python benchmarks/run_assay.py          # amplicon prediction
python benchmarks/reproducibility.py    # provenance seal audit
python benchmarks/diff_anchored_sites.py  # matcher equivalence, fixed seed
python benchmarks/chance_tir.py         # chance inverted-repeat ceiling, fixed seed
python benchmarks/testsuite.py          # test suite, counts recorded like any other result
python benchmarks/make_table1.py        # tables
python benchmarks/make_table2.py
python benchmarks/make_table3.py
python benchmarks/make_tables.py
python benchmarks/fig_calibration.py    # Figure 3
python benchmarks/fig_workflow.py       # Figure 1
```

## What is deliberately not here

**Per-case raw output** — the fetched GenBank records and the comparator's intermediate files, about
3.3 GB, of which a single human chromosome accounts for 326 MB. Every byte is regenerable from the
accessions in `corpora/` by the scripts above. The input SHA-256 of each analysed sequence is recorded in
the scored output that *is* here, so a regenerated fetch can be confirmed as the sequence we analysed.

**The benchmark specification** (pre-registration, written against version 3.4.0 before these panels were
executed) is available from the corresponding author on request.

## Reading the numbers

Three things are worth knowing before comparing a figure here against the paper.

**Accuracy is on an answered-only denominator, at all three ranks.** Abstentions are counted in their own
column and never scored as errors, so accuracy and abstention have to be read together — a tool that
answers rarely can post a high accuracy on a shrinking denominator. Each file states its denominators.

**The three ranks do not share a denominator and are not nested.** Class and order are scored over the
same set; superfamily over a smaller one, which excludes rows whose superfamily this classifier cannot
represent and rows carrying no superfamily label. `scores.json` records the exclusions by reason.

**Containing records are stratified out.** A deposit that contains the element among other sequence is not
a single-element input. `corpus.tsv` records the distinction per row as a curated field; rows narrowed to a
curated coordinate span before analysis are retained, because the engine then sees one element.
`CORPUS_SCOPE.md` gives the rule and the counts.

## Licence

AGPL-3.0-or-later, as the repository. The corpora are derived from public GenBank records; each row cites
the publication its label came from.
