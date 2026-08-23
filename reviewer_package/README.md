# TEagle — reviewer package

Everything needed to check a number in the manuscript, in one place. This is a curated copy of material
that also lives in `benchmarks/` and is assembled so that a reviewer does not have to navigate the working
tree. Nothing here is unique to this folder; it is a convenience, not a separate deposit.

TEagle 3.8.0. This is the current local reviewer evidence set after independent reconciliation. It is a
reviewer convenience set, not a release asset or independent deposit. The
benchmark evidence is deposited at https://doi.org/10.5281/zenodo.21848814. The exact reviewed application
source is identified by tag `v3.8.0` and commit `642b271e2cd4d496a5277f92dd4abe9d6e8578b3` at
https://github.com/tunabirgun/TEagle/releases/tag/v3.8.0. No repository-independent software-source DOI
has been assigned.

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
| `corpus.tsv` | the development corpus: 135 rows and 122 distinct GenBank accessions. Its `organism` field has 56 distinct literal annotations and 50 normalized primary labels after parenthetical qualifiers are removed; these are record-label counts, **not a taxonomic species count**. One row per case carries the expected class, order and superfamily, the publication each label was taken from, and the curated `record_scope` field. Repeated inputs remain a caveat: cases are not independent observations. |
| `corpus_holdout.tsv` | the development-example-excluded retrospective literature-derived corpus: 43 records carrying 36 distinct literal organism annotations, disjoint from `corpus.tsv` at accession level. Labels were transcribed from primary literature or other recorded authorities. The exclusion does not establish a prospective or preregistered holdout. `X01005` was removed because the specification lists it among bundled development examples; the superseded row and raw record are retained outside this package |
| `assay_corpus.json` | 17 published PCR assays from 8 sources, each with the primers as printed, the template accession, the published product size, and the passage in the source that supplies them |
| `CORPUS_SCOPE.md` | what the corpora do and do not support, including which panels are stratified out of the primary figures and why |

Every row records where its label came from. The corpora are built entirely from public GenBank
accessions, so they can be reconstructed independently of this repository.

### results/

Scored output, one file per experiment, plus the generated tables as CSV.

| file | experiment |
|---|---|
| `scores.json` | classification over the development corpus |
| `scores_holdout.json` | classification over the development-example-excluded retrospective literature-derived corpus |
| `head_to_head.json` | paired comparison against TEsorter 1.5.1, both database modes |
| `assay_scores.json` | amplicon prediction against published product sizes |
| `sim_divergence.json` | 600-run LTR divergence gradient |
| `reproducibility.json` | provenance seal: stability, sensitivity, threshold coverage |
| `diff_anchored_sites.json` | 6,000-trial matcher equivalence test |
| `chance_tir.json` | 4,500-trial chance inverted-repeat measurement |
| `testsuite.json` | test-suite counts |
| `provenance_index.json`, `provenance_summary.json` | compact, deterministic, sequence-free per-case provenance evidence and its derived accounting summary; source-manifest and input hashes permit audit without publishing raw manifests |
| `tool_capability_matrix.json` | the 16-tool, 16-axis capability survey, every cell sourced |
| `table1..4`, `tableS1` `.csv` | the manuscript tables, generated from the files above |

### scripts/

Figure, table and consumer scripts run from the supplied public aggregates. Commands that refetch accessions, rerun the tools or rebuild the compact provenance index require a matching environment and locally regenerated, ignored per-case outputs. Run commands from the repository root, not from this folder; they resolve paths relative to the project.

```
python benchmarks/sim_divergence.py     # divergence gradient
python benchmarks/figures.py            # Figure 2
python benchmarks/run_teagle.py         # classification corpus
python benchmarks/score.py              # scoring: class, order and superfamily
python benchmarks/run_tesorter.py       # the comparator, both database passes
python benchmarks/head_to_head.py       # paired comparison, McNemar exact
python benchmarks/run_assay.py          # amplicon prediction
python benchmarks/reproducibility.py    # provenance seal audit
python benchmarks/make_provenance_index.py  # index from locally regenerated ignored per-case manifests
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

**Per-case raw output** — the fetched GenBank records and the comparator's intermediate files, about 3.3 GB, of which a single human chromosome accounts for 326 MB. The cited accessions and runners permit regeneration, but future upstream retrieval can change text representation or record version. The compact index therefore records two deliberately different input hashes: `input_sha256` covers the fetch-level FASTA text, including its header and text representation, while `sealed_input_sha256` covers the analysed residue string in TEagle's content-addressed manifest. They are expected to differ and must not be equated. A regenerated fetch can be checked at the appropriate level.

**Published primer oligos** — the assay corpus and scored assay summary retain the short primer sequences transcribed from the cited publications because they are the inputs required to audit the reported in-silico PCR reconstruction. They do not contain fetched genomic records or per-case detector payloads.

**Sequence-bearing manifest content** â€” `provenance_index.json` intentionally excludes raw nucleotide
and amino-acid sequences, previews, motifs, primers, ORFs, domain payloads, structural coordinates,
free-text explanations, and citations. It is an audit extract, not a replacement for raw reruns or for
comparator intermediates. Its structural fields are bounded LTR/TIR/TSD measurements only; the summary
marks the historical orphan-file count as non-reconstructable from the current manifests.

**The benchmark specification and its amendments** are tracked at
`benchmarks/protocol/BENCHMARK_SPEC_v1.md` (mirrored at `report/BENCHMARK_SPEC_v1.md`). They are not
duplicated in this convenience directory.

## Evidence status and limitations

This is a single-locus evaluation on curated, non-random corpora; its case-level estimates are not
population-level performance estimates. Answered-case accuracy must always be read with the abstention
denominator because unanswered cases are excluded from that accuracy denominator. The assay panel is a
retrospective reconstruction of published primer assays, not prospective wet-lab validation. TEagle is
distributed under AGPL-3.0-or-later. The benchmark evidence has persistent public access; the application
is identified by its exact public tag and commit but does not yet have a repository-independent software-source DOI.

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
