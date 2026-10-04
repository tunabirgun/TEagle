# TEagle changes from 4.0.0 to 4.1.0: inventory for benchmark design

This report lists every behaviour change between tag v4.0.0 and version 4.1.0 (branch audit/4.1-usability), as found by reading every hunk of the source diff and checking the code. The same content is available as `CHANGES_4.0.0_to_4.1.0.json` (validated against `CHANGES_4.0.0_to_4.1.0.schema.json`) and, one row per change, as `CHANGES_4.0.0_to_4.1.0.tsv`. It contains 129 change records: 67 result changing, 28 output format, 18 display only, 15 infrastructure, 1 performance.

A record is *result changing* when a scientific output value or label can differ for the same input, *output format* when the science is the same but field names, columns or files differ, *display only* when only the screen changes, and *infrastructure* or *performance* otherwise. Each record gives the 4.0.0 and 4.1.0 behaviour, the inputs that trigger the difference, the affected outputs, the parameters involved, the tests that pin it, any measurement made during development, residual limitations, and a benchmark probe: a positive, a negative and a boundary case that would detect the change.

The DOI-bound 4.0.0 evaluation describes 4.0.0 only; none of its values transfers to 4.1.0. Four readers walked every hunk of git diff v4.0.0..HEAD by subsystem (classification and domains; structural detection; primers, in-silico PCR, genome scan, fetch and backend; outputs and interface) and checked the code; a catalogue reader mapped the three existing benchmark sets. The outputs/interface records were written from the CHANGELOG and the code. Pending at the time of writing: the measured 4.0.0-vs-4.1.0 behaviour diff over the benchmark caches and the Dfam consensus set, the endpoint-level impact mapping, and two independent checks of this inventory.

## What did not change

| Parameter | Value | Defined in |
|---|---|---|
| MIN_LTR_IDENTITY | 80.0 | `teagle_core/structural.py` |
| MIN_TIR_IDENTITY | 80 | `teagle_core/structural.py` |
| MIN_TIR_STANDALONE | 20 | `teagle_core/structural.py` |
| MIN_TSD / MAX_TSD | 2..12 bp | `teagle_core/structural.py` |
| domain E-value acceptance | < 0.001 | `teagle_core/domains.py` |
| MAX_ORFS_SCANNED | 12 | `teagle_core/domains.py` |
| Pfam panel | 30 models; te_domains.hmm byte-identical to 4.0.0 | `app/backend/data/te_domains.hmm` |
| LARGE_INPUT_BP (large-input warning) | 200000 | `app/backend/engine.py` |
| Primer3 Standard preset | 18-27 nt, Tm 57-63 C, GC 40-60%, product 150-500 bp, 5 candidates | `teagle_core/primers.py and the panel 04 preset` |

Schemas: analysis result 1.0 -> 1.1 (1.0 documents still validate); manifest 0.2; source context 1.0. The bundled HMM file is byte-identical, so the 30 model names and Pfam accessions are unchanged; what changed is how one of them (RVT_3) is interpreted.

Behaviours the readers checked and found unchanged:

- Bundled profile file te_domains.hmm is byte-identical to v4.0.0 (no diff); panel size 30; scan E-value cut-off 1e-3; MAX_ORFS_SCANNED 12; minimum ORF length 40 aa. Only the code and label of the RVT_3 model changed.
- classify._GAG_CORE {Gag_p24, Gag_p24_C, Gag_p10, PEG10_N-capsid}, _CODE_ORDER, _LINE_CORROBORATION {ORF1} and _RETRO_CORROBORATION {RNaseH, ENV, PR} are unchanged (though RNaseH now also receives RVT_3 hits).
- Copia and Gypsy canonical domain orders (GAG-PR-INT-RT-RNaseH-ENV and GAG-PR-RT-RNaseH-INT-CHR-ENV) and the rule INT before RT = Copia, INT after RT = Gypsy for a single same-strand non-overlapping integrase are unchanged; High still needs a Copia or Gypsy call with >= 2 domain codes and a resolvable order.
- The class (klass) vocabulary (Class I retrotransposon, Class II DNA transposon, unclassified), the confidence vocabulary (High, Moderate, Candidate) and the outcome vocabulary (complete, withheld) are unchanged; no te_class string was added or removed (DNA/unclassified is still emitted by the structure-aware 'DNA transposon' fallback).
- The five bundled examples keep their calls: M11240 LTR/Copia High, M12927 LTR/Gypsy High (env-like ORF is TLV_coat-type evidence, still not a withhold), M80343 LINE Moderate (tier 'intact / autonomous-consistent'), X01005 DNA/Tc1-Mariner Moderate (TIR 54, no TSD bare), X05424 DNA/hAT Moderate. Cached public elements AJ289709 (HERV-H, TLV_coat), X01472, M18706, X04456, K01964, U19263, X13291 keep te_class and confidence; AJ289709 (HERV-H) remains LTR/Gypsy High.
- HERV-K113 (AY037928) completeness (intact / autonomous-consistent, ERV env-bearing LTR retroelement), is_erv True and the retroviral transcript architecture (leader 968-1141, env 6792-8499, polyprotein intron) are unchanged; retroviral.py has no diff.
- MIN_TIR_STANDALONE 20 bp, MIN_TSD/MAX_TSD 2-12, LTR and TIR identity floors 80%, LTR seed window 6000 bp, and the hAT expected TSD length 8 bp are unchanged. CACTA, MULE, IS4-like and generic DDE calls still have no expected TSD length ('not diagnostic').
- No Bel-Pao, Retrovirus/ERV or LINE-superfamily (L1, CR1, Jockey) label exists in either version; the classifier still stops at order for those. Helitron remains domain-only: 'Helitron-family helicase present — superfamily not called', te_class unclassified, outcome complete.
- The lone-RT refusal (RT without LTR/INT and without GAG/PR/RNaseH/ENV/ORF1 or an EN-RT arrangement) is unchanged: 'reverse transcriptase, class unassigned', Candidate, withheld (but RVT_3 hits now count as RNase H corroboration, see CLS-003).
- The 'order' field is retained with its 4.0.0 semantics (None for LINE and unresolved calls); the new fields are additive.
- Calls on bare reference elements for tier names: every completeness tier and kind string of 4.0.0 is still emitted except where listed; two tiers and two kinds were added, none removed.
- Primer design, in-silico PCR, genome scan, splice and fetch behaviours are outside this scope and not assessed here.
- Frozen thresholds with values (v4.0.0 = 4.1.0, AST-identical): MIN_LTR_IDENTITY 80.0; MIN_TIR_IDENTITY 80; MIN_TIR_STANDALONE 20; MIN_TSD 2 and MAX_TSD 12; _LTR_SEED_WINDOW 6000; find_ltr k 13, min_ltr 80, min_anchors 4; _extend x_drop 30 and +1/-2 scoring; find_tir k 11, min_tir 10, max_tir 1000, min_anchors 3; _TIR_SEED_WINDOW value 1200 (was an inline literal, now named); find_polya min_run 8; find_pbs search 44, min_ident 0.55, strong 0.72; find_ppt window 30, min_len 9, min_purine 0.82, max_defects 2; find_polya_signal lo 20, hi 60, min_gt_frac 0.65, min_t 4. The 23 keys sealed in v4.0.0 have identical values in 4.1.0.
- Function bodies identical by AST comparison: _termini_motif, _CANONICAL_MOTIF, _NONCANONICAL_MOTIFS, _motif_mismatches, _PRIMER_TRNA (tRNA-Lys3 panel only; no initiator tRNA-Met added), find_ppt, find_polya_signal, _dse_state, _POLYA_HEXAMERS.
- LTR and TIR remain mutually exclusive in detect_all: a TIR is searched only when no LTR was accepted. The near-miss row type string 'Sub-threshold terminal direct repeat (advisory)' and its 80.0 threshold are unchanged (only the set of pairs that fall into it changes, STR-001).
- Bare specimens give identical structural rows in both versions: M11240.1 LTR 276 bp at 100.0 (span [0,5146]), M12927.1 LTR 482 bp at 100.0 (span [0,7469]), X01005.1 TIR 54 bp at 100.0, X05424.1 TIR 11 bp at 90.9, M80343.1 poly-A tail 27 bp at [6023,6050]; engine superfamily and confidence for these are unchanged. Forward PPT and poly(A)-signal coordinates on the copia and gypsy specimens are unchanged.
- Substitution-only LTR pairs are measured as before (300-bp pairs: 20/20 accepted at 5%, 19/20 at 20% in both; accepted counts differ by 1-2 of 60 only between 20% and 22% divergence).
- The anchored terminal-path TIR read of bare records is unchanged: synthetic bare TIRs (11/1 mismatch, 13/0, 20/2, 30/3, 15/3 mismatches; 20 seeds each) are found 20/20 with identical lengths and identities in both versions.
- Indels longer than 32 bp (deletion or insertion) still halve the LTR to the in-phase part, as in v4.0.0; LTRs of 6 kb or longer remain unseedable; the short/imperfect TIR joint window is 400 bp and a 1200-bp window governs 13-bp exact cores, neither of which widened.
- The random-sequence chance-TIR ceiling is unchanged in magnitude (17 bp across 22,500 sequences; 17 bp over 1,800 re-measured), and the 20-bp standalone floor is unchanged.
- PBS random-leader rate is unchanged: 23 of 200 random 300-nt leaders return a PBS row in both versions; the reporting floor and confident threshold did not move.
- Poly-A tail detection of a clean homopolymer at the record end, and poly-T on the reverse complement, is unchanged (L1.2 27 bp both; A8 and A20 bare tails 20/20 both). A tail followed by a flank with an 8-bp or 6-bp TSD is not called in either version; an uncorroborated tract in a flank (A14 in random flanks) is not a tail in either version.
- Lower-case input through the engine gives identical structural rows (the engine already upper-cases); only direct detector callers see the new upper-casing.
- TSD search order and range are unchanged (expected length first, then 12 down to 2); a TSD that is a (TA)n flank run is still emitted by detect_all, now with a chance_probability that reads as non-chance.
- The domain-panel, ORF and classification logic are not in this scope; the structural rows they consume change as listed, but no classification threshold in structural.py moved.
- Primer3 design defaults other than the product range are identical: opt_size 20, min_size 18, max_size 27, opt_tm 60, min_tm 57, max_tm 63, GC 40-60, num_return 5, max_poly_x 4, gc_clamp 0. The four GUI presets are identical (Standard 150-500 bp; qPCR 70-150, Tm 58-62, clamp 1, max size 24; High-specificity 150-500, Tm 60-64, size 20-26, GC 45-60, poly-X 3, clamp 2, 8 pairs; Permissive 100-1000, Tm 52-65, size 18-30, GC 30-70, poly-X 5, 10 pairs).
- Untargeted Standard-preset design on a record without ambiguity codes returns the same pairs, products, Tm and penalties (measured on M11240); Primer3's Tm model and condition values are unchanged; primer3-py, pyhmmer, build and backend requirement pins and the backend tool pins (RepeatMasker 4.2.4, Dfam 4.0 file list) show no diff.
- oligoqc thresholds (hairpin -2 caution / -9 warn, dimer -5 / -9, 3'-end -6 / -9 kcal/mol), QC conditions (50 mM Na+, 0 Mg2+, 0 dNTP, 250 nM, 37 C), the ViennaRNA cross-check, flag tiers and the hairpin, self-dimer and hetero-dimer values are unchanged (measured hetero-dimer -4.39 Primer3 / -3.14 ViennaRNA in both).
- In-silico PCR parameters and rules not listed as changes: max_mm 2, tp 5, min_anneal 18, product window 70-1000, MAX_SITES and MAX_AMPS 4000, intended endpoint tolerance 0, the IUPAC primer-base table, both-strand search, single-primer (F+F, R+R) products and amplicon ranking. For a single sequence with no background, coordinates, mismatches, tails and anneal lengths are identical; anchored_sites without a core argument returns the same sites (tests/test_anchored_equivalence.py).
- The specimen's lone-search rule is unchanged: per primer and strand the longest core wins, so a shorter imperfect copy is hidden by a longer perfect site (measured: the background alone as specimen still gives 2 products in both versions).
- Genome-scan engine: isPcr min_perfect 15, min_good 15, default product window 4000 (min 0), query row names pair/fwdonly/revonly, isPcr headers 1-based inclusive (confirmed live: the designed span maps onto the real product endpoints), the sealed content (assembly accession, source checksum, isPcr version and parameters) and its independence of any locus; the product set for the same pair and assembly is identical.
- Genome-scan primer input still goes through strip().upper() and a DNA/IUPAC regex, not normalize_primer: internal whitespace, U and '-' are rejected in both versions (measured) while lower case is accepted, which is inconsistent with run_pcr after PCR-004.
- RepeatMasker invocation (-qq, -pa, -species only when an organism is named, curated library by default, uncurated only when requested) and its hits are identical in both versions (measured on M11240 with a blank organism and with Drosophila melanogaster); Dfam curated-only family counts are unchanged (Homo sapiens 1439, Drosophila melanogaster 399 measured; the rest of the coverage table is not in the diff).
- minimap2 splice invocation (-x splice --secondary=no), the canonical splice-site set (GT-AG, GC-AG, AT-AC) and the exon and intron blocks are unchanged; only the strand logic and new fields change (GEN-010).
- Fetch: E-utilities endpoints, ENA fallback, coordinate semantics (1-based inclusive efetch seq_start/seq_stop), the 17 curated organisms and 26 bundled assemblies with their pinned accessions, acceptance of ordinary accessions (X04456, M11240.1, NC_003075.7) and the on-disk cache location are unchanged.
- Provenance canonicalisation (teagle-json-utf8-v1), manifest schema 0.2 and source-context schema 1.0 are unchanged, and 4.0.0 manifests still verify under 4.1.0 code; every primer, PCR (local-input), annotate and splice manifest hash nevertheless changes (INP-010, INP-011, PRM-001, GEN-006).
- The right-click 'Design primers here' on a domain row is unchanged (routed: the domain is Primer3's target and the product flanks and spans it); the Custom range goal over the whole record equals the old untargeted design; the engine API still accepts primer3_target and included_region exactly as in 4.0.0.
- benchmarks/score.py remains the frozen 4.0.0 scorer and benchmarks/diff_anchored_sites.py (matcher equivalence harness) is untouched; no 4.1.0 label in this scope (default-lineage species string, on-target tier names, template_ambiguity, annealing_core) is understood by the 4.0.0 scorer.

## Result-changing records

### Classification and domain panel

#### CLS-001 — Retrovirus-specific capsid/matrix or envelope profiles withhold the Gypsy call (RETROVIRAL_UNDETERMINED)

**Before (4.0.0).** In the RT + DDE integrase + accepted paired LTR branch, an integrase C-terminal to the RT on the same strand (non-overlapping) always gave superfamily 'Gypsy (Ty3)' (suffix ' · chromovirus' when a chromodomain was present), te_class LTR/Gypsy, and confidence High when >=2 domain codes and the INT/RT order was resolvable. The identity of the gag or env profile played no role, so HERV-K113 (HML-2) was LTR/Gypsy High and its explanation said 'endogenous retrovirus (ERV) / errantivirus lineage'.

**After (4.1.0).** If the same branch is reached with INT after RT and any detected hit has hmm name in {Gag_p24, Gag_p24_C, Gag_p10, HERV-K_env_2, GP41}, superfamily is the constant RETROVIRAL_UNDETERMINED, te_class LTR/unclassified, confidence Moderate (the Copia/Gypsy High rule needs a Copia or Gypsy superfamily), outcome stays 'complete', is_erv is unchanged (True when an ENV hit and LTRs are present), 'order' keeps the Gypsy-type canonical string (e.g. GAG–PR–RT–RNaseH–INT–ENV), wicker_order is LTR, completeness tier and kind are unchanged ('intact / autonomous-consistent', 'ERV (env-bearing LTR retroelement)' for HERV-K113), and the retroviral transcript architecture (env splice model) is unchanged. An evidence sentence says the Gypsy-type pol order is shared by Gypsy, Bel-Pao and Retrovirus/ERV so the superfamily is not called. Because the superfamily key no longer starts with Copia/Gypsy, no TSD-congruence sentence about LTR target sites is produced. TLV_coat, zf-CCHC_5, PEG10_N-capsid and Chromo do not trigger the withhold. A Copia-order call (INT before RT) is unaffected by a retroviral profile.

**Trigger.** Records with an accepted paired LTR row, an RT hit and a same-strand non-overlapping INT hit lying after the RT (translation order), plus at least one hit from the five profile names. In practice HML-2-like betaretroviral gag (Gag_p24/p24_C/p10) or a betaretroviral/lentiviral/gamma-like env hit by HERV-K_env_2 or GP41. Not triggered when INT precedes RT, when no LTR was accepted, when the only gag/env hits are TLV_coat, zf-CCHC_5 or PEG10_N-capsid, or when the INT is on another strand/overlaps the RT.

**Affected outputs.** classification.superfamily: 'Gypsy (Ty3)' -> RETROVIRAL_UNDETERMINED string for affected inputs; classification.te_class: LTR/Gypsy -> LTR/unclassified; classification.confidence: High -> Moderate; classification.explanation and evidence (new retrovirus-specific sentence; old Gypsy/Ty3 order sentence removed); classification.evidenceLedger statements; classification.wicker_order (new field, LTR); classification TSV / RepeatMasker class string / Wicker code and GFF3 TE parent attributes derived from te_class and superfamily; any scorer keyed on the superfamily string (the string is new to the frozen 4.0.0 scorer, which raises on it).

**Labels.** added: LTR retrotransposon (superfamily undetermined; retroviral-type gag/env — Gypsy, Bel-Pao and Retrovirus/ERV are not separated by the bundled panel); renamed: Gypsy (Ty3) -> RETROVIRAL_UNDETERMINED string (conditional re-label for inputs meeting the trigger; 'Gypsy (Ty3)' is still emitted for other inputs).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _RETROVIRAL_PROFILES (HMM names that withhold Gypsy) | absent | {Gag_p24, Gag_p24_C, Gag_p10, HERV-K_env_2, GP41} | no | Wicker et al. 2007 (Gypsy and Retrovirus/ERV are sibling superfamilies sharing the PR-RT-RH-INT pol order); TLV_coat excluded on purpose because it also annotates the env of insect Ty3/gypsy elements (errantiviruses) |
| _GAG_CORE (capsid/matrix profiles counted as core gag) | {Gag_p24, Gag_p24_C, Gag_p10, PEG10_N-capsid} | {Gag_p24, Gag_p24_C, Gag_p10, PEG10_N-capsid} | yes | unchanged set; the withhold set is a separate constant |

**Measured.** Public GenBank record AY037928 (HERV-K113, 9,472 bp; cached fetch), run through engine.analyze in separate processes, both strands: v4.0.0 LTR/Gypsy 'Gypsy (Ty3)' High; 4.1.0 LTR/unclassified, RETROVIRAL_UNDETERMINED, Moderate, outcome complete, is_erv True, order GAG–PR–RT–RNaseH–INT–ENV, completeness 'intact / autonomous-consistent', retroviral architecture identical (leader 968-1141, env 6792-8499, intron 1141-6792). Domains there: Gag_p10, Gag_p24, zf-CCHC_5, RVP, RVT_1, RNase_H, rve, HERV-K_env_2, GP41. Hand-built rows (LTR 0-968/8504-9472; GAG 1141-1390, PR, RT 4130-4613, RNaseH, INT 5804-6095, ENV 6792-8499) and each of the five profile names substituted into an RT+INT+LTR layout give the same switch (8 profile cases: 5 withheld, 3 control profiles keep LTR/Gypsy High). Unchanged between versions: M12927 (mdg4/gypsy) and X01472 (Drosophila gypsy) LTR/Gypsy High, AJ289709 (HERV-H/env62; env hit TLV_coat) LTR/Gypsy High, M11240, M18706 and X04456 LTR/Copia High.

**Residual limitations.** HERV-H (AJ289709, TLV_coat env) and any retrovirus/ERV whose gag/env is not hit by the five named profiles is still called Gypsy (Ty3) High; the list is five profile names, not a retroviral/non-retroviral discriminator, and the panel has no retroviral RT or integrase profile. The withheld label sits at confidence Moderate with outcome 'complete' (not Candidate/withheld), so consumers that treat outcome 'withheld' as the abstention signal will not see this as one. Bel-Pao and Retrovirus/ERV are still never named; the call stops at LTR order. An insect Ty3/gypsy element that happens to hit Gag_p10 or Gag_p24 would also be withheld (not measured).

**Benchmark probe.** Positive: AY037928 (HERV-K113) and, from the 4.0.0 corpus ERV rows, other HML-2 proviruses such as AY037929 and M14123: 4.0.0 expected LTR/Gypsy High, 4.1.0 expected LTR/unclassified Moderate (verified only for AY037928). Negative/control: M12927 and X01472 (Drosophila gypsy): LTR/Gypsy High in both versions; M11240 and X04456 (Copia) High in both. Boundary: AJ289709 (HERV-H, TLV_coat env) stays LTR/Gypsy High in 4.1.0 (a known residual, expected superfamily Retrovirus/ERV per Wicker). Synthetic: paired 300 bp direct repeats around back-translated RVT_1 consensus then rve consensus (Gypsy order) gives Gypsy High in both versions; adding a back-translated Gag_p24 consensus ORF 5' of the RT is the 4.1.0 trigger (hand-built-row equivalent verified; synthetic with the gag ORF not run).

**Pinned by.** `tests/test_41_classify.py::test_hervk_layout_withholds_the_superfamily_instead_of_naming_gypsy`, `tests/test_41_classify.py::test_each_retrovirus_specific_profile_withholds_gypsy`, `tests/test_41_classify.py::test_gypsy_type_profiles_keep_the_gypsy_call_and_errantivirus_wording`, `tests/test_41_classify.py::test_copia_order_is_unaffected_by_a_retroviral_profile`, `tests/test_41_classify.py::test_withheld_retroviral_label_with_tsd_does_not_print_an_ltr_target_site_claim`, `tests/test_41_classify.py::test_withheld_retroviral_call_keeps_the_architecture_and_the_erv_transcript_model`, `tests/test_41_classify.py::test_regression_anchors_keep_their_calls`.

Files: `app/backend/teagle_core/classify.py`, `tests/test_41_classify.py`.

#### CLS-003 — Pfam RVT_3 (PF13456) recoded from reverse transcriptase to an RNase H-like domain

**Before (4.0.0).** DOMAIN_INFO['RVT_3'] = ('RT', 'reverse transcriptase', 'retro', 'PF13456'). A hit from the RVT_3 model was reported as domain code RT with label 'reverse transcriptase' and counted as reverse-transcriptase evidence; a back-translated RNase_H or RVT_3 consensus returned an RT hit and a Class I call (retro/unassigned, Moderate).

**After (4.1.0).** DOMAIN_INFO['RVT_3'] = ('RNaseH', 'RNase H-like (Pfam RVT_3)', 'retro', 'PF13456'). The hit is reported as domain code RNaseH with that label (hmm name and pfam accession unchanged) and merges with an overlapping same-strand RNase_H hit under the existing same-code dedup (higher score kept). An RNase H-like region alone no longer licenses a Class I call. The bundled profile file, the 30-model panel size and the E-value are unchanged.

**Trigger.** Any ORF in which RVT_3 scores below the E-value cut-off (1e-3), notably RNase H regions and records whose only 'RT' evidence was RVT_3. Records with RVT_1 or RVT_2 hits are unaffected except where an RVT_3 hit was the highest-scoring 'RT' (RT representative used for strand/order tests).

**Affected outputs.** domains[].domain: RT -> RNaseH for hmm RVT_3; domains[].label: 'reverse transcriptase' -> 'RNase H-like (Pfam RVT_3)'; processing.domains.count (merged same-code hits); classification te_class/superfamily/confidence/tier/kind for RVT_3-only records; completeness.present/missing lists (RNaseH now includes RVT_3 hits); domain architecture strings (RT -> RNaseH); domain table exports and GFF3 domain features.

**Labels.** added: RNase H-like (Pfam RVT_3) (domain label); renamed: domain code RT -> RNaseH for hmm RVT_3, domain label 'reverse transcriptase' -> 'RNase H-like (Pfam RVT_3)'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| DOMAIN_INFO['RVT_3'] domain code | RT | RNaseH | no | InterPro/Pfam: PF13456 'Reverse transcriptase-like' sits in the RNase H clan (CL0219) and is integrated under IPR002156 (Ribonuclease H domain); an RNase H-like fold, not the polymerase core |
| bundled panel size (te_domains.hmm models) | 30 | 30 | yes | profile file is byte-identical to 4.0.0 (no diff) |
| domain E-value cut-off (scan_domains evalue) | 1e-3 | 1e-3 | yes | signature default unchanged |

**Measured.** Synthetic single-ORF back-translation (recipe R) via engine.analyze. RVT_3 consensus (375 bp): v4.0.0 domains [RNaseH(RNase_H, score 34.3), RT(RVT_3, 183.4)], retro/unassigned, Moderate, tier 'coding evidence present; order unassigned'; 4.1.0 domains [RNaseH(RVT_3, 183.4)], unclassified, Candidate, superfamily 'RNaseH domain(s) detected; no transposable-element assignment', tier 'coding fragment (RNaseH detected; no structural context)'. RNase_H consensus: same switch, domain RNase_H kept in 4.1.0. Hand-built control: RVT_1 + a separate RVT_3 hit: 4.0.0 withheld (lone RT, 'reverse transcriptase, class unassigned', Candidate); 4.1.0 retro/unassigned Moderate complete (see limitations). RVT_1 alone: withheld/Candidate in both.

**Residual limitations.** Side effect measured on hand-built rows: an RVT_1 RT plus a separate RVT_3 hit was withheld as a lone RT in 4.0.0 but is called 'retrotransposon order unassigned' Moderate in 4.1.0, because RNaseH is in the corroboration set (_RETRO_CORROBORATION); the non-specific RVT_3 hit now corroborates a lone RT. Elements whose only RT-like hit was RVT_3 lose their Class I call. No real accession carrying an RVT_3 hit was found among the cached public records, so only synthetic and hand-built cases were run.

**Benchmark probe.** Positive: back-translated RVT_3 consensus (recipe R): 4.0.0 retro/unassigned Moderate, 4.1.0 unclassified Candidate. Control: back-translated RVT_1 consensus: unclassified/withheld/Candidate in both. Boundary: an ORF with an RVT_1 region plus an RVT_3 region on the same strand (4.0.0 withheld, 4.1.0 retro/unassigned Moderate).

**Pinned by.** `tests/test_41_classify.py::test_rvt3_is_coded_as_an_rnase_h_like_domain`, `tests/test_41_classify.py::test_an_rnase_h_fold_alone_is_never_reported_as_reverse_transcriptase_or_class_i`.

Files: `app/backend/teagle_core/domains.py`, `app/backend/teagle_core/classify.py`, `tests/test_41_classify.py`.

#### CLS-004 — Overlapping integrase and transposase hits on one strand keep only the stronger

**Before (4.0.0).** _dedup_domains kept the best hit per (domain code, strand, overlapping nt region) only. The rve integrase model and DDE transposase models (DDE_3, etc.) cross-score the same DDE fold, so an integrase region also returned a TPase hit; a retroviral integrase consensus gave INT + TPase(DDE_3), a Class II call (DNA/Tc1-Mariner, Moderate), and an LTR element with a cross-scoring TPase hit was flagged composite, capping High at Moderate.

**After (4.1.0).** Where an INT hit and a TPase hit lie on the same strand and overlap by at least half of the shorter span, only the higher-scoring one is kept (symmetric: a stronger TPase also removes a weaker INT). Other code pairs (for example RNaseH with TPase, including the RNase-H-fold hAT models) and opposite strands are untouched.

**Trigger.** Records with an rve (INT) region that also scores a TPase profile above the E-value cut-off, or the converse, on the same strand with >=50% overlap of the shorter hit.

**Affected outputs.** domains[] (a hit disappears) and processing.domains.count; classification class (Class I vs II), te_class, superfamily, confidence, tier, kind for such records; composite flag effect on confidence (High no longer capped when the cross-hit was the only TPase).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| INT/TPase shadowing overlap rule | none | same strand and overlap*2 >= shorter span, INT vs TPase only | no | INT and TPase are both DDE folds whose models cross-score; the stronger model identifies the region |

**Measured.** Synthetic (recipe R) rve consensus alone (312 bp): v4.0.0 domains [INT, TPase DDE_3 score 17.7], DNA/Tc1-Mariner Moderate, tier 'partial (transposase present)'; 4.1.0 [INT], retro/partial Candidate, 'coding fragment (INT; RT/transposase not detected)'. Synthetic LTR element (300 bp direct repeats, spacer, rve ORF then RVT_1 ORF, or the reverse): v4.0.0 carried a cross-hit TPase:DDE_3 and returned Gypsy/Copia Moderate; 4.1.0 returns Gypsy/Copia High with no TPase hit. No cached public real element showed an INT/TPase overlap (13 accessions checked), so no real-sequence change was measured.

**Residual limitations.** A real DDE transposase region where rve outscores the transposase model would now lose its TPase hit and become an INT-only coding fragment (not observed in the cached real elements; untested on a corpus). The 2x overlap threshold is a heuristic, not a literature value.

**Benchmark probe.** Positive: back-translated rve consensus (recipe R): 4.0.0 DNA/Tc1-Mariner Moderate, 4.1.0 retro/partial Candidate. Negative: X01005 (Tc1) and X05424 (Ac): TPase hits intact, calls unchanged. Boundary: a hand-built INT and TPase pair with overlap exactly half the shorter span (kept = stronger one) and a pair on opposite strands (both kept).

**Pinned by.** `tests/test_41_classify.py::test_integrase_consensus_does_not_also_name_a_transposase`, `tests/test_41_classify.py::test_int_and_tpase_overlap_keeps_the_stronger_only`, `tests/test_41_classify.py::test_int_tpase_shadowing_is_narrow`.

Files: `app/backend/teagle_core/domains.py`, `tests/test_41_classify.py`.

#### CLS-005 — Withheld or conflicting classifications always report Candidate confidence

**Before (4.0.0).** Confidence was forced to Candidate only for conflicting transposase classes. A lone RT beside a TIR, an unattached EN hit or a Helitron helicase (outcome 'withheld') could still show Moderate through the 'rt and (has_ltr or has_tir)' and 'ndom >= 2' rules.

**After (4.1.0).** After the rule ladder, confidence is set to Candidate whenever tpase_conflict, order_conflict, rt_class_withheld or classification_withheld applies. Outcome vocabulary (complete/withheld) and the three confidence values are unchanged. The retroviral-undetermined and nested-integrase conflict labels are NOT in this set: they stay Moderate with outcome complete.

**Trigger.** A lone RT with one of: a TIR of any length (including >=20 bp), an EN hit outside an EN-RT arrangement, a Helitron helicase hit; a withheld terminal direct-repeat candidate; the new LINE-vs-LTR order conflict.

**Affected outputs.** classification.confidence: Moderate -> Candidate; classification explanation prefix '(candidate confidence)'; confidence column in tables and exports.

**Measured.** Hand-built rows: lone RT + 12 bp TIR, + 25 bp TIR, + EN in another ORF, + HEL: v4.0.0 Moderate/withheld, 4.1.0 Candidate/withheld (lone RT alone: Candidate in both). End-to-end synthetic (recipe R): RVT_1 consensus ORF between 25 bp inverted arms: v4.0.0 Moderate, 4.1.0 Candidate, both outcome withheld; same ORF without arms Candidate in both.

**Residual limitations.** The CHANGELOG states 'withheld outright always reports Candidate', which is true; but the RETROVIRAL_UNDETERMINED and nested-integrase labels are not withheld and report Moderate.

**Benchmark probe.** Positive: synthetic RVT_1 consensus ORF with 25 bp terminal inverted repeats (4.0.0 Moderate, 4.1.0 Candidate). Negative: same ORF without arms (Candidate in both). Boundary: arms of 12 bp versus 25 bp (both flip; the TIR length does not matter).

**Pinned by.** `tests/test_41_classify.py::test_withheld_lone_rt_is_never_moderate`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-006 — Terminal direct repeats no longer raise a transposase-only call to Moderate

**Before (4.0.0).** The confidence rule '(tpase and (has_ltr or tir_ok))' gave Moderate to a transposase call that had paired terminal DIRECT repeats (LTR row) but no enclosing TIR pair.

**After (4.1.0).** The rule is '(tpase and tir_ok)': only a TIR pair that encloses (or is not contradicted by) the transposase lifts the call from Candidate (when only one domain code) to Moderate. A transposase call with LTR-type direct repeats stays Candidate.

**Trigger.** A transposase (TPase) call with an accepted LTR row and no usable TIR; mainly DNA-transposase ORFs inside a segmental/tandem duplication or a mis-scanned flanked window.

**Affected outputs.** classification.confidence: Moderate -> Candidate for the trigger; wicker_order TIR remains (new field).

**Measured.** Hand-built LTR row + one hAT TPase hit: v4.0.0 DNA/hAT Moderate, 4.1.0 Candidate. With a 30 bp TIR instead: Moderate in both.

**Residual limitations.** Rare input (a transposase together with LTR-type direct repeats and no TIR); no real-sequence case was found in the cached records.

**Benchmark probe.** Positive: a hand-built or synthetic hAT transposase ORF flanked by a 300 bp direct-repeat pair (no inverted repeat): 4.0.0 Moderate, 4.1.0 Candidate. Negative: X05424 (Ac, TIR present): Moderate in both. Boundary: 30 bp TIR present together with direct repeats (Moderate in both).

**Pinned by.** `tests/test_41_classify.py::test_terminal_direct_repeats_do_not_lift_a_transposase_call`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-007 — Nested or adjacent integrases: every same-strand non-overlapping integrase is read against the RT

**Before (4.0.0).** Copia versus Gypsy used only the highest-scoring INT hit (_rep) and compared its translation position with the RT's. With several INT hits (nested or adjacent loci) the call followed the best bit score; an INT in another ORF could decide it.

**After (4.1.0).** All INT hits on the RT's strand that do not overlap it are collected. If any lie in the RT's own ORF only those are used. If INT hits remain on both sides of the RT, superfamily is 'LTR retrotransposon (superfamily undetermined)', te_class LTR/unclassified, confidence Moderate and the evidence says the integrase-vs-RT order is conflicting. With one side only, INT before RT gives Copia and INT after RT gives Gypsy as before. A single-INT locus behaves exactly as in 4.0.0. 'order' is None for the conflict case.

**Trigger.** Records with an accepted LTR pair, an RT and two or more INT hits on the RT's strand, at least one on each side of the RT, or an extra INT in the RT's own ORF with a stronger INT elsewhere.

**Affected outputs.** classification.te_class/superfamily (Copia or Gypsy -> LTR/unclassified, or Copia -> Gypsy); classification.confidence (High -> Moderate for the conflict); classification.order (INT–RT -> None for conflict); evidence sentence ('integrase hits lie on both sides of the RT ... conflicting').

**Measured.** Hand-built rows (RT 3000-3500 orf0 score 90): own-ORF INT 4700-5000 (23) plus upstream INT 1000-1300 orf2 (34): v4.0.0 LTR/Copia High (stronger upstream INT), 4.1.0 LTR/Gypsy. INTs in other ORFs on both sides: 4.0.0 LTR/Copia High, 4.1.0 LTR/unclassified Moderate. Single INT: High Gypsy in both. End-to-end synthetic (recipe R; 300 bp direct repeats around rve consensus, RVT_1 ORF, mutated rve): 4.0.0 Copia Moderate (an additional cross-hit TPase), 4.1.0 LTR/unclassified Moderate with the conflict sentence.

**Residual limitations.** The own-ORF preference treats an INT in the RT's ORF as decisive even when a second INT lies on the other side. Only INT hits on the RT's strand count; hits on the other strand are ignored here (a different rule, CLS-008 for YR).

**Benchmark probe.** Positive: synthetic LTR construct (recipe R) with INT-ORF, RT-ORF, INT-ORF in that order on one strand: 4.0.0 Copia/Gypsy by best score, 4.1.0 LTR/unclassified Moderate. Negative: single INT after RT (Gypsy High both) and single INT before RT (Copia High both). Boundary: two INTs on the same side of the RT (call unchanged).

**Pinned by.** `tests/test_41_classify.py::test_integrase_in_the_rts_own_orf_beats_a_stronger_integrase_elsewhere`, `tests/test_41_classify.py::test_bracketing_integrases_in_other_orfs_are_conflicting_not_a_confident_copia_or_gypsy`, `tests/test_41_classify.py::test_single_integrase_stays_high`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-008 — DIRS-group call requires the RT and tyrosine recombinase on the same strand

**Before (4.0.0).** DIRS-group = an RT hit and a YR hit and no INT, with no strand test. The same test gated the DIRS completeness ledger, the structure-free DIRS call and the ERV exclusion (is_erv false when 'yr and not intg').

**After (4.1.0).** yr now means an RT hit and a YR hit that share a strand (_same_strand). Opposite strands give no DIRS call: the lone-RT state applies ('reverse transcriptase, class unassigned', Candidate, outcome withheld) and, in the structure-free path, '{RT–YR} domain(s) detected; transposable-element class unassigned'. The same predicate now gates the DIRS ledger and the ERV exclusion, so an RT(+) with YR(-), an ENV and LTRs can be an ERV (LTR/unclassified) instead of DIRS.

**Trigger.** Records with both RT and YR hits on opposite strands (YR on the other strand).

**Affected outputs.** classification te_class/superfamily/confidence/outcome/tier/kind (DIRS -> unclassified/withheld); is_erv (False -> True where ENV and LTRs present); order (RT–YR -> None); retroviral stage (now runs for the ERV case).

**Measured.** Hand-built rows RT + YR(-): v4.0.0 DIRS Moderate complete, 4.1.0 unclassified Candidate withheld (superfamily 'reverse transcriptase, class unassigned'); same YR on the plus strand: DIRS Moderate in both. Structure-free path: DIRS-group (coding architecture; structure unassessed) -> 'RT–YR domain(s) detected; transposable-element class unassigned'. Case with LTR, RT(+), YR(-), ENV: v4.0.0 DIRS (is_erv False), 4.1.0 LTR/unclassified is_erv True, kind 'ERV (env-bearing LTR retroelement)'. End-to-end synthetic (recipe R; RVT_1 ORF and a reverse-complemented Phage_integrase consensus ORF 1.5 kb apart): 4.0.0 DIRS Moderate, 4.1.0 unclassified Candidate withheld.

**Residual limitations.** Strand is the only added test; no distance or ORF-proximity criterion, so a same-strand RT and YR kilobases apart is still called DIRS.

**Benchmark probe.** Positive: synthetic construct RT ORF (+) and YR ORF (-) (recipe R): 4.0.0 DIRS Moderate, 4.1.0 withheld Candidate. Negative: RT and YR ORFs on the same strand (DIRS Moderate in both). Boundary: LTR + RT(+) + YR(-) + ENV (4.0.0 DIRS, 4.1.0 ERV LTR/unclassified).

**Pinned by.** `tests/test_41_classify.py::test_dirs_needs_rt_and_recombinase_on_one_strand`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-009 — A LINE's tail is credited only on the strand it belongs to; the poly(A)-signal motif row is not a tail

**Before (4.0.0).** Any structural row whose type starts with 'poly' counted as a tail (has_polya), regardless of RT strand, and the advisory 'polyA-signal (motif, advisory)' row also started with 'poly'. A plus-strand LINE with a 5' poly-T (or minus-strand LINE with a 3' poly-A) reached 'intact / autonomous-consistent' and Moderate.

**After (4.1.0).** Tail types are matched exactly. A 3' poly-A tail is credited only when the RT hit is on the plus strand, a 5' poly-T only when it is on the minus strand (tail_ok). The opposite tract is reported in the evidence as present but not credited. The tail gates LINE confidence (Moderate if tail_ok else Candidate) and the LINE completeness tier ('intact / autonomous-consistent' needs the credited tail; otherwise 'near-complete (ORF1 + ORF2 endonuclease and RT recovered)'). The advisory polyA-signal motif is never a tail.

**Trigger.** A LINE call (te_class LINE) with a poly-A/poly-T row on the strand opposite to the RT, or a polyA-signal row alone.

**Affected outputs.** classification.confidence (Moderate -> Candidate); completeness.tier ('intact / autonomous-consistent' -> 'near-complete (...)'); completeness.missing/present tail entries; evidence sentence.

**Measured.** Hand-built (ORF1-EN-RT LINE, 4 strand/tail combinations): matching strand Moderate/'intact' in both; wrong strand v4.0.0 Moderate/'intact', 4.1.0 Candidate/'near-complete'. End-to-end synthetic: M80343 (L1.2, plus strand, 27 bp 3' poly-A) Moderate/'intact' in both; its reverse complement Moderate/'intact' in both; M80343 with the poly-A removed and 30 T prepended: v4.0.0 Moderate/'intact', 4.1.0 Candidate/'near-complete'. polyA-signal motif row alone with ORF1+RT: 4.0.0 Moderate, 4.1.0 Candidate.

**Residual limitations.** Applies to real inputs only when a tail of the wrong strand is detected, which is uncommon; the common LINE case (tail on the RT strand) is unchanged. Flanked LINE records change for a different reason (CLS-035).

**Benchmark probe.** Positive: synthetic M80343 with the 3' poly-A replaced by a 5' poly-T on the plus strand: 4.0.0 Moderate/intact, 4.1.0 Candidate/near-complete. Negative: M80343 and its reverse complement (unchanged Moderate/intact). Boundary: M80343 without any tail (Candidate in both).

**Pinned by.** `tests/test_41_classify.py::test_line_tail_is_credited_only_on_the_strand_it_belongs_to`, `tests/test_41_classify.py::test_polya_signal_motif_row_is_not_counted_as_a_tail`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-010 — LINE-specific coding evidence plus paired terminal direct repeats withholds the order (retro/unassigned)

**Before (4.0.0).** With an RT, accepted paired LTR-type direct repeats and no integrase, the call was 'LTR retrotransposon (superfamily undetermined)', LTR/unclassified, Moderate, tier 'fragment', kind 'LTR retroelement', even when ORF1p or an accepted EN-RT arrangement was present.

**After (4.1.0).** If ORF1p or an EN-RT arrangement is present with LTRs and no integrase the result is superfamily 'retrotransposon order unassigned', te_class retro/unassigned, Candidate (order_conflict), tier 'coding evidence present; order unassigned', kind 'retroelement (order unassigned)', and the ledger lists ORF1/EN. Evidence names both lines of evidence. wicker_order is 'unassigned'.

**Trigger.** Records with accepted LTR rows (including direct repeats from segmental or tandem duplication), an RT, no INT, and ORF1p or an EN-RT pair on one ORF/strand.

**Affected outputs.** classification te_class/superfamily/confidence/tier/kind; completeness.present (ORF1/EN); wicker_order (unassigned).

**Labels.** renamed: LTR/unclassified + 'LTR retrotransposon (superfamily undetermined)' -> retro/unassigned + 'retrotransposon order unassigned' (for the trigger).

**Measured.** Hand-built rows (ORF1+EN+RT, ORF1+RT, EN+RT, each with LTR 0-968/8504-9472): v4.0.0 LTR/unclassified Moderate tier 'fragment'; 4.1.0 retro/unassigned Candidate. End-to-end synthetic (300 bp direct repeats around M80343): v4.0.0 LTR/unclassified Moderate, 4.1.0 retro/unassigned Candidate, tier 'coding evidence present; order unassigned'.

**Residual limitations.** 'retrotransposon order unassigned' is not in the frozen 4.0.0 scorer vocabulary (it raises).

**Benchmark probe.** Positive: M80343 flanked on both sides by one 300 bp repeat copy (synthetic): 4.0.0 LTR/unclassified Moderate, 4.1.0 retro/unassigned Candidate. Negative: bare M80343 (LINE Moderate, unchanged). Boundary: ORF1-only plus direct repeats (also flips), RT-only plus direct repeats (stays LTR/unclassified, unchanged).

**Pinned by.** `tests/test_41_classify.py::test_line_coding_evidence_plus_direct_repeats_is_a_conflict_not_an_ltr_call`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-011 — RT + chromodomain calls use the LTR completeness ledger

**Before (4.0.0).** A call placed in the LTR order by RT + chromodomain (no LTR, no integrase; te_class LTR/unclassified) fell into the LINE ledger: kind 'non-LTR retroelement (LINE-like)', tier 'partial (RT only)', with ORF1, EN and the 3' tail listed as missing.

**After (4.1.0).** The LTR ledger gate also accepts te_class starting with LTR (and excludes retro/unassigned); CHR is added to expected/present; tier is 'partial (RT + chromodomain; integrase not confirmed)' and kind 'retroelement (LTR not confirmed)'. The call (LTR/unclassified, confidence) is unchanged.

**Trigger.** RT + chromodomain (CHR) with no accepted LTR and no integrase.

**Affected outputs.** completeness.tier; completeness.kind; completeness.present/missing/expected.

**Labels.** added: partial (RT + chromodomain; integrase not confirmed); renamed: tier 'partial (RT only)' -> 'partial (RT + chromodomain; integrase not confirmed)' (for the trigger), kind 'non-LTR retroelement (LINE-like)' -> 'retroelement (LTR not confirmed)' (for the trigger).

**Measured.** Hand-built RT 100-700 + Chromo 800-950: v4.0.0 LTR/unclassified, tier 'partial (RT only)', kind 'non-LTR retroelement (LINE-like)'; 4.1.0 same call, tier 'partial (RT + chromodomain; integrase not confirmed)', kind 'retroelement (LTR not confirmed)', CHR in present. te_class, confidence and outcome unchanged. (The 4.0.0 benchmark row AF050455 is the corpus case named in the wave record; not re-run.)

**Residual limitations.** The strand/ORF co-location of CHR with RT is still not required (documented 3.6.0 rule, left for an owner decision).

**Benchmark probe.** Positive: a synthetic ORF with RVT_1 and Chromo consensus segments (recipe R) and no LTR/INT: tier differs between versions. Negative: RT + ORF1 (LINE ledger unchanged). Boundary: RT + CHR with a detected INT (full LTR ledger in both).

**Pinned by.** `tests/test_41_classify.py::test_rt_plus_chromodomain_uses_the_ltr_ledger`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-012 — A long terminal inverted repeat beside only host-like domains is again a non-autonomous Class II candidate

**Before (4.0.0).** With a TIR of any length and only host-like coding domains (not RT/TPase and none of GAG/ENV/INT/ORF1), the call was 'unclassified', superfamily '{codes} domain(s) detected; no transposable-element assignment', and the TIR was not mentioned. Completeness kind read 'DNA transposon'.

**After (4.1.0).** If the TIR is >= MIN_TIR_STANDALONE (20 bp) the result is Class II, superfamily 'non-autonomous TIR element (superfamily undetermined)', te_class DNA/non-autonomous, Candidate, evidence naming the TIR length, the chance ceiling and the domain (possible captured gene fragment, Pack-MULE). A shorter TIR leaves the call 'unclassified' and adds a sentence that it is too short to attribute. Completeness kind follows the call (see CLS-013).

**Trigger.** TIR >= 20 bp and the detected domains are a subset of {CHR, PR, RNaseH, EN (not in EN-RT arrangement), YR}; no RT, TPase, GAG, ENV, INT, ORF1 or HEL. (RNaseH now includes RVT_3 hits, CLS-003.)

**Affected outputs.** classification class/te_class/superfamily (unclassified -> DNA/non-autonomous); wicker_order (TIR); completeness.kind; evidence.

**Labels.** renamed: '{codes} domain(s) detected; no transposable-element assignment'/unclassified -> 'non-autonomous TIR element (superfamily undetermined)'/DNA/non-autonomous (for the trigger).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| MIN_TIR_STANDALONE (bp, TIR length that licenses a Class II call without transposase) | 20 | 20 | yes | above the longest chance inverted repeat measured (17 bp over 22,500 random sequences, 35-65% GC) |

**Measured.** Hand-built (31 bp TIR + Chromo): v4.0.0 unclassified, 4.1.0 DNA/non-autonomous Candidate. 12 bp TIR + Chromo: unclassified in both, evidence adds 'too short to attribute'. End-to-end synthetic (recipe R): 31 bp inverted arms around a Chromo consensus ORF: v4.0.0 unclassified (kind 'DNA transposon'), 4.1.0 DNA/non-autonomous Candidate; 12 bp arms: unclassified in both, kind 'terminal inverted repeat (class unassigned)' in 4.1.0.

**Residual limitations.** A long inverted repeat flanking an unrelated host-like domain (for example a gene) is now a Class II candidate; the composition caution (CLS-021) is the only brake. The 4.0.0 label 'DNA/non-autonomous' with a plain long TIR and no domain at all is unchanged.

**Benchmark probe.** Positive: synthetic 31 bp inverted arms around a back-translated Chromo consensus (recipe R): 4.0.0 unclassified, 4.1.0 DNA/non-autonomous. Negative: 12 bp arms (unclassified both). Boundary: 19 bp versus 20 bp arms (the 20 bp floor is unchanged).

**Pinned by.** `tests/test_41_classify.py::test_long_tir_with_a_host_like_domain_keeps_the_class_ii_candidate_and_mentions_the_tir`, `tests/test_41_classify.py::test_working_domain_scan_still_calls_a_long_tir_non_autonomous`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-013 — Completeness kind follows the call (no 'DNA transposon' for a solo LTR or a short TIR)

**Before (4.0.0).** For records with an LTR or TIR row but no RT/transposase, completeness.kind was 'LTR retroelement' if an LTR row existed, else 'DNA transposon' regardless of the call. A solo copia LTR (276 bp, TIR scan trips), a short TIR beside GAG or EN, or a Helitron with a TIR all read 'DNA transposon'.

**After (4.1.0).** kind is 'LTR retroelement' if an LTR row exists; 'DNA transposon' if te_class starts with DNA/; 'retroelement fragment (order unassigned)' if te_class starts with retro; otherwise 'terminal inverted repeat (class unassigned)'. Tier strings are unchanged.

**Trigger.** Records with a TIR row (no LTR row) whose te_class is not DNA/*: solo LTRs, TIRs shorter than 20 bp, TIR + non-transposase domains, Helitron + TIR.

**Affected outputs.** completeness.kind.

**Labels.** added: retroelement fragment (order unassigned), terminal inverted repeat (class unassigned); renamed: kind 'DNA transposon' -> 'terminal inverted repeat (class unassigned)' or 'retroelement fragment (order unassigned)' (for the trigger).

**Measured.** M11240[:276] (solo copia 5' LTR) end-to-end: v4.0.0 te_class repeat/structural-only, kind 'DNA transposon'; 4.1.0 same class, kind 'terminal inverted repeat (class unassigned)'. Hand-built: 25 bp TIR + GAG: kind 'DNA transposon' -> 'retroelement fragment (order unassigned)'; 12 bp TIR + EN or + Chromo or YR: -> 'terminal inverted repeat (class unassigned)'; 12 bp TIR alone with no domain: -> 'terminal inverted repeat (class unassigned)'.

**Benchmark probe.** Positive: M11240[:276] (public GenBank copia 5' LTR): kind changes. Negative: X01005 (kind 'DNA transposon' in both). Boundary: TIR of 19 vs 20 bp with no domain (kind 'terminal inverted repeat (class unassigned)' vs 'DNA transposon' in 4.1.0).

**Pinned by.** `tests/test_41_classify.py::test_solo_copia_ltr_is_not_filed_under_dna_transposon_in_the_completeness_line`, `tests/test_41_classify.py::test_tir_plus_non_transposase_domain_kind_follows_the_call`, `tests/test_41_classify.py::test_helitron_evidence_does_not_claim_tir_absence_when_a_tir_was_detected`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-014 — Failed protein-domain scan withholds every structure-only call

**Before (4.0.0).** With domains_ok False and no RT/TPase evidence, a TIR record was called 'non-autonomous TIR element' / DNA/non-autonomous (outcome complete, 'no coding domain detected' in the ledger), a poly-A-only or empty record 'no clear TE signature' (te_class none), asserting an absence from a scan that never ran. A no-signal record with credited rows was 'complete'.

**After (4.1.0).** With domains_ok False: a TIR record gives superfamily 'terminal inverted repeat, class unassigned', te_class repeat/structural-only, outcome withheld, Candidate, completeness tier 'structural-only (protein-domain scan unavailable; coding content unassessed)' with 'protein domains' unassessed; any other structure-only or empty record gives 'classification withheld (protein-domain scan unavailable)', te_class unclassified, outcome withheld, completeness null. The explanation starts 'Classification withheld'. domains_ok True behaviour is unchanged.

**Trigger.** domains_ok False (the domain stage failed, e.g. pyhmmer missing or erroring) and no RT/INT/TPase evidence.

**Affected outputs.** classification superfamily/te_class/outcome/confidence/tier/kind; processing.classification.status (complete -> withheld); processing.retroviral (complete -> not_run); wicker_order (withheld), domain_architecture ('unavailable').

**Labels.** added: classification withheld (protein-domain scan unavailable); renamed: 'no clear TE signature' -> 'classification withheld (protein-domain scan unavailable)' (domains_ok False), 'non-autonomous TIR element (superfamily undetermined)' -> 'terminal inverted repeat, class unassigned' (domains_ok False).

**Measured.** Hand-built domains_ok False: 31 bp TIR: v4.0.0 DNA/non-autonomous complete, 4.1.0 repeat/structural-only withheld; poly-A only and empty rows: 'no clear TE signature' (te_class none) -> 'classification withheld (protein-domain scan unavailable)' (unclassified, withheld). End-to-end with pyhmmer forced unavailable (domains.pyhmmer=None): X01005 v4.0.0 DNA/non-autonomous complete, 4.1.0 terminal inverted repeat class unassigned, withheld, processing.classification 'withheld', retroviral 'not_run'; M11240 stays 'terminal direct-repeat candidate' withheld in both.

**Residual limitations.** Only the failed-stage case is handled; a domain stage that 'ran' on zero ORFs (status not_run, CLS-029) is still classified as if searched.

**Benchmark probe.** Positive: any record analysed with the domain scan disabled (fault injection, not an input property): X01005 gives DNA/non-autonomous (4.0.0) vs withheld TIR candidate (4.1.0). Negative: the same record with the scan working (Tc1/Mariner both). Boundary: LTR-only record (withheld candidate in both).

**Pinned by.** `tests/test_41_classify.py::test_failed_domain_scan_withholds_every_structure_only_call`, `tests/test_41_classify.py::test_working_domain_scan_still_calls_a_long_tir_non_autonomous`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-016 — Structure-free transposase call carries the advisory strongest-profile family

**Before (4.0.0).** When structural detection failed (structural_ok False) a transposase hit always gave superfamily 'DNA transposon (transposase detected; termini unassessed)' and te_class DNA/unclassified; the family of the strongest hit was ignored and conflicts were not noted.

**After (4.1.0).** The strongest-hit family mapping is a shared helper (_tpase_family). The structure-free call is '{family} (advisory: strongest transposase profile; termini unassessed)' with te_class DNA/{hAT|Tc1-Mariner|CACTA|MULE|IS4-like|DDE|unclassified}, Candidate, termini listed as unassessed, and a ledger sentence when several transposase classes are detected. The structure-aware path uses the same helper with unchanged output. wicker_order is TIR.

**Trigger.** structural_ok False (structural stage failed) with a TPase hit and no RT.

**Affected outputs.** classification superfamily/te_class; evidenceLedger (ambiguity sentence); wicker_order.

**Labels.** added: {family} (advisory: strongest transposase profile; termini unassessed) with family in: hAT, Tc1/Mariner, CACTA (En/Spm), MULE (Mutator), IS4-like DDE (piggyBac group), DDE transposon, DNA transposon; removed: DNA transposon (transposase detected; termini unassessed); renamed: 'DNA transposon (transposase detected; termini unassessed)' -> '{family} (advisory: strongest transposase profile; termini unassessed)', te_class DNA/unclassified -> DNA/{family token} on the structure-free path.

**Measured.** Hand-built structural_ok False with one hAT TPase: v4.0.0 DNA/unclassified 'DNA transposon (transposase detected; termini unassessed)'; 4.1.0 DNA/hAT 'hAT (advisory: ...)'. CACTA: DNA/CACTA 'CACTA (En/Spm) (advisory: ...)'. Two classes (hAT strongest + Tc1/Mariner): DNA/hAT with the ambiguity sentence.

**Residual limitations.** Fault-injection path only (structural stage failure); the frozen 4.0.0 scorer raises on all seven advisory labels.

**Benchmark probe.** Positive: any hAT record with the structural stage disabled. Negative: same record with structure working (DNA/hAT 'hAT', Moderate with TIR). Boundary: two transposase classes of different score.

**Pinned by.** `tests/test_41_classify.py::test_transposase_without_structure_keeps_the_advisory_strongest_profile_family`.

Files: `app/backend/teagle_core/classify.py`.

#### CLS-033 — Knock-on: flanked-record structural changes alter classification confidence, tier, TSD verdicts and evidence

**Before (4.0.0).** On records carrying genomic flanks, detect_all missed the TIR of flanked hAT/Ac elements (TSD 0/30), mislocated the TSD of LTR elements, and found no poly-A tail or TSD for a flanked LINE (confidence Candidate, tier 'near-complete (ORF1 + ORF2 endonuclease and RT recovered)').

**After (4.1.0).** Structural scope changes (flank-aware TIR search, boundary snapping with the TSD, poly-A tail with chance-controlled flank TSD) feed the unchanged classify rules: a flanked LINE with a credited 3′ tail reaches Moderate and 'intact / autonomous-consistent'; flanked hAT records get a TSD row and the 'captured' sentence; flanked LTR records get congruent 4-6 bp TSDs; a 49.7 kb human window changes from 'no clear TE signature' to a withheld 'terminal direct-repeat candidate (LTR origin unconfirmed)' because a sub-threshold pair became an accepted LTR row.

**Trigger.** Records with flanking genomic sequence (coordinate fetches, genome windows); bare reference elements are unchanged. Detailed detector changes belong to the structural-scope records.

**Affected outputs.** classification.confidence and completeness.tier for flanked LINEs; TSD congruence/evidence for flanked LTR and hAT; te_class/outcome for windows where an LTR row appears or disappears.

**Labels.** renamed: 'no clear TE signature' (none, complete) -> 'terminal direct-repeat candidate (LTR origin unconfirmed)' (repeat/structural-only, withheld) for one cached 49.7 kb human coordinate window.

**Measured.** engine.analyze, 30 seeded flanked records per cell (300 bp random flanks around bundled public records): M80343 with a 12 bp TSD: v4.0.0 Candidate in 30/30 (no tail/TSD), 4.1.0 Moderate in 30/30 (tail 27 bp and 12 bp TSD found, tier 'intact / autonomous-consistent'); M80343 with no TSD: Candidate in both (30/30). X05424 with 8 bp TSD: DNA/hAT Moderate in both, TSD row 0/30 -> 28/30. M11240/M12927 calls unchanged (LTR/Copia, LTR/Gypsy High, 30/30 each). Cached public coordinate windows: five of seven unchanged, one (human, 49,721 bp) changed class as above, one yeast 12,001 bp window kept LTR/Copia High with LTR 325 -> 334 bp and a 5 bp TSD found.

**Residual limitations.** The classification differences are consequences of detector changes in other scopes; bare reference elements (the five bundled examples and the cached element accessions) keep their te_class, confidence and tier.

**Benchmark probe.** Positive: flanked M80343 (300 bp random flanks, 12 bp TSD): Candidate (4.0.0) vs Moderate/intact (4.1.0). Negative: bare M80343 (Moderate/intact both) and flanked M80343 without TSD (Candidate both). Boundary: flanked X05424 with a 7-10 bp TSD.

**Pinned by.** `tests/test_41_structural.py::test_flanked_l1_tail_and_its_target_site_duplication_are_found`, `tests/test_41_classify.py::test_bare_specimens_keep_their_anchors_through_the_tc1_wiring`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

### Structural detection: LTR, TIR, TSD, poly-A, cis-elements

#### CLS-017 — TSD congruence uses superfamily length ranges and a required motif (Tc1/Mariner TA)

**Before (4.0.0).** TSD_EXPECT held a single exact expected length per superfamily (Copia 5, Gypsy 5, hAT 8, Tc1/Mariner 2); tsd_congruence(length, superfamily) judged congruent only when length == expected and took no motif. A 4 or 6 bp TSD on a Copia/Gypsy call was 'incongruent'; any 2 bp repeat on Tc1/Mariner was congruent; the LTR incongruence sentence said 'likely coincidental' and quoted 'duplicate 5 bp'.

**After (4.1.0).** TSD_EXPECT entries are (preferred length, accepted range, basis, required motif): Copia and Gypsy 5 (4-6), hAT 8 (8-8), Tc1/Mariner 2 (2-2) with motif TA. tsd_congruence(length, superfamily, motif=None) is congruent when length is in range and (if a motif is supplied and one is required) the motif matches; otherwise incongruent. classify passes the TSD row's motif; the sentences quote the accepted range and motif ('duplicate 4-6 bp', '2 bp (GC) where Tc1/Mariner elements duplicate 2 bp (TA)') and say the repeat may be coincidental or the boundary/superfamily call may be wrong. The preferred length is still what _refine_tsd asks find_tsd for. Superfamilies absent from the table (CACTA, MULE, IS4-like, DDE, unclassified) remain 'not diagnostic'.

**Trigger.** Any LTR (Copia/Gypsy) or DNA (hAT, Tc1/Mariner) call with a TSD row.

**Affected outputs.** TSD row tsd_congruence verdict (congruent/incongruent) for 4 and 6 bp LTR TSDs and for non-TA 2 bp Tc1/Mariner TSDs; evidence sentences; GFF3 tsd_congruence attribute (reads the verdict); sealed parameters (tsd_expectation_by_superfamily, CLS-022).

**Labels.** renamed: evidence 'likely coincidental' -> 'the repeat may be coincidental, or the boundary or superfamily call may be wrong'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| TSD_EXPECT['Copia'] | length 5, no range | length 5, range 4-6 | no | Wicker et al. 2007 Nat Rev Genet 8:973-982 (4-6 bp for LTR retrotransposons), checked against the HML-2 survey; 5 bp preferred per Ou et al. 2019 Genome Biol 20:275 |
| TSD_EXPECT['Gypsy'] | length 5, no range | length 5, range 4-6 | no | same as Copia |
| TSD_EXPECT['hAT'] | length 8 | length 8, range 8-8 | yes | Ou et al. 2019 |
| TSD_EXPECT['Tc1/Mariner'] | length 2, no motif | length 2, range 2-2, motif TA | no | Wicker et al. 2007 (TA target site defines the superfamily) |
| MIN_TSD / MAX_TSD | 2 / 12 | 2 / 12 | yes | search range unchanged |

**Measured.** Hand-built rows: Copia with a 4 bp or 6 bp TSD: v4.0.0 incongruent ('duplicate 5 bp'), 4.1.0 congruent; 2 bp TSD on Copia: incongruent in both, now quoting '4-6 bp'; Tc1/Mariner with 2 bp TA: congruent in both; with 2 bp GC: congruent in 4.0.0, incongruent in 4.1.0 (text names 'GC' versus 'TA'). End-to-end, 30 flanked records per cell (300 bp random flanks, X01005/M11240/M12927 with planted TSDs): Copia 5 bp TSD congruent 15/30 (4.0.0) vs 29/30 (4.1.0); Copia 6 bp 9/30 vs 29/30; Gypsy 5 bp 20/30 vs 30/30; Copia 4 bp 0/30 vs 29/30. These end-to-end shifts also reflect LTR boundary snapping (CLS-035/structural scope), not only the range.

**Residual limitations.** The 4-6 bp LTR range is citation-based, not fitted to a corpus. A 5 bp TSD on an LTR element is flagged chance_level True in 4.1.0 (search-adjusted probability), though the LTR congruent sentence does not use that flag.

**Benchmark probe.** Positive: flanked M11240 with a planted 4 bp or 6 bp TSD (synthetic flanks around the public GenBank copia record): 4.0.0 incongruent, 4.1.0 congruent. Negative: planted 5 bp (congruent both once the boundary is placed). Boundary: planted 3 bp (incongruent both), and a Tc1/Mariner call with a GC 2 bp TSD (4.0.0 congruent, 4.1.0 incongruent).

**Pinned by.** `tests/test_41_structural.py::test_ltr_tsd_congruence_accepts_the_4_to_6_bp_range`, `tests/test_41_structural.py::test_tc1_tsd_congruence_requires_the_ta_motif_when_a_motif_is_given`, `tests/test_41_classify.py::test_the_tsd_motif_reaches_the_congruence_test`, `tests/test_41_classify.py::test_tsd_expectation_is_quoted_as_the_range_the_congruence_test_uses`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### CLS-019 — Tc1/Mariner TIR trimmed to the true arms and TA target site recovered inside classify

**Before (4.0.0).** On a flanked Tc1/Mariner record the TIR arms run 2 bp through the palindromic TA target site (tir_len 56-58 instead of 54) and either no TSD row exists or a chance 3 bp repeat (for example GTA) is flagged incongruent. Nothing in classify trimmed the arms.

**After (4.1.0).** After an unambiguous Tc1/Mariner call (not when transposase classes conflict) and a credited TIR row with the record sequence, classify calls structural.palindromic_target_site with the TA motif and updates the TIR row in place (arms trimmed up to 6 bp until TA flanks both ends, boundary_adjustment recorded) and sets or inserts a single TSD row (length 2, motif TA, tsd_congruence congruent). The TSD row carries selected_for_motif True when the arms were trimmed to expose the motif, chance_level True (2 bp), and the evidence says the TA match is by construction and not independent evidence. has_tsd, tir_encloses_tpase and tir_ok are recomputed. Gated on the call: hAT or other superfamilies are never trimmed. Bare records (no flanks) return no TSD and are unchanged.

**Trigger.** Flanked records (genomic windows) that classify as Tc1/Mariner with a detected TIR and a TA at the trimmed element boundaries.

**Affected outputs.** structural TIR row: tir_len (56-58 -> 54), five_prime/three_prime/element_span, boundary_adjustment, identity; TSD row: added or replaced (length 2, motif TA, tsd_congruence congruent, selected_for_motif, chance_level True); evidence sentences (hedged); GFF3 TE feature span and TSD rows; structural list length (see CLS-032).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _TIR_TRIM_MAX / palindromic_target_site max_trim (bp) | absent (function not wired) | 6 | no | arms found by extension may run through the target site by at most a few bp |

**Measured.** engine.analyze on X01005 (Tc1, GenBank) with 300 bp random flanks and a planted TA TSD, 30 seeds: v4.0.0 TIR 54 bp in 0/30 (TIR 56-57, TA TSD row 0/30, one 3 bp GTA row), 4.1.0 TIR exactly 54 bp 30/30 and a 2 bp TA TSD row 30/30 (all congruent, chance_level True, selected_for_motif True, evidence hedged). Controls (30 each): no TSD, random flanks: no TA TSD fabricated in 30/30 (TIR 54 in 24/30 in both versions); planted GC TSD: no trim in either version (TIR 54 0/30 both); planted random 3 bp TSD: TIR 54 in 25/30 both, incongruent 25/30 both. Bare X01005: TIR 54, no TSD, unchanged.

**Residual limitations.** A 2 bp TA is always chance level (p about 0.06-0.29), so the recovered TSD is congruent by construction and is not credited as insertion-site evidence; completeness tier is unchanged. The trim is skipped when transposase classes conflict, and applies only to a Tc1/Mariner call; flanked mariner-like elements called another superfamily keep arms 2 bp too long. The recovered TA TSD is reported as congruent although the search required TA; selected_for_motif is the only flag. Inserting the TSD row after classification desynchronises processing.structural.count (CLS-032).

**Benchmark probe.** Positive: public GenBank X01005 with 300 bp random flanks and TA TSD (synthetic flanks): 4.0.0 TIR 56-57 no TA TSD, 4.1.0 TIR 54 and 2 bp TA TSD. Negative: the same element with no TSD between random flanks (no TSD fabricated; 24/30 TIR 54 in both), and flanked X05424 (hAT; never trimmed). Boundary: GC 2 bp TSD (unchanged, no trim) and a 3 bp random TSD (incongruent both).

**Pinned by.** `tests/test_41_classify.py::test_flanked_tc1_reports_the_trimmed_tir_and_its_ta_target_site`, `tests/test_41_classify.py::test_target_site_trim_is_gated_on_the_tc1_call`, `tests/test_41_classify.py::test_bare_specimens_keep_their_anchors_through_the_tc1_wiring`.

Files: `app/backend/teagle_core/classify.py`, `app/backend/teagle_core/structural.py`.

#### CLS-020 — Refined TSD keeps the boundary-placement count of the span it flanks

**Before (4.0.0).** _refine_tsd re-ran find_tsd(seq, up[1], dn[0], expect=exp) and tsd.update(refined) in place; TSD rows had no chance_probability, n_tests or chance_level at all.

**After (4.1.0).** find_tsd takes placements and the refined row inherits the LTR or TIR row's boundary_tests (36 for a snapped LTR), so the refined TSD's chance_probability, n_tests (lengths x placements) and chance_level reflect the multiple testing of the boundary. The re-detection preference for the superfamily's target-site length (shorten a coincidental longer TSD, never lengthen) is unchanged.

**Trigger.** Flanked LTR or TIR records with a TSD row and a Copia, Gypsy, hAT or Tc1/Mariner call.

**Affected outputs.** TSD row chance_probability, n_tests, chance_level (new fields; value after refinement).

**Measured.** End-to-end flanked records (30 per cell): Copia with a 5 bp planted TSD carries chance_level True in 29/30 (4.0.0: key absent); Gypsy 5 bp True 30/30; hAT 8 bp False 28/30; LINE 12 bp False 30/30. Copia 6 bp: False 17/30, True 12/30.

**Residual limitations.** A correctly detected 5 bp LTR TSD is flagged chance-level (search-adjusted), which a UI badge could read as weak evidence; the LTR evidence text does not use the flag.

**Benchmark probe.** Positive: flanked M11240 with a 5 bp TSD: TSD row has chance_level True in 4.1.0 and no chance fields in 4.0.0. Negative: flanked M80343 with a 12 bp TSD (chance_level False). Boundary: 6 bp versus 5 bp TSD (a 6 bp TSD is often below the chance-level threshold).

**Pinned by.** `tests/test_41_structural.py::test_tsd_rows_carry_composition_adjusted_chance`.

Files: `app/backend/teagle_core/classify.py`, `app/backend/teagle_core/structural.py`.

#### CLS-030 — LTR cis-element orientation taken from the RT/INT/TPase strand when the structure cannot orient the element

**Before (4.0.0).** structural.detect_all ran once without an orientation; an LTR row with cis_orientation 'undetermined' kept no PBS/PPT/poly(A)-signal rows on the right strand.

**After (4.1.0).** After the domain scan, if every RT, INT and TPase hit shares one strand and the LTR row's cis_orientation is 'undetermined', detect_all is re-run with orientation=that strand; the new LTR row carries cis_orientation_source 'domain-strand'; a cis-element-determined orientation is never overridden; a failure keeps the first result and adds a note. Shared with the structural/cis-element scope.

**Trigger.** Records with an accepted LTR whose cis-elements do not fix the orientation and whose RT/INT/TPase hits are unanimous in strand (typically a reverse-complemented element).

**Affected outputs.** structural rows PBS, PPT, polyA-signal and their strand; LTR row cis_orientation_source; processing.structural.count; classification structural summary text.

**Labels.** added: cis_orientation_source: domain-strand.

**Measured.** Reverse complements of cached public records, engine.analyze: AY037928_rc gains PBS and PPT rows (LTR 968 bp unchanged); AJ289709_rc gains PBS and PPT; X04456_rc (copia) gains PPT; X01472_rc: polyA-signal row replaced by PPT. Forward-strand records: structural rows unchanged for these accessions. Classification calls unchanged in all of them.

**Residual limitations.** Mixed-strand domain hits give no vote; a nested element on the other strand blocks the re-orientation.

**Benchmark probe.** Positive: reverse complement of X04456 (PPT appears in 4.1.0 only). Negative: forward X04456 (unchanged) and a record whose PBS fixes the orientation. Boundary: RT and INT hits on opposite strands (no vote).

**Pinned by.** `tests/test_41_engine.py::test_undetermined_ltr_orientation_is_taken_from_the_domain_strand`, `tests/test_41_engine.py::test_orientation_is_not_taken_from_a_split_or_irrelevant_domain_vote`, `tests/test_41_engine.py::test_a_cis_element_orientation_is_not_overridden_by_the_domain_strand`.

Files: `app/backend/engine.py`.

#### STR-001 — LTR pair measured with a gapped alignment; identity redefined; lower-bound cause named

**Before (4.0.0).** v4.0.0: find_ltr seeded 13-mers on the best diagonal (+/-3), then extended that single diagonal with an ungapped X-drop (+1 match, -2 mismatch, x_drop 30). Identity was matches / min(len) of the two ungapped copies (N==N never a match, but N counted in the denominator). One indel put everything downstream out of phase: a copia LTR pair with a 1-bp deletion read 75.0% (advisory 'Sub-threshold terminal direct repeat'), a 10-bp deletion truncated the LTR to the in-phase part (166 bp, 100%, span shifted by 110 bp). length_is_lower_bound was set when the repeat touched a record end or an extension stopped at its limit; bound_reason was always the single 'still matched where the record ends' sentence. Method string was always 'k-mer seed (k=13, W bp window) + X-drop diagonal extension'.

**After (4.1.0).** 4.1.0: the ungapped measurement is kept as a candidate. A gapped candidate is also built, and attempted only when there is no ungapped candidate, the ungapped identity is below 97.0, or seed anchors lie 4-32 diagonals from the best diagonal. The gapped candidate chains collinear anchors across diagonals up to 32 away (_chain), extends both ends with a banded affine-gap X-drop (gap open 4, extension 1, band = indel margin 8 + 4), and measures identity from one banded global alignment of the final extents (identity = matches / (matches + mismatches + gap columns); ambiguous columns excluded; skipped, leaving the ungapped result, when the band exceeds 250,000 cells). The candidate with the higher score (matches - 2 mismatches - 4 per gap event - 1 per gap column) is reported. Ungapped identity is now matches / columns where both bases are A/C/G/T. Rows from the gapped path carry indels (gap events) and a method string 'banded gapped alignment'. The lower-bound flag still fires on a record edge, and bound_reason now names the cause ('the repeat still matched where the record ends', 'the two copies met each other', or both joined by ' and'). Accepted rows then pass through _snap_ltr (STR-002).

**Trigger.** Any record in which the two LTR copies differ by an indel of up to 32 bp, or whose ungapped identity is below 97%. Bare specimens whose copies are identical or substitution-only stay on the ungapped path and give the same length and identity. Records with anchors only on the best diagonal and ungapped identity at or above 97% are unaffected.

**Affected outputs.** LTR row ltr_len and identity: higher and full length where an indel previously truncated or deflated the pair; LTR row new key indels (gap events), present only when non-zero; LTR row method: 'banded gapped alignment' on the gapped path; LTR row bound_reason: text depends on the cause; Near-miss row 'Sub-threshold terminal direct repeat (advisory)' becomes 'LTR (terminal direct repeat)' where the gapped identity clears 80.0 (and classification follows); element_span for pairs previously truncated (10-bp deletion case: [110,5136] -> [0,5136]).

**Labels.** added: method value 'k-mer seed (k=13, W bp window) + banded gapped alignment', bound_reason 'the two copies met each other, so the copies may extend beyond it ...', bound_reason 'the repeat still matched where the record ends and the two copies met each other, ...'; renamed: 'Sub-threshold terminal direct repeat (advisory)' -> 'LTR (terminal direct repeat)' for gapped-rescued pairs (row type, not a renamed string).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _GAP_OPEN | absent | 4 | no | new; affine gap open cost on the +1/-2 scale (matches the extension scoring) |
| _GAP_EXT | absent | 1 | no | new; per gap column |
| _LTR_MAX_INDEL | absent | 32 | no | new; anchors up to this many diagonals from the best are chained across an indel; measured limit: 32 bp recovered, 33 bp not |
| _GAP_BAND_MARGIN | absent | 8 | no | new; band slack around the diagonals visited |
| _GAPPED_BELOW_IDENTITY | absent | 97.0 | no | new; gapped attempt only below this ungapped identity (keeps bare specimens on the 4.0.0 path) |
| _MAX_BAND_CELLS | absent | 250000 | no | new; larger banded alignments are not attempted and the ungapped measurement stands |
| MIN_LTR_IDENTITY | 80.0 | 80.0 | yes | frozen identity floor; what changed is the identity measure it is applied to |
| find_ltr k / min_ltr / min_anchors | 13 / 80 / 4 | 13 / 80 / 4 | yes | signature defaults, unchanged |
| _extend x_drop and scoring | 30; +1/-2 | 30; +1/-2 | yes | unchanged; now sealed as ltr_extend_x_drop |

**Measured.** Own run, both engines in separate processes, same seeds. M11240.1 (copia) bare: 1-bp deletion 100 bp into the 3' LTR: v4.0.0 advisory row, 276 bp, 75.0%; 4.1.0 LTR row, 276 bp, 99.6%, indels 1. 2-bp deletion: 75.7% -> 99.3%. 10-bp deletion: v4.0.0 LTR 166 bp at 100% span [110,5136]; 4.1.0 LTR 276 bp at 96.4% span [0,5136]. Synthetic pairs (300-bp random LTR, 2 kb random interior, 3' copy mutated, 20 seeds 0-19): 5% substitution + 1% indel (indels up to 3 bp): accepted 4/20 -> 20/20, length within 5% of truth 2/20 -> 20/20; 10% sub + 2% indel: 2/20 -> 20/20; 15% sub + 2% indel: 1/20 -> 19/20 (18 full length). Single deletion in a 400-bp LTR (20 seeds): full 400 bp recovered at deletion 6, 20, 31, 32 bp in 20/20 (4.1.0) versus 0/20 (v4.0.0, reports ~200 bp at 100%); at 33 and 36 bp both versions report ~200 bp. Insertion in the 3' copy: 20 and 32 bp recovered 20/20, 40 bp not. Substitution-only controls are unchanged: 300-bp LTR at 5% sub 20/20 both; 20% sub 19/20 both. Realised substitution sweep (60 seeds each, 300-bp LTR): accepted at 20% divergence 38 -> 39, at 22% 26 -> 28, 25% 8 -> 8, 28% 2 -> 2, 30% 0 -> 0 (the gapped measure moves a few pairs across the 80% floor). Bare copia M11240.1 and gypsy M12927.1: LTR 276 and 482 bp at 100.0%, unchanged. Wave record: 40 diverged pairs (4% sub + 3 indels, 500 bp) full length 4/100 -> 40/40.

**Residual limitations.** Indels longer than 32 bp (deletion or insertion) still halve the LTR to the in-phase part, as in v4.0.0 (measured at 33, 36, 40, 60 bp). The near-miss advisory row built from a gapped candidate still carries the method text 'X-drop diagonal extension' because that string is hard-coded in the rejection branch. Gapped identity is a score-optimal alignment and is slightly optimistic for heavily diverged pairs: a copia pair with 120 bp of 3'-LTR randomised read 79.7% ungapped (v4.0.0 advisory) and 82.5% gapped (accepted, classification moves from superfamily-undetermined to Copia). LTRs of 6 kb or more cannot be seeded (seed window 6000 bp); gapped alignment is skipped above 250,000 band cells so very long diverged pairs fall back to the ungapped measurement.

**Benchmark probe.** Positive: take M11240.1, delete 1 bp at offset +100 inside the 3' LTR (3' LTR starts at 4870); run the analysis. v4.0.0: advisory 'Sub-threshold terminal direct repeat', 75.0%, structural evidence 'LTR' absent; 4.1.0: 'LTR (terminal direct repeat)', ltr_len 276, identity 99.6, indels 1. Also a 400-bp random LTR with a 20-bp deletion at offset 200 in the 3' copy flanked by 2 kb random sequence: v4.0.0 ~200 bp at 100%; 4.1.0 400 bp at ~95%, indels 1. Negative/control: the same construct with substitutions only (5%), or the unmodified bare specimens, must be identical in both versions (copia 276 @100, gypsy 482 @100). Boundary: sweep the deletion length 31, 32, 33 bp in the 400-bp LTR: 4.1.0 returns 400 bp at 31 and 32 and ~200 bp at 33; v4.0.0 returns ~200 bp throughout.

**Pinned by.** `tests/test_41_structural.py::test_single_indel_in_a_ltr_copy_reads_as_one_indel_not_divergence`, `tests/test_41_structural.py::test_ltr_with_an_indel_is_classified_as_ltr_not_advisory`, `tests/test_41_structural.py::test_diverged_ltr_pairs_with_several_indels_keep_their_full_length`, `tests/test_41_structural.py::test_lower_bound_names_the_actual_cause`, `tests/test_structural_bounds.py::test_diverged_ltr_pair_is_detected`, `tests/test_structural_bounds.py::test_repeat_running_off_the_record_is_flagged_as_a_lower_bound`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-002 — LTR outer boundaries re-chosen jointly with the flanking TSD (boundary snapping)

**Before (4.0.0).** v4.0.0: the LTR outer ends were wherever the X-drop extension stopped. An extension past the true end keeps a chance-matching base about a quarter of the time, so a flanked element's boundary landed 1-4 bp too wide and its real TSD was searched at the wrong place: a 5-bp TSD was typically reported as 4, 3 or 2 bp (flagged incongruent downstream) and ltr_len came out 277 or 483 instead of 276 or 482.

**After (4.1.0).** 4.1.0: every accepted LTR pair is passed through _snap_ltr. Over 36 shifts (5' end -1..+4, 3' end -4..+1) a shifted boundary is taken only if it gives a TSD of at least 4 bp, total evidence at least 5, and a gain of at least 1.5 over the unshifted boundary after a cost of 0.5 per base moved. Evidence is the TSD length capped at 6, plus 1 when it equals the Copia modal length (5), plus 3 when both LTR ends read canonical TG..CA. On a move the row gets new ltr_len, five_prime, three_prime, element_span, termini, boundary_adjustment [a, b], boundary_tests 36; termini.selected_for_motif is true when the canonical-motif bonus decided the move. Ungapped rows get identity recomputed; gapped rows keep the alignment identity. detect_all then searches the TSD at the snapped span with placements = boundary_tests, so the TSD chance is adjusted (STR-011).

**Trigger.** Flanked records (at least 4 bp of flank on both sides) holding an accepted LTR pair whose unshifted boundary has no TSD or a shorter one than a shifted boundary offers. Bare records and records with the TSD already exact are unchanged.

**Affected outputs.** LTR row five_prime/three_prime/element_span/ltr_len (shift of up to 4 bp); LTR row termini (canonical flag may change); TSD row length, motif, upstream/downstream coordinates; LTR row new keys boundary_adjustment, boundary_tests, termini.selected_for_motif; PPT, PBS and polyA-signal coordinates follow the snapped boundary; tsd_state and downstream congruence verdicts.

**Labels.** added: boundary_adjustment (LTR row), boundary_tests (LTR row, value 36), termini.selected_for_motif.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _SNAP_START | absent | range(-1,5) | no | new; 5' end shift window, positive moves inward |
| _SNAP_END | absent | range(-4,2) | no | new; 3' end shift window |
| _SNAP_MIN_BP | absent | 5 | no | new; minimum evidence to move a boundary |
| _SNAP_MIN_TSD | absent | 4 | no | new; bottom of the 4-6 bp LTR TSD range |
| _SNAP_CAP | absent | 6 | no | new; TSD length counted up to the LTR range |
| _CANONICAL_BP | absent | 3 | no | new; bonus for TG..CA at both ends (1 in 256 event) |
| _SHIFT_COST | absent | 0.5 | no | new; per base moved |
| _SNAP_MARGIN | absent | 1.5 | no | new; gain required beyond the cost |

**Measured.** Own run: M11240.1 copia and M12927.1 gypsy embedded between 300-bp random flanks with a random 5-bp TSD, 30 seeds (0-29): exact element span plus TSD at the planted position copia 11/30 -> 29/30, gypsy 20/30 -> 29/30. Flank sweep, 10 seeds per cell (flanks 50, 300, 800, 1500, 3000 bp): copia exact 3/10 -> 9/10 in every cell; gypsy 7/10 -> 9-10/10. Engine, copia seed 4: v4.0.0 LTR 277 bp, TSD 4 bp; 4.1.0 LTR 276 bp, TSD 5 bp; gypsy seed 3: LTR 483 + TSD 4 bp -> 482 + 5 bp. Wave records (60 seeds, 300-bp flanks): copia exact 34/60 -> 58/60 (0.97), gypsy 25/60 -> 0.82; the gypsy ceiling is a genuine tie between the two ends (one overshooting end leaves a matching TSD prefix at the other).

**Residual limitations.** The snap prefers the LTR modal 5 bp; a real 4- or 6-bp TSD is retained only when it is the longest TSD on offer. Gypsy-type records tie between ends in about one in six flanked cases and are left at an unsnapped boundary. A TSD shorter than 4 bp never moves a boundary. The snap applies to every superfamily's LTR pair, including ones with a different true TSD length (HERV 4-6 bp is the documented range; others are untested).

**Benchmark probe.** Positive: embed M11240.1 between two random 300-bp flanks with a random 5-bp TSD (copy the 5 bp to both sides), 30 seeds; score element_span == [305, 5451]-type exact span and TSD 5 bp. v4.0.0 ~11/30, 4.1.0 ~29/30; v4.0.0 often gives ltr_len 277 and TSD 4 bp. Negative: the same without a planted TSD must give no 12-bp TSD and no boundary_adjustment on most seeds (test_no_tsd_control). Boundary: planted TSD of 3, 4, 6 bp: 3 bp never moves a boundary; 4 and 6 bp move it only if the shifted boundary yields the longer match.

**Pinned by.** `tests/test_41_structural.py::test_flanked_element_boundaries_and_tsd[M11240-LTR-0.9]`, `tests/test_41_structural.py::test_ltr_boundary_snap_records_the_shifts_searched_and_the_motif_it_used`, `tests/test_41_structural.py::test_no_tsd_control_does_not_fabricate_long_target_sites`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-003 — Terminal repeats longer than about 5.1 kb no longer rejected as simple repeats

**Before (4.0.0).** v4.0.0: _is_simple rejected a sequence when its number of distinct 4-mers was below max(4, 0.05 * (len - 3)). The proportional term exceeds the 256 possible 4-mers at about 5,123 bp, so any LTR candidate longer than that was rejected as a simple repeat although the seed window is 6000 bp.

**After (4.1.0).** 4.1.0: the floor is max(4, min(192, 0.05 * (len - 3))). Behaviour is identical below about 3,840 bp (where 0.05*(len-3) < 192) and the change matters above about 5.1 kb. Exact microsatellites remain simple (e.g. (ATTA)n and (AT)n at 6 kb).

**Trigger.** LTR candidates (and TIR arms checked with _is_simple) longer than about 3.84 kb; decisive above 5.12 kb.

**Affected outputs.** find_ltr returns a row (instead of None) for random-sequence LTR pairs of 5.2-5.9 kb; downstream LTR classification and element_span for such elements.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _is_simple distinct-4-mer cap | none (floor 0.05*(len-3), reaching 256 at 5,123 bp) | 192 (0.75 x 256) | no | new cap; chosen so a long random repeat stays judgeable; the literal 192 is not sealed |
| _is_simple proportional term and absolute floor | 0.05 and 4 | 0.05 and 4 | yes | unchanged |
| _LTR_SEED_WINDOW | 6000 | 6000 | yes | unchanged; an LTR of 6 kb or more still cannot be seeded |

**Measured.** Own run: random 5,500-bp LTR pair, 200-bp flanks, 1.5 kb interior (seed 11): v4.0.0 no LTR row; 4.1.0 LTR 5,500 bp, 100.0%. A 4,500-bp pair (seed 12) is identical in both (4,500 bp, 100.0%). Wave records: pairs of 5.2, 5.5 and 5.9 kb at 0% and 2% divergence all found within 3 bp of truth in 4.1.0, all six rejected in v4.0.0; 3.0 and 4.5 kb unchanged.

**Residual limitations.** Candidates between roughly 3.84 and 5.1 kb get a slightly lower floor than v4.0.0 (still met by random sequence, so practically unchanged). LTRs of 6 kb and longer remain undetectable by the seed window.

**Benchmark probe.** Positive: a random 5.5-kb sequence used as both LTRs, 1.5 kb random interior, 200-bp random flanks (one construct per seed): v4.0.0 no LTR row; 4.1.0 LTR 5,500 +/- 3 bp at 100%. Negative: (AT)n or (ATTA)n at 6 kb used as both LTRs must still be rejected in both versions. Boundary: 4,500-bp pair identical in both; 6,100-bp pair not detected in either.

**Pinned by.** `tests/test_structural_bounds.py::test_a_long_random_ltr_is_not_rejected_as_a_simple_repeat`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-005 — TIR search in flanked genomic records: seeded arms, joint trimming with the TSD, boundary_uncertain

**Before (4.0.0).** v4.0.0: find_tir first read the record ends (_terminal_tir, anchored at positions 0 and n); if that failed it seeded 11-mers within 1200 bp of each end (>= 3 anchors on one anti-diagonal), took the single densest diagonal and reported the extent of the anchor cluster (identity >= 80, not _is_simple) with no TSD test. In a flanked record the winner was often a 13-bp subterminal repeat inside the element, or an arm anchored in the flank: flanked Ac gave tir_len 13 at [element+113, element end-82], TSD absent. Short or imperfect TIRs (hAT 11 bp at 90.9%) had no 11-mer seed and were not found.

**After (4.1.0).** 4.1.0: the anchored read of the record ends is tried first and returns as in v4.0.0 (a chance repeat at the ends of a flanked record can stand in for the real one, about 1% of the time). Otherwise two seeded searches run, one candidate per diagonal: (a) exact 11-mers (>= 3 clustered on a diagonal, i.e. a ~13-bp perfect core) within 1200 bp of each end, arms measured up to 1000 bp; (b) exact 6-mers (>= 2 on a diagonal) within 400 bp of each end, arms measured up to 60 bp. Each diagonal is walked outward (X-drop 6), cut of any tandem run (STR-009), then scanned with an identity-aware X-drop of 60. For every candidate the arm pair is re-placed jointly with the TSD over shifts from -2 to +6 (trim at most 6 inward, restore at most 2 outward) using _snap_choice with a 5-bp minimum. A candidate is accepted when it is 'supported' (a flanking TSD of at least 6 bp that is not inside a microsatellite) or 'strong' (an exact core of at least 13 bp). Supported candidates have their arms re-measured at 90% identity (in place of the 80% floor); unsupported ones are reported as their exact core. Candidates are ranked by (supported and core at least 20 bp, core at least 20 bp, supported, TSD length, strong, length), so a corroborated terminus outranks a longer subterminal repeat. Rows carry boundary_adjustment [a,-a] when trimmed, boundary_tests when TSD-supported, and boundary_uncertain=true when not. The method string gains the window ('k-mer seed vs reverse-complement (k=11, 1200 bp window)'). length_is_lower_bound is set only when the arm reached the scan limit.

**Trigger.** Records that have flank on both sides of a TIR element (any genomic window), where the anchored read of the record ends finds nothing; TIR elements with a 6-bp or longer TSD benefit most. Bare elements (the record starts and ends on the element) are unchanged.

**Affected outputs.** TIR row element_span, tir_len, five_prime, three_prime, identity; TIR row new keys boundary_uncertain, boundary_adjustment, boundary_tests; TIR row method text; TSD row presence and coordinates for TIR elements; tsd_state on TIR rows; classification of flanked TIR elements through tir_len versus the 20-bp standalone floor.

**Labels.** added: boundary_uncertain (TIR row), boundary_tests (TIR row), boundary_adjustment (TIR row), method 'k-mer seed vs reverse-complement (k=11, 1200 bp window)'; renamed: method 'k-mer seed vs reverse-complement (k=11)' -> 'k-mer seed vs reverse-complement (k=11, 1200 bp window)'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _TIR_SEED_WINDOW | 1200 (inline literal) | 1200 (named constant) | yes | value unchanged; now named and sealed |
| _TIR_JOINT_WINDOW | absent | 400 | no | new; the 6-mer joint search window |
| _TIR_JOINT_K | absent | 6 | no | new |
| _TIR_JOINT_ANCHORS | absent | 2 | no | new; 6-mer seeds needed on one diagonal |
| _TIR_JOINT_MAX | absent | 60 | no | new; longest arm measured for a 6-mer-seeded candidate |
| _TIR_WEAK_TSD | absent | 6 | no | new; shortest flanking TSD that corroborates a short or imperfect candidate |
| _TIR_WEAK_IDENTITY | absent | 90.0 | no | new; arm identity for TSD-corroborated candidates (80% floor let chance TIR + chance TSD reach 19 bp) |
| _TIR_OUT_MAX | absent | 2 | no | new; bases the TSD may restore outward |
| _TIR_TRIM_MAX | absent | 6 | no | new; bases the TSD may trim inward |
| _TIR_XDROP | absent | 60 | no | new; scan stop for a seeded candidate |
| _TIR_MAX_SEEDS | absent | 64 | no | new; diagonals kept per seed search |
| MIN_TIR_LEN | 10 (find_tir min_tir default) | 10 (named constant) | yes | value unchanged |
| find_tir k / max_tir / min_anchors | 11 / 1000 / 3 | 11 / 1000 / 3 | yes | signature defaults unchanged |
| MIN_TIR_IDENTITY | 80 | 80 | yes | frozen; TSD-corroborated candidates are additionally re-measured at 90.0 |
| MIN_TIR_STANDALONE | 20 | 20 | yes | frozen; chance ceiling re-measured under the flanked search (see STR-020) |

**Measured.** Own run: X05424.1 (Ac) with 300-bp random flanks and a random 8-bp TSD, 30 seeds: exact span plus TSD 0/30 (v4.0.0, reports 13-bp subterminal TIR) -> 25/30. 10 seeds per flank length, TSD 8: exact 10/10 at flank 300, 350, 380 and 0/10 at 400, 450, 600 (joint window); v4.0.0 0/10 everywhere. Perfect 20-bp TIR with no TSD (10 seeds): exact 7/10 -> 6-7/10 up to flank 1150, 0/10 at 1250 in both versions. Synthetic 13-bp perfect TIR, 8-bp TSD, 2.5 kb interior, 300-bp flanks, 20 seeds: found 19 -> 20, exact span 12 -> 19. 11-bp arm with 1 mismatch and 8-bp TSD: found 1/20 -> 15/20 (exact 1 -> 13); with no TSD 1 -> 2 found, exact 0 -> 0. Engine, Ac seed 3: v4.0.0 TIR 13 bp [421,4791], tsd_state absent; 4.1.0 TIR 11 bp [308,4873] with an 8-bp TSD, tsd_state found. 200 random 4-kb records at 50% GC (seeds 1000-1199): TIR rows 1 -> 3, longest 17 bp in both, 1 boundary_uncertain; at 35% GC: 10 -> 11 rows, longest 14 in both, 8 boundary_uncertain. Bare specimens X01005.1 (TIR 54 @100) and X05424.1 (TIR 11 @90.9) are identical. Bare terminal-path controls (synthetic 11/1 mm, 13/0, 20/2, 30/3, 15/3, 20 seeds): found 20/20 with identical lengths in both. Wave record: chance-TIR ceilings 15, 17, 16, 17, 15 bp across five seeds; 17 bp over 1,800 sequences (44 hits).

**Residual limitations.** Short or imperfect TIRs are searched only within 400 bp of each record end; a 13-bp exact core within 1200 bp. A flank beyond these is not searched (measured: Ac with 8-bp TSD exact up to flank 380, none at 400). A TSD shorter than 6 bp can never corroborate: families with 2-4 bp TSDs (Tc1/Mariner TA, piggyBac TTAA, and 3-bp CACTA or PIF/Harbinger in typical descriptions) are reported only as an uncorroborated exact core (STR-006) unless classification recovers the target site (STR-007). With no TSD, the 13-bp subterminal repeat still outranks the 11-bp Ac terminus: Ac with 300-bp flanks and no TSD gives tir_len 13 at element offset +113 in both versions, now flagged boundary_uncertain. A chance inverted repeat at the very ends of a flanked record is read by the anchored path and can stand in for the real TIR (about 1%). A 12-bp (TA)6 'TSD' inside a flanking (TA)n run is still emitted by detect_all with chance_level false (the microsatellite-TSD guard acts only on corroboration inside the TIR ranking). Chance TIR rows on random 4-kb records rose slightly (1 -> 3 per 200 at 50% GC), all below 20 bp and flagged boundary_uncertain.

**Benchmark probe.** Positive: X05424.1 (Ac) flanked by 300-bp random sequence with a random 8-bp TSD on each side, 30 seeds; score element_span and TSD. v4.0.0: TIR 13 bp inside the element, tsd_state 'absent', 0/30 exact; 4.1.0: TIR 11 bp, TSD 8 bp, ~25/30 exact. Synthetic: 11-bp arm with one mismatch + 8-bp TSD + 2.5 kb random interior: v4.0.0 found ~1/20, 4.1.0 ~15/20. Negative: random 4-kb sequence (50% GC, 200 seeds): no TIR longer than 17 bp in either version; 4.1.0 rows carry boundary_uncertain. Same Ac with no TSD: v4.0.0 and 4.1.0 both give the 13-bp subterminal row, 4.1.0 adds boundary_uncertain=true and tsd_state 'not assessed ...'. Boundary: Ac with 8-bp TSD at flank 380 (exact) versus 400 (not found).

**Pinned by.** `tests/test_41_structural.py::test_flanked_element_boundaries_and_tsd[X05424-TIR-0.85]`, `tests/test_41_structural.py::test_imperfect_tir_with_a_central_mismatch_is_found_when_a_tsd_corroborates_it`, `tests/test_41_structural.py::test_weak_tir_without_a_tsd_is_not_called_in_flanked_records`, `tests/test_41_structural.py::test_subterminal_inverted_repeat_does_not_outrank_the_corroborated_termini`, `tests/test_41_structural.py::test_tsd_corroborated_tir_row_counts_its_boundary_placements`, `tests/test_41_structural.py::test_chance_tir_ceiling_with_the_flanked_search_stays_below_the_standalone_floor`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-006 — Uncorroborated flanked TIR candidates reported as their exact core: long imperfect arms are shortened

**Before (4.0.0).** v4.0.0: in a flanked record a seeded TIR candidate was reported with the extent of its anchor cluster whether or not a TSD supported it. A flanked 60-bp arm with 10% mismatches came back near 53 bp (median) and a 200-bp arm with 10% mismatches near 186 bp, with start offsets of about 1 bp from the true end in most records. tsd_state read 'absent (flanks present, no direct repeat found)'.

**After (4.1.0).** 4.1.0: a candidate with no TSD of 6 bp or longer is accepted only if it has an exact core of at least 13 bp, and is then reported as that exact core, never as an extended arm, so that a seeded search over a long flank cannot reach 20+ bp by chance (MIN_TIR_STANDALONE). The row has boundary_uncertain=true and the TSD search is not run against its boundary: tsd_state is 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'. Cores shorter than min_tir (10) are not reported. With a TSD of 6 bp or more the arms are extended and measured at 90% identity.

**Trigger.** Flanked records whose TIR element has no flanking TSD of at least 6 bp, i.e. (a) no TSD planted or lost, and (b) superfamilies with TSDs shorter than 6 bp, notably Tc1/Mariner TA, piggyBac TTAA, CACTA and PIF/Harbinger (3 bp in typical descriptions), unless a superfamily call later recovers the site (STR-007). Long imperfect arms (MULE/Foldback-like) are affected most.

**Affected outputs.** TIR row tir_len (lower, the exact core); TIR row element_span (can sit tens of bp inside the element: median start offset 56 bp for 200-bp arms); TIR row boundary_uncertain=true and tsd_state 'not assessed ...'; classification label where the core falls below 20 bp: 'non-autonomous TIR element (superfamily undetermined)' -> 'terminal inverted repeat, class unassigned'; identity of the reported row.

**Labels.** added: tsd_state 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'; renamed: tsd_state 'absent (flanks present, no direct repeat found)' -> 'not assessed (element boundary uncertain ...)' on uncorroborated TIR rows.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| exact-core acceptance without TSD | not applicable (cluster extent reported) | core >= min_anchors + k - 1 = 13 bp | no | derived from find_tir defaults (3 + 11 - 1); chance ~1% (benchmarks chance_tir) |
| _TIR_WEAK_TSD | absent | 6 | no | new; see STR-005 |

**Measured.** Own run, synthetic arm + 2.5 kb random interior + reverse-complemented arm with 10% substituted bases, 300-bp random flanks, detect_all, 20 seeds. No TSD: 30-bp arms median tir_len 21 -> 18, boundary_uncertain 0 -> 17/17; 60-bp arms median 53 -> 23, exact span 8 -> 2 of 20, start offset median 1 -> 10; 200-bp arms median 186 -> 43, exact span 4 -> 0 of 20, median start offset 1 -> 56. With an 8-bp TSD the result reverses: 60-bp arms median 60 and 200-bp arms median 200, exact in nearly all records (boundary_uncertain 1/20 each). Through engine.analyze (6 seeds, no TSD): 60-bp arms 'non-autonomous TIR element (superfamily undetermined)' 6/6 in v4.0.0 -> 4/6, with 2/6 'terminal inverted repeat, class unassigned' (cores 19 bp); 30-bp arms: one record moved from 'non-autonomous TIR element' to 'class unassigned' (core 16 bp); 200-bp arms stay 'non-autonomous TIR element' (cores 32-70 bp, all boundary_uncertain, lengths 175-201 -> 32-70). Confidence stays Candidate in all.

**Residual limitations.** Long imperfect TIRs in records with flank and no TSD are under-measured by up to a factor of four and their boundary sits inside the element; v4.0.0 measured them close to truth in these constructs (and equally by luck in others). The 6-bp TSD requirement cannot be met by 2-4 bp TSD families. The exact-core length is a lower bound on the arm length; the row does not say so beyond boundary_uncertain.

**Benchmark probe.** Positive: build element = arm(200 bp random) + 2.5 kb random + reverse complement with 10% substitutions, flank with 300-bp random sequence on both sides and no TSD, 20 seeds. v4.0.0: tir_len median ~186, start offset ~0-1 bp; 4.1.0: tir_len median ~43, start offset median ~56 bp, boundary_uncertain true, tsd_state 'not assessed ...'. Negative: same element with an 8-bp random TSD on both sides: 4.1.0 tir_len ~200, exact span, no boundary_uncertain. Boundary: 60-bp arms at 10% mismatches, no TSD: 4.1.0 cores ~14-36 bp, some below 20 bp so the label changes to 'terminal inverted repeat, class unassigned'; with a 5-bp TSD (below the 6-bp rule) the uncorroborated path applies.

**Pinned by.** `tests/test_41_structural.py::test_uncorroborated_arms_are_reported_as_uncertain_and_the_tsd_state_says_so`, `tests/test_41_structural.py::test_tir_lengths_are_not_extended_by_chance_when_the_arms_have_no_corroboration`, `tests/test_41_structural.py::test_fast_live_chance_tir_check_at_the_recorded_seed`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### STR-007 — Superfamily-gated recovery of the palindromic TA target site on Tc1/Mariner elements

**Before (4.0.0).** v4.0.0: on a flanked Tc1/Mariner insertion the TIR arms ran 2 bp through the palindromic TA target site, so tir_len came out 54 + 2..4 (56-58) with element_span widened by 2 on each side, and the TSD row was absent or a chance 3-bp match (e.g. GTA/TGT) judged incongruent; no helper existed.

**After (4.1.0).** 4.1.0: structural.palindromic_target_site(seq, tir, motif='TA') trims the smallest number of outer bp (0 to 6, never below 10 bp of arm) that exposes the motif at both ends and returns (trimmed TIR row, TSD row via find_tsd with expect=len(motif) and placements = trim + 1). It is not called inside detect_all because trimming palindromic outer bases would manufacture a TSD on about one TIR in four; classify calls it only after an unambiguous Tc1/Mariner call and updates the TIR row (tir_len, five/three_prime, element_span, boundary_adjustment, identity) and inserts or overwrites the TSD row. detect_all alone still reports TIR 56 and no corrected TSD.

**Trigger.** Flanked records classified as Tc1/Mariner (TIR plus a Tc1/Mariner transposase profile); bare records and other superfamilies unchanged.

**Affected outputs.** TIR row tir_len (56-58 -> 54 on X01005.1); TIR row element_span and boundary_adjustment; TSD row: 2-bp TA with motif TA; TSD row chance_level (true for a 2-bp match); evidence wording about the insertion site (hedged because a 2-bp TSD is at chance level); confidence/ledger recomputed from the trimmed arms.

**Labels.** added: boundary_adjustment on a recovered TIR row.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| palindromic_target_site max_trim | absent | 6 (= _TIR_TRIM_MAX) | no | new; bound on bases the arms may have run through the site |
| TSD_EXPECT Tc1/Mariner | (2, basis) | (2, (2,2), basis, motif 'TA') | no | shape changed, preferred length 2 and TA motif unchanged in meaning |

**Measured.** Own run: X01005.1 (Tc1) with 300-bp random flanks and a TA TSD, 30 seeds: detect_all alone exact 0/30 in both versions (TIR 56 at offsets -2/+2). palindromic_target_site on seed 3: TIR 56 [300,1914] -> 54 [302,1912], TSD 'TA'. Through engine.analyze (seeds 3 and 4): v4.0.0 TIR 56 with a 3-bp chance TSD (or absent); 4.1.0 TIR 54 at [302,1912] with a 2-bp TA TSD (chance_probability 0.18, chance_level true). Wave records: 59/60 recovered exactly by the helper; 30/30 through classify on the seeded flanked Tc1 test.

**Residual limitations.** Only effective after a Tc1/Mariner call; an unflagged Tc1-like element keeps arms 2 bp too long. The recovered 2-bp TA TSD is always chance level (p about 0.18), so it is hedged rather than credited as insertion-site evidence. tsd_state on the corrected TIR row is not refreshed: engine seed 4 shows tsd_state 'not assessed (element boundary uncertain...)' beside a 2-bp TSD row, and a row can carry boundary_uncertain=true with tsd_state 'found' (Tc1 in a (TA)n run). A (TA)n stretch shorter than about 10 bp beside the element is not cut and can lengthen the arms by a few bp.

**Benchmark probe.** Positive: X01005.1 flanked by 300-bp random sequence with TA duplicated at both ends, run the full analysis (not detect_all). v4.0.0: tir_len 56, TSD 3 bp or absent; 4.1.0: tir_len 54, element_span [302,1912]-type exact, TSD 2 bp TA. Negative: bare X01005.1 (identical, TIR 54, tsd_state 'not assessable'); a hAT (X05424.1) flanked record must receive no TA trim. Boundary: Tc1 flanked with a non-TA 2-bp site (e.g. GC GC): the motif gate gives no trim and the congruence verdict is incongruent (STR-012).

**Pinned by.** `tests/test_41_structural.py::test_flanked_tc1_target_site_recovered_for_the_superfamily_that_owns_it`, `tests/test_41_structural.py::test_palindromic_target_site_is_not_invented_for_other_elements`, `tests/test_41_classify.py::test_flanked_tc1_reports_the_trimmed_tir_and_its_ta_target_site`, `tests/test_41_classify.py::test_target_site_trim_is_gated_on_the_tc1_call`, `tests/test_41_classify.py::test_bare_specimens_keep_their_anchors_through_the_tc1_wiring`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### STR-008 — Microsatellite arm and microsatellite-TSD guards for TIR candidates (periods 1-6, self-paired units)

**Before (4.0.0).** v4.0.0: a TIR candidate was rejected only if its 5' arm was _is_simple (fewer distinct 4-mers than the floor). A diverged tandem repeat passes that test: a 35-bp (TA)n arm with one substitution, a 10-bp CTGCAGCAGC arm, and (AATT)n, (GATC)n, (ACGT)n arrays all pair with their own reverse complement, so they were reported as TIRs (arms up to the 1-kb limit) and classified 'non-autonomous TIR element'; a TSD formed by the run also counted.

**After (4.1.0).** 4.1.0: _simple_arms rejects a candidate when the 5' arm is _is_simple, or when any of six windows (5' arm, 3' arm, each arm with 16 bp outside it, each arm with 16 bp inside it) is a tandem run by _period_run. _period_run tests periods 1-2 for windows under 20 bp and periods 1-6 from 20 bp, allowing at most 15% broken comparisons (30% from 20 bp); a period-4-6 reading counts only when its unit is its own reverse complement up to rotation (GATC, AATT, AAATTT) or repeats a unit of period 1-3, so a unit that cannot pair with itself, such as the Ginger1-type (TCTA)n terminal repeat, is kept. _tsd_in_microsatellite zeroes a candidate's TSD when a TSD of 6 bp or more lies in a tandem run (any unit, period 1-6) that continues 8 bp outside it, so it cannot corroborate. The same guard applies on the anchored path (_scan_tir) and on every seeded path.

**Trigger.** Records holding or flanked by microsatellite or minisatellite arrays (including diverged ones up to about 10%), and elements whose terminal repeat begins with a self-paired run.

**Affected outputs.** TIR rows disappear for such arrays and flanks; classification: 'non-autonomous TIR element' 31 -> 0-1 on diverged arrays; TIR row arms no longer extended through a flanking microsatellite; tsd_state/TSD corroboration for candidates whose TSD is a run.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _period_run tolerance and window | absent | 0.15 (0.30 from 20 bp); periods 1-6 from 20 bp | no | new; inline literals, not sealed; false rejection of random arms measured (see measured) |
| _simple_arms reach (outside/inside context) | absent | 16 bp | no | new; sealed as structural_simple_arms_reach (outside); the 16-bp inside context is an inline literal |
| _tsd_in_microsatellite reach | absent | 8 bp | no | new |
| _scan_tir min_ident | 80 (literal in _terminal_tir) | 80 | yes | MIN_TIR_IDENTITY unchanged |

**Measured.** Own run, 3-kb arrays, 5 seeds (0-4) per cell, detect_all, TIR rows: AT 10% 3 -> 0, CAG 10% 1 -> 0, GATC 5% 5 -> 0, AATT 5% 1 -> 0, AATT 10% 4 -> 0; zero rows at 2% in both. Wave records: 160 diverged arrays through the engine: TIR rows 51 (intermediate pre-guard build) -> 6 (all at 20% divergence, arms 10-19 bp); 'non-autonomous TIR element' 31 -> 0; 60 palindromic arrays 52 -> 0; five specimens between 2-kb diverged (AT)n flanks: TIR rows 60 -> 0 (arms no longer lie in a flank). False rejection of random windows by _period_run, GC 0.5 / 0.35 / 0.2: 10 bp 0.04% / 0.12% / 0.67%; 20 bp 0.08% / 0.17% / 1.56%; of a real random arm by _simple_arms with 16 bp random sides: 10 bp 0.25% (GC 0.5), 0.73% (GC 0.35), 4.31% (GC 0.2). Controls: (TCTA)6 terminal repeat kept (24 bp at 100%, both versions); a 10-bp CAGCTGCAGC terminal arm on a bare record kept (11 bp, both).

**Residual limitations.** Arrays at 15-20% divergence still yield short (10-19 bp) TIR rows reported as class unassigned (6 of 160 arrays; AT 5, CAG 1). The guard rejects random arms in AT-rich sequence more often than v4.0.0 (at 20% GC, about 3.3% of 10-bp and 1.0% of 20-bp arms against 2.9% and 0.8% in the preceding build). A terminal repeat whose own outer end is a self-paired run cannot be told from an insertion into such a run; the run is cut (STR-009). A microsatellite arm at the far ends of a flanking region of two different non-pairing runs, e.g. (CA)n on one side and (TG)n on the other, still extends the arms into the runs. Tolerances and the 16-bp inside window are inline literals and are not sealed.

**Benchmark probe.** Positive: 3,000-bp (AT)n, (CAG)n, (GATC)n and (AATT)n arrays with 5% and 10% random substitutions, 5 seeds each. v4.0.0: TIR rows in several cells (e.g. GATC 5%: 5/5, AATT 10%: 4/5, AT 10%: 3/5) with class 'non-autonomous TIR element'; 4.1.0: no TIR row. Specimen X01005.1 between two 2-kb (AT)n flanks at 5%: v4.0.0 TIR arm inside the flank; 4.1.0 none from the flank. Negative: a real random 20-bp TIR with 2.5 kb interior must still be found; a (TCTA)6 terminal repeat must be kept. Boundary: arrays at 20% divergence retain 10-19 bp rows in 4.1.0 (documented residual); false-rejection rate of random 10-bp arms at GC 0.2 is about 4%.

**Pinned by.** `tests/test_structural_bounds.py::test_diverged_microsatellite_arrays_yield_no_terminal_inverted_repeat`, `tests/test_structural_bounds.py::test_specimens_between_diverged_microsatellite_flanks_carry_no_arm_inside_a_flank`, `tests/test_structural_bounds.py::test_arms_that_are_a_microsatellite_with_a_few_substitutions_are_not_a_terminal_repeat`, `tests/test_structural_bounds.py::test_a_tsd_inside_a_microsatellite_run_does_not_corroborate_a_short_terminal_repeat`, `tests/test_structural_bounds.py::test_palindromic_period_4_to_6_arrays_yield_no_terminal_inverted_repeat`, `tests/test_structural_bounds.py::test_specimens_between_diverged_palindromic_flanks_carry_no_arm_inside_a_flank`, `tests/test_structural_bounds.py::test_short_arm_at_the_edge_of_a_record_cut_from_a_microsatellite_array_is_not_a_terminal_repeat`, `tests/test_structural_bounds.py::test_the_tandem_run_test_rarely_rejects_a_random_window`, `tests/test_structural_bounds.py::test_a_terminal_repeat_made_of_a_unit_that_cannot_pair_with_itself_is_kept`, `tests/test_structural_bounds.py::test_a_tsd_inside_any_tandem_run_is_uninformative`, `tests/test_structural_bounds.py::test_simple_repeats_are_not_called_terminal_repeats`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-009 — Tandem run cut from TIR arms instead of discarding the candidate (transposon inserted next to a (TA)n run)

**Before (4.0.0).** v4.0.0: the seeded TIR arms ran outward through the target site into an adjacent (TA)n run, which pairs with the other half of a split run, so a Tc1 inserted at TA next to a run reported arms of 75, 115 or 154 bp for runs of 20, 60 and 100 bp, element_span widened to include part of each run, and no TSD.

**After (4.1.0).** 4.1.0: after the outward walk each seeded candidate has its outer ends read as a tandem run of period 1-6 (+1/-2 comparison walk, X-drop 14, started up to 40 bp outside the arm, run allowed to start within 10 bp of the arm end, smallest period within 90% of the best score, net score 10 or more). If the run's unit is its own reverse complement up to rotation (self-paired: TA, GATC, AATT) it is cut from both ends (the longer of the two sides), arms are then measured on the trimmed span and checked with reach 0. A run that fills 1 kb (max_tir) or half the span drops the candidate. A run of an unpaired unit (A(n), (CA)n) at both arm ends is treated as the element's own.

**Trigger.** Flanked records in which the element sits next to a self-paired tandem run of net score 10 or more at the arm end: Tc1/Mariner next to (TA)n is the realistic case.

**Affected outputs.** TIR row tir_len, element_span (back to the true arm, e.g. 54 bp for Tc1); TSD row presence (TA or the run-derived TSD), tsd_state; classification confidence for Tc1 next to a run (Candidate returns to Moderate).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _RUN_TRIM_SCORE | absent | 10 | no | new; net matching comparisons that define a run |
| _RUN_TRIM_SLACK | absent | 10 | no | new; run may start this far in |
| _RUN_TRIM_OUT | absent | 40 | no | new; bases outside the arm read with it |
| _run_extent x_drop | absent | 14 | no | new |

**Measured.** Own run, X01005.1 with 100-bp random flanks outside a (TA)n run of 20, 60 and 100 bp on each side (detect_all): v4.0.0 TIR 75, 115, 154 bp with spans widened into the runs and tsd_state 'absent'; 4.1.0 TIR 54 bp at the element, with a TSD row. Engine, run 60: v4.0.0 'Tc1/Mariner' Moderate, TIR 114 bp, no TSD; 4.1.0 'Tc1/Mariner' Moderate, TIR 54 bp at [160,1770], 2-bp TA TSD. Wave record: runs of 20, 40, 60, 100, 200 bp at 0% divergence: 30/30 exactly 54 bp in every cell; at 5% divergence 28, 30, 30, 29, 27 of 30 exact. In the intermediate build before this trim the same inputs lost the TIR entirely (27/30 Candidate); 4.1.0 vs v4.0.0 is the comparison shown here. 905-record corpus and 200 random-flank specimens: no differences from the trim.

**Residual limitations.** When the TIR's own outer end is a self-paired run ((TA)n, (GATC)n) the run is cut from the reported arms: the element's own run cannot be told apart from an insertion into a run. When a run of another unit begins the element's terminal repeat and sits inside a (TA)n run, the longer run decides: a longer (TA)n run is cut together with the element's own, a longer element run leaves (TA)n in the arms or loses the TIR. A (TA)n stretch shorter than about 10 bp is not cut. detect_all can still emit a 12-bp (TA)6 'TSD' inside the flanking run with chance_level false, and the TIR row carries boundary_uncertain=true with tsd_state 'found' (Tc1 run 60 through the engine). A tandem run of 1 kb or more beside an element removes its candidate.

**Benchmark probe.** Positive: X01005.1 with TA at both ends flanked by (TA)30 (60 bp) on each side and 100-bp random sequence outside. v4.0.0: TIR 114-115 bp, span [99,1831]-type (through the run), no TSD; 4.1.0: TIR 54 bp, span [160,1770]. Repeat with run lengths 20, 100, 200: 4.1.0 54 bp throughout. Negative: Ac X05424.1 between the same runs (no TIR is expected from the run; 4.1.0 reports the uncorroborated 13-bp subterminal row). Boundary: a (TA)n stretch of 6-8 bp beside the element is not trimmed, so tir_len may exceed 54 by a few bp.

**Pinned by.** `tests/test_structural_bounds.py::test_a_tc1_next_to_a_ta_run_keeps_its_own_terminal_repeat`, `tests/test_structural_bounds.py::test_a_terminal_repeat_that_begins_with_a_simple_run_keeps_its_full_length`, `tests/test_structural_bounds.py::test_only_a_run_whose_unit_is_its_own_reverse_complement_is_trimmed`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-011 — TSD rows carry a composition-adjusted chance probability; chance_level flag; ACGT-only TSD

**Before (4.0.0).** v4.0.0: a TSD row had type, length, motif, upstream, downstream and matched_expected only. A 2-bp match and a 12-bp match were indistinguishable by any reported number; any equal flank strings without N were accepted.

**After (4.1.0).** 4.1.0: find_tsd(placements=1) adds chance_probability, n_tests and chance_level. chance_probability = 1 - [prod over lengths l from L to max_tsd of (1 - (sum of squared base frequencies of the whole record)^l)]^placements, with a single length when the expected length matched first; n_tests = lengths x placements; chance_level = chance_probability > CHANCE_LEVEL_P (1e-3). placements is 1 for a fixed span, boundary_tests for a snapped LTR or a TSD-corroborated TIR, and trim + 1 for the palindromic target-site recovery. At uniform composition chance_level is true for TSDs of 2-5 bp (e.g. 0.085, 0.023, 0.0061, 0.0014 for 2, 3, 4, 5 bp when searched without an expected length) and false from 6 bp (0.00037); AT-rich records read higher. Classification no longer says the insertion site is captured unless chance_level is false (a row without the field is treated as chance level).

**Trigger.** Every TSD row. The chance_level flag affects wording and scoring wherever a consumer reads it.

**Affected outputs.** TSD row new keys chance_probability, n_tests, chance_level; 'insertion site captured' wording in the TIR evidence (classify); TSD row absence for IUPAC-only matches.

**Labels.** added: chance_probability (TSD row), n_tests (TSD row), chance_level (TSD row).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| CHANCE_LEVEL_P | absent | 1e-3 | no | new; a single pre-specified 5-bp match is 9.8e-4 at uniform composition; unspecified-length matches and snapped boundaries need more |
| MIN_TSD / MAX_TSD | 2 / 12 | 2 / 12 | yes | frozen search range |

**Measured.** Own run, find_tsd on random flanks with an exact planted TSD (seed = length; 60-bp random pads; ~uniform composition): chance_probability / chance_level / n_tests for L = 2: 0.0852 / true / 11; 3: 0.0226 / true / 10; 4: 0.0061 / true / 9; 5: 0.00142 / true / 8 (expect=5: 0.00106 / true); 6: 0.000368 / false / 7; 8: 2.1e-05 / false / 5; 12: 6.5e-08 / false / 1. Engine on flanked specimens (300-bp random flanks): copia and gypsy 5-bp TSD chance_level true (p 0.00212 in detect_all); Ac 8-bp false (5.7e-4); L1.2 12-bp false (5.4e-4); Tc1 2-bp true (0.18). Random 4-kb records, seeds 1000-1199: 1 TSD row at 50% GC and 1 at 35% GC, none non-chance.

**Residual limitations.** Every genuine 4-6 bp LTR TSD and every 2-bp Tc1/Mariner TA reads chance_level true; only 6 bp or longer reads false. A scorer that counts 'TSD found' as chance_level == false will count true positives as failures. The chance model uses whole-record base composition and assumes independent bases, so a TSD inside a tandem run is scored as extremely unlikely: a 12-bp (TA)6 flank inside a (TA)n run gets chance 1e-7 and chance_level false (copia between (TA)30 runs: 2.4e-07). The 1e-3 level is applied per search, not corrected across the records of a run.

**Benchmark probe.** Positive: M11240.1 flanked with a random 5-bp TSD (300-bp random flanks): v4.0.0 TSD row without chance fields; 4.1.0 chance_probability ~2e-3, chance_level true, n_tests 8. Same with a 12-bp TSD (L1.2 M80343.1 construct): chance ~5e-4, false. Negative: random 4-kb record at 50% GC: no TSD row with chance_level false. Boundary: TSD length 5 (true) versus 6 (false) at uniform composition; (TA)n flank run of 12 bp (chance_level false although the repeat is uninformative).

**Pinned by.** `tests/test_41_structural.py::test_tsd_rows_carry_composition_adjusted_chance`, `tests/test_41_structural.py::test_tsd_chance_is_adjusted_for_the_lengths_and_boundaries_tested`, `tests/test_41_classify.py::test_insertion_site_captured_requires_a_tsd_that_is_not_chance_level`, `tests/test_41_classify.py::test_a_tsd_row_without_a_chance_estimate_is_not_credited_as_the_insertion_site`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### STR-012 — TSD congruence uses superfamily length ranges and a required motif

**Before (4.0.0).** v4.0.0: TSD_EXPECT mapped a superfamily to one length (Copia 5, Gypsy 5, hAT 8, Tc1/Mariner 2) and a TSD was congruent only at exactly that length; any 2-bp match was congruent for Tc1/Mariner; tsd_congruence took no motif.

**After (4.1.0).** 4.1.0: TSD_EXPECT holds (preferred length, accepted range, basis, required motif): Copia and Gypsy 5 (accepted 4-6 bp; Wicker et al. 2007, Ou et al. 2019), hAT 8 (8-8), Tc1/Mariner 2 (2-2, motif TA). tsd_congruence(length, superfamily, motif=None) returns congruent when the length is inside the range, adds expected_range, and when a motif is supplied and a motif is required it returns incongruent unless the motif matches (expected_motif/observed_motif). The preferred length 'expected' is unchanged and still drives find_tsd(expect=) in classification. 'Not diagnostic' superfamilies and 'not assessed' (length None) are unchanged. classify passes the motif and quotes range and motif in the wording; its incongruent sentence no longer says 'likely coincidental'.

**Trigger.** Any LTR or Tc1/Mariner/hAT classification with a TSD.

**Affected outputs.** tsd_congruence verdict for 4 and 6 bp on Copia/Gypsy: incongruent -> congruent; verdict for non-TA 2-bp TSD on Tc1/Mariner: congruent -> incongruent (when the motif is passed); new keys expected_range, expected_motif, observed_motif; evidence sentences quoting the accepted range.

**Labels.** added: expected_range, expected_motif, observed_motif.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| TSD_EXPECT Copia/Gypsy | 5 bp, exact | 5 preferred, 4-6 accepted | no | Wicker et al. 2007 Nat Rev Genet 8:973-982; Ou et al. 2019 Genome Biol 20:275 |
| TSD_EXPECT hAT | 8 bp | 8 (8-8) | yes | unchanged value, Ou et al. 2019 |
| TSD_EXPECT Tc1/Mariner | 2 bp | 2 (2-2) with motif TA | no | motif requirement is new; length unchanged |

**Measured.** Own run, tsd_congruence: (4, Copia) v4.0.0 incongruent -> 4.1.0 congruent; (6, Gypsy) incongruent -> congruent; (7, Copia) incongruent both; (5, Copia) and (8, hAT) congruent both; (3, hAT) incongruent both; (3, CACTA) not diagnostic both; (2, Tc1/Mariner, motif GC): v4.0.0 TypeError (no motif argument), 4.1.0 incongruent with expected_motif TA, observed_motif GC.

**Residual limitations.** Without a motif only length is judged (backward-compatible path). Only four superfamilies have an expectation; CACTA, PIF/Harbinger, piggyBac etc. remain 'not diagnostic'. hAT's range is a single value; real hAT 8-bp TSDs are the documented norm but variants exist.

**Benchmark probe.** Positive: tsd_congruence(4, 'Copia (Ty1)') and (6, 'Gypsy (Ty3)'): v4.0.0 'incongruent', 4.1.0 'congruent' with expected_range [4, 6]. tsd_congruence(2, 'Tc1/Mariner', motif='GC'): v4.0.0 raises TypeError, 4.1.0 'incongruent'. Negative: (5, Copia) and (8, hAT) congruent in both. Boundary: 7 bp on Copia stays incongruent in both; (None, Copia) 'not assessed' with expected_range added.

**Pinned by.** `tests/test_41_structural.py::test_ltr_tsd_congruence_accepts_the_4_to_6_bp_range`, `tests/test_41_structural.py::test_tc1_tsd_congruence_requires_the_ta_motif_when_a_motif_is_given`, `tests/test_41_classify.py::test_the_tsd_motif_reaches_the_congruence_test`, `tests/test_41_classify.py::test_tsd_expectation_is_quoted_as_the_range_the_congruence_test_uses`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### STR-013 — tsd_state gains 'partially assessable' and 'not assessed (boundary uncertain)'

**Before (4.0.0).** v4.0.0: _tsd_state returned 'found', 'not assessable (no flanking sequence in the record)' (flank shorter than 2 bp on either side) or 'absent (flanks present, no direct repeat found)' (any flank of 2 bp or more with no TSD, including a flank of 3 bp).

**After (4.1.0).** 4.1.0: room = min(left flank, right flank). Below MIN_TSD (2): 'not assessable (no flanking sequence in the record)'. When the TIR row has boundary_uncertain: 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'. Between 2 and 11 bp: 'partially assessable (flanks of N bp allow only TSDs of up to N bp)'. From 12 bp (MAX_TSD): 'absent (flanks present, no direct repeat found)'. 'found' takes precedence over every other state when a TSD row exists.

**Trigger.** LTR and TIR rows in records with 2-11 bp of flank on either side, and TIR rows with boundary_uncertain.

**Affected outputs.** tsd_state string on LTR and TIR rows; GFF3 tsd_state attribute on long_terminal_repeat and terminal_inverted_repeat features (export scope).

**Labels.** added: 'partially assessable (flanks of N bp allow only TSDs of up to N bp)', 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'; renamed: 'absent (flanks present, no direct repeat found)' -> 'partially assessable (...)' when the smaller flank is 2-11 bp, 'absent (flanks present, ...)' -> 'not assessed (element boundary uncertain ...)' on uncorroborated TIR rows.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _tsd_state min_tsd | 2 | 2 (MIN_TSD) | yes | unchanged value |
| partial-assessability limit | absent | MAX_TSD = 12 | no | derived from MAX_TSD (the longest TSD searched) |

**Measured.** Own run, _tsd_state called directly on a 100-bp string with no TSD: flanks 1/1 and 0/0 'not assessable' in both versions; 3/3: 'absent' -> 'partially assessable (flanks of 3 bp ...)'; 8/8: 'absent' -> 'partially assessable (flanks of 8 bp ...)'; 12/12 and 30/30 'absent' in both. uncertain=True on 30/30: 'not assessed (element boundary uncertain ...)' (parameter absent in v4.0.0). Engine, Ac with no TSD at 300-bp flanks: 'absent' -> 'not assessed'.

**Residual limitations.** After classification recovers a Tc1/Mariner target site (STR-007) the TIR row's tsd_state is not refreshed: seed 4 shows 'not assessed' beside a 2-bp TSD row. A TIR row can carry boundary_uncertain=true and tsd_state 'found' at once (Tc1 inside a (TA)n run) because 'found' takes precedence.

**Benchmark probe.** Positive: copia with random flanks of exactly 3 bp each side and no TSD (or call the state function directly): v4.0.0 'absent (flanks present, no direct repeat found)'; 4.1.0 'partially assessable (flanks of 3 bp allow only TSDs of up to 3 bp)'. Negative: 30-bp flanks with no TSD stay 'absent' in both. Boundary: 11 bp -> partial, 12 bp -> absent; 1 bp -> not assessable in both.

**Pinned by.** `tests/test_41_structural.py::test_tsd_state_distinguishes_missing_partial_and_absent_flanks`, `tests/test_41_structural.py::test_uncorroborated_arms_are_reported_as_uncertain_and_the_tsd_state_says_so`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-014 — Poly-A tail detection: A-rich tract with interruptions and a 3-bp end slack

**Before (4.0.0).** v4.0.0: find_polya counted the run of A's that reached the last base of the record (>= 8) as the 3' tail, and the run of T's that began at the first base (>= 8) as the 5' poly-T. A single non-A base, or any 1-3 bases after the tail, ended the run and the tail was missed or truncated (L1.2 tail with one G: 13 bp instead of 27).

**After (4.1.0).** 4.1.0: the tail is the best-scoring A-rich tract (A +1, any other base -9, about one interruption per ten bases at the floor; the tract starts and ends on an A) in the last 600 bp, with at least 8 A's; it counts when it ends within 3 bases of the record end. The row length is the tract span (interruptions included), not the A count. Only the 3'-most qualifying tail is returned. The 5' poly-T of a reverse-oriented element is found by the same search on the reverse complement. find_polya still returns the row list; min_run (8) now counts A's in the tract.

**Trigger.** Records whose tail has a non-A interruption or 1-3 bases after it; bare tails that are a clean homopolymer at the end give the same row.

**Affected outputs.** poly-A tail (3') row length and pos; poly-T (5') row length and pos; classification credit for a tail (LINE confidence).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _TAIL_PENALTY | absent | 9 | no | new; about one interruption per ten bases at the purity floor |
| _TAIL_SLACK | absent (tail had to reach the last base) | 3 | no | new |
| _TAIL_WINDOW | absent | 600 | no | new; only the last (first, for poly-T) 600 bases are searched |
| POLYA_MIN_RUN / find_polya min_run | 8 | 8 | yes | unchanged value; semantics now A count within a tract |

**Measured.** Own run: M80343.1 (L1.2) bare tail 27 bp at [6023,6050] in both versions; with one A->G at the 14th tail position: v4.0.0 13 bp [6037,6050], 4.1.0 27 bp [6023,6050]. Reverse complement of L1.2: poly-T 27 bp at [0,27] in both. Synthetic: 20-A tail followed by 2 random bases (20 seeds): tail called 0/20 -> 20/20; followed by 4 bases: 0/20 -> 2/20 (slack is 3); 'ACGT...A20CGT' (3 trailing bases): none -> 21-bp tail at [79,100]. Clean A8 and A20 bare tails: 20/20 in both.

**Residual limitations.** A tail more than 3 bases from the record end counts only when corroborated (STR-015). The reported length includes interruptions, so it is not the A count. Only the 3'-most tail is returned.

**Benchmark probe.** Positive: M80343.1 with one A replaced by G in the tail: v4.0.0 poly-A 13 bp; 4.1.0 27 bp. A 20-A tail followed by 2 random bases: v4.0.0 none, 4.1.0 tail. Negative: bare M80343.1 identical (27 bp at [6023,6050]); a random record with an internal A14 run in flanks must give no tail in either. Boundary: trailing bases 3 (called) versus 4 (mostly not).

**Pinned by.** `tests/test_41_structural.py::test_bare_l1_tail_anchor`, `tests/test_41_structural.py::test_tail_tolerates_an_interruption_and_keeps_its_length_otherwise`, `tests/test_41_structural.py::test_reverse_oriented_element_has_a_poly_t_tract_and_its_tsd`.

Files: `app/backend/teagle_core/structural.py`.

#### STR-015 — Poly-A tail corroborated by a flanking TSD anywhere in the last 600 bp, with family-wise chance control

**Before (4.0.0).** v4.0.0: a tail followed by flank (a genomic window) was never found, so a flanked LINE had no tail row and no target-site duplication; classification fell to Candidate (flanked L1.2 + 12-bp TSD + 300-bp flanks: 0/60).

**After (4.1.0).** 4.1.0: a tail that does not end within 3 bases of the record end counts when it is followed by a direct repeat of at least 8 bp of the sequence before the element (the TPRT target-site duplication). The repeat is searched from 12 bp down to 8 bp, left copy within the first 600 bp of the record, the tail end tried at each of the last few A positions (longest duplication wins). It is accepted only if its family-wise chance, 1 - (1 - p)^n_tests with p from the record's composition and n_tests = tail ends tried x left-copy positions x flanked tracts, is at most 1e-3. The row gets tsd_corroborated=true, and the TSD is emitted as a TSD row (context 'flanks a poly-A tail (target-primed reverse transcription)', matched_expected false, chance_probability, n_tests, chance_level) only when no LTR or TIR row exists. The mirror search on the reverse complement gives the poly-T row and its TSD.

**Trigger.** Flanked records of non-LTR, non-TIR elements (LINEs) with a TSD of about 11-12 bp or longer and flanks of up to roughly 590 bp on the left and 600 bp on the right.

**Affected outputs.** poly-A tail (3') / poly-T (5') row presence and pos in flanked records; TSD row (new, tail-framing) with context; tsd_corroborated key; LINE classification confidence (Candidate -> Moderate on the flanked L1.2 construct).

**Labels.** added: tsd_corroborated (tail row), TSD row context 'flanks a poly-A tail (target-primed reverse transcription)'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| _TAIL_TSD_MIN | absent | 8 | no | new; seeds the search; in practice only 11-12 bp passes the family-wise gate |
| CHANCE_LEVEL_P (tail gate) | absent | 1e-3 | no | new; applied to the family-wise probability |
| MAX_TSD | 12 | 12 | yes | caps the tail TSD search |

**Measured.** Own run: M80343.1 + 12-bp random TSD + 300-bp random flanks, 30 seeds: tail row 0 -> 30, TSD row 0 -> 30, exact (tail end and 12-bp TSD at the true positions) 0 -> 29. Engine, seeds 3 and 4: 'LINE (non-LTR)' Candidate -> Moderate with poly-A 27 bp at [6335,6362] and a 12-bp TSD (chance 5.4e-4, chance_level false). 10-seed flank sweep: recovered 10/10 at left flank 100, 300, 590 (right 300) and at right flank 500; 0/10 at left flank 700 and at right flank 700, in both orientations (poly-T mirrored). TSD length (20 seeds, 20-A tail, 300-bp flanks, 1.5-kb random body): 12 bp 20/20, 8 bp 0/20, 6 bp 0/20 in both versions. Flank composition: 12-bp TSD with flank GC 0.50 / 0.35 / 0.30: 20/20, 20/20, 19/20 (the chance term uses the whole-record composition, here dominated by the L1 body). Controls: an internal A14 run in random flanks: no tail in either version; wave record: credited tails in 400 records without a planted TSD, with an internal A14 run, 16/400 and 43/400 (intermediate build) -> 0/400; recall for planted TSD of 8/10/11/12/15 bp at 50% GC: 3/59/99/100/100%.

**Residual limitations.** Only a TSD of about 11-12 bp or longer corroborates; the search is capped at MAX_TSD 12, so an 8-10 bp LINE TSD is not credited (0/20 at 8 bp). The left copy must start within the first 600 bases of the record and the tail within the last 600: a left flank over about 590 bp or a right flank over 600 bp hides the tail. The family-wise test uses whole-record composition: a record that is AT-rich overall (not only its flanks, about 30% GC or lower) cannot corroborate even a perfect 12-mer. The TSD framing a tail is emitted only when there is no LTR or TIR row. The test is per strand; a record searched on both strands carries up to about twice the per-strand false rate (0/400 measured per record).

**Benchmark probe.** Positive: M80343.1 with a random 12-bp TSD on both sides and 300-bp random flanks, 30 seeds. v4.0.0: no poly-A row, no TSD row, classification Candidate; 4.1.0: poly-A 27 bp ending at the TSD, TSD 12 bp with chance_level false, Moderate. Negative: random flanks of 300 bp around a body with an internal A14 run and no TSD: no tail in either version; the same with an 8-bp TSD: no tail in either. Boundary: left flank 590 bp (found) versus 700 bp (not found); TSD 12 bp (found) versus 10-11 bp (marginal) versus 8 bp (not found).

**Pinned by.** `tests/test_41_structural.py::test_flanked_l1_tail_and_its_target_site_duplication_are_found`, `tests/test_41_structural.py::test_uncorroborated_tract_inside_a_flank_is_not_a_tail`, `tests/test_41_structural.py::test_chance_direct_repeat_does_not_corroborate_a_tail`, `tests/test_41_structural.py::test_tail_tsd_row_carries_its_family_wise_chance`, `tests/test_41_structural.py::test_reverse_oriented_element_has_a_poly_t_tract_and_its_tsd`, `tests/test_41_classify.py::test_polya_signal_motif_row_is_not_counted_as_a_tail`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/teagle_core/classify.py`.

#### STR-016 — LTR cis-elements (PBS, PPT, poly(A) signal) read on the element's own strand; orientation from cis-elements or the domain strand

**Before (4.0.0).** v4.0.0: PBS, PPT and the poly(A) signal were always searched on the forward strand of the record, with no strand field and no statement of orientation. A reverse-complemented copia lost its PPT and reported an antisense poly(A) signal in the wrong LTR; forward searches on a record of unknown orientation produced chance rows.

**After (4.1.0).** 4.1.0: detect_all(seq, orientation=None) reads the cis-elements per orientation. A caller-supplied '+' or '-' is used as given; otherwise both strands are read and the one with more evidence wins, evidence being a PPT row, a confident PBS (>= 72%), or a poly(A) signal with its downstream element found; on a tie (including no evidence) the three rows are withheld and the LTR row carries cis_orientation 'undetermined' and cis_note. On '-' the detectors run on the reverse complement and are mapped back to record coordinates. Rows carry strand '+' or '-'; the accepted LTR row carries cis_orientation, cis_orientation_source ('caller' or 'cis-elements') and pbs_scope. In engine.analyze, when the LTR row is 'undetermined' after the domain scan and every RT/INT/TPase hit shares a strand, detect_all is re-run with that strand and the LTR row gets cis_orientation_source 'domain-strand'; an orientation already set by cis-elements is never overridden. find_ppt and find_polya_signal are unchanged (AST-identical) but now run on strand-specific sequence.

**Trigger.** Every record with an accepted LTR row. Reverse-oriented elements and elements with no clear PPT/PBS/PAS evidence change most.

**Affected outputs.** PBS, PPT, polyA-signal rows: presence, pos, and new strand key; LTR row new keys cis_orientation, cis_orientation_source, cis_note, pbs_scope; PPT/PAS rows on reverse-oriented elements (previously missing or on the wrong strand); PBS/PPT/PAS rows withheld when the orientation is undetermined; rec['processing']['structural'].count after a domain-strand re-run.

**Labels.** added: cis_orientation values '+', '-', 'undetermined', cis_orientation_source values 'caller', 'cis-elements', 'domain-strand', strand on PBS/PPT/polyA-signal rows, cis_note, pbs_scope.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| find_pbs strong threshold / min_ident / search | 0.72 / 0.55 / 44 | 0.72 / 0.55 / 44 | yes | unchanged; the 0.72 strong threshold is an inline literal |
| find_ppt window / min_len / min_purine / max_defects | 30 / 9 / 0.82 / 2 | 30 / 9 / 0.82 / 2 | yes | unchanged |
| find_polya_signal lo / hi / min_gt_frac / min_t | 20 / 60 / 0.65 / 4 | 20 / 60 / 0.65 / 4 | yes | unchanged; now sealed (polya_signal_dse_*) |

**Measured.** Own run: forward M11240.1 and M12927.1: PPT [4861,4870] / [6976,6987] and poly(A) signal [4954,4960] / [7138,7144] identical to v4.0.0 plus strand '+', LTR row cis_orientation '+', source 'cis-elements'. Reverse complement of M11240.1: v4.0.0 only an antisense polyA-signal at [4902,4908], no PPT; 4.1.0 PPT [276,285] strand '-' and polyA-signal [186,192] strand '-'. Reverse complement of M12927.1: v4.0.0 polyA-signal [7447,7453] only; 4.1.0 PPT [482,493] and polyA-signal [325,331], strand '-'. 200 random-interior LTR pairs (300-bp random LTR, 3 kb random interior, seeds 0-199): records with any PBS/PPT/PAS row 38/200 (PBS 12, PPT 6, PAS 23) in v4.0.0 versus 14/200 in 4.1.0 (PBS 2, PPT 7, PAS 8); cis_orientation 'undetermined' 186, '+' 9, '-' 5. Domain-strand path (engine.analyze): copia M11240.1 with the PPT and poly(A)-signal segments shuffled gave cis_orientation '+', source 'domain-strand' in 4.1.0; the PPT row was still reported there because the shuffled tract is still purine-rich.

**Residual limitations.** A chance poly(A)-signal with a downstream element on the opposite strand can set the orientation wrongly and block the domain strand: copia with 50-70 bp around the PPT and signal replaced by random C/T or C/A sequence read cis_orientation '-' (source 'cis-elements') with a '-' strand AATAAA row, although the element is forward; the domain strand was not consulted because the cis rows had already decided. The cis_orientation of '+' on bare, correctly oriented specimens rests on the PPT rather than the PBS; the PBS panel (tRNA-Lys3 only) is not an orientation witness unless confident. Elements without domain hits and without cis evidence keep 'undetermined' and no PBS/PPT/PAS rows, where v4.0.0 would have reported forward-strand rows (often chance rows). PBS rows on '-' elements are matched on the reverse complement; motifs are reported as read on the element's strand.

**Benchmark probe.** Positive: take M11240.1 and reverse-complement it. v4.0.0: LTR + polyA-signal only (antisense), no PPT; 4.1.0: LTR + PPT [276,285] strand '-' + polyA-signal [186,192] strand '-', LTR row cis_orientation '-'. Negative: forward M11240.1 gives the same PPT and signal coordinates in both versions; random-interior LTR pairs (200 seeds) give cis rows in about 19% (v4.0.0) versus 7% (4.1.0). Boundary: an LTR element whose PPT and signal are destroyed but with an RT/INT hit: 'undetermined' first, then 'domain-strand' source and restored rows in 4.1.0; if a chance AATAAA+DSE exists on the wrong strand the orientation is set wrongly by cis-elements (documented limitation).

**Pinned by.** `tests/test_41_structural.py::test_reverse_oriented_ltr_element_reports_mapped_cis_elements_on_the_right_strand`, `tests/test_41_structural.py::test_caller_orientation_wins`, `tests/test_41_structural.py::test_unknown_orientation_withholds_cis_elements_and_says_why`, `tests/test_41_engine.py::test_undetermined_ltr_orientation_is_taken_from_the_domain_strand`, `tests/test_41_engine.py::test_orientation_is_not_taken_from_a_split_or_irrelevant_domain_vote`, `tests/test_41_engine.py::test_a_cis_element_orientation_is_not_overridden_by_the_domain_strand`, `tests/test_41_engine.py::test_a_failed_orientation_rerun_keeps_the_first_structural_result`.

Files: `app/backend/teagle_core/structural.py`, `app/backend/engine.py`.

### Primer design and oligo QC

#### PRM-002 — IUPAC ambiguity codes in the template are masked to N for Primer3 (design no longer fails) and counted

**Before (4.0.0).** A template containing R, Y, S, W, K, M, B, D, H or V reached Primer3 unchanged and Primer3 raised 'Rejecting invalid sequence character (only ACGTN allowed)', which run_primers reported as BadRequest 'primer parameters rejected: OSError: ...'. A '-' or other non-nucleotide character produced the same misleading parameter error.

**After (4.1.0).** design_primers replaces RYSWKMBDHV by N in a copy of the template (same length, same coordinates), sets PRIMER_MAX_NS_ACCEPTED 0 explicitly so no primer covers a masked base, and returns template_ambiguity {masked_to_n, note}. run_primers adds a display cautions entry when bases were masked and seals max_ns_accepted and template_ambiguity_masked. Characters outside the IUPAC set raise primers.TemplateError, reported as a template-content BadRequest ('the template contains characters that are not nucleotide codes (...)') instead of a parameter error.

**Trigger.** Any primer-design template with at least one IUPAC ambiguity code (common in consensus sequences and in genome windows with N-gaps); any template with a non-nucleotide character (message change only).

**Affected outputs.** run_primers success vs BadRequest for ambiguous templates; result.template_ambiguity: absent -> {masked_to_n, note}; result.cautions: absent -> list (non-empty when bases were masked); candidate list: primers never cover a masked base, so pairs whose primers overlapped an ambiguity code are replaced; provenance.parameters.template_ambiguity_masked and max_ns_accepted (sealed).

**Labels.** added: template_ambiguity, cautions, parameters.max_ns_accepted, parameters.template_ambiguity_masked; renamed: BadRequest 'primer parameters rejected: OSError ...' -> 'the template contains characters that are not nucleotide codes (...)' for non-IUPAC characters.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| PRIMER_MAX_NS_ACCEPTED | not passed (Primer3 default 0) | 0 (explicit, sealed as max_ns_accepted) | yes | Primer3 default; made explicit so a masked base is never covered |
| ambiguity handling in the template | rejected with an error | RYSWKMBDHV -> N (count sealed) | no | changed behaviour; length and coordinates preserved |

**Measured.** Run in this audit on M11240, one R at 0-based position 300: 4.0.0 BadRequest 'primer parameters rejected: OSError: Rejecting invalid sequence character (only ACGTN allowed): R'; 4.1.0 returns 5 pairs, masked_to_n 1 and the caution '1 ambiguous template base(s) were treated as N for Primer3; no primer covers them.' With R at 4600 (inside the best clean pair's left primer, 4591-4611) 4.1.0 drops that pair (first pairs become 1125-1305 and 271-673) and no returned primer covers 4600. A '-' at 300 fails in both versions with different messages.

**Residual limitations.** A masked base can still lie between the primers inside the product. In-silico PCR still counts a template ambiguity code as a mismatch (PCR-002) and lists 'IUPAC-ambiguous template bases' under notSearched.

**Benchmark probe.** Positive: M11240 with one R at position 300 (or 4600): 4.0.0 BadRequest, 4.1.0 designs and reports masked_to_n 1. Negative: the unmodified record: masked_to_n 0, cautions [], identical candidates in both versions. Boundary: R at 4600, inside the top pair's left primer: 4.1.0's first pair changes from 4591-4611/4746-4766 to 1125-1145/1285-1305 and no primer overlaps position 4600; also a template with '-' fails in both versions but the 4.1.0 message names the character.

**Pinned by.** `tests/test_41_primers.py::test_iupac_template_designs_and_primers_avoid_the_ambiguous_bases`, `tests/test_41_primers.py::test_single_ambiguity_code_in_a_full_element_still_designs`, `tests/test_41_primers.py::test_clean_template_reports_no_ambiguity_masking`, `tests/test_41_primers.py::test_template_rejection_is_not_blamed_on_parameters`.

Files: `app/backend/teagle_core/primers.py`, `app/backend/engine.py`.

#### PRM-003 — API default product-size range changes from 120-700 to 150-500 bp

**Before (4.0.0).** DESIGN_DEFAULTS prod_min 120, prod_max 700, so run_primers with no params designed 120-700 bp products.

**After (4.1.0).** DESIGN_DEFAULTS prod_min 150, prod_max 500 (the Standard preset). The native panel always sends explicit values, so the change reaches only callers that omit params; the seal records the values used.

**Trigger.** run_primers or resolved_design_parameters called with no product-size params (scripts, benchmarks, tests); the GUI is unaffected.

**Affected outputs.** candidate list for a params-less run_primers call; provenance.parameters.prod_min/prod_max (120/700 -> 150/500).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| prod_min (default) | 120 | 150 | no | documented Standard preset |
| prod_max (default) | 700 | 500 | no | documented Standard preset |

**Measured.** Run in this audit: run_primers on M11240 with params {} returns pairs of 606, 607, 643, 644 and 175 bp in 4.0.0 (range 120-700) and 175, 180, 402, 403 and 167 bp in 4.1.0 (range 150-500, identical to the explicit Standard run).

**Residual limitations.** Any harness that relied on the implicit 4.0.0 range must now pass prod_min and prod_max explicitly to reproduce 4.0.0 designs.

**Benchmark probe.** Positive: run_primers(M11240, params={}) in both versions; 4.0.0 yields products up to 644 bp, 4.1.0 none above 500 bp. Negative: the same call with explicit prod_min 120 and prod_max 700 gives the 4.0.0 list in both versions. Boundary: prod_min 150/prod_max 500 explicit is identical in both.

**Pinned by.** `tests/test_41_primers.py::test_design_defaults_equal_the_standard_preset`.

Files: `app/backend/teagle_core/primers.py`.

#### PRM-004 — Minimum product size must be at least twice the maximum primer size; the primer form mirrors it and adds guards

**Before (4.0.0).** resolved_design_parameters accepted any prod_min >= 1 up to prod_max, so prod_min 30 with 18-27 nt primers returned pairs whose primers overlapped (products of 30-34 bp); each such pair later failed sealed in-silico PCR with 'sealed primer spans contradict oligo lengths or orientation'.

**After (4.1.0).** prod_min < 2 x max_size raises ValueError ('the minimum product size (N bp) must be at least twice the maximum primer size (M bp), otherwise the two primers can overlap and no PCR product can be called'); run_primers reports it as BadRequest. The primer form mirrors the rule before submitting and adds UI guards: primer sizes 18-36 nt (Primer3's limit), pairs 1-50, Tm 30-90 C, GC clamp not above the minimum primer size, plain messages with units, and unit-bearing field labels and tooltips (display-only; the 36-nt, 50-pair and Tm bounds are UI guards, not scientific limits).

**Trigger.** Any design with prod_min < 2 x max_size (default max_size 27 -> prod_min < 54); all four presets satisfy the rule (qPCR 70 >= 48, Permissive 100 >= 60).

**Affected outputs.** run_primers: candidates -> BadRequest for the affected range; primer form: error text before submission.

**Labels.** renamed: Min size -> Min size (nt), Max size -> Max size (nt), Return -> Pairs to return, Opt size -> Opt size (nt), Tm min -> Tm min (°C), Tm max -> Tm max (°C), Max poly-X -> Max poly-X (bases), GC clamp -> GC clamp (bases).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| minimum prod_min | 1 | 2 x max_size (54 at max_size 27) | no | guarantees non-overlapping primers (sealed PCR geometry rule) |
| max primer size accepted by the form (nt) | no cap | 36 | no | Primer3's built-in maximum |
| pairs returned accepted by the form | no cap | 1-50 | no | UI guard (each pair is folded and tabulated) |
| Tm accepted by the form (C) | no bounds | 30-90 | no | UI guard, not a scientific limit |

**Measured.** Run in this audit on M11240: prod 30-40 returns 5 overlapping pairs (30-34 bp) in 4.0.0 and BadRequest in 4.1.0; prod_min 40 with prod_max 200 returns pairs of 56-84 bp in 4.0.0 and BadRequest in 4.1.0 ('40 bp must be at least twice ... 54 bp').

**Residual limitations.** The 54-bp floor follows max_size, not the actual primer lengths, so a user with max_size 27 cannot request 40-bp products even with 18-nt primers unless max_size is lowered too.

**Benchmark probe.** Positive: prod_min 40/prod_max 200 on M11240: 4.0.0 returns 56-84 bp pairs, 4.1.0 rejects. Negative: prod_min 150/prod_max 500 identical. Boundary: prod_min 54 with max_size 27 accepted, 53 rejected; the same pair of values with max_size 24 accepts prod_min 48.

**Pinned by.** `tests/test_41_primers.py::test_product_range_that_forces_overlapping_primers_is_rejected_at_design`, `tests/test_41_primers.py::test_every_designed_pair_passes_sealed_pcr`, `tests/test_41_ui_M1c.py::test_primer_form_rejects_values_the_designer_would_refuse_with_a_plain_message`, `tests/test_41_ui_M1c.py::test_every_preset_passes_the_primer_form_rules`.

Files: `app/backend/teagle_core/primers.py`, `app/native/panel_primers.py`.

#### PRM-005 — The Design primers button now runs the selected assay goal; the default is a design inside the strongest domain

**Before (4.0.0).** designBtn was connected to _design: an untargeted Primer3 design over the whole analysed record (no target, no included region). A domain-specific design was reachable only from a domain row's right-click 'Design primers here', which sets the domain as Primer3's required target so the product flanks and spans it with the primers outside the domain.

**After (4.1.0).** designBtn is connected to _design_for_goal. A new Assay goal chooser (Detect or quantify this element / Genotype this insertion / Custom range) is shown above the preset. For a record that has a domain hit the default goal is Detect: the strongest HMMER-scoring domain (coords.domain_choices sorted by score) becomes Primer3's included region (primers and product both inside it, no target), clipped to the interior between an accepted terminal-repeat pair so the primers avoid both repeat copies. A record without a domain defaults to Custom range over the whole record, which is the old untargeted _design. The right-click 'Design primers here' is unchanged (routed, flank-and-span).

**Trigger.** GUI path only (the engine API still accepts the same primer3_target / included_region arguments). Affects every analysed record that has at least one domain hit: all five bundled examples.

**Affected outputs.** designed pairs: whole-record pairs -> pairs lying inside the default domain; design_assessment of a downstream in-silico PCR: the design carries included_region, no primer3_target, so target coverage is reported unassessed and the intended product is checked; provenance.parameters.included_region/primer3_target (sealed): null/null -> domain span/null; result note 'Designed inside the <domain> (N bp)...'.

**Labels.** added: Assay goal, Detect or quantify this element, Domain, design note 'Designed inside the ...'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| Primer3 included_region sent by the Design button | null (whole record) | [domain start, domain length] when the goal is Detect | no | UI behaviour |
| Primer3 target sent by the Design button | null | null | yes | Detect uses an included region, not a target |

**Measured.** Run in this audit, Standard preset, 5 pairs, bundled examples; first three products (start-end) from the 4.0.0 button (untargeted) vs the 4.1.0 default goal: M11240 4591-4766, 1125-1305, 271-673 vs RT 3188-3923 (interior 276-4870) -> 3303-3564, 3303-3564, 3303-3563; M12927 4266-4450, 4265-4450, 4266-4449 vs RT 3031-3532 -> 3247-3437, 3247-3425, 3247-3440; M80343 1684-2144, 2124-2593, 452-790 vs ORF1 1372-1666 (score 215.6) -> 1471-1627, 1470-1627, 1472-1627; X01005 524-893, 524-1012, 873-1140 vs TPase 525-744 -> 527-677, 527-678, 527-682; X05424 787-985, 340-807, 340-787 vs TPase 3169-3418 (score 109.0) -> 3228-3382, 3231-3382, 3237-3387. Default domains: RT (253.1), RT (134.8), ORF1 (215.6), TPase (82.5), TPase (109.0); all at least 219 bp, so none hits the refusal below.

**Residual limitations.** If the strongest domain is shorter than Prod min the GUI refuses ('The <domain> is N bp, shorter than Minimum product size'); 4.0.0 designed a whole-record pair. Strongest by HMMER score is not the longest or most informative (M80343 defaults to the second of three ORF1 hits, not RT); a domain wholly inside a terminal repeat is refused. Detect sets no Primer3 target, so domain coverage is by construction. Engine-driven benchmarks are unaffected; GUI-driven ones see the new region.

**Benchmark probe.** Positive: load M11240 in the GUI (or call run_primers with included_region [3188, 735]) and press Design with the Standard preset: 4.1.0 pairs lie inside RT 3188-3923 (first product 3303-3564, 261 bp) while 4.0.0 returned 4591-4766. Negative: choose Custom range 1..N (whole record): 4.1.0 reproduces the 4.0.0 untargeted list. Boundary: a record whose strongest domain is under 150 bp: 4.1.0 refuses, 4.0.0 designs.

**Pinned by.** `tests/test_41_ui_M3.py::test_goal_defaults_to_detect_with_the_strongest_domain_and_designs_inside_it`, `tests/test_41_ui_M3.py::test_a_detect_design_returns_pairs_whose_primers_and_product_stay_in_the_domain`, `tests/test_41_ui_M3.py::test_the_domain_is_clipped_to_the_interior_between_the_terminal_repeats`, `tests/test_41_ui_M3.py::test_a_domain_shorter_than_prod_min_is_not_submitted_and_says_why`.

Files: `app/native/panel_primers.py`, `app/native/coords.py`, `app/native/panel_records.py`.

#### PRM-006 — Genotype assay goals: junction and whole-insertion Primer3 targets

**Before (4.0.0).** No genotyping assay existed. A junction design had to be assembled by hand from a flank selection.

**After (4.1.0).** Genotype this insertion offers three designs on the analysed record, all routed (Primer3 target = the stated bases, included region = whole record): 5' junction (target = the 2 bases straddling the element's 5' edge, [start-1, start+1)), 3' junction ([end-1, end+1)) and Whole insertion (target = the element span). The element edges come from the accepted terminal-repeat pair's element_span (coords.insertion_bounds). An option is disabled when its flank is shorter than the minimum primer size (18 nt); with no terminal-repeat pair the goal says the edges are unknown; a record that ends at the element says to fetch the insertion with flank by coordinate.

**Trigger.** Records with an accepted LTR or TIR pair (element_span present); GUI only.

**Affected outputs.** primer pairs for the junction: Primer3 target 2 bp; the pair must span the junction, so one primer sits upstream of the edge and the other downstream; Location text, e.g. '5′ flank → 5′ LTR/internal'; Whole insertion: target = element, product = element + 2 primers.

**Labels.** added: Genotype this insertion, 5' junction - upstream flank into the element, 3' junction - element into the downstream flank, Whole insertion - flank to flank (short insertions).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| junction target length (bp) | n/a | 2 | no | UI design: two bases straddling the edge |
| minimum flank for a junction option (nt) | n/a | primers.MIN_PRIMER_SIZE = 18 | no | a primer needs a flank at least as long as the shortest primer |

**Measured.** Run in this audit on a synthetic genomic window (400-bp random flank + 5-bp TSD + M11240 + 5-bp TSD + 400-bp random flank, seed 11; flanks synthetic): LTR pair element_span 405-5551 (4.0.0: 405-5552, TSD length 4; 4.1.0 TSD length 5), flanks 405 bp. 5' junction (target [404,2)): best pair 280-300 / 675-695, 415 bp, Location '280–695 · 5′ flank → 5′ LTR/internal'. 3' junction (target [5550,2)): best 5341-5361 / 5792-5812, 471 bp, '3′ LTR → 3′ flank'. Whole insertion (target = 5,146-bp element): no pair at Prod max 500.

**Residual limitations.** The CHANGELOG wording 'one primer in each flank' is inaccurate: a junction pair has one primer in the flank and one inside the element (5' LTR or internal), and a flank-to-flank pair exists only for the Whole insertion goal, which needs Prod max of at least element length + 2 primers (about 5.2 kb here). A primer inside a terminal repeat can also prime the other copy and solo copies; the Notes column flags this but does not remove the pair. Junction coordinates follow the structural element_span, which moved by 1 bp between versions for the synthetic window; a benchmark must treat the junction position as dependent on the structural call.

**Benchmark probe.** Positive: build the synthetic flanked M11240 above and design the 5' junction goal: expect a pair with one primer upstream of position 405 and the other inside the element, product spanning 404-406. Negative: M11240 as bundled (no flank): the junction options are disabled ('record ends at the element'). Boundary: flank of 17 nt vs 18 nt on one side: the corresponding junction option is disabled below 18.

**Pinned by.** `tests/test_41_ui_M3.py::test_a_junction_design_targets_the_two_bases_either_side_of_the_edge`, `tests/test_41_ui_M3.py::test_flank_to_flank_over_a_long_insertion_is_refused_with_the_arithmetic`, `tests/test_41_ui_M3.py::test_genotype_on_a_record_that_ends_at_the_element_says_what_to_fetch`, `tests/test_41_ui_M3.py::test_genotype_without_a_terminal_repeat_pair_says_the_edges_are_unknown`.

Files: `app/native/panel_primers.py`, `app/native/coords.py`.

#### PRM-007 — Custom-range goal and pre-flight refusals of domain designs

**Before (4.0.0).** A routed domain design was submitted whatever the product window; Primer3 then returned no pair or a generic message. No sub-range design existed.

**After (4.1.0).** Custom range designs inside a 1-based inclusive range (included region, no target); the whole record is the old untargeted design. _design_for_domain refuses before submitting when (a) an inside/range design is shorter than Prod min, (b) a routed target is longer than Prod max - 2 x min primer size (length + 2 x min_size), or (c) the selected domain lies wholly inside a terminal repeat; each states the arithmetic and what to change. An empty result explains routed/within/masked-base causes instead of blaming parameters.

**Trigger.** GUI path; Custom range, Detect, junction and right-click domain designs.

**Affected outputs.** no request is sent in the refusal cases (4.0.0 sent it and returned no pair); no-pair explanation text.

**Labels.** added: Custom range, refusal notes with the product-size arithmetic.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| routed target feasibility check | none | target length + 2 x min_size <= Prod max | no | Primer3 places one primer on each side of a required target |
| inside-range feasibility check | none | range length >= Prod min | no | no product can fit in a shorter range |

**Measured.** From code and tests (not run as a GUI session in this audit): the refusal thresholds are computed in _design_for_domain from the form values.

**Residual limitations.** The arithmetic uses the minimum primer size, not the primer lengths Primer3 later picks, so a design the check allows can still return no pair.

**Benchmark probe.** Positive: Detect on a domain of 120 bp with Prod min 150: 4.1.0 refuses, 4.0.0 (whole-record button) designs. Negative: Custom range 1..N reproduces the untargeted result. Boundary: a routed 480-bp domain with Prod max 500 and min size 18 needs 516 bp: refused; with Prod max 516 accepted.

**Pinned by.** `tests/test_41_ui_M3.py::test_custom_range_over_the_whole_record_is_the_untargeted_design_and_a_subrange_is_inside_only`, `tests/test_41_ui_M1c.py::test_routed_domain_longer_than_the_product_window_is_not_submitted_and_says_what_is_needed`, `tests/test_41_ui_M1c.py::test_empty_design_names_the_domain_and_never_blames_the_parameters_for_the_template`, `tests/test_41_ui_M1c.py::test_a_failed_design_does_not_blame_the_parameters`.

Files: `app/native/panel_primers.py`.

#### PRM-012 — Oligo QC: a non-negative Primer3 dG is no structure; QC Tm is labelled; displayed thresholds are derived and labelled

**Before (4.0.0).** _p3_kcal returned (dg/1000, tm) whenever Primer3 set structure_found, including a dG of 0 or above; _p3_end_stability did the same. The 3'-end entry hard-coded agree 'single'.

**After (4.1.0).** A Primer3 dG of 0 or above returns (None, None), as ViennaRNA already did, so the metric renders as no structure; the end-stability agree field is 'none' when neither engine reports a structure. Warn and caution flags are unchanged (a non-negative dG never reached a threshold). Also: qc_oligo returns tm_basis 'qc_conditions' beside tm (QC conditions 50 mM Na+, 0 Mg2+, 0 dNTP, 250 nM, 37 C; the comment that called them sealed now says advisory); the primer-table legend, export notes and structure dialog read the dG thresholds from oligoqc, attribute only the -9 warn tier to IDT's rule of thumb and label the caution tiers (-2 hairpin, -5 dimer, -6 3'-end) 'TEagle heuristic'; the Struct column states the worst flag as text; 'dG checks unavailable' appears when no QC succeeded (display-only).

**Trigger.** Any pair or oligo for which Primer3 reports a thermodynamically unfavourable hairpin, dimer or 3'-end structure.

**Affected outputs.** qc.end_stability.p3 / qc.hairpin.p3 / self_dimer.p3 / hetero_dimer.p3: positive number -> null; qc.end_stability.agree: 'single' -> 'none'; table cells: a positive value -> '-'; qc.left/right.tm_basis: absent -> 'qc_conditions'; legend, export-note and Struct-column text.

**Labels.** added: agree value 'none', tm_basis, Struct (text flag), TEagle heuristic.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| structure accepted from Primer3 | structure_found (any dG) | structure_found and dG < 0 | no | an unfavourable structure is no structure |
| oligoqc caution/warn thresholds (kcal/mol) | hairpin -2/-9, dimer -5/-9, 3'-end -6/-9 | hairpin -2/-9, dimer -5/-9, 3'-end -6/-9 | yes | _HAIRPIN_*, _DIMER_*, _END_* constants have no diff |

**Measured.** Run in this audit: oligoqc.qc_pair('TGCTGCGAGATTTGTGGAGT', 'GTCACCACTGTGGCAGAGAA') (the left primers of the first two M11240 candidates, used here as a pair) gives end_stability {p3 0.4, vrna null, flag ok, agree single} in 4.0.0 and {p3 null, vrna null, flag ok, agree none} in 4.1.0; the heterodimer value (-4.39 Primer3, -3.14 ViennaRNA) is identical. Also: qc tm 56.4 in both versions with tm_basis null (4.0.0) vs 'qc_conditions' (4.1.0); oligoqc.CONDITIONS and the _HAIRPIN_/_DIMER_/_END_ thresholds have no diff.

**Residual limitations.** A dG that is only slightly negative still counts as a structure.

**Benchmark probe.** Positive: the pair above: 4.0.0 shows a 3'-end value of +0.4 kcal/mol, 4.1.0 none. Negative: a pair with a negative end dG gives the same value in both versions. Boundary: dG exactly 0.0 is no structure in 4.1.0.

**Pinned by.** `tests/test_41_primers.py::test_positive_primer3_dg_is_no_structure`, `tests/test_41_primers.py::test_real_pair_with_positive_end_stability_reports_none`, `tests/test_41_ui_M1c.py::test_positive_delta_g_is_shown_as_no_structure_not_as_a_value`, `tests/test_41_primers.py::test_qc_tm_is_labelled_with_its_conditions`.

Files: `app/backend/teagle_core/oligoqc.py`, `app/native/ui_common.py`, `app/native/panel_primers.py`.

### In-silico PCR and whole-genome scan

#### GEN-001 — Whole-genome scan: the on-target product is the one at the designed product's genome position

**Before (4.0.0).** The engine took design_locus entries {accession, start, stop} as a window and marked as the single on-target the pair product with the largest overlap with any window (ties to the lowest start); any product overlapping the window could win, including a long product spanning the designed one. The shipped panel built the locus from state['source'].regions, whose normalised entries carry 'end', not 'stop', and read 'stop', so the engine received stop = None (coerced to 0) and a window [0, start]: the true product at positions above start overlapped nothing and was called off-target (tier intended-locus-negative). The panel also sent every fetched region.

**After (4.1.0).** A locus entry is {accession, start, end, strand (1/2 or +/-), intended_product_span [a, b), intended_endpoint_tolerance_bp}: the fetched window (1-based inclusive) plus the designed product's 0-based half-open span inside the fetched record. _design_loci maps it to genome coordinates, strand-aware: forward lo = start + a, hi = start + b - 1; reverse lo = end - b + 1, hi = end - a. Exactly one product is marked on_target: the non-single-primer product on the same sequence whose two endpoints are within the tolerance (default 0 bp), smallest error then lowest start; a designed product that is absent promotes nothing (tier intended-locus-negative). An entry lacking accession, window ends, strand or span names no intended product and is ignored (has_locus False, neutral priming sites); a bad strand or an inconsistent span is a BadRequest. The panel sends the locus only for a pair whose own sealed design input is a one-region coordinate fetch of the assembly being scanned and whose record is still on screen (the region the pair was designed on in a multi-region fetch). Labelling is derived and outside the seal.

**Trigger.** Whole-genome scan of a pair designed on a specimen fetched by coordinate from the scanned assembly; any engine caller that sends design_locus.

**Affected outputs.** amplicons[].on_target; n_on_target; has_locus; summary.tier (intended-only / limited-off-target / many-off-target / intended-locus-negative / none); summary.n_on / n_off / verdict; gel and table on/off-target calls.

**Labels.** added: design_locus.end, .strand, .intended_product_span, .intended_endpoint_tolerance_bp; removed: design_locus.stop as the window end (a legacy 'stop' is still read for 'end' but a locus without span/strand is now ignored).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| on-target rule | largest overlap with the fetched window | both endpoints within tolerance of the mapped designed product | no | the designed product is the intended target |
| intended_endpoint_tolerance_bp | n/a for the scan | 0 (primers.INTENDED_ENDPOINT_TOLERANCE_BP, carried from the design) | no | design geometry tolerance |
| isPcr min_perfect / min_good | 15 | 15 | yes | isPcr defaults, unchanged |

**Measured.** Live (this audit; WSL backend with yeast R64, GCF_000146045.2, prepared; NCBI coordinate fetch): window NC_001136.10:100001-102000 fetched on the + and the - strand, first Standard-preset pair (150-500 bp); real isPcr returns exactly one product, NC_001136.10 100394-100812 (419 bp, 1-based inclusive), on both strands. 4.1.0: on_target True, tier intended-only ('1 on-target product; no additional priming site detected in this assembly under current settings'). 4.0.0 as the shipped panel sent the locus (emulated by posting {accession, start, stop: None}; not a GUI run): on_target False, tier intended-locus-negative ('1 off-target candidate priming site(s); no product at the intended locus under current settings'). Mocked isPcr output, engine called directly, window NC_1:1000-3000, designed span [500,700): products 1500-1699 (designed), 1101-1900 (spanning), NC_9:1500-1699. 4.0.0 with a well-formed {start 1000, stop 3000}: the spanning 1101-1900 product is on-target; 4.1.0: 1500-1699. Reverse-strand fetch (span maps to 2301-2500): 4.0.0 marks 1500-1699, 4.1.0 marks 2301-2500. Designed product absent but spanning product present: 4.0.0 intended-only, 4.1.0 intended-locus-negative. Off-by-one product 1501-1700: 4.0.0 on-target; 4.1.0 off-target at tolerance 0, on-target at tolerance 1. Window-only locus (no span/strand): 4.0.0 marks the spanning product, 4.1.0 neutral (has_locus False, tier limited-site). Two equal copies, designed second: 4.0.0 marks the left copy, 4.1.0 the second. The sealed manifest is identical with and without a locus in both versions.

**Residual limitations.** On-target needs exact endpoint agreement (tolerance 0 by default); a product shifted one base is off-target. If fetched sequence names (RefSeq NC_) differ from the scanned assembly's names (e.g. GenBank GCA names) nothing can be on-target; the panel then says the split is not assessed. Only single-region coordinate-fetched designs carry a locus; single-primer products are never on-target.

**Benchmark probe.** Positive: fetch S. cerevisiae chrIV:100,001-102,000 (assembly GCF_000146045.2, + or - strand), design with the Standard preset, take the first pair, scan against the prepared genome with the locus built from the fetch: 4.1.0 labels the single 419-bp product NC_001136.10:100394-100812 on-target; a 4.0.0-style request labels it off-target. Negative: same pair, specimen pasted without a fetch (no locus): both versions neutral (limited/single-site tiers, n_on_target 0). Boundary: shift the designed span by one base with tolerance 0 (off-target) and 1 (on-target); and a designed span that matches a second equal-length copy.

**Pinned by.** `tests/test_41_engine.py::test_genome_on_target_is_the_designed_product_not_the_longest`, `tests/test_41_engine.py::test_genome_on_target_beats_a_spanning_product`, `tests/test_41_engine.py::test_genome_tie_resolves_to_the_designed_copy_not_the_leftmost`, `tests/test_41_engine.py::test_genome_span_maps_through_a_minus_strand_fetch`.

Files: `app/backend/engine.py`, `app/native/panel_genome.py`, `app/native/main.py`, `app/backend/teagle_core/genomepcr.py`.

#### GEN-002 — Genome-scan verdicts state single-primer products and name an intended-locus miss

**Before (4.0.0).** With no forward+reverse product the verdict said there was no genome-wide product and that the primers do not amplify anywhere in the assembly, even when F+F or R+R products existed; with one on-target product it said no additional priming site was detected regardless of single-primer products.

**After (4.1.0).** n_single is counted and, when above zero, '; N single-primer (F+F/R+R) product(s) also predicted' is appended to every verdict; the sentences 'no additional priming site', 'do not amplify anywhere' and 'no genome-wide priming site' are not emitted while a single-primer product exists. The has_locus zero-product verdict becomes 'no forward+reverse product at the intended locus or elsewhere in this assembly under current settings'. Tier names and thresholds (none, single-site, limited-site (2-5), many-site (>5), intended-only, limited-off-target, many-off-target, intended-locus-negative) are unchanged and still pair-based.

**Trigger.** Scans whose product set contains single-primer rows or no pair product.

**Affected outputs.** summary.verdict text; summary.tier is unchanged.

**Labels.** added: '; N single-primer (F+F/R+R) product(s) also predicted'; renamed: 'no genome-wide product - the primers do not amplify anywhere in this assembly' -> 'no forward+reverse product at the intended locus or elsewhere in this assembly under current settings', '1 on-target product; no additional priming site detected ...' -> '... no additional forward+reverse product detected ... ; N single-primer ... predicted' (when single-primer products exist), 'no genome-wide priming site (...)' -> 'no forward+reverse priming-site pair in this assembly, at the current 3' match length; N single-primer ...' (when single-primer products exist).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| tier thresholds (limited vs many sites/off-targets) | <=5 / >5 | <=5 / >5 | yes | unchanged |

**Measured.** Run in this audit (genomepcr.summarize on synthetic amplicon rows): single-primer product only, has_locus True: 4.0.0 'no genome-wide product - the primers do not amplify anywhere in this assembly', 4.1.0 'no forward+reverse product at the intended locus or elsewhere in this assembly under current settings; 1 single-primer (F+F/R+R) product(s) also predicted' (tier none in both). One on-target plus one single-primer product: 4.0.0 '1 on-target product; no additional priming site detected in this assembly under current settings', 4.1.0 '... no additional forward+reverse product detected ...; 1 single-primer (F+F/R+R) product(s) also predicted' (tier intended-only in both). No products, no locus: identical in both.

**Residual limitations.** The tier still ignores single-primer products, so intended-only can coexist with predicted single-primer products.

**Benchmark probe.** Positive: a scan with one pair product plus one fwdonly row: verdict text differs, tier equal. Negative: no single-primer rows and a locus: intended-only verdict text unchanged. Boundary: exactly 5 vs 6 off-target products keep limited-off-target / many-off-target.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_verdict_with_locus_and_single_primer_product_does_not_claim_no_additional_site`.

Files: `app/backend/teagle_core/genomepcr.py`.

#### GEN-003 — Genome-scan product-size window widens for a designed product longer than 4,000 bp (sealed)

**Before (4.0.0).** The scan always sent prod_max 4000 (and prod_min 0), so a pair designed for a product above 4 kb could not report its own product.

**After (4.1.0).** prod_max stays 4000 unless the designed product_size exceeds it, then prod_max = int(1.25 x product_size); the window in force is stated on the panel, in the busy text and in exports, and is sealed as max_size.

**Trigger.** Pairs with product_size > 4000 (Permissive preset allows up to 1000 bp by default; only manually raised ranges reach this).

**Affected outputs.** isPcr -maxSize / sealed max_size; scan manifest seal for those designs; products between 4000 and 1.25 x size now reported.

**Labels.** added: product-size window statement.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| scan prod_max (bp) | 4000 | 4000, or int(1.25 x designed product) when the designed product is longer | no | UI margin with no scientific basis; sealed |
| scan prod_min | 0 | 0 | yes | unchanged |

**Measured.** From code and tests.

**Residual limitations.** The 1.25 factor is a margin, not a derived value; products longer than the window are still not reported.

**Benchmark probe.** Positive: a pair with a 5,000-bp designed product: 4.1.0 sends max_size 6250 and reports the product; 4.0.0 sends 4000 and does not. Negative: a 400-bp pair: max_size 4000 in both. Boundary: product 4000 (no widening) vs 4001 (max_size 5001).

**Pinned by.** `tests/test_41_ui_M1c.py::test_scan_window_covers_a_designed_product_longer_than_4000_bp_and_says_so`, `tests/test_41_ui_M1c.py::test_scan_result_discloses_the_size_cap_and_identifies_pair_organism_and_criteria`.

Files: `app/native/panel_genome.py`, `app/backend/engine.py`.

#### PCR-001 — A product must be longer than the longer annealed primer core (crossed 3' ends and primer-length self-primed products are no longer reported)

**Before (4.0.0).** A pair of binding sites formed a product whenever left < right and prod_min <= templated span + tails <= prod_max. With prod_min 1, crossed 3' ends gave 15-bp 'products' shorter than a primer, a tailed crossed pair was rescued by its tails, and a palindromic 20-mer gave a 20-bp F+F product.

**After (4.1.0).** The pair must satisfy right - left > max(forward anneal_len, reverse anneal_len) (strict) in addition to the unchanged size window, so the two annealed cores cannot cross or coincide and no product is shorter than a primer core. Valid overlapping non-crossed pairs (for example a 30-bp span with 20-nt cores) are still reported.

**Trigger.** Short product windows (prod_min below about 2 x primer length), primers that are palindromic or self-complementary at their 3' ends, and pairs that bind nearly the same position.

**Affected outputs.** result.amplicons: fewer products (those with right - left <= longer anneal_len are removed); gel lanes and verdicts that counted such a product.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| product acceptance rule | left < right and prod_min <= length <= prod_max | right - left > max(anneal_len F, anneal_len R) and prod_min <= length <= prod_max | no | a product shorter than a primer core cannot be a PCR product |
| prod_min / prod_max defaults (run_pcr) | 70 / 1000 | 70 / 1000 | yes | unchanged |

**Measured.** Run in this audit with primers.in_silico_pcr on M11240 (prod_min 1, prod_max 1000, max_mm 2, tp 5, min_anneal 18). Crossed ends (F = bases 100-120, R = reverse complement of bases 95-115): 4.0.0 reports two 15-bp F+R products (100-115 and a second site at 4970-4985 in the 3' LTR copy); 4.1.0 reports none. Palindromic 20-mer (ACGT)5 inserted at 500 as a primer: 4.0.0 F+F products of 20, 24, 24 and 28 bp; 4.1.0 24, 24 and 28 bp (the 20-bp product is dropped, the 24-bp ones carry 2 permitted mismatches). A normal 320-bp pair at 1000-1320 is identical in both.

**Residual limitations.** The 4.0.0 default prod_min is 70, so the change matters only for short windows or for self-primed palindromic primers; most benchmark runs at defaults see no difference.

**Benchmark probe.** Positive: F = M11240[100:120], R = revcomp(M11240[95:115]) at prod_min 1: 4.0.0 gives 15-bp products, 4.1.0 none. Negative: F = M11240[1000:1020], R = revcomp(M11240[1300:1320]): one 320-bp product in both. Boundary: palindromic (ACGT)5 as both primers: product lengths 20 (dropped in 4.1.0), 21 (kept) when anneal_len is 20.

**Pinned by.** `tests/test_41_primers.py::test_crossed_three_prime_ends_make_no_product`, `tests/test_41_primers.py::test_tailed_crossed_pair_is_not_rescued_by_its_tails`, `tests/test_41_primers.py::test_primer_length_self_primed_product_is_not_reported`, `tests/test_41_primers.py::test_overlapping_footprints_keep_product_length_and_start_with_the_forward_oligo`.

Files: `app/backend/teagle_core/primers.py`.

#### PCR-002 — A template ambiguity code always counts as a mismatch, including against an identical primer code

**Before (4.0.0).** _base_ok returned True when the primer base equalled the template base, so primer N matched template N and primer R matched template R at zero cost; a degenerate primer base otherwise matched any template base it covers.

**After (4.1.0).** A primer base matches only when its IUPAC set contains the template base and that base is A, C, G or T; an ambiguous or unknown template base (N, R, Y, ...) is a mismatch whatever the primer base. The degenerate-primer table and the rest of the matcher are unchanged.

**Trigger.** Templates with N or IUPAC codes at a primer footprint (assembly gaps, consensus sequences) searched with primers that carry the same code.

**Affected outputs.** fwd_mm / rev_mm of products across an ambiguous template base (+1 per such base); products that exceed max_mm after the extra mismatch disappear.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| max_mm default | 2 | 2 | yes | unchanged |
| tp (strict 3' window) | 5 | 5 | yes | unchanged |

**Measured.** Run in this audit on M11240 with a template N (and R) planted at 1010 inside the forward primer site and the primer carrying the same code: max_mm 0 -> 4.0.0 one product with fwd_mm 0, 4.1.0 none (also for R vs R); max_mm 2 -> 4.0.0 fwd_mm 0, 4.1.0 fwd_mm 1 for the same product.

**Residual limitations.** engine.run_pcr still lists 'IUPAC-ambiguous template bases' under notSearched, and an ambiguity inside the 3' strict window (tp) now removes the site at any max_mm.

**Benchmark probe.** Positive: template with N at a mid-primer position, primer with N there, max_mm 0: 4.0.0 finds the product, 4.1.0 does not. Negative: same template, primer with the real base: identical. Boundary: the N sits in the last 5 bases (strict 3' window): no site in either version for a perfect-base primer; 4.0.0 alone accepts primer N.

**Pinned by.** `tests/test_41_primers.py::test_primer_n_does_not_match_template_n`.

Files: `app/backend/teagle_core/primers.py`.

#### PCR-004 — Primer input is normalised (whitespace, case, U to T) and anything outside the IUPAC codes is rejected

**Before (4.0.0).** in_silico_pcr only upper-cased the primers; a space inside a primer silently gave zero products (status complete_zero_products), '-', '1' or letters outside IUPAC did not raise, and RNA 'U' primers were sealed unnormalised. run_pcr gated primers with the generic _require_nt check.

**After (4.1.0).** primers.normalize_primer strips all whitespace, upper-cases, reads U as T and raises ValueError naming the offending characters ('forward primer contains characters that are not IUPAC nucleotide codes: ...'); run_pcr applies it before matching, the design-provenance comparison and sealing (BadRequest on error); in_silico_pcr applies the same step. The genome-scan route does not use it (see unchanged).

**Trigger.** Pasted or typed primers containing spaces, lower case, U, digits or gap characters.

**Affected outputs.** products for a primer with an internal space: none -> normal product; error instead of silent zero products for illegal characters; sealed primer sequences are the normalised DNA text.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| accepted primer alphabet | run_pcr: nucleotide-ish by _require_nt; in_silico_pcr: any upper-cased text | ACGTRYSWKMBDHVN after stripping whitespace and U->T | no | IUPAC nucleotide codes |

**Measured.** Run in this audit on M11240 (primers.in_silico_pcr): a space inside the forward primer returns [] in 4.0.0 and the normal 320-bp product in 4.1.0; a '-' returns [] in 4.0.0 and raises 'forward primer contains characters that are not IUPAC nucleotide codes: '-'' in 4.1.0; lower-case and U primers give the same product in both.

**Residual limitations.** Primer minimum length (18) and the 5' tail rule are enforced elsewhere (engine and panel).

**Benchmark probe.** Positive: forward primer with an embedded space: 4.0.0 zero products, 4.1.0 one. Negative: clean upper-case primers identical. Boundary: a primer with 'U' is read as T in both matchers, but only 4.1.0 seals the T form.

**Pinned by.** `tests/test_41_primers.py::test_primer_whitespace_case_and_u_are_normalised`, `tests/test_41_primers.py::test_non_iupac_primer_is_rejected_with_the_offending_characters`, `tests/test_41_primers.py::test_in_silico_pcr_normalises_whitespace_too`.

Files: `app/backend/teagle_core/primers.py`, `app/backend/engine.py`.

#### PCR-005 — Annealing cores are taken from the specimen and applied to background sequences

**Before (4.0.0).** Every sequence (specimen and each custom background) was searched alone: per primer and strand, only the sites that reach the longest 3'-proximal core found in that sequence were kept, and the core was reported as anneal_len with the 5' remainder as tail. Results carried no record of the cores used.

**After (4.1.0).** When backgrounds are supplied, run_pcr derives each primer's per-strand core from the specimen alone (primers.specimen_cores) and passes it to the backgrounds: in a background, best = min(its own longest reach, the specimen core) and every site reaching at least that length is kept and called at it. The specimen (k = 0) is always searched as a lone search, so its result no longer depends on the backgrounds. A strand where the specimen has no site gives core None and the background keeps its own. criteria.annealing_core reports the cores used (the largest per primer and strand over the searched sequences).

**Trigger.** In-silico PCR with a custom background where the specimen's best site for a primer is shorter than a site in the background (an imperfect specimen copy against a perfect background copy).

**Affected outputs.** background amplicons: perfect sites re-called with shorter anneal_len and a tail (anneal 26 -> 20, tail 0 -> 6, start shifted by the tail length, product_length unchanged); background amplicons: additional products from sites between the specimen core and the background's own best reach; criteria.annealing_core {F/R: {forward_strand, reverse_strand}}: absent -> object; same sequence can give different product sets as specimen and as background.

**Labels.** added: criteria.annealing_core.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| core selection per sequence | longest core reached in that sequence | specimen: unchanged; background: min(own longest, specimen core) | no | makes an identical site get the same call in specimen and background |
| min_anneal | 18 | 18 | yes | primers.MIN_PRIMER_SIZE, unchanged |

**Measured.** Run in this audit through engine.run_pcr; 26-mer forward and 24-mer reverse primer, max_mm 2, tp 5, synthetic random sequences (seed 7). Specimen holds only a forward-site copy with 8 mismatches in its 5' end (core 20, tail 6); background holds the perfect site plus the mismatched copy, each followed by the reverse site. 4.0.0: specimen product (anneal 20, tail 6, 2 mm) plus two background products at the perfect site (anneal 26, tail 0; 200 and 470 bp). 4.1.0: same specimen product; the two background products are re-called at anneal 20 / tail 6 with the start shifted 200 -> 206 (lengths unchanged) AND a third background product from the mismatched copy (start 506, 170 bp, anneal 20, tail 6, 2 mm) that the background's own lone search drops. A perfect specimen with the same background gives the same three products at core 26 in both versions; a background with only the mismatched copy gives the same product in both; searched alone, the background still hides the mismatched copy when a perfect site exists (2 products in both).

**Residual limitations.** A perfect background site is reported with the specimen's shorter core and a tail (26-nt perfect site reads 20 + 6); tail bases are not scored. The same sequence now gives different product sets as specimen (lone rule hides the shorter imperfect copy) and as background (both copies reported): a new asymmetry. CHANGELOG says 'without dropping any product'; it does not say 4.1.0 can add products and re-call perfect sites.

**Benchmark probe.** Positive: specimen with only an imperfect primer site, background with perfect + imperfect sites: 4.0.0 reports the perfect site only; 4.1.0 reports both, re-called at the specimen core. Negative: single-sequence run (no background): identical. Boundary: specimen with a perfect site: background results identical in both versions (core 26).

**Pinned by.** `tests/test_41_primers.py::test_core_decision_is_the_same_for_a_site_in_specimen_and_background`, `tests/test_41_primers.py::test_specimen_result_does_not_depend_on_the_backgrounds`, `tests/test_41_primers.py::test_background_with_a_longer_site_is_still_reported_at_the_specimen_core`, `tests/test_41_primers.py::test_background_without_a_specimen_site_keeps_its_own_core`.

Files: `app/backend/teagle_core/primers.py`, `app/backend/engine.py`.

#### PCR-006 — Search-cap flags are set only when something was skipped

**Before (4.0.0).** amplicons_capped and sites_capped were set when the count reached MAX_AMPS / MAX_SITES (4000 each), so a search that found exactly 4000 products or sites was reported as truncated.

**After (4.1.0).** The flags are set only when a further site or product exists beyond the limit (checked before appending).

**Trigger.** Searches that find exactly MAX_AMPS (4000) products or MAX_SITES (4000) sites; very rare in practice.

**Affected outputs.** stats.amplicons_capped / stats.sites_capped; engine truncated flag and the 'stopped early' banner.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| MAX_SITES / MAX_AMPS | 4000 / 4000 | 4000 / 4000 | yes | unchanged limits |

**Measured.** Run in this audit with MAX_AMPS patched to N on a synthetic template built to hold exactly 3 (then 4) F+R products: with the limit equal to the product count 4.0.0 reports amplicons_capped True and 4.1.0 False; with one more product than the limit both report True and return the limit's worth of products.

**Residual limitations.** At the true limit of 4000 the tool still stops after 4000; only the flag is corrected.

**Benchmark probe.** Positive: template with exactly MAX_AMPS products (patch the constant to 3 for a fast test): 4.0.0 flags capped, 4.1.0 not. Negative: below the limit both false. Boundary: limit + 1 products: both true.

**Pinned by.** `tests/test_41_primers.py::test_site_cap_is_set_only_when_a_further_site_exists`, `tests/test_41_primers.py::test_amplicon_cap_is_set_only_when_a_further_product_exists`.

Files: `app/backend/teagle_core/primers.py`.

#### PCR-008 — In-silico PCR accepts a primer pair typed by the user

**Before (4.0.0).** Only pairs designed in panel 04 (or sent from the family panel) could be queued; a pair you already own could not be checked.

**After (4.1.0).** Panel 05 adds forward, reverse and an optional name; the pair is normalised (PCR-004), must have at least 18 bases including any 5' tail, and is queued as an unsealed custom pair (no design_provenance, so design_assessment is not_run and the intended-product check does not apply). The queue row reads 'your own pair (name)'.

**Trigger.** GUI use of the new fields.

**Affected outputs.** PCR queue entry; design_assessment 'not_run' for those pairs.

**Labels.** added: Or check a primer pair you already have, your own pair (<name>).

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| minimum primer length incl. tail (nt) | n/a | 18 | yes | primers.MIN_PRIMER_SIZE |

**Measured.** From code; the engine path is the same run_pcr call (PCR-004).

**Residual limitations.** A typed pair has no designed target, so no on/off-target-style design verdict is produced.

**Benchmark probe.** Positive: type a 6-nt EcoRI tail + 20-mer pair for M11240: the product is found with tail 6 and design_assessment not_run. Negative: a sealed designed pair on its own record gets design_assessment assess. Boundary: a 17-base primer is refused.

**Pinned by.** `tests/test_41_ui_m1a.py::test_a_pair_you_already_have_can_be_queued_and_run_with_a_tail`, `tests/test_41_ui_m1a.py::test_a_pair_with_a_bad_letter_or_too_short_is_refused_with_a_reason`.

Files: `app/native/panel_pcr.py`.

#### PCR-009 — The in-silico PCR product-size window is widened to hold every staged designed product

**Before (4.0.0).** The window prod_min/prod_max (defaults 70/1000) filtered products, so a designed product outside it vanished and the pair read as 'zero products' although it amplifies its own template.

**After (4.1.0).** _run_pcr widens the window to min(prod_min, smallest designed product) and max(prod_max, largest designed product), updates the two fields, shows a warning that names the pair and the old and new window, and repeats it in the result note and the table export. The widened values are sent to the engine and sealed as parameters. A pair inside the window changes nothing.

**Trigger.** GUI runs where a staged designed pair has product_size outside the window (for example a designed 1,100-bp product with Prod max 1000).

**Affected outputs.** engine params prod_min/prod_max (sealed); result note and export note; the previously hidden designed product now appears.

**Labels.** added: window-widening warning.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| PCR product window used for the run | as typed | widened to hold the designed products | no | UI decision; shown and sealed |

**Measured.** From code and tests.

**Residual limitations.** The widening is per run and also widens the window for background products.

**Benchmark probe.** Positive: design a pair with a 1,100-bp product (Permissive preset) and run PCR with Prod max 1000: 4.0.0 reports zero products, 4.1.0 widens to 1100 and reports the designed product. Negative: a 400-bp product: no change. Boundary: product exactly 1000: no widening.

**Pinned by.** `tests/test_41_ui_M1c.py::test_a_designed_product_outside_the_pcr_window_widens_it_visibly_instead_of_vanishing`, `tests/test_41_ui_M1c.py::test_a_pair_inside_the_window_leaves_it_alone`, `tests/test_41_ui_M1c.py::test_the_pcr_note_reports_the_widened_window`.

Files: `app/native/panel_pcr.py`.

### Family annotation and splice detection

#### GEN-006 — Family annotation with no organism is labelled and sealed as RepeatMasker's default lineage (Homo sapiens)

**Before (4.0.0).** With no organism chosen the result and the sealed parameters.species read '(all installed families)', although RepeatMasker (no -species) searches only its default lineage; resolve_species advised leaving the organism blank 'to search all installed families'.

**After (4.1.0).** The annotate script prints the line after 'Species/Taxa Search:' from RepeatMasker's log and the result reports species 'Homo sapiens (RepeatMasker default lineage)' (the log's lineage, or the documented default if the log has none) plus species_requested, species_defaulted, lineage_searched, lineage_taxid; the genome path does the same via SPECIES_SEARCHED. A named organism keeps its name and gains lineage_searched / lineage_taxid. The advice for an unknown name now says 'leave the organism blank to use RepeatMasker's default lineage (Homo sapiens)'. The RepeatMasker command and its hits are unchanged.

**Trigger.** Every family annotation and every whole-genome annotation run with the organism blank.

**Affected outputs.** result.species and provenance.parameters.species (sealed value); new result fields species_requested, species_defaulted, lineage_searched, lineage_taxid; every annotate and genome_annotate manifest hash for blank-organism runs.

**Labels.** added: species_requested, species_defaulted, lineage_searched, lineage_taxid; renamed: species '(all installed families)' -> 'Homo sapiens (RepeatMasker default lineage)'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| RepeatMasker -species for a blank organism | not passed (default lineage) | not passed (default lineage) | yes | search unchanged; only the label and seal change |
| sealed species value for a blank organism | (all installed families) | Homo sapiens (RepeatMasker default lineage) | no | label corrected to what is searched |

**Measured.** Live (this audit; WSL backend: RepeatMasker 4.2.4, Dfam 4.0 curated partition 0 only, no uncurated): wsl.annotate on bundled M11240 with no organism returns the same 2 hits in both versions (both 'A-rich' Low_complexity; no Copia family because the Homo sapiens default lineage is searched). Species field: 4.0.0 '(all installed families)'; 4.1.0 'Homo sapiens (RepeatMasker default lineage)', species_defaulted True, lineage_searched 'Homo sapiens', lineage_taxid 9606. With 'Drosophila melanogaster' both versions return the same 2 hits (Copia_LTR, Copia_I; LTR/Copia); 4.1.0 adds lineage_searched and lineage_taxid 7227.

**Residual limitations.** A blank-organism run still searches only the human lineage, so a non-mammalian element receives no family name; the label now says so but the result is the same empty or low-complexity list. benchmark seals of annotate runs from 4.0.0 do not match 4.1.0 for blank-organism runs.

**Benchmark probe.** Positive: wsl.annotate(M11240) with species None: identical hits in both versions, sealed species string differs. Negative: species 'Drosophila melanogaster': hits and label identical (plus new lineage fields). Boundary: a name that is ambiguous ('drosophila', 141 taxa) fails before the search in both versions.

**Pinned by.** `tests/test_41_wsl.py::test_no_organism_is_labelled_with_the_lineage_repeatmasker_defaults_to`, `tests/test_41_wsl.py::test_no_organism_falls_back_to_the_documented_default_when_the_log_has_no_taxon_line`, `tests/test_41_wsl.py::test_named_organism_keeps_its_name_and_the_resolved_taxon`, `tests/test_41_wsl.py::test_unknown_species_advice_no_longer_claims_a_search_of_all_families`.

Files: `app/backend/teagle_core/wsl.py`, `app/backend/engine.py`.

#### GEN-007 — Bundled examples set their own organism for family naming; an organism taken from a specimen is cleared when the specimen changes

**Before (4.0.0).** _set_species(org) ran after an accession or coordinate fetch but not after loading a bundled example, so Run family annotation on an example with the organism left blank searched the default (human) lineage; an organism carried over from an earlier fetch or example stayed selected for the next specimen.

**After (4.1.0).** _load_example selects the example's organism from examples.CATALOG (M11240 and M12927 Drosophila melanogaster, M80343 Homo sapiens, X01005 Caenorhabditis elegans, X05424 Zea mays; a name not in the dropdown goes to 'Other...' free text). A species set this way is flagged automatic and is reset to '- select organism -' on any specimen change or edit; a choice the user made stays.

**Trigger.** GUI: load a bundled example (or fetch) and press Run family annotation without touching the organism.

**Affected outputs.** the -species value sent to RepeatMasker, hence the family hits and the sealed species.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| organism after loading an example | unchanged (blank unless chosen) | the example's own organism | no | catalogue organism of the published element |

**Measured.** Effect of the organism, measured live at the wsl.annotate level with the argument the GUI sends in each case: M11240 with species None -> 2 hits, both Low_complexity 'A-rich' (4.0.0 GUI path after loading the example); M11240 with 'Drosophila melanogaster' -> Copia_LTR and Copia_I (4.1.0 GUI path). The GUI path itself was not driven in this audit.

**Residual limitations.** The example's organism is fixed by the catalogue, so a user who wants the human default must clear it; the Zea mays and C. elegans names were not run live here.

**Benchmark probe.** Positive: GUI, load the M11240 example, press Run family annotation with no other action: 4.0.0 returns A-rich low complexity only; 4.1.0 returns Copia_LTR and Copia_I. Negative: pick an organism manually: identical. Boundary: load an example, then fetch another accession: the example's organism is replaced by the fetch's organism.

**Pinned by.** `tests/test_41_ui_m1a.py::test_example_sets_its_own_organism_and_a_paste_clears_it`, `tests/test_41_ui_m1a.py::test_a_chosen_organism_survives_a_specimen_change`, `tests/test_41_ui_m1a.py::test_fetch_after_a_human_example_does_not_keep_human_for_a_record_without_organism`.

Files: `app/native/panel_family.py`, `app/native/panel_specimen.py`.

#### GEN-008 — Dfam family-count lookup fails closed for ambiguous names and dead sessions

**Before (4.0.0).** A famdb failure, a dead WSL session or an ambiguous organism name (famdb exits 0 and lists no families, for example 'drosophila' matching 141 taxa) returned ok True with 0 families.

**After (4.1.0).** The script captures famdb's exit status (TEAGLE_FAMDB_FAILED, exit 4) and the function fails on a non-zero rc, a missing count or a famdb failure ('the Dfam library could not list the families for <name> ...; no count is available'); a zero count is reported only after resolve_species confirms the name resolves to one taxon (otherwise the ambiguity or not-found error is returned with ambiguous_species).

**Trigger.** Annotation cost dialog and library-choice probes with an ambiguous, unknown or unreachable organism.

**Affected outputs.** run_annotate_budget.library_families_for_species (None instead of 0); error text.

**Labels.** added: ambiguous_species.

**Measured.** Live (this audit): dfam_lineage_families('drosophila') -> 4.0.0 ok True, families 0; 4.1.0 ok False with the ambiguity error and ambiguous_species True. 'Notaspecies xyz' -> 4.0.0 ok True, families 0; 4.1.0 ok False, 'was not found ... leave the organism blank to use RepeatMasker's default lineage (Homo sapiens)'. Valid names unchanged: Drosophila melanogaster 399 and Homo sapiens 1439 curated family models in both versions.

**Residual limitations.** Counts are for the installed partitions (curated partition 0 here); an uncurated install changes the numbers.

**Benchmark probe.** Positive: families for 'drosophila': 0 in 4.0.0, error in 4.1.0. Negative: 'Drosophila melanogaster': 399 in both. Boundary: 'human' (one taxon by common name) resolves.

**Pinned by.** `tests/test_41_wsl.py::test_lineage_family_count_fails_closed`, `tests/test_41_wsl.py::test_lineage_count_script_propagates_famdb_failure`, `tests/test_41_wsl.py::test_zero_family_count_for_an_ambiguous_name_is_not_reported_as_a_count`.

Files: `app/backend/teagle_core/wsl.py`.

#### GEN-010 — Splice detection takes the gene strand from minimap2's ts tag; minus-strand exon FASTA is numbered in transcript order

**Before (4.0.0).** The strand was the SAM FLAG-16 orientation of the read, so an antisense read (3' ESTs, first-strand cDNA, a reverse-complemented transcript) was reported as the opposite gene strand and its splice motifs were read on the wrong strand: CT..AC instead of GT..AG, 0 canonical introns.

**After (4.1.0).** The gene strand combines the alignment orientation with the ts:A tag (minimap2 -x splice reports the transcript strand relative to the read): strand = alignment strand, flipped when ts is '-'. The result also returns alignment_strand, transcript_strand and strand_determined (False when no ts tag is present; the strand then falls back to the alignment strand). The splice panel states 'strand undetermined' in that case and shows alignment orientation and transcript strand separately. Exon and intron blocks and the canonical set (GT-AG, GC-AG, AT-AC) are unchanged. Exports: splice exon FASTA is sorted in transcript order (exon_1 is the mRNA's 5' exon for a minus-strand gene) with transcript_rank= and coordinate_system=0-based-half-open in the header; ORF, domain and exon FASTA headers state the convention and write a missing Pfam as NA.

**Trigger.** Transcripts that align to the reverse strand or are antisense; every splice result gains the new fields.

**Affected outputs.** result.strand; introns[].canonical and canonical_introns (counted on the corrected strand); alignment_strand, transcript_strand, strand_determined (new); splice GFF/FASTA exports that carry strand; exon numbering and header tokens in the splice exon FASTA.

**Labels.** added: alignment_strand, transcript_strand, strand_determined, transcript_rank=, coordinate_system=0-based-half-open.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| minimap2 invocation | -x splice --secondary=no | -x splice --secondary=no | yes | unchanged; the ts tag was already emitted |
| canonical splice-site set | GT-AG, GC-AG, AT-AC | GT-AG, GC-AG, AT-AC | yes | _SPLICE_CANON unchanged |

**Measured.** Live (this audit, real minimap2 in the WSL backend; synthetic 3-exon gene, GT..AG introns, 150-bp random flanks, seed 5, all synthetic). Gene on + with a sense transcript: strand + in both versions, 2 of 2 canonical. Gene on + with the antisense (reverse-complemented) transcript: 4.0.0 strand '-' and 0 canonical; 4.1.0 strand '+', alignment_strand '-', transcript_strand '-', 2 canonical. Gene on - (genomic window reverse-complemented) with a sense transcript: strand '-', 2 canonical in both versions; with the antisense transcript: 4.0.0 strand '+', 0 canonical; 4.1.0 strand '-', 2 canonical. Pure parser probes on a SAM line with ts:A:- (FLAG 0): 4.0.0 '+', 4.1.0 '-'; without a ts tag strand_determined False.

**Residual limitations.** A transcript that is unspliced or has no ts tag cannot fix the strand; the alignment orientation is shown as such. For a minus-strand gene aligned with a sense read the alignment strand and gene strand coincide, so a read-orientation-only reading was right only in the common case.

**Benchmark probe.** Positive: align a reverse-complemented cDNA of a + strand gene: 4.0.0 reports the gene on '-' with 0 canonical introns, 4.1.0 reports '+', canonical introns equal to the intron count. Negative: sense transcript on a + gene: identical. Boundary: a single-exon transcript (no ts tag): strand_determined False in 4.1.0.

**Pinned by.** `tests/test_41_wsl.py::test_gene_strand_combines_alignment_orientation_with_the_ts_tag`, `tests/test_41_wsl.py::test_spliced_alignment_without_a_ts_tag_does_not_claim_a_determined_strand`, `tests/test_41_wsl.py::test_antisense_transcript_reads_splice_sites_on_the_gene_strand`, `tests/test_41_ui_m1b.py::test_splice_header_separates_alignment_orientation_from_transcript_strand`.

Files: `app/backend/teagle_core/wsl.py`, `app/native/panel_splice.py`.

### Input, fetch and record validity

#### CLS-015 — Missing pyhmmer is reported as a failed domain stage, not as 'no domains detected'

**Before (4.0.0).** scan_domains returned [] when pyhmmer could not be imported; the engine recorded domains status 'complete' (count 0, orfsScanned = ORF count) and the call read as a complete negative.

**After (4.1.0).** scan_domains raises RuntimeError('pyhmmer unavailable: ...'); the engine records processing.domains {'status':'failed','errorType':'RuntimeError'}, adds the note 'domain scan unavailable: RuntimeError', passes domains_ok False to classify (CLS-014) and the top-level status becomes complete_with_stage_failures.

**Trigger.** Environments where pyhmmer is absent or broken (frozen build with a packaging fault, hermetic CI without the package).

**Affected outputs.** processing.domains; result.status (complete -> complete_with_stage_failures); notes; classification.domains_unavailable (False -> True).

**Measured.** engine.analyze with domains.pyhmmer set to None on M11240: v4.0.0 status complete, domains complete (11 ORFs scanned, 0 hits); 4.1.0 status complete_with_stage_failures, domains failed, note present, domains_unavailable True. Same on X01005 (see CLS-014).

**Residual limitations.** Fault-injection only; not an input property, so a benchmark on bundled data cannot trigger it.

**Benchmark probe.** Positive: run any element with pyhmmer made unimportable. Negative: normal environment (identical). Boundary: a record with zero ORFs (domain status not_run, not failed; CLS-029).

**Pinned by.** `tests/test_41_classify.py::test_scan_domains_raises_when_pyhmmer_is_unavailable`, `tests/test_41_classify.py::test_engine_reports_the_domain_stage_as_failed_when_pyhmmer_is_unavailable`.

Files: `app/backend/teagle_core/domains.py`, `app/backend/engine.py`.

#### CLS-028 — Empty and all-N records are invalid with a stated reason and are not classified

**Before (4.0.0).** A header with no sequence, or a record of only N/ambiguity codes, was valid and classified 'no clear TE signature' (te_class none, outcome complete) with all stages 'complete' (count 0); the run status was 'complete'.

**After (4.1.0).** A record with no A/C/G/T is invalid with invalidReason ('empty record (a header with no sequence)' or 'no informative bases (only N and ambiguity codes)'), all stages not_run, classification null; with no other valid record the run status is no_valid_records and the warning names the skipped record. The result schema is 1.1; 1.0 documents (no invalidReason, no domains not_run) still validate under 4.0.0 rules.

**Trigger.** Empty-sequence FASTA entries and records whose sequence has no unambiguous base.

**Affected outputs.** record valid/invalidReason; processing.* (not_run); classification (null); result.status (complete -> no_valid_records); result.warning; resultSchemaVersion 1.0 -> 1.1.

**Labels.** added: invalidReason (record field); renamed: resultSchemaVersion '1.0' -> '1.1'.

**Measured.** engine.analyze on '>a' and '>a' + 300 N: v4.0.0 status complete, valid True, 'no clear TE signature'; 4.1.0 status no_valid_records, valid False with the reasons above.

**Benchmark probe.** Positive: '>a\n' and 300 N. Negative: any normal record (unchanged). Boundary: 300 N plus a single A (valid, classified).

**Pinned by.** `tests/test_41_engine.py::test_empty_record_is_invalid_with_a_reason_and_never_classified`, `tests/test_41_engine.py::test_ordinary_record_is_unchanged_by_the_empty_record_rule`, `tests/test_41_engine.py::test_schema_is_1_1_and_both_it_and_a_4_0_0_shaped_1_0_document_validate`.

Files: `app/backend/engine.py`.

#### GEN-009 — Backend failure handling: a dead WSL session fails annotation; probe timeouts give plain errors; 'backend not ready' names the failing condition

**Before (4.0.0).** annotate parsed whatever came back: a session that died (return code 4294967295, empty output) or returned exit 1 beside an 'EXIT 0' marker produced ok True with 0 hits, which could be sealed; the not-ready error read 'not ready (RepeatMasker=None, Dfam=False)'.

**After (4.1.0).** annotate requires the EXIT marker and a zero session return code (the script ends with exit $ec) and otherwise returns ok False ('the WSL session ended before RepeatMasker reported, so nothing was searched' or 'RepeatMasker reported success but the WSL session returned exit N, so the result was discarded'). _backend_not_ready names the condition (WSL not installed or running; RepeatMasker not installed; Dfam root and curated libraries not both installed; FamDB could not open the library) and points to the Backend installer; custom-library genome runs skip the Dfam conditions. genome_scan's cached-genome checksum probe now has a 1800 s ceiling (was 30 s; covers a 25 GB target at 15 MB/s) and genome_prepare's tool probe is guarded: a slow WSL gives 'checking the cached genome's checksum timed out after 1800 s ... try again' instead of a raw TimeoutExpired.

**Trigger.** WSL shut down mid-run, crashes, or an unready backend.

**Affected outputs.** annotate.ok / error text; no sealed zero-hit result from a dead session.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| checksum probe timeout (s) | 30 | 1800 | no | bounded by the slowest plausible read of a large target; WSL liveness is probed first |

**Measured.** Unit tests with a stubbed session; the live annotate on M11240 (GEN-006) returns ok True with the EXIT marker in both versions.

**Residual limitations.** Other marker-parsed probes already failed closed in 4.0.0.

**Benchmark probe.** Positive: kill wsl.exe during RepeatMasker on a large input: 4.0.0 returns ok True, 0 hits; 4.1.0 ok False. Negative: normal run: identical. Boundary: RepeatMasker non-zero exit already failed in 4.0.0.

**Pinned by.** `tests/test_41_wsl.py::test_annotate_dead_session_is_not_a_sealed_zero_hit_result`, `tests/test_41_wsl.py::test_annotate_nonzero_session_rc_fails_even_with_a_zero_exit_marker`, `tests/test_41_wsl.py::test_not_ready_error_names_the_failing_condition`, `tests/test_41_wsl.py::test_scan_checksum_probe_timeout_is_a_plain_error_with_room_for_a_large_genome`.

Files: `app/backend/teagle_core/wsl.py`.

#### INP-001 — Accession validation: NZ_ accessions accepted; protein and assembly accessions get specific refusals

**Before (4.0.0).** ACC_RE was ^[A-Za-z]{1,6}_?\d+(\.\d+)?$, so NZ_CP012345.1 was refused; protein (NP_, XP_, YP_, WP_, AP_, ZP_) and assembly (GCF_, GCA_) accessions matched the pattern and went to NCBI as nucleotide requests; the refusal example was 'M11240, NC_003075.7' (an 18.6-Mb chromosome).

**After (4.1.0).** ACC_RE accepts an optional upper-case NZ_ prefix; protein prefixes are refused with 'is a protein accession - TEagle fetches nucleotide records (e.g. X04456)' and GCF_/GCA_ with 'is an assembly accession, not a sequence record - use Fetch by coordinate'; the generic refusal example is 'X04456, AB087615'.

**Trigger.** Typed accessions with an NZ_ prefix, protein accessions, assembly accessions.

**Affected outputs.** validate_accession result / FetchError text; fetch-by-accession availability for NZ_ records.

**Labels.** added: protein-accession refusal, assembly-accession refusal; renamed: refusal example 'M11240, NC_003075.7' -> 'X04456, AB087615'.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| accepted accession grammar | [A-Za-z]{1,6}_?digits(.version) | (NZ_)?[A-Za-z]{1,6}_?digits(.version), protein and assembly prefixes refused | no | RefSeq WGS/plasmid accessions; protein and assembly records are not nucleotide records |

**Measured.** Run in this audit (fetch.validate_accession): 'NZ_CP012345.1' refused in 4.0.0, accepted in 4.1.0; 'GCF_000001405.40' accepted in 4.0.0 (would be sent as a nucleotide id), refused in 4.1.0 with the assembly message; 'NP_000509.1' accepted in 4.0.0, refused in 4.1.0; X04456, M11240.1 and NC_003075.7 accepted in both; bare 'NZ_' and 'not!an!acc' refused in both (new example text).

**Residual limitations.** The NZ_ prefix is matched in upper case only: 'nz_cp012345.1' is still refused (measured), whereas other accessions are accepted in any case. A syntactically valid accession can still be absent at NCBI.

**Benchmark probe.** Positive: fetch accession NZ_CP012345.1 (syntax only; a real NZ_ record such as a bacterial plasmid accession for a live test): 4.0.0 FetchError at validation, 4.1.0 proceeds. Negative: X04456 identical. Boundary: 'GCF_000001405.40': 4.0.0 reaches NCBI, 4.1.0 refuses locally.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_nz_and_lowercase_accessions_validate`, `tests/test_41_gff3_fetch.py::test_bare_nz_prefix_is_still_invalid`, `tests/test_41_gff3_fetch.py::test_protein_and_assembly_accessions_get_a_targeted_message`, `tests/test_41_gff3_fetch.py::test_error_example_is_a_te_sized_record_not_a_chromosome`.

Files: `app/backend/teagle_core/fetch.py`.

#### INP-002 — Fetch uses the resolved accession.version; ENA headers are rewritten; unversioned and ENA-served records keep fetched provenance

**Before (4.0.0).** retrieve() fetched the FASTA with the text the user typed and built endpoint and sourceUrl from it. A typed unversioned accession (X04456, M11240) produced a source context that did not match the verified accession.version, and an ENA-served record kept the header '>ENA|acc|acc.ver ...': in both cases normalize_source_context returned local_input, so the sealed manifest lost accession, taxon and source URL.

**After (4.1.0).** retrieve() fetches with the NCBI-resolved accession.version from NCBI or ENA alike, rewrites an ENA header to '>acc.ver ...' only when its third field equals that accession (otherwise untouched, so the record stays local_input), and builds endpoint and sourceUrl from the versioned accession (a cache written from the typed text is repaired on read). provenance accepts a source URL or ENA endpoint built from the typed unversioned form of the resolved accession and always emits the canonical versioned URL; URLs naming another record, version or host and headers naming another accession still downgrade to local_input.

**Trigger.** Accessions typed without a version; any record served by the ENA fallback (NCBI efetch failure).

**Affected outputs.** provenance.input.inputType: local_input -> fetched_accession; input.accession, taxid, source, sourceUrl and endpoint in the manifest; identitySource: local_input -> validated_fetch; FASTA header of an ENA-served record: '>ENA|X|X.1 ...' -> '>X.1 ...'.

**Labels.** renamed: inputType 'local_input' -> 'fetched_accession' for unversioned and ENA-served fetches.

**Measured.** Run in this audit with real engine/fetch code and a mocked network layer (canned esummary and FASTA bodies for M11240.1; refresh True; the shared on-disk cache was not touched). Typed 'M11240.1', NCBI: fetched_accession in both versions. Typed 'M11240', NCBI-served: 4.0.0 manifest inputType local_input, accession absent, identitySource local_input; 4.1.0 fetched_accession, accession M11240.1, validated_fetch, versioned sourceUrl. ENA-served (typed with or without version): 4.0.0 local_input (header '>ENA|M11240|M11240.1 desc'); 4.1.0 fetched_accession, header rewritten to '>M11240.1 desc', endpoint https://www.ebi.ac.uk/ena/browser/api/fasta/M11240.1.

**Residual limitations.** Both code versions share the same on-disk fetch cache; a 4.1.0 run rewrites legacy entries on read, so a 4.0.0 run after a 4.1.0 run can see a repaired cache entry. An ENA record whose third header field names a different accession stays local_input by design.

**Benchmark probe.** Positive: fetch M11240 (typed without a version; resolves to M11240.1) and run the analysis: 4.0.0 seals local_input, 4.1.0 seals fetched_accession with taxon 7227 and the versioned NCBI URL. Negative: type M11240.1: identical in both. Boundary: force the ENA fallback (NCBI blocked) for X01005: header rewritten and provenance fetched in 4.1.0, local_input in 4.0.0.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_ena_fallback_yields_fetched_accession_provenance`, `tests/test_41_gff3_fetch.py::test_ena_record_naming_a_different_accession_stays_local_input`, `tests/test_41_gff3_fetch.py::test_unversioned_typed_accession_still_seals_as_fetched_from_ncbi`, `tests/test_41_gff3_fetch.py::test_cache_written_from_the_typed_accession_is_repaired_on_read`.

Files: `app/backend/teagle_core/fetch.py`, `app/backend/teagle_core/provenance.py`.

#### INP-003 — Fetched gene models: absent codon_start means frame 1, introns only within one feature or gene, derived exons and introns flagged, stale caches rebuilt

**Before (4.0.0).** A CDS without codon_start got phase None (GFF3 export failed); introns were the gaps between all exons of a record, so separate CDS features or genes gave an intergenic 'intron'; derived introns had no derived flag; derived exons inherited cds_group/phase; complete_gene_model re-derived gaps globally (quadratic).

**After (4.1.0).** Absent codon_start gives phase 0 (running phase for later segments); a present invalid value (empty, 0, 4, x) still gives None and the GFF3 writer fails closed. Introns are derived only within one unit (one CDS/mRNA feature, or explicit exons grouped by gene/locus_tag/transcript_id plus strand) and carry derived=True (GFF3 derived=true, BED |derived=true, FASTA ' derived=true'); derived exons drop cds_group/phase; complete_gene_model splits derived introns around synthesised exons (bisect index). GENE_MODEL_VERSION 2 stamps the cache: older models are rebuilt once; a failed feature fetch is stored as featuresStatus 'failed' and retried on the next cache hit.

**Trigger.** Any fetched accession with a feature table (CDS, mRNA, exon, intron features); GFF3/BED export of a fetched gene model.

**Affected outputs.** gene_model.cds[].phase; gene_model.introns (count and derived flag); GFF3/BED/FASTA export of fetched models; meta.featuresStatus, meta.featuresModel (new).

**Labels.** added: derived=true on introns, featuresStatus, featuresModel; removed: cds_group/phase on derived exon rows.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| GENE_MODEL_VERSION | none | 2 | no | forces rebuild of cached models built by the old logic |
| codon_start default | None (unavailable) | 1 (phase 0) when the qualifier is absent | no | INSDC feature-table convention |

**Measured.** Run in this audit (fetch.build_gene_model on synthetic feature lists): one CDS feature joining 10-100 and 200-301, no codon_start: 4.0.0 phases [None, None], intron 100-200 without a derived flag; 4.1.0 phases [0, 0], the same intron with derived True. Two single-interval CDS features of different genes (10-100 and 400-600): 4.0.0 one intron 100-400; 4.1.0 none (counts introns 1 -> 0).

**Residual limitations.** Records that annotate separate CDS features of one transcript without a joined mRNA get no introns even where biology implies them. Exon and intron counts of fetched records change, so any benchmark that scores fetched gene models must be rebuilt.

**Benchmark probe.** Positive: fetch a multi-CDS record such as M80343 (L1.2, two ORF CDS features) and read the gene model: 4.0.0 reports a derived 'intron' between ORF1 and ORF2, 4.1.0 none. Negative: a record with a single joined CDS (for example the insulin example used for splice) keeps its intron in both, flagged derived in 4.1.0. Boundary: codon_start absent vs '2': phase 0 vs 1.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_absent_codon_start_means_frame_one_and_gff3_exports`, `tests/test_41_gff3_fetch.py::test_present_but_invalid_codon_start_is_still_unavailable`, `tests/test_41_gff3_fetch.py::test_separate_unspliced_cds_features_give_no_introns`, `tests/test_41_gff3_fetch.py::test_gypsy_style_orf_cds_features_give_no_introns`.

Files: `app/backend/teagle_core/fetch.py`, `app/backend/teagle_core/gff3.py`.

#### INP-004 — The record-transcript picker extracts accessions per field, so TSA and local ids no longer give wrong accessions

**Before (4.0.0).** transcript_id was searched with a regex ([A-Z]{2}_\d+ then [A-Z]{1,2}\d+) anywhere in the string, so GAAA01000001.1 yielded AA01000001.1 and gnl|WGS:JAABCD|FUN_000123-T1 yielded UN_000123.

**After (4.1.0).** The value is split on '|' and each whole field is tested against the accession grammar (protein and assembly prefixes excluded); gnl| and lcl| local identifiers and non-accession values give nothing.

**Trigger.** mRNA features whose /transcript_id is a TSA accession, a gnl|/lcl| id or a multi-field GenBank string.

**Affected outputs.** gene_model.transcripts[] offered by the splice picker.

**Measured.** Run in this audit: 'GAAA01000001.1' -> 4.0.0 'AA01000001.1', 4.1.0 'GAAA01000001.1'; 'gnl|WGS:JAABCD|FUN_000123-T1' -> 4.0.0 'UN_000123', 4.1.0 no transcript; 'ref|NM_001.2|' -> 'NM_001.2' in both.

**Residual limitations.** A transcript whose id is a local identifier is not offered at all.

**Benchmark probe.** Positive: mRNA feature with transcript_id gnl|WGS:JAABCD|FUN_000123-T1: 4.0.0 offers UN_000123, 4.1.0 none. Negative: ref|NM_001.2| equal. Boundary: a TSA accession with four letters.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_transcript_picker_accessions_are_anchored`.

Files: `app/backend/teagle_core/fetch.py`.

#### INP-005 — Chromosome resolution: mitochondrion in six bundled assembly maps, plastid synonyms, paginated custom maps, error lists valid names

**Before (4.0.0).** Sequence reports were read from the first page only (page_size 1000), so organelles listed after the unplaced scaffolds were missing: six bundled maps lacked the mitochondrion (D. melanogaster, D. rerio, A. aegypti, B. taurus, M. mulatta, T. aestivum) and 'chrM' failed for them; a cached custom map was reused whatever its age; plastid names were not aliased; the not-found message gave no valid names; an alias could shadow an exact name.

**After (4.1.0).** Both readers follow next_page_token (verified live against NCBI Datasets v2); the six bundled maps gained their MT molecule (NC_024511.2, NC_002333.2, NC_035159.1, NC_006853.1, NC_005943.1, NC_036024.1 respectively; 20 other maps are byte-identical); cached custom maps lacking a 'paginated' stamp are rebuilt; resolve_chrom matches the exact name first, then MT aliases (M, MT, chrM, chrMT) and plastid aliases (C, ChrC, Pt, Pltd, chloroplast); an unknown name lists up to 40 available names.

**Trigger.** Coordinate fetches on the mitochondrion of the six assemblies; plastid names in plant assemblies; custom assemblies.

**Affected outputs.** resolve_chrom result / CoordError text; bundled assembly maps (6 files, +1 molecule each).

**Labels.** added: MT molecules in six maps, '(available: ...)' in the not-found error.

| Parameter | 4.0.0 | 4.1.0 | Frozen and unchanged | Basis |
|---|---|---|---|---|
| assembly maps bundled | 26 (6 without the mitochondrion) | 26 (+1 molecule in 6) | no | generated under a diff gate: existing molecule lists identical as a prefix, additions only |

**Measured.** Run in this audit (fetch.resolve_chrom, bundled maps, offline): D. melanogaster GCF_000001215.4 'chrM', 'MT', 'M': 4.0.0 not found; 4.1.0 NC_024511.2 (19,524 bp). D. rerio GCF_000002035.6 chrM: not found vs NC_002333.2. Arabidopsis GCF_000001735.4 'Pt', 'chrC', 'C', 'chloroplast': not found vs NC_000932.1 ('Pltd' resolved in both). Human and yeast chrM unchanged (NC_012920.1, NC_001224.1); 'mitochondrion' unresolved in both, with a names list only in 4.1.0.

**Residual limitations.** The plastid aliases 'C' and 'Pt' apply only when no chromosome has that exact name. Aliases other than the listed ones are not guessed.

**Benchmark probe.** Positive: fetch chrM:1-300 of D. melanogaster (GCF_000001215.4): 4.0.0 'chromosome chrM not found', 4.1.0 NC_024511.2:1-300. Negative: human chrM unchanged. Boundary: Arabidopsis 'Pt' resolves to NC_000932.1.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_bundled_large_assembly_maps_include_the_mitochondrion`, `tests/test_41_gff3_fetch.py::test_live_assembly_map_follows_next_page_token`, `tests/test_41_gff3_fetch.py::test_build_script_sequence_report_follows_next_page_token`, `tests/test_41_gff3_fetch.py::test_stale_cached_custom_map_without_the_pagination_stamp_is_rebuilt`.

Files: `app/backend/teagle_core/fetch.py`, `app/backend/teagle_core/build_assembly_maps.py`, `app/backend/teagle_core/data/assemblies/GCF_000001215.4.json`, `app/backend/teagle_core/data/assemblies/GCF_000002035.6.json`, `app/backend/teagle_core/data/assemblies/GCF_002204515.2.json`, `app/backend/teagle_core/data/assemblies/GCF_002263795.3.json`.

#### INP-006 — Only one assembly per organism: adding a second assembly of a pinned organism is refused

**Before (4.0.0).** add_custom_assembly stored a second assembly for an already pinned organism and reported success; the entry could never be selected (the curated or first entry wins) and a second D. simulans add silently replaced the first.

**After (4.1.0).** Inside the lock, an organism that is curated or already stored with a different accession raises CoordError ('<org> is already pinned to <name> (<accession>); one assembly per organism is supported, so <new> was not added'); re-adding the identical accession is a no-op.

**Trigger.** Adding a custom assembly for an organism that already has one (including the 17 curated organisms).

**Affected outputs.** add_custom_assembly return / error; custom assembly store contents.

**Measured.** From tests (network assembly resolution not run).

**Residual limitations.** A different assembly of a pinned organism cannot be used at all.

**Benchmark probe.** Positive: add the T2T human assembly while hg38 (GCF_000001405.40) is pinned: 4.0.0 reports success with an unusable entry, 4.1.0 refuses. Negative: add a new organism. Boundary: re-add the pinned accession: no-op.

**Pinned by.** `tests/test_41_gff3_fetch.py::test_adding_an_alternative_assembly_of_a_curated_organism_is_refused`, `tests/test_41_gff3_fetch.py::test_second_assembly_of_a_custom_organism_is_refused_and_keeps_the_first`, `tests/test_41_gff3_fetch.py::test_re_adding_the_same_accession_is_idempotent_and_new_organisms_are_added`, `tests/test_41_gff3_fetch.py::test_re_adding_the_curated_accession_is_a_noop`.

Files: `app/backend/teagle_core/fetch.py`.

#### INP-009 — Multi-region coordinate fetches can be analysed and exported

**Before (4.0.0).** Analysing a two-region coordinate fetch produced an analysis manifest for which provenance.verify_manifest returned False (contradictory identity rules), so the result document could not be exported.

**After (4.1.0).** A fetched_coordinate input with more than one region requires top-level identitySource 'generated' and id 'multiple_records'; a single region keeps 'validated_fetch'. The two-region analysis verifies and exports; the specimen panel states that every region is analysed and sealed, the record to view is chosen in the summary table, and genome-anchored GFF3/BED export covers the first region only.

**Trigger.** Coordinate fetches with more than one region.

**Affected outputs.** manifest verification and analysis-result export; input.identitySource / id for multi-region fetches.

**Measured.** Run in this audit (real engine/fetch code, mocked efetch bodies of synthetic random sequence for chr1:1,000,001-1,000,600 and chr2:2,000,001-2,000,600 on the bundled human map, refresh True): one region: manifest verifies and the document validates in both versions (identitySource validated_fetch). Two regions: 4.0.0 manifest_verifies False and the export raises 'analysis result is not bound to its provenance manifest'; 4.1.0 manifest verifies, the document validates (identitySource generated, id multiple_records, 2 records, status complete).

**Residual limitations.** Only the first region can be exported as genome-anchored GFF3/BED.

**Benchmark probe.** Positive: fetch two regions of one assembly in one request and analyse: 4.0.0 cannot export, 4.1.0 can. Negative: single region identical. Boundary: two regions on different chromosomes.

**Pinned by.** `tests/test_41_engine.py::test_multi_region_coordinate_analysis_verifies_and_exports`, `tests/test_41_engine.py::test_multi_region_identity_rules_still_reject_contradictions`.

Files: `app/backend/teagle_core/provenance.py`, `app/backend/engine.py`, `app/native/panel_specimen.py`.

#### INP-012 — Results for a previous specimen are discarded; a slow fetch cannot overwrite a later sequence; a fetch can be cancelled

**Before (4.0.0).** A late family-annotation, splice, primer-design, in-silico-PCR or genome-scan result from a previous specimen was shown as current; a slow accession or coordinate fetch could overwrite a sequence pasted, uploaded or loaded after it; a failed fetch erased the loaded specimen's label and genome binding; there was no way to cancel a fetch.

**After (4.1.0).** Requests are tagged '<op>#n' with the analysis generation and settled in _retire; a result whose specimen changed is dropped and busy buttons reset; both fetch types share _start_fetch (key 'fetch#n'), any specimen change retires a pending fetch so its result is ignored, and a Cancel button appears while a fetch runs; a failed, cancelled or superseded fetch only clears the 'fetching...' indicator and restores the prior label.

**Trigger.** GUI use with slow networks or quick specimen changes.

**Affected outputs.** which result is displayed; no engine output changes.

**Labels.** added: Cancel (fetch).

**Measured.** From tests (GUI-only).

**Residual limitations.** A cancelled fetch still completes on the network; its result is ignored.

**Benchmark probe.** Positive (GUI): start a coordinate fetch, paste another sequence before it returns: 4.0.0 replaces the paste, 4.1.0 keeps it. Negative: no interleaving: identical. Boundary: cancel then fetch again.

**Pinned by.** `tests/test_41_ui_m1a.py::test_late_family_result_after_specimen_change_is_dropped`, `tests/test_41_ui_m1a.py::test_late_splice_result_after_specimen_change_is_dropped`, `tests/test_41_ui_m1a.py::test_late_primers_after_specimen_change_never_render_current`, `tests/test_41_ui_m1a.py::test_genome_scan_started_for_previous_specimen_is_discarded`.

Files: `app/native/panel_specimen.py`, `app/native/ui_common.py`, `app/native/panel_family.py`, `app/native/panel_splice.py`, `app/native/panel_primers.py`.

#### PCR-010 — In-silico PCR batches are bound to the record they were run on

**Before (4.0.0).** A batch still running when another record was selected was shown against the new record; a finished result stayed on screen after the record changed.

**After (4.1.0).** Each batch stores the record index and digest it ran on; switching records discards a batch still running for the previous record and clears a finished result that belongs to another record. Requests are tagged ('pcr#<gen>#<i>', '<op>#n') and a late result from a previous specimen is discarded; family and splice results returning during a same-specimen re-analysis are kept and shown afterwards.

**Trigger.** Multi-record input or re-analysis while a job runs.

**Affected outputs.** displayed PCR/gel/amplicon table.

**Measured.** From code and tests.

**Residual limitations.** GUI-only; engine runs are unaffected.

**Benchmark probe.** Positive (GUI): run a PCR batch on record 1 of a two-record FASTA, switch to record 2 before it finishes: 4.0.0 renders it under record 2, 4.1.0 discards it. Negative: no switch: same result. Boundary: switch after completion: cleared in 4.1.0.

**Pinned by.** `tests/test_41_ui_m1a.py::test_pcr_batch_for_another_record_is_dropped_and_the_button_restored`, `tests/test_41_ui_m1a.py::test_a_finished_pcr_for_another_record_does_not_stay_on_screen_as_current`, `tests/test_41_ui_m1a.py::test_pcr_batch_in_flight_is_dropped_on_specimen_change`, `tests/test_41_ui_m1a.py::test_a_result_arriving_while_the_same_specimen_is_reanalysed_is_kept_and_shown`.

Files: `app/native/panel_pcr.py`, `app/native/panel_records.py`, `app/native/ui_common.py`.

#### STR-004 — Ambiguity codes are neither matches nor mismatches; ACGT-only TSD; input upper-cased

**Before (4.0.0).** v4.0.0: LTR and TIR identity divided matches by the full arm length, so N and IUPAC columns counted as mismatches (only N and n excluded from the match count); the LTR X-drop extension penalised N as a mismatch and stopped at an N-run; find_tsd accepted any equal flank strings that held no 'N', so an IUPAC TSD such as ACRYT was reported as a 5-bp TSD; PBS identity counted only N/n as unknown. find_ltr, find_tir and detect_all did not upper-case their input.

**After (4.1.0).** 4.1.0: a column counts only when both bases are A, C, G or T. Identity is matches / known columns (LTR ungapped, TIR scan) or matches / (matches + mismatches + gap columns) with ambiguous columns excluded (gapped). The X-drop extension and the gapped alignment treat an ambiguous column as neutral, so an N-run no longer stops extension. LTR rows carry ambiguous_columns when non-zero. find_tsd (via _tsd_at) accepts only ACGT motifs. PBS identity counts a match only for A/C/G/T. find_ltr, find_tir and detect_all upper-case the sequence (the engine already normalises case, so engine results are the same for lower-case input).

**Trigger.** Records containing N or IUPAC codes inside or next to a terminal repeat, a TSD or a PBS leader window; direct callers passing lower-case sequence to the detector functions.

**Affected outputs.** LTR/TIR row identity (higher when ambiguity was present); LTR row new key ambiguous_columns; TSD row absent where the only match was an IUPAC string; LTR length where an N-run previously stopped the extension.

**Labels.** added: ambiguous_columns (LTR row).

**Measured.** Own run: M11240.1 copia with 40 N written over the 3' LTR: identity v4.0.0 85.5%, 4.1.0 100.0% with ambiguous_columns 40 (length 276 both). M01005... X01005.1 Tc1 with 3 N inside the 5' arm: TIR identity 94.4% -> 100.0% (54 bp both). Copia flanked by the IUPAC string ACRYT as TSD: v4.0.0 reports a 5-bp TSD 'ACRYT'; 4.1.0 reports no TSD. Lower-case copia and Tc1 through detect_all: identical rows in both versions.

**Residual limitations.** An IUPAC code that is compatible with its partner (R against A) is still not counted as a match, so identity is a lower bound on true identity in that column but the column is excluded from the denominator. No ambiguity handling for PPT or poly(A)-signal motifs (those functions are unchanged).

**Benchmark probe.** Positive: M11240.1 with 40 consecutive N inside the 3' LTR: v4.0.0 identity 85.5, 4.1.0 identity 100.0 and ambiguous_columns 40. A flanked element whose both flanks carry the 5-bp IUPAC string ACRYT: v4.0.0 TSD row of 5 bp, 4.1.0 none. Negative: the same record with no N must give identical identity in both versions. Boundary: N fraction at which the extension stops (4.0.0 stops at any N-run of about 15 bases; 4.1.0 does not).

**Pinned by.** `tests/test_41_structural.py::test_ambiguous_bases_are_neither_matches_nor_mismatches`, `tests/test_41_structural.py::test_ambiguity_codes_are_not_a_target_site_duplication`.

Files: `app/backend/teagle_core/structural.py`.

## Output-format records

These change what a parser reads, not the science. Each needs an adapter or parser update before a 4.1.0 run is scored.

| ID | Change | 4.1.0 behaviour | Affected outputs | Probe |
|---|---|---|---|---|
| CLS-022 | Classification constants and the TSD expectation table are sealed in the provenance parameters | The sealed parameter set has 88 keys: the 23 originals unchanged in value plus classify_chance_tir_max_bp (17), classify_chance_tir_sequences (22,500), structural_chance_level_p (1e-3), structural_chance_tir_gc_range ([0.35, 0.65]), tsd_expectation_by_superfamily (Copia/Gypsy length 5 range 4-6, Tc1/Mariner length 2 range 2-2 motif TA, hAT length 8 range 8-8) and the other numeric structural constants and signature defaults under derived keys. Values are read from the defining code. | provenance.parameters; analysis manifest seal / resultPayloadSha256 inputs | Positive: compare the parameters block of a 4.0.0 and a 4.1.0 manifest for the same record (23 vs 88 keys, hashes differ). Negative: none (always differs). Boundary: n/a. |
| CLS-026 | New classification fields wicker_order and domain_architecture | Two new fields: domain_architecture (canonical order for Copia/Gypsy and the withheld retroviral label; otherwise the detected codes in canonical order; 'none detected'; 'unavailable' when the domain scan failed) and wicker_order with the vocabulary LTR, TIR, LINE, DIRS, withheld, unassigned (te_class LTR* -> LTR; DNA/* including DNA/non-autonomous -> TIR; LINE; DIRS; outcome withheld -> withheld; everything else, including retro/*, unclassified, repeat/structural-only (when not withheld) and the Helitron call -> unassigned). 'order' is unchanged. GFF3/BED export both new attributes with classification_order kept as a deprecated alias of domain_architecture, and the BED name carries the Wicker order. | classification.domain_architecture (new); classification.wicker_order (new); GFF3 TE parent attributes wicker_order, domain_architecture, classification_order; BED name order= | Positive: any 4.1.0 result (field present; 4.0.0 absent). Negative: 4.0.0 results (KeyError on the new fields). Boundary: lone RT (withheld), RT+CHR (LTR) and CHR alone (unassigned) differ. |
| CLS-029 | Domain stage with no ORF to read is 'not_run'; classification still treats it as searched | processing.domains status is 'not_run' (count 0, orfsScanned 0) when no ORF was searched (schema 1.1; 1.0 validators refuse it); result.stageSummary counts records per stage state. The call is made as before: domains_ok stays True, so classify reads absence of domains as 'no coding domain detected'. | processing.domains.status; result.stageSummary (new) | Positive: a 60 bp element with two 30 bp inverted arms and no ORF >= 40 aa (4.0.0 domains 'complete', 4.1.0 'not_run'; call identical). Negative: any record with at least one ORF. Boundary: 40-aa versus 39-aa ORF. |
| CLS-031 | Analysis-result export of flanked Tc1/Mariner records whose TA target site is recovered in classification | Fixed in 4.1.0 before release: the engine refreshes processing.structural.count after classification, before the result is sealed, so the TSD row that classification adds is counted and the analysis-result JSON exports and validates. (Found on the audit branch by this inventory: the count was recorded before classification, so 4 of 6 flanked TA-TSD constructs could not be exported.) | analysis-result JSON export (fails); processing.structural.count (understated by 1); validate_analysis_result_document (False) | Positive: X01005 + 300 bp random flanks + TA TSD (synthetic, seeds 0, 2, 4, 5): document build raises in 4.1.0; 4.0.0 builds. Negative: bare X01005 and flanked X05424. Boundary: a flanked Tc1 whose detect_all output already had a TSD row (seed 1, 3) exports. |
| CLS-032 | Provenance reference names no longer state Dfam or RepeatMasker versions; Pfam name states the accession-version count | Pfam name 'Pfam-A profiles (via InterPro, 30 accession versions, retrieved 2026-07-17/23/28)'; Dfam name 'Dfam (transposable-element family database)'; RepeatMasker name 'RepeatMasker (RMBLAST engine)'. The module docstring no longer cites an internal deliverables file. data/README.md lists all 30 profiles with versioned accessions in file order. Citations, DOIs and licences are unchanged. | result.references[].name; manifest references (seal content) | Positive: compare references[].name in a 4.0.0 and a 4.1.0 result. Negative: citation DOIs (unchanged). Boundary: n/a. |
| GEN-004 | Genome-scan table, gel menu and exports: capped table, numeric Start/End, descriptive notes, unmatched-locus message | First 500 rows are drawn (on-target first) with a 'Showing the first N of M' note while sort and exports use every row; Start and End are separate numeric 1-based columns; export notes name both primers, organism, assembly accession and name, isPcr version and thresholds, size window and coordinate convention; default file names carry pair and assembly; the gel menu offers 'Copy locus, 1-based' and '1 of N in this band'; an unmatched locus says the on/off-target split is not assessed. | scan table columns and notes; gel context menu | Positive: scan a repetitive TE pair with more than 500 products: 4.1.0 shows 500 and exports all; 4.0.0 draws all. Negative: a unique pair: identical rows. Boundary: exactly 500 products: no capped note. |
| GEN-012 | Whole-genome annotation manifests name the annotated assembly and refuse to seal a mismatched accession | databases gains {'name': 'RefSeq genome assembly', 'version': <accession>, 'sha256': <genome checksum>} (library stays databases[0]); the source={'assemblyAccession'} argument is dropped; a backend accession that differs from the request, or an accession that is not public, returns ok False / partial True unsealed. The sealed species follows GEN-006. | provenance.databases (+assembly entry); provenance input source fields (assemblyAccession removed); manifest hash | Positive: annotate a prepared yeast genome: 4.1.0 manifest lists the RefSeq assembly with its checksum; 4.0.0 does not. Negative: custom-library run: library entry is databases[0] in both. Boundary: backend returns a different accession: 4.1.0 refuses. |
| INP-010 | Primer, PCR, annotation and splice manifests for bare pasted sequence carry the analysis record's generated identity | _operation_record_identity maps 'input_sequence' to a generated identity (id 'input_sequence', identitySource generated), matching analyze(); the GUI sends a local record under its original header (accepted only when the public identity of that header reproduces the selected record's id and hash, otherwise the bare sequence) so operation manifests join the analysis record by id; headers are never written to a result; PCR tables show a fetched accession readably. | provenance.input.id / idSha256 / identitySource of primer, PCR, annotate, splice manifests; amplicon source labels; manifest hashes | Positive: run_primers on a bare M11240 sequence: 4.0.0 id differs from the analysis record id, 4.1.0 equals it. Negative: headed input equal in both. Boundary: multi-FASTA with the GUI: each record's operations carry that record's id in 4.1.0. |
| INP-011 | Seals are independent of how numbers were typed and of fetch time; the top-level input id is tied to the sealed record bindings | build_manifest writes integral floats in sealed parameters as ints and strips retrievedUtc, source, endpoint and sourceUrl from the input bindings nested in parameters (the top-level input still records them); verify_manifest ties the unsealed top-level id, idSha256 and identitySource to records[0] for one record and requires 'multiple_records' / generated for several. Canonicalisation 'teagle-json-utf8-v1', manifest schema 0.2 and the verification of 4.0.0 manifests are unchanged. | manifestSha256 for runs with integral float parameters; PCR seal invariance to retrievedUtc and serving database | Positive: design twice with opt_tm 60 and 60.0: equal seals in 4.1.0, different in 4.0.0. Negative: different Tm values give different seals in both. Boundary: a manifest re-serialised through jq/JavaScript still verifies in 4.1.0. |
| OUT-001 | Analysis-result schema 1.0 -> 1.1: domains 'not_run' state, invalidReason, stageSummary | resultSchemaVersion 1.1. Invalid records carry invalidReason and every stage 'not_run'. processing.domains may be 'not_run' (no ORF to read). A top-level stageSummary counts stage states over records. Documents at schema 1.0 still validate under the 4.0.0 rules (accepted schemas 1.0 and 1.1). | resultSchemaVersion; records[].invalidReason; records[].processing.domains.status ('not_run'); stageSummary | Schema contract test: validate a 4.0.0 schema-1.0 document and a 4.1.0 schema-1.1 document with the 4.1.0 validator (both pass); an empty and an all-N record must give invalidReason with every stage not_run; a record with no ORF of the minimum length must give domains not_run. A parser written for 1.0 must accept the new keys. |
| OUT-002 | Record-summary table: columns renamed and reordered | Columns lead with Classification, Call confidence, completeness and 'Call status' (renamed from 'Outcome'); 'Wicker order' and 'Domain architecture' are separate columns; the per-domain column is 'Hit confidence'. CSV/TSV/XLSX export headers follow. | record summary CSV/TSV/XLSX headers; domain table header 'Hit confidence' | Any benchmark that reads the app's table exports must map the new headers; parse a 4.1.0 export and assert the header set. |
| OUT-003 | GFF3/BED: element span, neutral direct-repeat terms, candidate PBS, new attributes | The TE feature spans the element (from its terminal repeats or tail), not the record; the neutral terminal direct-repeat candidate is exported with direct_repeat terms; sub-confident PBS matches are exported as candidates; attributes add wicker_order, domain_architecture and tsd_state (classification_order kept as an alias). | TE feature start/end in GFF3 and BED; GFF3 attributes wicker_order, domain_architecture, tsd_state | Export GFF3 for M11240 between 300-bp random flanks and assert the TE feature equals the LTR outer span, not the record; 4.0.0 gives the record span. |
| OUT-004 | Per-record classification TSV in TEsorter .cls.tsv column order | A per-record classification TSV in TEsorter's .cls.tsv column order plus the RepeatMasker-style class string and the Wicker three-letter code. | new export file | Join the 4.1.0 TSV with TEsorter 1.5.1 output on the same FASTA; the join must be one-to-one on record id. |
| OUT-005 | Assay report (PDF/HTML) and notebook entry | An assay report (PDF or HTML) and a copyable notebook entry combine the call, its cautions, the element map, the chosen pairs with locations and QC, the predicted gel, the reaction conditions and the run seals. | new export files | Not a scientific endpoint; export-parity checks (R1 type) should include it. |
| OUT-007 | No author or creator field in end-user exports | PNG, PDF, SVG, XLSX and the assay report carry no author or creator field; the metadata checks take author=None. | export file metadata | Export checks that compare metadata must expect no author or creator. |
| OUT-008 | ID map export only from the record row menu; 'Your label' column on screen only | The ID map export is offered only from a record row's right-click menu after a confirmation; an on-screen 'Your label' column shows the first word of each FASTA header; exports keep the neutral record identifiers. | ID map file; record table (screen only) | Benchmarks keyed on record identifiers are unaffected; exported identifiers stay neutral. |
| OUT-010 | CSV/TSV caveat lines readable by comment-skipping readers | Caveat lines start with '# ' so pandas (comment='#') and R (comment.char='#') skip them; Ctrl+C copies clean values; export failures are reported. | CSV/TSV preamble lines | Readers of CSV/TSV exports must skip '#' lines. |
| PCR-003 | Amplicon sequence carries the primer's own bases over its whole footprint | seq = left oligo + template between the two annealed footprints + reverse complement of the right oligo (overlapping footprints: the left oligo covers the overlap and the leading characters of the reverse-complemented right oligo are dropped), so len(seq) == product_length and the product begins with the forward oligo and ends with the reverse complement of the reverse oligo. Coordinates, mismatch counts, tail and anneal lengths are unchanged. Degenerate primers now put IUPAC letters into seq. | amplicon.seq / amplicon FASTA; templated_seq is unchanged | Positive: tailed F on M11240 1000-1020: 4.1.0 seq starts with the full tailed oligo; 4.0.0 does not. Negative: a perfectly matching untailed pair gives identical seq. Boundary: a pair whose footprints overlap by a few bases: len(seq) == product_length in both. |
| PCR-007 | Effective 3' strict window and cautions are reported | criteria.three_prime_strict_effective = min(tp, min_anneal) is added; a display 'cautions' entry explains the case tp > min_anneal; result.cautions is always present (an empty list when nothing applies). | result.cautions: absent -> list; criteria.three_prime_strict_effective: absent -> int | Positive: run_pcr with tp 20: caution present only in 4.1.0. Negative: tp 5: cautions [] in 4.1.0, key absent in 4.0.0, identical products. Boundary: tp = min_anneal = 18: effective 18, no caution. |
| PCR-011 | PCR output presentation: amplicon table gains tail and anneal columns and per-primer mismatches; gel marks failed and capped lanes | Adds 'Tail F/R' and 'Anneal F/R' beside 'Mism F/R'; Mism F/R maps per-site values to the primer at that site (left_primer/right_primer), so F/R means forward/reverse primer for R+F, F+F and R+R products; tails are flagged in the caution colour; export headers "Mism F/R (annealed 3' core)", "Tail F/R (unmatched 5' nt)", 'Anneal F/R (nt)'; a note states coordinates are 0-based half-open; the FASTA header carries F= and R=; export notes list every lane's primer sequences. Gel (figures.svg_gel): failed lanes are hatched and labelled 'failed', capped searches carry a 'capped' badge, a lane with no designed target prints no red 'no on-target' (the tooltip says why), and the legend lists only the band colours and states actually drawn, including mixed designed-target / target-not-assessed gels. | amplicon table / export columns and headers | Positive: pair with a 5' tail: Tail F/R shows the tail length and Anneal F/R the core, in 4.1.0 only. Negative: perfect pair: tails 0/0 and anneal equal to primer lengths. Boundary: an R+F product: mismatches appear under the correct primer in 4.1.0. |
| PRM-001 | Primer3 reaction conditions are passed explicitly, sealed and reported; Tm values are unchanged | PRIMER_SALT_MONOVALENT 50, PRIMER_SALT_DIVALENT 1.5, PRIMER_DNTP_CONC 0.6 and PRIMER_DNA_CONC 50 are passed explicitly (primers.TM_CONDITIONS) and sealed as parameters.tm_conditions. The result carries design_tm_conditions (sealed) and oligoqc_conditions (QC set, advisory, unsealed). The manifest seals 19 parameters; software entries are 'primer3-py (Primer3 bindings)' and 'libprimer3 (Primer3 core)' (read at run time). Tm headers in the primer table and order sheet name the design conditions. | result.design_tm_conditions: absent -> {mv_conc 50, dv_conc 1.5, dntp_conc 0.6, dna_conc 50}; result.oligoqc_conditions: absent -> {mv_conc 50, dv_conc 0, dntp_conc 0, dna_conc 250, temp_c 37}; provenance.parameters: +tm_conditions, +max_ns_accepted, +template_ambiguity_masked (16 -> 19 keys); provenance.software: 'primer3 (primer3-py)' -> 'primer3-py (Primer3 bindings)', +'libprimer3 (Primer3 core)'; every primer manifestSha256 and every candidate design_provenance seal changes; left_tm / right_tm / penalty / candidate list: unchanged for the same template and parameters | Positive: run_primers on M11240 (Standard preset) in both versions; expect identical candidate tuples (left_seq, right_seq, product_size, left_tm, right_tm, penalty) but parameters.tm_conditions present only in 4.1.0 and 19 vs 16 sealed parameter keys. Negative/control: compare candidate lists only (they must be equal, 0 differences). Boundary: pass a non-default opt_tm (for example 62) and confirm the Tm columns still equal the 4.0.0 values at the same conditions, since the conditions equal the old implicit defaults. |
| PRM-008 | Primer pairs keep one name everywhere: design tag plus Primer3 number; oligos are named by record and pair | Each design gets a batch tag (D1, D2, ... in design order, once per design, not on a theme re-render); the pair label is '<tag>·<Primer3 id>' (for example D2·P3; a pair the user typed keeps its own name) and appears in the primer table ID column, scan picker, PCR queue, gel lane captions (shortened on the gel only), amplicon table, secondary-structure dialog, genome-scan header and every export. Oligo names are '<record id>[_<flank label>]_<pair>_<F\|R>' (vendor-safe), stamped once at submit time so a record switch cannot rename a design's oligos. | ID column: P1 -> D1·P1; FASTA header and filenames: pair token D1-P1 (filename/FASTA token); oligo names in copy-FASTA and the order sheet; amplicon export includes primer sequences and a FASTA header with F= and R= | Positive: press Design twice in the GUI: the second batch's first pair is D2·P1 in the table, queue, gel lane and exports, where 4.0.0 shows P1 twice. Negative: the engine result ids remain P1..Pn in both versions. Boundary: a custom typed pair keeps the user's name (or custom1) with no batch tag. |
| PRM-009 | Primer table gains Location and Notes columns and is reordered; exports gain the matching columns | PRIMER_COLS order: ID, Struct, Penalty, Location, Notes, Forward, Reverse, Product, Fwd pos, Rev pos, Tm F/R, GC% F/R, Hairpin, Self-dim, Hetero, 3'-end. Location = product span, then the segment each primer overlaps (5' flank, 5' LTR/TIR, internal, 3' repeat, 3' flank) from the accepted terminal-repeat pair of widest extent, plus '(DOMAIN)' when a primer lies wholly inside a domain hit; a pair designed on another template reads 'a-b on <label>'. Notes flags (a) a primer whose 3' 18-nt core lies inside one terminal-repeat copy ('may also bind the other copy and solo copies'; for a TIR 'may prime from the inverted ... arm'), (b) pairs whose forward and reverse sites both overlap ('≈P2'), (c) a primer shared by two pairs ('same F as ...'). It flags; it never removes a pair. Exports carry 'Location (0-based, half-open)' and Notes. | primer table and XLSX/CSV/TSV export columns; export header 'Location (0-based, half-open)'; Fwd/Rev pos headers carry the coordinate base | Positive: synthetic flanked M11240, 3' junction goal: pair with F inside the 3' LTR gets the 'may also bind 5' LTR / solo LTRs' note. Negative: M11240 designed untargeted with primers in the internal region: Location 'internal', Notes '-' or only pair-overlap flags. Boundary: a forward primer whose 3' 18-mer ends exactly at the LTR/internal boundary vs one base inside the LTR. |
| PRM-010 | Oligo order sheet export and copy-order-rows buttons | 'Copy order rows' (name TAB sequence per unique oligo) and 'Export order sheet' (XLSX/CSV/TSV) in panels 04 and 05: one row per unique oligo sequence (a primer shared by two pairs is listed once and names both under 'Used in'); columns Name, Sequence (5'-3'), Scale, Purification (both empty), Length (nt), Tm (header names the sealed Primer3 conditions, or says they differ between pairs), GC (%), Used in. Tm and GC are the Primer3 design values, never the QC Tm. Rows follow the on-screen sort; no author field in the file. | new export file; clipboard text | Positive: design three pairs where two share a reverse primer: the sheet has 5 rows, not 6, and the shared oligo's 'Used in' names both pairs. Negative: no pairs: buttons disabled. Boundary: a case-different duplicate of the same sequence merges (comparison is on the upper-cased sequence). |
| STR-017 | PBS rows report a composition-adjusted chance probability; PBS panel scope stated on LTR rows | 4.1.0: every PBS row carries chance_probability (Poisson-binomial over the 18 positions from the record's composition, over all leader windows tested), the panel list, and, when the match is sub-confident, a note quoting the chance ('a match this good arises by chance in about N% of leaders of this composition'). The accepted LTR row carries pbs_scope: 'PBS search covers tRNA-Lys3 only; an element primed by another tRNA (e.g. initiator tRNA-Met in Ty1/copia) is not assessed, so no PBS row is not no PBS'. Match counting uses A/C/G/T only. Reporting floor (0.55), confident threshold (0.72), search window (44 nt) and the panel are unchanged. | PBS row new keys chance_probability, panel; PBS row note text for sub-confident matches; LTR row new key pbs_scope | Positive: random 300-nt sequences (200 seeds) run through find_pbs at ltr end 0: both versions return a row in 23/200; only 4.1.0 reports chance_probability (mean ~0.14) and a note quoting it. Negative: bare M11240.1 and M12927.1 give no PBS row in either version. Boundary: a PBS match at exactly 72% (confident) versus 66.7% (sub-confident, hedged note). |
| STR-018 | scan_scope() and the engine record note for long records with no terminal repeat | 4.1.0: scan_scope(seq) returns record_bp, ltr_seed_window_bp (min(n//2, 6000)), tir_seed_window_bp (min(n//2, 1200)), tir_joint_window_bp (min(n//2, 400)), tail_window_bp (min(n, 600)) and covers_whole_record. After structural detection (and any orientation re-run), when no LTR or TIR row exists and 2 x a seed window is below the record length, engine.analyze appends a record note quoting the record length, the LTR and/or TIR seed windows, the 600-bp tail window and, when 2 x the TIR window is below the length, a sentence that more than 1200 bp of flank on a side hides the element's ends from the TIR scan. classify appends the same windows to its 'no terminal repeat' sentences for records longer than the windows. | record notes (new sentence); classification evidence wording for absent terminal repeats; scan_scope return value for callers | Positive: a 30-kb random record or bare M80343.1: the analysis carries the new note in 4.1.0 and none in v4.0.0. Negative: bare M11240.1 (LTR found) carries no such note. Boundary: records of 2,300 bp (no TIR note) versus 2,500 bp (TIR note) with no terminal repeat. |
| STR-019 | Provenance seal records every structural threshold (23 -> 88 keys) | 4.1.0: the seal is table-driven. The 23 v4.0.0 keys keep their names and values. New keys: ltr_extend_x_drop (30), polya_signal_dse_lo/hi/min_gt_frac/min_t (20, 60, 0.65, 4), classify_chance_tir_max_bp (17) and classify_chance_tir_sequences (22500), tsd_expectation_by_superfamily (nested {length, range, motif}), and 58 derived structural_* keys: every numeric module constant and every numeric signature default in structural.py not already keyed (ranges as [start, stop]). Examples: structural_tir_seed_window 1200, structural_tir_joint_window 400, structural_seed_max_occ 20, structural_chance_level_p 0.001, structural_snap_start [-1,5], structural_snap_end [-4,2], structural_tail_window 600, structural_gapped_below_identity 97.0, structural_scan_tir_min_ident 80. A threshold added to the detectors is sealed automatically (tested). | analysis manifest thresholds block (88 keys); manifest hash and every dependent seal | Positive: compare the thresholds block of the manifest of any run in both versions: 4.1.0 holds 88 keys and v4.0.0 holds 23, the 23 with identical values. Negative: ltr_min_identity_pct, tir_min_identity_pct, tir_min_standalone, tsd_min_len, tsd_max_len, ltr_seed_window must be unchanged. Boundary: a benchmark that hashes the manifest cannot compare hashes across versions (every hash changes). |
| UI-008 | Genome-viewer actions use the feature's strand; figure export fixes | Right-click actions use the feature's own strand; figure exports render correctly. | copied DNA for minus-strand features; figure files | Copy a minus-strand domain and compare it with the reverse complement of its coordinates. |

## Display-only, performance and infrastructure records

| ID | Relevance | Change | 4.1.0 behaviour |
|---|---|---|---|
| STR-010 | performance | Low-complexity seed masking and occurrence caps for LTR and TIR seed search (speed) | 4.1.0: seed words seen more than 20 times in one search window are dropped (both LTR windows and both TIR windows, query and index side). LTR seeds are additionally masked where a DUST score (64-bp window, triplet score above 2.0, dustmasker's default level 20) flags the sequence, and a word masked in either window seeds nothing. TIR 11-mer/6-mer seeds use the occurrence cap only (DUST would mask 9-20% of the Ac and Tc1 specimens and change real calls). At most 64 densest diagonals are kept per TIR seed search; the per-diagonal seed density uses a bisect over sorted anchors instead of a quadratic sum. Seed-pair loop order is preserved so diagonal ties break as before. |
| CLS-002 | display only | ERV evidence sentence names 'errantivirus' only for Gypsy-called elements; Copia/Gypsy evidence wording | The ERV sentence reads 'env-bearing LTR retroelement (ERV architecture)' and appends '; errantivirus lineage for a Gypsy element' only when the superfamily starts with Gypsy. The Copia/Gypsy sentences read 'integrase N-terminal to RT (INT before RT, the Copia/Ty1-type domain order) + paired LTRs → Copia superfamily' and the Gypsy equivalent, so a superfamily is no longer presented as a Wicker order. The is_erv flag and its rule (RT and LTR and ENV, and not a co-oriented DIRS call) are unchanged. |
| CLS-018 | display only | 'Insertion site captured / both termini present' wording requires a non-chance-level TSD | That sentence is emitted only when the TSD row has chance_level False. A chance-level row, or a row with no chance estimate, gets 'a target-site duplication flanks the inverted repeats (<congruence note>), but a flanking direct repeat of N bp is short enough to arise by chance (p = ...) [or: its chance probability was not assessed], so it is not credited as the insertion site and does not confirm the element ends'. The incongruent wording is unchanged apart from range/motif quoting (CLS-017). The LTR branch never claimed a captured site and has no chance clause. |
| CLS-021 | display only | Standalone-TIR evidence quotes the 22,500-sequence ceiling and adds a composition caution | The sentence quotes CHANCE_TIR_MAX_BP (17) 'across 22,500 random sequences (five seeds, checked at 35-65% GC)'. For a record outside 35-65% GC, or where the record's own composition makes an exact repeat of that length likely (probability > 1e-3), a caution is appended giving the GC and the composition-specific chance probability and asking for a flanking TSD or a transposase before relying on the call. The call, confidence and class are unchanged. |
| CLS-023 | display only | Scan-scope statements when no terminal repeat is found in a record longer than the scan windows | classify appends '(the record is N bp and terminal repeats were searched only within W bp of each end)' (TIR: 'short or imperfect terminal repeats were searched only within 400 bp' when the 1,200 bp seed window already covers the record) to the TIR-absent sentence and the two LTR-absent sentences (RT+chromodomain; RT+DDE integrase). engine.analyze adds a record note naming the LTR/TIR seed windows and the 600 bp tail window when no LTR/TIR row was found and the windows do not cover the record. |
| CLS-024 | display only | Helitron evidence no longer claims TIR absence when a TIR was detected | When a credited TIR exists the sentence reads '; a terminal inverted repeat was detected in this record, which is not a Helitron feature (it may belong to another element or be coincidental)'. The Helitron call itself (domain-only; 'Helitron-family helicase present — superfamily not called', te_class unclassified) is unchanged. |
| CLS-025 | display only | Coding-axis note counts ORFs and recognises the canonical LINE ORF1/ORF2 layout | The note reads 'N separate ORFs (frames +2, ...)' listing frames, and for two ORFs on one strand where ORF1p precedes an endonuclease-RT ORF it says 'consistent with the canonical LINE ORF1/ORF2 organisation (ORF1p, then the endonuclease–RT ORF)'. The completeness.coding.state string is unchanged. |
| CLS-027 | display only | Truncation note only for an open-ended ORF that carries a detected domain, naming the end | A note is added only when an ORF that hosts at least one detected domain is open-ended, naming the end (forward ORF open at the 3′ end, reverse-strand ORF open at the 5′ end): 'an ORF carrying a detected domain runs off the 3′ [and/or 5′] end(s) of the supplied sequence — the element may be truncated there ...'. |
| GEN-024 | display only | Family panel partitions RepeatMasker hits into TE, tandem/low-complexity, unclassified and other (RNA gene) and keys its blank-result caveats on the defaulted lineage | _partition_hits uses the backend classifier wsl.repeat_kind so the panel and the whole-genome landscape agree: TE, tandem (low-complexity and satellites), unclassified ('Unknown') and other (RNA genes, artefacts, listed apart and not counted as TE family calls); the panel header names the strongest TE class; the blank-result caveats key on species_defaulted / lineage_searched, so RepeatMasker's default lineage is no longer read as an organism choice; the family and splice controls are disarmed and the library choice cleared when the backend is lost. |
| INP-008 | display only | Fetch reports NCBI record status (replaced / replacedBy); bundled copia and L1.2 examples note their replacement accessions | resolve() adds status and replacedBy (display only, excluded from the seal); the M11240 example note says NCBI lists it as replaced by X04456 (an identical 5,146-bp sequence, compared live), the M80343 note says replaced by AH005269, and the M12927 note says its env-like ORF is not treated as retroviral-envelope evidence (the call stays Gypsy). |
| OUT-006 | display only | Provenance view, run history and export of all runs | A sectioned provenance view with units and full copyable SHA-256 values; a run history covering analysis, primer design, in-silico PCR, genome scan, family naming and splice alignment; an export of all runs. |
| OUT-009 | display only | Coordinate display toggle and convention labels | A display toggle (0-based half-open or 1-based inclusive) that never changes stored, exported or design coordinates; convention labels on coordinate columns and a coordinate_system token in exported FASTA headers. |
| UI-001 | display only | Start here strip, Help menu with quick start and glossary, Next step links | A Start here strip with one-click published examples, Next step links after each result, and a Help menu (F1) with a quick start and a glossary. |
| UI-002 | display only | Scope-caution strip and confirmation above 200 kb (display only) | A caution strip for large inputs, probable multiple copies and ORFs that run off the record, and a confirmation before analysing more than 200 kb (LARGE_INPUT_BP = 200,000; threshold unchanged); the call and its confidence are unchanged. |
| UI-003 | display only | Wicker rank path, RepeatMasker-style class string and Published-as check | Beside each call: the Wicker rank path with three-letter code (X = rank not called) and a RepeatMasker-style class string (a display mapping, not another classification); a Published-as check for the bundled examples. |
| UI-004 | display only | Visible record selection for multi-record input | Click or arrow keys select a record; the bound record is named in the detail title, the specimen panel and panels 04-05. |
| UI-005 | display only | Header labels, status chip, panel 03b name, accessibility | Hide input / Show input, Text size, Dark mode, Help; the status chip separates the built-in engine from the Linux backend; panel '03b Splice structure (transcript alignment)'; WCAG AA contrast, 2-px focus ring, figure text of at least 9 px; row-action buttons beside context menus; keyboard access to canvases. |
| UI-006 | display only | Fetch Cancel, close confirmation while jobs run, Show all N ORFs | A Cancel button for a running fetch (its late result is ignored); a 'Close TEagle?' confirmation while a job runs; the ORF table lists the ORFs searched for domains with a 'Show all N ORFs' link. |
| UI-007 | display only | Gel lane marks and legend | No 'no on-target' label for undesigned pairs; failed and capped lanes are marked; the legend explains every band colour on a mixed gel. |
| GEN-013 | infrastructure | Whole-genome annotation sizes disk from the source FASTA, checks space before the download, numbers resumed progress and marks completion | genome_prepare records fasta_bytes (genome_list returns it); annotate_budget(genome_bytes, fasta_bytes, fasta_cached) uses ceil(k x FASTA / GiB) + 2 with k = 2 (FASTA must be re-downloaded) or 1 (cached; the panel sends no fasta_cached, so it stays conservative); the script checks 2 x FASTA + 1 GiB free (smaller of WSL and Windows backing volume) before downloading (4 x the 2bit size for caches without fasta_bytes); progress starts at the chunks already done; a .complete marker makes a settings clash read 'a completed annotation ... exists with different settings' and genome_annotate returns annotation_state complete or part_finished; the cost dialog shows the error and a 'Check again' button when the budget probe fails. |
| GEN-014 | infrastructure | The backend skips container-runtime distributions and distributions without bash or python3, and remembers its distro | The distro table is parsed once; docker-desktop(-data), rancher-desktop(-data) and podman-machine-default are never candidates; candidates are tried as remembered distro (wsl_distro.txt), Windows default, then the rest, and one is chosen only if it has bash and python3 (a timeout counts as unknown); set_distro pins a choice; envcheck.probe_backends reports 'available' only when a non-runtime distro exists, else 'no distro'. |
| GEN-015 | infrastructure | WSL failures are diagnosed (still starting, no bash, start failed, restart pending) and unregister is advised only after a confirmed start failure | available() returns reason timeout \| starting \| no_bash \| start_failed (only on a WSL error code) \| unknown \| no_suitable_distro; components_status gives per-reason guidance and offers wsl --unregister only for start_failed, with a warning that it deletes the home directory, Dfam libraries and cached genomes. The Windows install log state is none \| running \| stale (no marker and older than 2 h or than the last boot) \| done \| failed; components_status returns installing False plus wsl2_installing, wsl2_install_stale and reboot_pending; a reopened dialog attaches to a running Windows install, releases it when stale, and retries a failed first probe up to 3 times before asking for Refresh. |
| GEN-016 | infrastructure | 'Install / update all' continues past a failing component and skips only its dependents | Each step is a function run in its own subshell (a failure ends only that step); a step whose prerequisite failed or was skipped is skipped (_STEP_NEEDS: repeatmasker, minimap2, miniprot, genomescan, viennarna need micromamba; famdb_conf needs repeatmasker); the run ends with one SUMMARY line and a single terminal '[teagle] FAILED: not installed: ...' (or DONE). Mid-run '[teagle] STEP x FAILED' lines are not terminal markers. |
| GEN-017 | infrastructure | The shared 'te' conda environment is adopted and removed only when TEagle created it; a different micromamba is backed up | 'te' is adopted only with ~/micromamba/envs/te/.teagle-owned or the pinned RepeatMasker from an earlier TEagle install (then marked); new environments are marked at creation; any other 'te' is left untouched and the step fails with a message. A different ~/bin/micromamba is copied to micromamba.pre-teagle first. The Inno clean install removes envs/te only by the same rule (marker, or repeatmasker-4.2.4-*.json in conda-meta, pin asserted equal to wsl._CONDA_PINS) and still removes teagle-vrna, genome caches and logs. |
| GEN-018 | infrastructure | The elevated WSL install passes its command directly (no staged script) and warns about a different approving account | ShellExecuteW 'runas' launches %SystemRoot%\System32\cmd.exe with a fixed '/v:on /d /s /c' command line that calls the system wsl.exe, appends to the Windows-side log and writes the DONE-WSL marker; the only install call is 'wsl --install -d Ubuntu --no-launch'. A standard user (no Administrators SID in whoami /groups) gets a NOTE in the log and a warning that the distro registers to the approving account. |
| GEN-019 | infrastructure | Backend robustness and status text: FamDB integrity graded on exit status and header; packaged build status; non-blocking exit; uncurated-row coverage text | The probe records FAMDB_RC=$? and integrity_check passes the check only when the exit status is 0 and parse_famdb_info finds the Database and Version header; the detail shows 'exit N:' plus the output when it fails. envcheck.check in a frozen build returns needs_install = not all packages present, first_run False, signature_matches True and bundled True (was needs_install True and first_run True with no env_state.json). engine_worker.Engine runs jobs on daemon threads (semaphore-gated) instead of a QThreadPool whose destructor waited for every job, so closing the app or the backend dialog never waits for a running job. The optional uncurated-root installer row quotes the measured per-lineage curated coverage (wsl.curated_coverage_sentence) instead of 'very few families'. |
| GEN-021 | infrastructure | Installer, launcher and build changes | MinVersion 10.0.17763 (Windows 10 1809); [InstallDelete] clears {app}\_internal before an upgrade; LICENSE, THIRD-PARTY-NOTICES.md and LGPL-3.0, GPL-3.0, HMMER and Easel texts install to {app}\licenses; clean install also deletes HKCU\Software\TEagle\TEagle and treats a missing wsl.exe or unregistered distro as nothing to clean; the launcher logs to %LOCALAPPDATA%\TEagle\logs\teagle.log when std streams are None; a second launch raises the running window or says 'TEagle is already running.'; the job also sets BREAKAWAY_OK; the version resource goes to PyInstaller's workpath; copy_metadata failures fail the build; hidden imports cover the new panel modules and drop 'server'; the backend dialog binds its deferred refreshes to the dialog; clean-install and uninstall prompts say the Linux backend and every downloaded genome are deleted. |
| INF-001 | infrastructure | Windowed app logs to %LOCALAPPDATA%\TEagle\logs | Logs go to %LOCALAPPDATA%\TEagle\logs, and the interface recovers after an engine fault. |
| INF-002 | infrastructure | Deferred installer-dialog refreshes bound to the dialog | Deferred refreshes use the dialog as context object and never run after it closes. |
| INF-003 | infrastructure | main.py split into ui_common.py and twelve panel modules (no behaviour change) | main.py (1,009 lines) plus ui_common.py and twelve panel mixins; every definition byte-identical; screenshots equivalent. |
| INF-004 | infrastructure | Test infrastructure: hermetic WSL, in-memory settings, window disposal | Tests stub the WSL distribution, use in-memory settings and delete leftover windows; the hermetic suite runs in about 5 minutes. |
| INP-007 | infrastructure | NCBI/ENA requests retry rate limits and server errors, are paced, and report failures in plain language; feature fetch failures are not cached | HTTP 429 and 5xx are retried at most twice (Retry-After honoured and capped at 10 s, else 1 s then 2 s); other 4xx, timeouts and connection failures are not retried; NCBI requests are spaced 0.34 s apart (3 per second without a key); messages explain rate limits, temporary unavailability and connectivity (cached accessions still load); the User-Agent carries the TEagle version; a cache hit with featuresStatus 'failed' re-requests features only; a coordinate cache hit keeps source.retrievedUtc as first retrieved and adds a top-level servedUtc. |
| OUT-011 | infrastructure | Legacy browser interface and --server flag removed | Removed; only the native app remains. |
| STR-020 | infrastructure | chance_tir() and the live seeded chance-TIR guard for the 20-bp standalone floor | 4.1.0: structural.chance_tir(seq, length) returns the chance, from the record's own base composition, that random sequence holds a perfect inverted repeat of that length between the two outer seed windows (p = 1 - exp(-W^2 q^length), q = 2(fA fT + fC fG)), with gc, placements, per_placement, chance_probability and in_measured_gc_range (35-65% GC). It is used by classify for a composition-aware chance caution. A live seeded measurement of find_tir on random sequence (seed 41, 450 sequences, and 1,800 sequences) runs in the test suite and fails if the longest chance TIR approaches the 20-bp floor. MIN_TIR_STANDALONE is unchanged. |

## Label space

The `classification.te_class` value set is identical in 4.0.0 and 4.1.0. The labels below were added, removed or renamed across all outputs; the JSON lists them per record and per output field.

Added:

- '(available: ...)' in the not-found error
- '; N single-primer (F+F/R+R) product(s) also predicted'
- 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'
- 'partially assessable (flanks of N bp allow only TSDs of up to N bp)'
- .teagle-owned marker
- 3' junction - element into the downstream flank
- 5' junction - upstream flank into the element
- 58 structural_* keys, ltr_extend_x_drop, polya_signal_dse_lo, polya_signal_dse_hi, polya_signal_dse_min_gt_frac, polya_signal_dse_min_t, classify_chance_tir_max_bp, classify_chance_tir_sequences, tsd_expectation_by_superfamily
- <record id>_<pair>_F/R oligo names
- Anneal F/R
- Assay goal
- Cancel (fetch)
- Copy locus, 1-based
- Copy order rows
- Custom range
- D<n>·P<m> pair labels
- Detect or quantify this element
- Domain
- Export order sheet
- FAMDB_RC
- Genotype this insertion
- LTR retrotransposon (superfamily undetermined; retroviral-type gag/env — Gypsy, Bel-Pao and Retrovirus/ERV are not separated by the bundled panel)
- Location
- MT molecules in six maps
- Notes
- Or check a primer pair you already have
- Penalty (first)
- RNase H-like (Pfam RVT_3) (domain label)
- Showing the first N of M products
- Struct (first)
- Struct (text flag)
- TEagle heuristic
- TSD row context 'flanks a poly-A tail (target-primed reverse transcription)'
- Tail F/R
- Whole insertion - flank to flank (short insertions)
- [teagle] STEP <key> FAILED / SKIPPED
- [teagle] SUMMARY installed/failed/skipped
- agree value 'none'
- alignment_strand
- ambiguous_columns (LTR row)
- ambiguous_species
- annotation_state
- assembly-accession refusal
- bound_reason 'the repeat still matched where the record ends and the two copies met each other, ...'
- bound_reason 'the two copies met each other, so the copies may extend beyond it ...'
- boundary_adjustment (LTR row)
- boundary_adjustment (TIR row)
- boundary_adjustment on a recovered TIR row
- boundary_tests (LTR row, value 36)
- boundary_tests (TIR row)
- boundary_uncertain (TIR row)
- bundled
- capped
- cautions
- chance_level (TSD row)
- chance_probability (PBS row)
- chance_probability (TSD row)
- chance_tir() return keys gc, placements, per_placement, chance_probability, in_measured_gc_range
- cis_note
- cis_orientation values '+', '-', 'undetermined'
- cis_orientation_source values 'caller', 'cis-elements', 'domain-strand'
- cis_orientation_source: domain-strand
- classification withheld (protein-domain scan unavailable)
- classify_chance_tir_max_bp
- classify_chance_tir_sequences
- coordinate_system=0-based-half-open
- criteria.annealing_core
- criteria.three_prime_strict_effective
- databases: RefSeq genome assembly
- derived=true on introns
- design note 'Designed inside the ...'
- design_locus.end, .strand, .intended_product_span, .intended_endpoint_tolerance_bp
- design_tm_conditions
- domain_architecture
- domain_architecture: none detected
- domain_architecture: unavailable
- expected_motif
- expected_range
- failed
- fasta_bytes
- featuresModel
- featuresStatus
- found
- invalidReason (record field)
- lineage_searched
- lineage_taxid
- method 'k-mer seed vs reverse-complement (k=11, 1200 bp window)'
- method value 'k-mer seed (k=13, W bp window) + banded gapped alignment'
- micromamba.pre-teagle
- n_tests (TSD row)
- observed_motif
- oligoqc_conditions
- order sheet columns Name, Sequence (5'-3'), Scale, Purification, Length (nt), Tm (C; Primer3: <conditions>), GC (%), Used in
- other non-TE annotations table (TEagle_other_annotations)
- panel (PBS row)
- parameters.max_ns_accepted
- parameters.template_ambiguity_masked
- parameters.tm_conditions
- partial (RT + chromodomain; integrase not confirmed)
- pbs_scope
- pbs_scope (LTR row)
- processing.domains.status: not_run
- product-size window statement
- protein-accession refusal
- reason
- reason: no_suitable_distro
- reboot_pending
- record note 'the record is N bp and no terminal repeat was found, but the terminal-repeat scans read only ...'
- refusal notes with the product-size arithmetic
- replacedBy
- retroelement fragment (order unassigned)
- servedUtc
- software: libprimer3 (Primer3 core)
- species_defaulted
- species_requested
- stageSummary (result field)
- status
- strand on PBS/PPT/polyA-signal rows
- strand_determined
- structural_chance_level_p
- structural_chance_tir_gc_range
- template_ambiguity
- terminal inverted repeat (class unassigned)
- termini.selected_for_motif
- tm_basis
- transcript_rank=
- transcript_strand
- tsd_corroborated (tail row)
- tsd_expectation_by_superfamily
- tsd_state
- tsd_state 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)'
- warning
- wicker_order
- wicker_order: DIRS
- wicker_order: LINE
- wicker_order: LTR
- wicker_order: TIR
- wicker_order: unassigned
- wicker_order: withheld
- window-widening warning
- wsl2_install_stale
- wsl2_installing
- your own pair (<name>)
- {app}\licenses
- {family} (advisory: strongest transposase profile; termini unassessed) with family in: hAT, Tc1/Mariner, CACTA (En/Spm), MULE (Mutator), IS4-like DDE (piggyBac group), DDE transposon, DNA transposon

Removed: --server; DNA transposon (transposase detected; termini unassessed); cds_group/phase on derived exon rows; design_locus.stop as the window end (a legacy 'stop' is still read for 'end' but a locus without span/strand is now ignored); input.source.assemblyAccession; install_wsl.bat; installer/_version_info.txt (tracked); red 'no on-target' under an untargeted lane; server module in the bundle.

Renamed: '1 on-target product; no additional priming site detected ...' -> '... no additional forward+reverse product detected ... ; N single-primer ... predicted' (when single-primer products exist); 'DNA transposon (transposase detected; termini unassessed)' -> '{family} (advisory: strongest transposase profile; termini unassessed)'; 'Dfam 4.0 (curated)' -> 'Dfam (transposable-element family database)'; 'Pfam-A profiles (via InterPro, retrieved 2026-07-17)' -> 'Pfam-A profiles (via InterPro, 30 accession versions, retrieved 2026-07-17/23/28)'; 'RepeatMasker 4.2.4 (RMBLAST engine)' -> 'RepeatMasker (RMBLAST engine)'; 'Sub-threshold terminal direct repeat (advisory)' -> 'LTR (terminal direct repeat)' for gapped-rescued pairs (row type, not a renamed string); 'absent (flanks present, ...)' -> 'not assessed (element boundary uncertain ...)' on uncorroborated TIR rows; 'absent (flanks present, no direct repeat found)' -> 'partially assessable (...)' when the smaller flank is 2-11 bp; 'no clear TE signature' (none, complete) -> 'terminal direct-repeat candidate (LTR origin unconfirmed)' (repeat/structural-only, withheld) for one cached 49.7 kb human coordinate window; 'no clear TE signature' -> 'classification withheld (protein-domain scan unavailable)' (domains_ok False); 'no genome-wide priming site (...)' -> 'no forward+reverse priming-site pair in this assembly, at the current 3' match length; N single-primer ...' (when single-primer products exist); 'no genome-wide product - the primers do not amplify anywhere in this assembly' -> 'no forward+reverse product at the intended locus or elsewhere in this assembly under current settings'; 'non-autonomous TIR element (superfamily undetermined)' -> 'terminal inverted repeat, class unassigned' (domains_ok False); '{codes} domain(s) detected; no transposable-element assignment'/unclassified -> 'non-autonomous TIR element (superfamily undetermined)'/DNA/non-autonomous (for the trigger); BadRequest 'primer parameters rejected: OSError ...' -> 'the template contains characters that are not nucleotide codes (...)' for non-IUPAC characters; Confidence (domain table) -> Hit confidence; Fwd pos -> Fwd pos (0-based) / (1-based); export 'Fwd pos (0-based, half-open)'; GC clamp -> GC clamp (bases); Gypsy (Ty3) -> RETROVIRAL_UNDETERMINED string (conditional re-label for inputs meeting the trigger; 'Gypsy (Ty3)' is still emitted for other inputs); LTR/unclassified + 'LTR retrotransposon (superfamily undetermined)' -> retro/unassigned + 'retrotransposon order unassigned' (for the trigger); Max poly-X -> Max poly-X (bases); Max size -> Max size (nt); Min size -> Min size (nt); Opt size -> Opt size (nt); Outcome -> Call status; Return -> Pairs to return; Tm max -> Tm max (°C); Tm min -> Tm min (°C); combined coordinates -> separate Start and End numeric columns; domain code RT -> RNaseH for hmm RVT_3; domain label 'reverse transcriptase' -> 'RNase H-like (Pfam RVT_3)'; evidence '17 bp over 4,500 random sequences' -> '17 bp across 22,500 random sequences (five seeds, checked at 35-65% GC)'; evidence 'insertion site itself is captured ... both element termini are present' -> chance-hedged sentence (for chance-level or unassessed TSD rows); evidence 'likely coincidental' -> 'the repeat may be coincidental, or the boundary or superfamily call may be wrong'; evidence text 'Copia/Ty1 order' -> 'Copia superfamily' (and Gypsy/Ty3 equivalent); evidence text 'ERV / errantivirus lineage' -> 'ERV architecture' (+ errantivirus clause only for Gypsy); export 'Mism F/R' -> "Mism F/R (annealed 3' core)"; input id 'record_<hash of input_sequence>' -> 'input_sequence' (bare input); amplicon source 'record_<hash> (input)' -> 'input_sequence (input)'; inputType 'local_input' -> 'fetched_accession' for unversioned and ENA-served fetches; installing (folded Windows install flag) -> separate wsl2_installing; kind 'DNA transposon' -> 'terminal inverted repeat (class unassigned)' or 'retroelement fragment (order unassigned)' (for the trigger); kind 'non-LTR retroelement (LINE-like)' -> 'retroelement (LTR not confirmed)' (for the trigger); method 'k-mer seed vs reverse-complement (k=11)' -> 'k-mer seed vs reverse-complement (k=11, 1200 bp window)'; note 'sequence appears 3′-truncated — some ORFs are open-ended ...' -> 'an ORF carrying a detected domain runs off the 3′/5′ end ...'; pair id 'P3' -> 'D2·P3' (display and exports); refusal example 'M11240, NC_003075.7' -> 'X04456, AB087615'; resultSchemaVersion '1.0' -> '1.1'; software 'primer3 (primer3-py)' -> 'primer3-py (Primer3 bindings)'; species '(all installed families)' -> 'Homo sapiens (RepeatMasker default lineage)'; te_class DNA/unclassified -> DNA/{family token} on the structure-free path; tier 'partial (RT only)' -> 'partial (RT + chromodomain; integrase not confirmed)' (for the trigger); tsd_state 'absent (flanks present, no direct repeat found)' -> 'not assessed (element boundary uncertain ...)' on uncorroborated TIR rows.

## Existing benchmark sets and what 4.1.0 does to them

### SET1 DOI-bound TEagle 4.0.0 prospective benchmark (protocol 1.0.0 to 1.4.0; results record)

FREEZE: each protocol version is a zero-results or post-closure archive with FREEZE.json (freeze_id chained through parent_freeze_id), config locks, validators run by a primary and an independent implementation, and a public-gate receipt after Zenodo publication; no stage may run until the exact archive is public, re-downloaded and hash-verified. Executed stage order: stage4 A1, stage8 original closure, stage9 A2, stage10 I2, stage11 D1-stress, stage12 P2, stage13 R1, stage14 curation (DEV14/15), stage16 domain completion (DEV16/17), stage17 reporting. Results freeze: reporting/REPORT_LOCK.json (12 source tables), RESULTS_FREEZE.json, ARCHIVE_MANIFEST, PACKAGE_ALLOWLIST, validators/validate_primary.py and validate_independent.py, deterministic rebuild; 152-member archive. Analysis families are never pooled: ORIGINAL_1_0_0, AMENDED_A1, AMENDED_A2, AMENDED_I2_INTEGRATION, AMENDED_D1_STRESS, AMENDED_P2_SYNTHETIC_SPECIFICITY, AMENDED_R1_DEPLOYMENT, AMENDED_FUNCTION_COMPLETION (deviations DEV-01..DEV-18). Of 10 original endpoints 9 are not_run and 1 (D1 admitted subset) is incomplete non-primary. IMPORTANT LIMIT: the stage harness code that actually invoked TEagle 4.0.0 and parsed its outputs for A1/A2/I2/D1-stress/P2/R1 is NOT in the results package, NOT in the protocol archives, and could not be located on disk (searched Desktop, Documents, Downloads, D:, local_archive, outdated; no stage* directory). TEagle fields read below are therefore taken from the protocol contracts (parameter_locks.json teagle_output_rules, a1_crosswalk_lock + apply_a1_crosswalk.py, r1_lock, i2 bundle schema, P2 templates) and are labelled 'protocol contract, adapter not inspected' where applicable. Native TEagle condition per INTEGRATION_COMPARATOR.md: operator uses documented GUI paths and native exports; scripts may validate but not repair outputs. Version binding: TEagle 4.0.0 only. Protocol 1.0.0 binds the v4.0.0 release by release commit, installer name/size and installer SHA-256 (config/benchmark.json) ; R1 lock re-asserts the installer hash and installed version 4.0.0; the A1 crosswalk lock binds by SHA-256 to the frozen 1.0.0 parameter lock and to the v4.0.0 classify.py source; D1/DEV16 bind the bundled te_domains.hmm by SHA-256 (its bytes are IDENTICAL in 4.1.0, checked by git blob comparison). Protocol DOIs: concept 10.5281/zenodo.22129307; 1.0.0 .22129308; 1.1.0 .22141494; 1.1.1 .22146620; 1.2.0 .22150368; 1.3.0 .22153554; 1.4.0 .22154316. Results DOI 10.5281/zenodo.22159977 (concept .22159976).

| Endpoint | Question | Units | TEagle outputs read |
|---|---|---|---|
| A1_classification (AMENDED_A1, DEV-01/02/07/12) | Selective classification: how often does TEagle (vs TEsorter 1.5.1, RepeatClassifier 2.0.9, PASTEClassifier 2.0) commit to class/order/superfamily labels that agree with curated-library header labels, and how often does it withhold, fail or fall outside the label space? Descriptive, post-raw/pre-score crosswalk, no superiority claim. | 180 records x 4 methods x 3 ranks = 2160 rank rows (TEagle 540). Sources: arabidopsis_edta_2_2_2 55, maize_mtec_2019_11_12 40, rice_edta_7_0_0 85. Reported TEagle 4.0.0 overall: class and order each 135 answered (77 correct, 58 incorrect), 18 withheld, 26 unsupported, 1 unresolved, 0 failed; superfamily 130 truth-scored (50 unscored), 30 answered/30 correct, 16 withheld, 19 unsupported, 65 unresolved. All strata and sources tabled in reporting/classification_metrics.tsv (312 rows, DEV-12). All inputs in stratum NO_NEAR_TARGET_RESOURCE_SIMILARITY. | Harness-level (protocol contract, adapter not inspected; shape taken from toy/a1_crosswalk_input.tsv and apply_a1_crosswalk.py): raw_execution_state in {completed, withheld, failed, not_run} derived from classification.outcome == 'withheld' (outcome overrides labels), process/parse failure, or not attempted; classification.te_class (stored as raw_class; the ONLY TEagle field the crosswalk maps: exact_te_class table + completed_unclassified_label 'unclassified'); classification.order (stored raw_order; domain-architecture string such as 'INT-RT', NOT read by the TEagle mapping); classification.superfamily (stored raw_superfamily; NOT read by the TEagle mapping); classification.class is not read |
| A2_structure (AMENDED_A2, DEV-03/13; panel protocol 1.1.1) | Within-record terminal direct-repeat (DTR) and terminal inverted-repeat (ITR) state retrieval and interval localization on controlled synthetic constructs with exact truth. | 48 records x 2 feature types. Reported: direct_terminal_repeat truth_positive 6/truth_negative 42, TP 6, FP 0, sensitivity 1.0 (Wilson 0.610-1.0), exact interval matches 6, endpoint MAE 0 bp; inverted_terminal_repeat positive 6, TP 6, FP 0, exact interval matches 5, endpoint MAE 0.083 bp, length MAE 0.167, IoU mean 0.997. Cluster = source template. | Protocol contract only, adapter not inspected: record-level state per feature type from machine-readable structural events (record.structural[] entries whose type starts 'LTR (terminal direct repeat)' for direct and 'TIR (terminal inverted repeat)' for inverted) and, for localization, their five_prime and three_prime [start,end) intervals; S1 common_task.state_contract (positive/negative/withheld/failed/unsupported/not_run). Inference: the 4.0.0 event keys are type, five_prime, three_prime, identity, ltr_len/tir_len, tsd_state. |
| I2_integration (AMENDED_I2_INTEGRATION, DEV-08) | Do workflow bundles detect eight injected defect classes (and accept pristine bundles), for TEagle vs three manual chains? Operational integrity only, not biological accuracy. | 48 base cases x 4 conditions (teagle, manual_chain_a/b/c) = 192 bundles; reported per condition 48/48 pristine accepted, 48/48 defects detected, 6 collapsed-unresolved and 6 removed-identity detected; 36 rows in reporting/integration_metrics.tsv. | Schema i2-canonical-task-bundle-1 (protocol contract): input.sequence_sha256; operation_states{classification,domains,primers,specificity,structure}; artifacts.classification{artifact_sha256,case_id}, artifacts.domains.artifact_sha256, artifacts.primers{artifact_sha256,table_sha256}, artifacts.specificity{artifact_sha256,cap,reported_hit_count,search_complete,truncated}, artifacts.structure{artifact_sha256,direct_intervals,inverted_intervals}; provenance.tools[], provenance.databases[]; bundle_sha256. Mapped from TEagle native exports (result JSON resultPayloadSha256, provenance manifest tool/database identities); producing adapter not inspected. |
| D1_stress amended (AMENDED_D1_STRESS, DEV-09) | Technical, model-derived domain recovery: does each tool assign the registered model on shuffled (negative) and 69% fragment (positive) constructs derived from each bundled model's HMMER consensus? | 30 models x 3 roles = 90 cases per condition; 60 scored (30 shuffled + 30 fragment), 30 full controls; TEagle TP 30, TN 30, FP 0, FN 0, 30/30 full controls pass; target_hmmer_assigned_call_disagreement_count 0; rpsblast_cdd unsupported on all 90 (PSSMs absent from frozen CDD archive); normalized_row_count 360. | Protocol contract, adapter not inspected: per case, which bundled model was assigned and its significance; per parameter_locks.json domain_thresholds the TEagle domain engine is pyhmmer on the 30-model te_domains.hmm with acceptance field independent_domain_evalue < 0.001 (TEagle record.domains[].evalue = domain i-Evalue), reporting filter full_sequence_evalue <= 0.001, ORF >= 40 aa, longest 12 ORFs, overlapping same-code same-strand hits deduplicated by highest bit score, opposite strands kept. Model identity = record.domains[].hmm (the model NAME); record.domains[].domain is a coarser TEagle code (RT, INT, TPase...). |
| P2_primer_specificity (AMENDED_P2_SYNTHETIC_SPECIFICITY, DEV-10) | Controlled in-silico product reconstruction on synthetic backgrounds with fixed primer pairs, for TEagle vs Primer-BLAST dated snapshot, MFEprimer 3.2.7 and UCSC isPcr. | 96 cases x 4 conditions; TEagle 4.0.0 exact recall 1.0 and precision 1.0 in every archetype (96 cases, 0 withheld/failed/unsupported); Primer-BLAST reverse_orientation unsupported. | Protocol contract (templates/p2_case_states.tsv and p2_product_calls.tsv), adapter not inspected: per-case operation_state, search_complete, truncated, tool_version, raw_output hash; per product call start_0_based, end_0_based_exclusive, strand, raw_call_label. Inference for TEagle: engine.run_pcr / teagle_core.primers.in_silico_pcr products (start, end, product_length, single_primer, on_target, fwd_mm, rev_mm). |
| R1_deployment amended (AMENDED_R1_DEPLOYMENT, DEV-11) | Is the installed GUI deployment and every native export semantically identical across four clean Windows guests? | 4 guests x 12 cases = 48 observations x 10 formats = 480 export observations; 48 verified result documents and manifests; 144 cross-guest identity sets; all passed. | resultPayloadSha256 and manifestSha256 (primary identities); manifest excluded fields: top-level createdUtc and environment; input.retrievedUtc, id, assemblyName, displayLocus, chromName, source, endpoint, sourceUrl; 10 native exports per case: CSV, TSV, XLSX (tables compared by headers, row order, values, missing states, units, notes after OOXML metadata removal), FASTA (identifiers and sequences), GFF3 and BED (records, coordinates, strands, attributes, caveats), SVG (parsed XML), PNG (decoded pixels, dimensions, mode), PDF (page count, geometry, rendered pixels), JSON manifest; installed-version identity 4.0.0, installer SHA-256, uninstall state (declared_residual_present false) |
| DEV14_tetrimmer_curation / DEV15_mchelper_curation (AMENDED_FUNCTION_COMPLETION, protocol 1.3.0) | Genome-aware consensus curation by TEtrimmer 1.7.2 and MCHelper 1.7.1 on a synthetic genome; family recovery, redundancy, absent-candidate retention. | 4 families per tool; TEtrimmer 4/4 recovered; MCHelper NR primary 3/4 and redundant R secondary 3/4 (byte-identical, never pooled); 2 absent negatives, 0 retained; preserved failures: 24 TE-Aid, 24 BLAST DB errors, final classification and proof clustering (TEtrimmer); first fatal MCHelper attempt. | none - comparator-only (TEagle is not run; no TEagle output is read) |
| DEV16_interproscan_domains / DEV17_rpsblast_domains (protocol 1.4.0 controller, scorer 1.3.0) | Source-model self-recovery of planted domain segments by InterProScan 5.78-109.0 (Pfam 38.2) and RPS-BLAST (CDD_NCBI.v3.21). | 30 cases per tool; both 30/30 planted segments recovered (positive hit requires registered id, E <= 0.001, >= 80% overlap); 2 RPS positives reached the 500-hit cap (censored); shuffles descriptive, 0 hits. | none - comparator-only. The only TEagle artefact consumed is the bundled te_domains.hmm (30 model NAMEs and Pfam accessions via the crosswalk); TEagle result fields are not read. The 4.1.0 HMM bytes are identical. |
| D1_admitted_subset (ORIGINAL_1_0_0, incomplete non-primary) | Per-model positive detection and biological near-neighbour false detection for TEagle, HMMER, InterProScan on a curated-biological corpus (13 admitted models). | 13 models x 8 = 104 constructs scored (completed, gate PASS_INCOMPLETE_NON_PRIMARY); full 30-model macro endpoint (240) not_run; registered RPS-BLAST assigned-model endpoint unsupported. | record.domains[] model identity (hmm, or the crosswalk teagle_model_id) and evalue per parameter_locks domain_thresholds (independent_domain_evalue < 0.001; ORF >= 40 aa; max 12 ORFs); crosswalk column ontology_target (TEagle code;Pfam accession) per model - see label_maps |
| Original endpoints not_run: I1_integration, E1_classification, S1_structure_synthetic, S1_structure_biological, D1_domain_full_biological, D1_domain_stress, P1_design, P1_specificity, R1_deployment_consistency | Registered 1.0.0 primary/secondary endpoints: I1 four co-primary integration estimands (48 base cases, 96 tasks; TEagle minus mean of 3 manual chains, bootstrap seed 410005); E1 180-record untouched external classification panel (6 strata x 30, >= 12 organism clusters); S1 two 48-case terminal-repeat panels (synthetic and biological, all-or-none comparator admission); D1 full 240-construct panel; D1 stress; P1 72 real loci design (Primer3, Primer-BLAST, PrimerServer2) and specificity; R1 four-VM deployment. | planned I1 48, E1 180, S1 48+48, D1 240, D1 stress 60, P1 72+72, R1 4 guests; observed 0 (D1 admitted subset 104 reported separately). | none observed - endpoint not run. Contracts that a future run must read: parameter_locks.json teagle_output_rules/domain_thresholds/primer_thermodynamics (Primer3 parameters, common acceptance fields forward_sequence, reverse_sequence, forward_tm_c, reverse_tm_c, forward_gc_percent, reverse_gc_percent, product_size_bp) and schemas/required_evidence_fields.json (31 evidence fields, 17 critical) |

### SET2 Historical retrospective benchmark (benchmarks/ and reviewer_package/)

Retrospective, literature-labelled, descriptive evaluation; not a holdout or population estimate. FREEZE MECHANISM: none in the preregistration sense - corpus TSVs tracked in git; run ledger raw/teagle/_run.json records corpus path, corpus text SHA-256, row count, analysed/failed counts, TEagle version, panel SHA-256 and detector parameters; score.py audits raw-case universe against the corpus (unique key accession+coords+expected_superfamily, missing/extra/duplicate raw files abort, skipped_existing must be 0); record_scope is a curated column and nothing defaults; make_provenance_index recomputes the manifest seal for every case and fails closed. Primary unit = distinct analysed input_sha256. Spec exclusion of development examples: M11240, M12927, M80343, X01005, X05424. Sections 5.2 (structural boundaries), 5.4 (domain architecture precision/recall) and 6 (robustness) of the spec are planned panels with no executing script in this folder. Version binding: Results scored on TEagle 3.8.0 (raw teagle_version fields in assay_scores, sim_divergence, reproducibility, chance_tir, diff_anchored_sites; reviewer_package README; Amendment 3 pin); spec written against 3.4.0 with pin amendments 3.5.0, 3.7.0, 3.8.0; deposit DOI 10.5281/zenodo.21848814 (immutable v3.8.0 evidence). Scorer score.py is the frozen 4.0.0 file. TEsorter 1.5.1 comparator. Not preregistered (spec explicitly makes no timing claim). Only the 43-record holdout was re-run as a set after the X01005 correction; per-case raw outputs (about 3.3 GB) are local and mostly not retained (only benchmarks/raw/teagle_holdout exists in the main folder), so regeneration needs refetching.

| Endpoint | Question | Units | TEagle outputs read |
|---|---|---|---|
| C1_broader_corpus_classification (corpus.tsv) | Class, order and superfamily answered-case accuracy and abstention of TEagle on 135 curated literature cases across 8 panels. | 129 analysed (6 failed fetches recorded in _run.json); 95 scored cases from 94 distinct inputs; 34 stratified out as whole containing records; class/order 94 gradable rows (93 unique inputs), 67 answered, 67 correct, abstention 0.2872; superfamily 59 gradable inputs (36 not gradable: 20 Retrovirus/ERV, 16 no label), 33 answered, 31 correct (0.9394, Wilson 0.8039-0.9832), abstention 0.4407; negative controls 11 with 0 false positives. | result.records[0].classification.te_class; classification.class; classification.superfamily; classification.confidence (High/Moderate/Candidate); classification.completeness (dict -> tier, or string); classification.orfs_unscanned; classification.n_domains; result.records[0].seq (reanalyse.py, run_tesorter.py read the stored analysed sequence); raw record fields corpus_row, accession, input_sha256, input_length, fetch_source{accession,organism,title,length,coords_applied} |
| C2_holdout_corpus_classification (corpus_holdout.tsv) | Same classification endpoints on the development-example-excluded retrospective literature corpus. | 43 scored from 43 distinct inputs; class/order 34 answered/34 correct of 43 (abstention 0.2093, Wilson 0.8985-1.0); superfamily 29 gradable (14 not gradable: 6 ERV, 8 no label), 19 answered/19 correct (abstention 0.3448); negatives 8. | same as C1 (classification.te_class, class, superfamily, confidence, completeness, n_domains, orfs_unscanned) |
| C3_stratum_containing_record_and_negative_controls_and_calibration | What does TEagle do on whole containing records, on negative controls, and does confidence tier track correctness (monotone gate in spec section 6)? | by_confidence_tier keys Candidate, High, Moderate; by_panel 8 panels; spec Amendment 2 records that the monotonicity gate cannot be evaluated at order rank (no errors) and fails at superfamily rank (both errors High); no tier was refitted. | classification.confidence; classification.te_class; classification.superfamily; classification.completeness.tier |
| C4_head_to_head_TEsorter | Paired order-level comparison of TEagle with TEsorter 1.5.1 (default rexdb pass and lineage-matched pass) on byte-identical input; abstention cost. | 82 paired per mode (default, matched); abstention_cost_teagle_silent_tesorter_right 9; McNemar exact on discordant pairs; historical TEsorter run had 136 corpus rows vs 135 current (mismatch flagged in output). | via score.py observed(): classification.te_class, classification.class; classification.confidence; abstained_at_order; raw records[0].seq residue hash vs the FASTA given to TEsorter; raw input_sha256 |
| C5_assay_validation (assay_corpus.json) | Does the in-silico PCR engine reproduce published amplicon sizes? | 17 analysed, 17 found a product; 16 single-locus cases all concordant (concordance 1.0), 1 multi-locus reported apart; tolerance max(20 bp, 10%); blind vs preworked and stated vs derived strata; min_anneal sensitivity with bound verification. | engine.run_pcr({'sequence','fwd','rev','params':{'prod_min','prod_max'}}).amplicons[] with keys product_length (fallback length, size), single_primer, fwd_tail5, rev_tail5, fwd_anneal_len, rev_anneal_len, fwd_mm, rev_mm; primers.in_silico_pcr(fwd, rev, seq, seq_id, prod_min, prod_max, min_anneal=); primers.MIN_PRIMER_SIZE (module constant) |
| C6_divergence_gradient (sim_divergence.py) | Where does LTR-LTR divergence make TEagle stop detecting a terminal repeat pair (acceptance floor MIN_LTR_IDENTITY 80, seeding limit)? | 600 records | structural.find_ltr(seq) return keys five_prime, three_prime, identity; structural._identity(r5, r3); structural.MIN_LTR_IDENTITY; engine.run_analyze records[0].structural[] entries: type startswith 'LTR', advisory flag, identity, ltr_len; classification.te_class, superfamily, confidence |
| C7_chance_TIR_and_matcher_equivalence | (a) Does MIN_TIR_STANDALONE (20 bp) clear chance inverted repeats? (b) Is the single-pass 3'-anchored primer matcher identical to the descending-scan semantics? | (a) 53 hits, longest 15 bp for the first seed, ceilings 15-17 over 5 seeds (22,500 sequences), margin 5 bp; (b) 6000/6000 agreed, 0 disagreements. | structural.find_tir(seq) return keys tir_len, identity; find_tir defaults k, min_tir, max_tir, min_anchors (inspect.signature); structural.MIN_TIR_STANDALONE, MIN_TIR_IDENTITY; primers.anchored_sites(primer, seq, max_mm, tp, min_anneal, on_reverse) site fields position, mismatches, positions, tail length, core length; primers.MAX_SITES, MIN_PRIMER_SIZE; sequtil.reverse_complement |
| C8_provenance_reproducibility (reproducibility.py) | Stability, sensitivity and fidelity of the provenance seal, including the assay path. | 3 inputs x 10 replicates; 23 sealed parameters, 17 detector thresholds enumerated across 7 detectors; reported stable, all distinct, coverage complete, assay coverage complete. | engine.run_analyze(...)['provenance'] keys createdUtc (volatile, excluded), teagleVersion, input.sha256, parameters (dict), software[].name, databases[].name/sha256, references[].doi, top-level keys; engine._detector_parameters() (live derivation via inspect.signature); defaults of sequtil.find_orfs(min_aa), domains.scan_domains(evalue,max_orfs), structural.find_ltr(k,min_ltr,min_anchors), find_tir(k,min_tir,max_tir,min_anchors), find_polya(min_run), find_pbs(search,min_ident), find_ppt(window,min_len,min_purine,max_defects); engine.run_pcr(...)['provenance']['parameters'] vs numeric defaults of primers.in_silico_pcr |
| C9_provenance_index_and_testsuite (make_provenance_index.py, testsuite.py) | Compact sequence-free per-case audit extract with recomputed seals; recorded test-suite counts. | 172 cases, 172/172 seals verified. | provenance: schemaVersion, runType, teagleVersion, input.sha256/length, parameters (flat numeric scalars only, <= 64 entries, keys snake_case), software[{name,version}], databases[], internal manifest seal (top-level keys manifestSha256/createdUtc/environment and input keys retrievedUtc,id,assemblyName,displayLocus,chromName,source,endpoint,sourceUrl excluded); records[0].classification: class, order, superfamily, confidence, completeness; records[0].structural: type 'LTR (terminal direct repeat)' (ltr_len, identity, length_is_lower_bound, threshold, tsd_state), 'TIR' (tir_len), 'TSD (target-site duplication)' (length, matched_expected, tsd_congruence) |

### SET3 Prospective framework benchmark_publication_v2 (benchmark 2.0.0, protocol 2.0.0-dev.2)

Four evidence levels never pooled: A implementation, B controlled synthetic, C retrospective biological, D prospective external biological. Freeze = clean-checkout create_freeze() hashing 26 defining files plus globs (benchpub/*.py, comparators/*.py, cache/sequences/*.fasta, scripts/*.py, runners/*.py, tests/test_phase*.py), corpus manifest hashes, comparator registry, seeds (synthetic_structure 3812001, robustness 3812002, off_target 3812003, bootstrap 3812004), then a preregistration Zenodo archive before any run; results record only after Windows 10 and Windows 11 clean-machine evidence and a PASS audit. Corpus state: holdout.tsv 109 gradable Level D rows (all curation_status approved, development_exposure none, record_scope element; class CLASS_I 75/CLASS_II 34; order LTR 66, HELITRON 20, TIR 14, LINE 9; superfamily COPIA 16, GYPSY 31, RETROVIRUS 19, TC1_MARINER 12, blank 31; taxon_group plant 47, animal 41, other 14, fungal 7; element_condition fragment 70, other 32, intact 5, internally_deleted 2; 51 organisms; 109 unique sequence hashes; E01 removed 3 duplicates of one 4-accession group; 75 exclusions logged in exclusions.tsv; discovery 184 candidates); negative_controls.tsv 102 rows (100 Level D from 85 accessions + 2 Level B generated); structural_boundaries 1 curated case (STR001, 2 LTR rows); domain_truth 24 rows (12 models x 1 PDB positive + 1 PDB negative); erv_gypsy_challenge 2 rows (ERV001 RETROVIRUS, GYPSY001); de_novo_primer_loci 3; retrospective_assays, wetlab_validation and clean_machine_evidence empty. Holdout breadth gates: >= 3 taxon groups, >= 10 organisms, plant/fungal/animal, orders LTR/LINE/TIR, superfamilies COPIA/GYPSY/RETROVIRUS, element conditions intact/fragment/internally_deleted. TEagle is called IN-PROCESS through engine.analyze (not the GUI) with network disabled. Version binding: TEagle v3.8.1 (release tag and release commit in config/publication.yaml; the tag resolves in the repository). STATUS: RECORD_FREEZE.json = UNFROZEN; the framework has NEVER executed a scientific stage; no derived_results/raw_results from a real run. require_frozen() refuses to run unless (i) RECORD_FREEZE status FROZEN with valid freeze_id, (ii) the frozen git commit is an ancestor of HEAD, (iii) git diff <v3.8.1 release commit> -- app is EMPTY (so ANY 4.0.0 or 4.1.0 app tree is refused), (iv) Python version/executable/OS/dependency versions (primer3-py, pyhmmer, PySide6, openpyxl, certifi, pypdf) and bundled HMM hash match, (v) every defining file/glob hash matches. It therefore cannot be pointed at 4.1.0 without a new benchmark version.

| Endpoint | Question | Units | TEagle outputs read |
|---|---|---|---|
| V2-1 holdout classification (classify stage) | Class, order and superfamily answered-case accuracy, coverage, abstention, committed-error rate and correct-call yield on prospective Level D sequences; per-order table; risk-coverage. | 109 gradable unique sequence hashes; ranks scored only where an expected label is supplied (superfamily 78 of 109 rows have labels). | engine.analyze(text, source={'accession','sha256'}) -> result.records[0].classification (stages._single_record: records list must have exactly 1 record); classification['class'] (substring 'class ii' -> CLASS_II, 'class i' -> CLASS_I, else None); classification['te_class'] head token before '/' (LTR/LINE/SINE -> same; DNA/TIR/MITE -> TIR; DIRS -> DIRS; otherwise None); classification['superfamily'] exact string looked up in label_space.teagle_superfamily_map (unmapped string -> None = withheld, silently); evidence tier: first present of classification['evidence_tier','support_tier','confidence','tier'] (4.x supplies 'confidence'); result['provenance']['manifestSha256']; status field of the run wrapper (ok/failure) |
| V2-2 biological negatives | False committed class/order/superfamily rate on 100 Level D biological negatives and detected evidence that did not propagate to a call; 2 Level B generated controls reported separately. | 102 rows (100 D + 2 B), cache FASTA per row. | same classification fields as V2-1 (expected None; any non-null call is a committed false positive); records[*].structural and records[*].domains non-empty with no call -> counted in detected_evidence_without_classification |
| V2-3 ERV / Gypsy challenge (reference-space analysis) | How does TEagle behave on endogenous-retrovirus and Gypsy challenge sequences within TEagle's representable space versus the full Wicker reference taxonomy (RETROVIRUS is not representable)? | 2 challenge rows (+19 retrovirus holdout rows in V2-1) | classification['superfamily'] -> teagle_superfamily_map (a committed Gypsy call on an ERV is an error in the full-taxonomy view; None = appropriate withholding); labels.reference_space_outcome(expected, observed): representable, correct_within_representable_space, correct_against_full_taxonomy, appropriate_withholding_for_unrepresentable, committed_error_for_unrepresentable |
| V2-4 structural truth (synthetic Level B + curated Level D boundaries) | Interval detection, boundary errors, length error, IoU and base-pair precision/recall/F1 for LTR, TIR, TSD, PBS, PPT, poly-A. | 259 synthetic cases (87 LTR, 86 TIR, 86 LINE) + 1 curated record; each case scored per truth interval. | records[0].structural[] (stages._structural_events); event['type'] prefix: 'LTR' -> feature LTR using five_prime [s,e) and three_prime; 'TIR' -> TIR likewise; 'TSD' -> upstream (copy left) and downstream (copy right); 'PBS' with pos; 'PPT' with pos; 'poly-A' with pos (a 'polyA-signal' event is NOT matched by this prefix); intervals assumed 0-based half-open |
| V2-5 domains (12 of 30 models, 24 constructs) | Construct-level domain recovery and false detection for 12 prespecified models; 18 models explicitly unscored. | 24 constructs (12 positive, 12 negative); descriptive only. | records[*].domains[] items: model id = first non-empty of hit['code'], hit['name'], hit['hmm'] - in 4.0.0 and 4.1.0 only 'hmm' (HMM NAME, e.g. RVT_1) is populated, 'code' is absent (the engine's own coarse code is the key 'domain'); hit['evalue'] (i-Evalue); HMM names from app/backend/data/te_domains.hmm NAME lines (parse_hmm_names) |
| V2-6 robustness curves | Order-level call rate and correctness of TEagle under controlled perturbations. | 259 observations; linked to source and condition, not independent. | classification['te_class'] head token -> order via teagle_rank_calls; evidence tier field |
| V2-7 comparator (TEsorter 1.5.1) paired analysis | Paired rank-level correctness of TEagle vs the sole frozen comparator on identical inputs. | <= 109 paired cases x 3 fair ranks. | same TEagle classification fields as V2-1 (teagle_rank_calls); TEagle record wrapper input_sha256, runtime_seconds, teagle_result_sha256 for assert_input_identity |
| V2-8 evidence-layer ablations | How much does each evidence channel contribute to the production classifier (full, structure_only, domain_only, forced_answer, exclusion_control)? | <= 109 cases x 5 modes x 3 ranks. | records[0] keys structural, domains, orfs, orfs_unscanned, classification; teagle_core.classify.classify(structural, domains, seq=, domains_ok=True, orfs_unscanned=, orfs=) - signature unchanged in 4.1.0; classification keys class, te_class, order (domain-architecture string), superfamily; withheld tokens in superfamily text: undetermined, unassigned, partial, fragment, no clear, not called, no coding; superfamily 'LINE (non-LTR)' counted as no call |
| V2-9 de novo primer design and controlled off-target | Design outcome, candidate on-target inclusion, product accounting; recovery of planted exact, paralog, 1-3 mismatch, tailed, degenerate, reverse-strand and single-primer products. | 3 loci; 11 controlled cases | teagle_core.primers.design_primers(template, params, target=[start,length]) -> keys candidates[] (left_seq, right_seq, left_tm, right_tm), explain_left, explain_right, explain_pair; teagle_core.primers.in_silico_pcr(left, right, background, background_id, target_span=, max_mm=, tp=, min_anneal=, prod_min=, prod_max=) - NOTE 'target_span' is not a parameter in 4.0.0 or 4.1.0 (they expose primer3_target and intended_product_span), so this call would raise TypeError on any app tree that has the 4.x signature; product keys start, end, product_length\|length, on_target, single_primer, fwd_mm, rev_mm, fwd_anneal_len, rev_anneal_len |
| V2-10 runtime performance | Operation cost vs input length after warm-up. | 56 operation x length cells x 5 trials | runners/operation_worker.py calls sequtil.parse_fasta, structural.detect_all(sequence), sequtil.find_orfs, domains.scan_domains, classify.classify(structures, domain_hits, seq=, domains_ok=True, orfs_unscanned=0, orfs=), engine.analyze(fasta), primers.design_primers(sequence, {'num_return':5}), primers.in_silico_pcr(left, right, sequence, prod_min=40, prod_max=...) |
| V2-11 reproducibility, invariance and export parity | Is the result invariant to FASTA presentation and repetition; do CSV/TSV/XLSX/GFF3/BED round-trip the feature rows; is clean-machine evidence consistent? | 1 canonical case; 8 summary checks | engine.analyze(...) full result with non-semantic keys excluded from semantic hash: generatedAt, generated_at, timestamp, retrievedAt, retrieved_at, endpoint, sourceUrl plus {id, sequence_id}; result['provenance']['manifestSha256']; records[0].classification['te_class']; teagle_core.gff3.build_features(record, 'benchmark') unpacked as 9-tuples (seqid, source, type, start, end, score, strand, phase, attributes) - arity unchanged in 4.1.0; gff3.to_gff3(record, 'benchmark', sequence), gff3.to_bed(record, 'benchmark'); app/native widgets.write_table(headers, values, path) for csv/tsv/xlsx and teagle_core.artifact_metadata.ooxml_metadata_violations(path) |
| V2-12 statistics, tables, figures, report, audit | Rebuild every table and figure from indexed row-level results; recompute denominators; verify freeze identity and paired-input identity. | all stages | no new TEagle fields; consumes stage artifacts |

### Breakpoints for a renewed benchmark

ozen score.py raises ValueError on any superfamily string outside its table - already true for several 4.0.0 labels and for all three new 4.1.0 labels; it cannot score a 4.1.0 run that touches ERVs or structure-free coding calls; v2 teagle_rank_calls instead silently maps unmapped strings to withheld (different failure mode, should be a deliberate choice); (ii) make_provenance_index.py fails closed on 4.1.0 (parameter count/shape, two new tsd_state strings); (iii) v2 primers_stage calls in_silico_pcr(..., target_span=) which does not exist in 4.0.0 or 4.1.0, and v2 require_frozen rejects any app tree that differs from v3.8.1, so SET3 is unusable on 4.x without a new benchmark version and adapter update; (iv) SET1 A1 mapping reads only te_class and the withheld outcome: ERV-type calls move from LTR/Gypsy (superfamily gypsy) to LTR/unclassified (superfamily unresolved), advisory DNA calls are still te_class DNA/<family> so they keep their canonical mapping while the superfamily text changes (the crosswalk does not read superfamily), and te_class values DNA/DDE, DNA/unclassified, LTR/partial, retro/*, repeat/structural-only remain 'unsupported' under the frozen crosswalk; (v) SET1 R1 semantic identity and SET1 I2/provenance fields will differ for every 4.1.0 output (schema 1.1, 88 sealed parameters, software labels) and the installer hash binding is to the 4.0.0 installer; (vi) domain-model identity read by v2 and the D1 harness (hmm NAME) is stable, but the code field now merges RVT_3 into RNaseH, so any metric keyed on the TEagle code or on the crosswalk ontology_target for RVT_3 must be re-keyed; (vii) reproducibility.py coverage gate checks only 17 listed detector defaults and passed 4.0.0 by that definition - it does not cover the 65 new sealed keys (not re-executed); (viii) chance_tir.py/sim_divergence.py/diff_anchored_sites.py are structurally compatible with 4.1.0 (names, defaults, constants present) but their deposited numbers describe 3.8.0 detectors and were not re-run (4.1.0 changed find_ltr/find_tir internals: gapped extension, microsatellite and low-complexity guards, seed masking).

### Drift checks run on the existing harnesses

(4.0.0 vs 4.1.0, by importing both trees in separate subprocesses with bytecode writing disabled, then replaying synthetic public constructs): (1) benchmarks/ and reviewer_package/ are byte-identical between v4.0.0 and 4.1.0; (2) te_domains.hmm blob identical; (3) engine API names and signatures unchanged for analyze, run_analyze, run_fetch, run_pcr, classify.classify, find_ltr/find_tir/find_polya/find_pbs/find_ppt, scan_domains, gff3 functions; additions: detect_all(orientation=), find_tsd(placements=), in_silico_pcr(cores=), anchored_sites(core=), design_primers result key template_ambiguity, ooxml_metadata_violations(author=); (4) te_class value set identical, superfamily vocabulary gained three strings (see label_maps); (5) resultSchemaVersion 1.0 -> 1.1, new top-level stageSummary, classification.wicker_order and domain_architecture, structural-event key additions, provenance parameters 23 -> 88 (the 23 original values unchanged), provenance software label for pyhmmer changed; (6) 18 A2-style synthetic constructs (6 direct 120 bp, 6 inverted 60 bp, 6 random negatives, 6 kb) gave identical terminal-repeat calls and intervals in both versions; (7) P2-style replay (5 archetypes x 2 orientations) gave identical product calls and keys; (8) D1-stress-style replay (all 30 consensus models, full and 69% fragment constructs back-translated) recovered the assigned model 60/60 in both versions, but off-model hits differ: RNase_H<->RVT_3 cross-hits are now suppressed (shared code RNaseH) and the DDE_3 co-hit on the rve full control disappears (new INT/TPase dedup); (9) v2 synthetic LTR case boundaries differ by 1 bp per end (4.0.0 off by one, 4.1.0 equal to truth), TSD 3 bp -> 5 bp with chance_level flag, TIR tsd_state wording changed.

### Parsers and adapters that read TEagle output

- SET1 protocol_22150368/x/scripts/apply_a1_crosswalk.py :: map_rank/truth_labels :: TEagle: classification.te_class (as raw_class) via exact_te_class + label 'unclassified'; raw_execution_state in {completed, withheld, failed, not_run}; other methods use raw_order/raw_superfamily. 4.1.0: te_class value set unchanged (every old literal still emitted; ERV-type LTR calls now carry te_class LTR/unclassified instead of LTR/Gypsy); the TEagle outcome 'withheld' path unchanged.
- SET1 (harness, not inspected) A1/A2/I2/D1-stress/P2/R1 TEagle adapters :: unavailable; contracts in parameter_locks.json, r1_lock.json, i2_bundle.schema.json, p2 templates. 4.1.0 changes visible at contract level: resultSchemaVersion 1.0 -> 1.1, new top-level stageSummary, provenance parameters 23 -> 88 (4 non-scalar), provenance software name 'HMMER (pyhmmer)' -> 'pyhmmer (HMMER 3 bindings)', classification.wicker_order and classification.domain_architecture added; resultPayloadSha256 and manifestSha256 necessarily differ for every input.
- SET2 benchmarks/run_teagle.py :: engine.run_fetch({'accession'}), engine.run_analyze({'sequence'}), engine._detector_parameters(), teagle_core.__version__, domains.DOMAIN_INFO, domains.HMM_SHA256, provenance.SCHEMA_VERSION. 4.1.0: all callables present (run_analyze(body), run_fetch(body), run_pcr(body), analyze(seq_text, source)); server.py removed but not imported; _detector_parameters now returns 88 keys.
- SET2 benchmarks/score.py :: observed(), expected_class/order, observed_superfamily() (exact-literal SF_OBSERVED, composed rules, chromovirus suffix), expected_superfamily(). Raises ValueError on any superfamily string not in its table. Verified by running it on both trees' labels: already RAISES at 4.0.0 for 'LTR retroelement (coding architecture; termini unassessed)', 'DIRS-group (coding architecture; structure unassessed)', 'DNA transposon (transposase detected; termini unassessed)', 'retrotransposon order unassigned', 'retrotransposon order unassigned (structure unavailable)', 'terminal direct-repeat candidate (LTR origin unconfirmed)', 'classification withheld (structural scan unavailable)', 'classification withheld (structural and domain scans unavailable)' (paths absent from the corpus); at 4.1.0 it additionally RAISES on the new RETROVIRAL_UNDETERMINED label (retroviral-type gag/env), the advisory '<family> (advisory: strongest transposase profile; termini unassessed)' label and 'classification withheld (protein-domain scan unavailable)'. Therefore the frozen scorer cannot score a 4.1.0 run containing any ERV, any structure-free coding call or a domain-scan failure; a 4.1 benchmark needs its own scorer.
- SET2 benchmarks/head_to_head.py :: imports score.expected_order and score.observed; reads raw/teagle/*.json (records[0].seq, input_sha256) and raw/tesorter; tesorter_order().
- SET2 benchmarks/run_assay.py :: engine.run_pcr(...).amplicons[] keys product_length/length/size, single_primer, fwd_tail5, rev_tail5, fwd_anneal_len, rev_anneal_len, fwd_mm, rev_mm; primers.in_silico_pcr(... min_anneal=). 4.1.0: P2-like replay of 5 archetypes in both orientations returned identical product calls and identical dict keys in 4.0.0 and 4.1.0; in_silico_pcr gained keyword 'cores'; anchored_sites gained 'core=None'.
- SET2 benchmarks/sim_divergence.py and chance_tir.py :: structural.find_ltr/find_tir/_identity, MIN_LTR_IDENTITY (80.0), MIN_TIR_IDENTITY (80), MIN_TIR_STANDALONE (20), event keys type/advisory/identity/ltr_len/tir_len. 4.1.0: all names and default signatures present and equal; constants unchanged; classify.CHANCE_TIR_MAX_BP=17 over 22,500 sequences now recorded as a constant (same seed-sensitivity result as the deposited chance_tir.json ceiling 15-17).
- SET2 benchmarks/reproducibility.py :: engine.run_analyze provenance, engine._detector_parameters, SEAL_KEY_FOR (17 detector defaults), inspect.signature of find_orfs, scan_domains, find_ltr, find_tir, find_polya, find_pbs, find_ppt and in_silico_pcr. 4.1.0: those signatures' defaults are unchanged and the original 23 sealed values are unchanged; 65 new sealed keys exist that this script's SEAL_KEY_FOR neither lists nor checks (it enumerates only the listed detectors). Not re-executed.
- SET2 benchmarks/make_provenance_index.py :: provenance and record fields listed in C9. FAILS CLOSED at 4.1.0: (1) 88 sealed parameters exceed MAX_PARAMETER_COUNT 64; (2) structural_chance_tir_gc_range, structural_snap_start, structural_snap_end and tsd_expectation_by_superfamily are non-scalar and rejected by the numeric-scalar rule (and one is a dict); (3) tsd_state strings 'not assessed (element boundary uncertain: the repeat arms were placed by repeat strength alone)' and 'partially assessable (flanks of N bp allow only TSDs of up to N bp)' are not in TSD_STATES (three strings); TSD congruence verdict set is unchanged; new TSD fields chance_probability, n_tests, chance_level are ignored by the extractor but present.
- SET3 benchpub/teagle_runner.py :: engine.analyze(text, source={'accession','sha256'}) via sys.path of app/backend; reads result['provenance']['manifestSha256']. 4.1.0: engine.analyze(seq_text, source=None) present and unchanged in signature.
- SET3 benchpub/labels.py teagle_rank_calls :: classification['class'], ['te_class'], ['superfamily'] (exact string) with label_space['teagle_superfamily_map']; unmapped superfamily string silently becomes None (withheld). 4.1.0: all ten mapped strings plus the chromovirus variant are still emitted; the new RETROVIRAL_UNDETERMINED, advisory '<family> (advisory...)' and 'classification withheld (protein-domain scan unavailable)' strings are unmapped, so they are scored as withholding rather than raising (asymmetry with score.py). ERV rows (expected RETROVIRUS, 19 in holdout.tsv, 1 in the challenge file) will move from a committed GYPSY/COPIA call (error in the full-taxonomy view) to withholding.
- SET3 benchpub/stages.py _structural_events + structure.teagle_intervals :: records[0].structural[].type prefixes LTR/TIR/TSD/PBS/PPT/poly-A; keys five_prime, three_prime, upstream, downstream, pos. 4.1.0: type strings and these keys unchanged; new event keys boundary_adjustment, boundary_tests, boundary_uncertain, cis_orientation, cis_orientation_source, pbs_scope, strand (PPT, polyA-signal), chance_probability, n_tests, chance_level (TSD). Boundary behaviour changed on the v2 synthetic canonical LTR case: 4.0.0 reported LTR [84,246]/[1644,1806] and a 3-bp TSD, 4.1.0 reports [85,245]/[1645,1805] (equal to the generator truth) and a 5-bp TSD [80,85]/[1805,1810] flagged chance_level; the canonical TIR case reports the same arms [79,112]/[1312,1345] (truth starts at 82) but tsd_state changes from 'absent (flanks present, no direct repeat found)' to 'not assessed (element boundary uncertain...)'.
- SET3 benchpub/domains.py teagle_domain_predictions :: record.domains[].code|name|hmm and .evalue. 4.1.0: hit keys unchanged (domain, label, class, hmm, pfam, score, evalue, confidence, orf, strand, aa, nt, dna, protein); only 'hmm' is used; te_domains.hmm bytes identical; DOMAIN_INFO code for RVT_3 changed RT -> RNaseH so _dedup_domains (keyed on code, strand, nt overlap) now merges RVT_3 and RNase_H hits, and INT/TPase same-strand cross-hits are dropped (new rule); scan_domains now RAISES RuntimeError when pyhmmer is unavailable instead of returning [].
- SET3 benchpub/ablation.py + runners/operation_worker.py :: classify.classify(structural, domains, seq=, domains_ok=, orfs_unscanned=, orfs=) (4.1.0 signature identical, still has structural_ok=True kw); structural.detect_all(sequence) (4.1.0 adds optional orientation=None); sequtil.find_orfs, domains.scan_domains, primers.design_primers (4.1.0 result adds key template_ambiguity), primers.in_silico_pcr.
- SET3 benchpub/stages.py primers_stage + assays.design_and_evaluate_locus :: in_silico_pcr(..., target_span=...) - keyword absent in BOTH 4.0.0 and 4.1.0 (signature has primer3_target, intended_product_span, intended_endpoint_tolerance_bp, assess_design_geometry, stats, min_anneal, and in 4.1.0 cores): the stage would TypeError on 4.x; intended only for 3.8.1.
- SET3 benchpub/stages.py reproducibility_stage :: gff3.build_features/to_gff3/to_bed (signatures unchanged; 9-tuple arity unchanged), widgets.write_table (present), artifact_metadata.ooxml_metadata_violations(path) (4.1.0 adds keyword author defaulting to the developer name; XLSX export now writes author=None).
- SET3 benchpub/corpus.py application_fetcher :: teagle_core.fetch.resolve(acc) and fetch.fetch_fasta(acc) (both present at 4.1.0).
- SET3 benchpub/freeze.py require_frozen/create_freeze :: reads teagle_core.__version__ by exec of app/backend/teagle_core/__init__.py; git diff <release commit> -- app must be empty.

### Label maps

- A1 canonical label space (protocol 1.2.0 a1_crosswalk_lock.json): class [class_i_retrotransposon, class_ii_dna_transposon, non_te, unresolved]; order [ltr_retrotransposon, non_ltr_retrotransposon, dirs, tir_dna_transposon, rolling_circle, non_te, unresolved]; superfamily [copia, gypsy, line, dirs, hat, tc1_mariner, cacta, mule, is4_like_dde, helitron, unresolved].
- A1 teagle exact_te_class (TEagle te_class -> class/order/superfamily): LTR/Copia -> class_i/ltr/copia; LTR/Gypsy -> class_i/ltr/gypsy; LTR/unclassified -> class_i/ltr/unresolved; LINE -> class_i/non_ltr/line; DIRS -> class_i/dirs/dirs; DNA/hAT -> class_ii/tir/hat; DNA/Tc1-Mariner -> class_ii/tir/tc1_mariner; DNA/CACTA -> class_ii/tir/cacta; DNA/MULE -> class_ii/tir/mule; DNA/IS4-like -> class_ii/tir/is4_like_dde; DNA/non-autonomous -> class_ii/tir/unresolved; none -> non_te/non_te/unresolved; 'unclassified' (completed) -> unresolved at all ranks; any other te_class (observed TEagle values not in the table: DNA/DDE, DNA/unclassified, LTR/partial, retro/partial, retro/unassigned, repeat/structural-only) -> unsupported. Outcome withheld overrides labels. In 4.1.0 an ERV-type retroviral gag/env LTR call is emitted as te_class LTR/unclassified (superfamily unresolved) rather than LTR/Gypsy.
- A1 truth_rules (library header regex -> canonical): ^biological_non_te$ -> non_te/non_te/unresolved; ^LTR/Copia$ -> copia; ^LTR/(Gypsy|Ty3)$ -> gypsy; ^LTR/ -> ltr, superfamily unresolved; ^LINE/ -> non_ltr, superfamily line; ^SINE/ -> non_ltr, unresolved; ^(DNA|DNAauto|DNAnona)/Helitron$ -> rolling_circle/helitron; /hAT$ -> hat; /CACTA$ -> cacta; /(MuDR|Mutator|MULE|MULEtir)$ -> mule; ^(DNAauto|DNAnona)/MLE$ -> tc1_mariner; ^(DNA|DNAauto|DNAnona|MITE)/ -> tir, unresolved. Corpus header->stratum: rice LTR/Gypsy|Copia|unknown -> LTR, LINE/* non-LTR, DNAauto/* (not Helitron) TIR, DNAauto|DNAnona/Helitron rolling circle, DNAnona/* (not Helitron), MITE/*, SINE/*, LTR/TRIM, LTR/Solo -> non-autonomous; maize LTR/*, LINE/*, DNA/* (not Helitron); Arabidopsis LTR/*, LINE/*, named DNA/* (not Helitron, not DNA/unknown), DNA/Helitron rolling circle.
- A1 comparator maps (for completeness; no TEagle output): TEsorter order_map LTR, LINE, SINE, TIR, Helitron, RC and superfamily_map Copia, Gypsy, hAT, Tc1-Mariner/Tc1_Mariner, CACTA, MuDR/Mutator/MuDR_Mutator/MULE, Helitron, unknown -> unresolved, absent row -> withheld; RepeatClassifier prefixes LTR/ LINE/ SINE/ DNA/ RC/ and superfamily patterns ^LTR/Copia, ^LTR/Gypsy, ^LINE/ (line), ^DNA/hAT, ^DNA/TcMar, ^DNA/CMC-EnSpm (cacta), ^DNA/MULE-MuDR, ^RC/Helitron, Unknown -> withheld; PASTEClassifier class I/II, order LTR/LARD/TRIM -> LTR, LINE/SINE -> non-LTR, DIRS, TIR/MITE -> TIR, Helitron -> rolling circle, PLE/Crypton/Maverick unsupported, Unclassified -> withheld at rank, NA unsupported, pipe-delimited ambiguous unsupported.
- D1 / DEV16 crosswalk (templates/domain_model_crosswalk.tsv, 30 rows, teagle_model_id = ontology_target): RVT_1=RT;PF00078, RVT_2=RT;PF07727, RVT_3=RT;PF13456, rve=INT;PF00665, RNase_H=RNaseH;PF00075, RVP=PR;PF00077, PEG10_N-capsid=GAG;PF03732, Gag_p24=GAG;PF00607, Gag_p24_C=GAG;PF19317, Gag_p10=GAG;PF02337, zf-CCHC_5=GAG;PF14787, HERV-K_env_2=ENV;PF13804, GP41=ENV;PF00517, TLV_coat=ENV;PF00429, Chromo=CHR;PF00385, HTH_Tnp_Tc3_2=TPase_Tc1-Mariner;PF01498, DDE_1=TPase_DDE;PF03184, DDE_3=TPase_Tc1-Mariner;PF13358, Transposase_1=TPase_Tc1-Mariner;PF01359, Dimer_Tnp_hAT=TPase_hAT;PF05699, hAT-like_RNase-H=TPase_hAT;PF14372, Transposase_22=ORF1_LINE;PF02994, Tnp_22_trimer=ORF1_LINE;PF17489, Tnp_22_dsRBD=ORF1_LINE;PF17490, Exo_endo_phos=EN_LINE_APE;PF03372, Phage_integrase=YR;PF00589, Helitron_like_N=HEL;PF14214, Transposase_24=TPase_CACTA;PF03004, MULE=TPase_MULE;PF10551, DDE_Tnp_1_7=TPase_IS4_like;PF13843. Columns also carry HMMER, InterProScan signature, RPS PSSM/CDD ids, independence_class PROGRAM_TRIO_SHARED_PFAM_LINEAGE;HMMER_ENGINE_PARITY. In 4.1.0 the TEagle code for RVT_3 (PF13456) is RNaseH (was RT) - the crosswalk ontology_target RVT_3=RT no longer matches TEagle's own code; the HMM and Pfam accessions are unchanged.
- TEagle DOMAIN_INFO code per HMM NAME at 4.1.0 (verified by import): Chromo=CHR; DDE_1, DDE_3, DDE_Tnp_1_7, Dimer_Tnp_hAT, HTH_Tnp_Tc3_2, MULE, Transposase_1, Transposase_24, hAT-like_RNase-H=TPase; Exo_endo_phos=EN; GP41, HERV-K_env_2, TLV_coat=ENV; Gag_p10, Gag_p24, Gag_p24_C, PEG10_N-capsid, zf-CCHC_5=GAG; Helitron_like_N=HEL; Phage_integrase=YR; RNase_H=RNaseH; RVP=PR; RVT_1, RVT_2=RT; RVT_3=RNaseH (4.0.0: RT); Tnp_22_dsRBD, Tnp_22_trimer, Transposase_22=ORF1; rve=INT. 4.1.0 retroviral-profile set for ERV withholding: Gag_p24, Gag_p24_C, Gag_p10, HERV-K_env_2, GP41 (TLV_coat excluded).
- SET3 config/label_space.json: classes [CLASS_I, CLASS_II]; orders [LTR, DIRS, PENELOPE, LINE, SINE, TIR, CRYPTON, HELITRON, MAVERICK]; superfamilies [COPIA, GYPSY, HAT, TC1_MARINER, CACTA, MULE, IS4_LIKE, DDE, DIRS, RETROVIRUS]; teagle_superfamily_map (10 entries, keyed on classification.superfamily): 'Copia (Ty1)'->COPIA; 'Gypsy (Ty3)'->GYPSY; 'Gypsy (Ty3) · chromovirus'->GYPSY; 'hAT'->HAT; 'Tc1/Mariner'->TC1_MARINER; 'CACTA (En/Spm)'->CACTA; 'MULE (Mutator)'->MULE; 'IS4-like DDE (piggyBac group)'->IS4_LIKE; 'DDE transposon'->DDE; 'DIRS-group (tyrosine-recombinase retroelement)'->DIRS; unrepresentable_reference_labels [RETROVIRUS, PENELOPE, CRYPTON, HELITRON, MAVERICK]. Not mapped (so withheld): 'DNA transposon', 'LINE (non-LTR)', every 'undetermined'/advisory/fragment/withheld label. Order is derived from te_class head (LTR/LINE/SINE; DNA/TIR/MITE->TIR; DIRS), class from classification.class text. Gradable holdout superfamily expectations present: COPIA, GYPSY, RETROVIRUS, TC1_MARINER.
- SET3 expected-label aliases (labels.py): class I/CLASS I -> CLASS_I, II -> CLASS_II; order DNA/MITE/TIR -> TIR, PENEOLOPE/PLE -> PENELOPE, POLINTON -> MAVERICK, others identity; superfamily TY1/COPIA -> COPIA, TY3/GYPSY -> GYPSY, TC1/MARINER -> TC1_MARINER, EN/SPM -> CACTA, MUTATOR -> MULE, IS4-LIKE -> IS4_LIKE, ENDOGENOUS RETROVIRUS -> RETROVIRUS; unknown value after normalisation raises (freeze blocker).
- SET3 comparators/tesorter.py: _ORDERS LTR, DIRS, PENELOPE|PLE, LINE, SINE, DNA|TIR -> TIR, HELITRON|RC -> HELITRON, MAVERICK|POLINTON -> MAVERICK; _SUPERFAMILIES COPIA|TY1_COPIA, GYPSY|TY3_GYPSY, HAT, TC1_MARINER, ENSPM_CACTA|EN_SPM_CACTA|CACTA, MUDR_MUTATOR|MULE|MUTATOR, DIRS, RETROVIRUS; class I orders {LTR, DIRS, PENELOPE, LINE, SINE}; pararetrovirus -> class only; unknown/mixture/total -> withheld.
- SET2 score.py SF_OBSERVED (TEagle superfamily string -> token): 'Copia (Ty1)'->COPIA; 'Gypsy (Ty3)'->GYPSY (chromovirus suffix ' · chromovirus' stripped); 'hAT'->HAT; 'Tc1/Mariner'->TCMAR; 'CACTA (En/Spm)'->CACTA; 'MULE (Mutator)'->MULE; 'IS4-like DDE (piggyBac group)'->PIGGYBAC; 'DDE transposon'->DDE_UNPLACED; 'DNA transposon'->DNA_UNPLACED; 'DIRS-group (tyrosine-recombinase retroelement)'->DIRS; ABSTAIN for 'LTR retrotransposon (superfamily undetermined)', 'LTR retrotransposon (no coding domains detected)', 'retrotransposon (partial)', 'reverse transcriptase, class unassigned', 'LINE (non-LTR)', 'non-autonomous TIR element (superfamily undetermined)', 'terminal inverted repeat, class unassigned', 'no clear TE signature'; composed ABSTAIN rules: 'LTR retroelement fragment (' ... '; RT/pol not detected)', 'coding fragment (' ... '; RT/transposase not detected)', '' ... ' domain(s) detected; no transposable-element assignment', 'Helitron-family helicase present' ... 'superfamily not called'. REPRESENTABLE tokens derived: COPIA, GYPSY, HAT, TCMAR, CACTA, MULE, PIGGYBAC, DDE_UNPLACED, DNA_UNPLACED, DIRS.
- SET2 score.py SF_EXPECTED (curated label -> token) about 60 literal keys: Copia variants (10) -> COPIA; Gypsy/Ty3 variants (9) -> GYPSY; Pao/BEL -> BELPAO; L1, CR1, Jockey (+ non-LTR LINE variant), Tad1, R1, R2, I, RTE, RTE/BovB -> own tokens; 7SL-RNA-derived Alu and B1 -> SINE_7SL; hAT; Tc1/Mariner variants and Tc1/Mariner-related MITE -> TCMAR; CACTA, CACTA (En/Spm); MULE (Mutator) variants and MULE-related MITE -> MULE; IS4-like/piggyBac -> PIGGYBAC; PIF/Harbinger -> PIFHARB; generic DDE (P-element family), P-element family -> P_ELEMENT; blank and 'MusD element, ERVL-related (superfamily left unresolved, see note)' -> not gradable (no label); 20 Retrovirus/ERV strings (HML-2, HML-2 (LTR5Hs), HERV-K/H/W, IAP, Gamma/Betaretrovirus-like ERV1/ERVK lineages, Retroviridae genera, ZFERV, PyERV) -> NG_RETROVIRUS (not gradable). Any expected label outside the table raises.
- SET2 score.py class/order: observed() as in C1; expected_class(): startswith 'not a te' or contains 'falls outside the wicker' -> NOT_TE; 'class ii' -> CLASS_II; 'class i' -> CLASS_I else OTHER; expected_order(): sine -> SINE; line or non-ltr -> LINE; ltr -> LTR; tir or class ii -> TIR else OTHER.
- SET2 head_to_head tesorter_order(): order column lower-case; empty/unknown/unclassified/mixture/na -> NO_CALL; ltr* -> LTR; line* -> LINE; sine* -> SINE; tir*/dna*/mite*/helitron*/maverick* -> TIR; else OTHER. run_tesorter organism->database: PLANT regex -> rexdb-plant; METAZOA regex -> rexdb-metazoa; unknown -> rexdb (matched pass only; default pass always rexdb).
- SET2 make_provenance_index closed vocabularies: TSD_STATES {found, 'not assessable (no flanking sequence in the record)', 'absent (flanks present, no direct repeat found)'}; TSD_CONGRUENCE {congruent, incongruent, not assessed, not diagnostic, not recorded}; PARAMETER_KEY_RE snake_case, parameters numeric scalar, max 64.
- TEagle te_class vocabulary actually emitted (classify.py, derived by static parse plus runtime): 4.0.0 and 4.1.0 identical: DIRS, DNA/non-autonomous, LINE, LTR/Copia, LTR/Gypsy, LTR/partial, LTR/unclassified, none, repeat/structural-only, retro/partial, retro/unassigned, unclassified, and (from the transposase branch, dynamic 'DNA/'+token) DNA/hAT, DNA/Tc1-Mariner, DNA/CACTA, DNA/MULE, DNA/IS4-like, DNA/DDE, DNA/unclassified. TEagle superfamily vocabulary: 4.0.0-only 'DNA transposon (transposase detected; termini unassessed)'; 4.1.0-only RETROVIRAL_UNDETERMINED ('LTR retrotransposon (superfamily undetermined; retroviral-type gag/env — Gypsy, Bel-Pao and Retrovirus/ERV are not separated by the bundled panel)'), '<family> (advisory: strongest transposase profile; termini unassessed)', 'classification withheld (protein-domain scan unavailable)'; all others common.
- SET1 R1 canonical case archetypes (12 cases, labels used only for workflow): terminal_direct_repeat (120-bp exact), terminal_inverted_repeat (60-bp exact), terminal_negative (<= 60% identity), partial_direct_repeat (120 bp at 80%); required export formats CSV, TSV, XLSX, FASTA, GFF3, BED, SVG, PNG, PDF, JSON manifest.
- SET1 I2 defect catalogue mapped to archetypes (i2_lock.json): complete_terminal_evidence=input_sequence_substitution_after_analysis; partial_terminal_evidence=stale_classification_artifact; conflicting_terminal_evidence=one_base_coordinate_shift; synthetic_non_te_control=primer_table_swap_between_cases; unsupported_classification_rank=collapsed_unresolved_state; declared_module_failure=removed_tool_or_database_identity; no_admissible_primer_pair=truncated_specificity_presented_complete; capped_specificity_search=manifest_digest_mismatch.

## Development exposure

benchmark_publication_v2/corpora/development_exposure_accessions.tsv lists 170 accessions (122 broad historical corpus, 43 historical holdout, 5 bundled examples) under rule E02_DEVELOPMENT_EXPOSURE.

Records and data added to the development exposure by the 4.1.0 work:

- X04456 (same sequence as bundled M11240)
- K01964
- J00265
- NM_000207
- AH005269
- M11241
- AJ289709
- A84748
- AY037928 (already in the 170-list)
- X13291 (already in the 170-list)
- GCF_000146045.2 with NC_001136.10 (yeast chrIV 100,001-102,000, live genome-scan test)
- GCF_000001215.4
- GCF_000002035.6

Constructed corpora used to measure or tune 4.1.0 behaviour:

- All 4,946 curated Dfam 4.0 DNA-transposon consensi (root, curated partition), placed between seeded random flanks with a random 5-bp TSD (TA for TcMar), were used to measure and choose the 4.1.0 TIR run-trim and microsatellite rules; any benchmark drawn from curated Dfam DNA consensi is development-exposed.
- Synthetic probes: Tc1 (X01005) inserted at TA between (TA)n runs of 8-200 bp at 0-10% divergence; elements whose TIR begins with A(n), (CA)n, (CAG)n, (GATC)n, (TA)n or (TCTA)n heads; palindromic and non-palindromic microsatellite arrays of period 1-6 at 0-20% divergence; random sequence at GC 0.2-0.65; long LTRs of 3.0-5.9 kb; the 905-record structural corpus of the test suite.
- Test fixtures in tests/fixtures and the five bundled examples (M11240, M12927, M80343, X01005, X05424).

A prospective holdout for 4.1.0 must exclude these records and any sequence identical to them, or report them as a separate development-exposed stratum.

## Known limitations of 4.1.0

- Microsatellite arrays with about 15–20% substitutions can still yield a short (10–19 bp) terminal inverted repeat reported as class unassigned.
- A terminal inverted repeat whose own outer end is a run that pairs with itself, such as (TA)n or (GATC)n, cannot be told apart from an insertion into such a run, and the run is cut from the reported arms.
- When a transposon whose terminal repeat begins with a run of another unit, such as A(n) or (CA)n, sits in a (TA)n run, the longer of the two runs decides the reported arms: a longer (TA)n run is cut together with the element's own run, and a longer element run leaves the (TA)n run in the arms or loses the terminal repeat. A (TA)n stretch shorter than about 10 bp beside a Tc1/Mariner element is not cut and can lengthen the arms by a few base pairs.
- An element flanked by an inverted pair of runs of a unit that cannot pair with itself, such as (CA)n on one side and (TG)n on the other, has its terminal-repeat arms extended into those runs.
- A tandem run of 1 kb or more beside an element removes its terminal-repeat candidate, and terminal inverted repeats are searched within about 1.2 kb of each record end.
- The microsatellite guard rejects random arms in AT-rich sequence slightly more often than before (at 20% GC, 3.3% of 10-bp and 1.0% of 20-bp arms, against 2.9% and 0.8%).
- Analysis-result JSON export of flanked Tc1/Mariner records was fixed during this inventory (CLS-031); any 4.1.0 build made before that fix cannot export those records.

## Still to come

The measured behaviour diff (4.0.0 and 4.1.0 run side by side on the benchmark-v2 sequence cache, the test fixtures, the curated Dfam DNA-transposon consensi between flanks and the synthetic probes), the endpoint-by-endpoint impact mapping with proposed new cases, and two independent checks of this inventory are running; they will be added as `behaviour_diff.jsonl`, `behaviour_diff_summary.json` and a second revision of this report.
