# te_domains.hmm — bundled TE protein-domain profiles

**File:** `te_domains.hmm` — 14 concatenated Pfam-A HMMER3 profiles used by `teagle_core.domains`
for TE protein-domain detection (reverse transcriptase, integrase, RNase H, protease, gag,
chromodomain, transposases).

**Source:** InterPro (`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<ACC>/?annotation=hmm`)
**Retrieved:** 2026-07-17
**License:** Pfam data is **CC0** (public domain) — redistributable; see `Deliverables/01_resource_matrix.md`.
**Cite:** Mistry J. et al. (2021) *Pfam: The protein families database in 2021.* NAR 49:D412–D419. doi:10.1093/nar/gkaa913

| Profile | Pfam | Domain (TEagle code) | Used for |
|---|---|---|---|
| RVT_1 | PF00078 | RT | retro reverse transcriptase |
| RVT_2 | PF07727 | RT | retro reverse transcriptase (copia-type) |
| RVT_3 | PF13456 | RT | retro reverse transcriptase |
| rve | PF00665 | INT | integrase (Copia/Gypsy discriminator) |
| RNase_H | PF00075 | RNaseH | RNase H |
| RVP | PF00077 | PR | aspartic protease |
| PEG10_N-capsid | PF03732 | GAG | gag / capsid |
| Chromo | PF00385 | CHR | chromodomain (chromoviruses) |
| HTH_Tnp_Tc3_2 | PF01498 | TPase | Tc1/mariner transposase |
| DDE_1 | PF03184 | TPase | DDE transposase |
| DDE_3 | PF13358 | TPase | DDE transposase (Tc1/mariner) |
| Transposase_1 | PF01359 | TPase | mariner-type transposase |
| Dimer_Tnp_hAT | PF05699 | TPase | hAT transposase |
| hAT-like_RNase-H | PF14372 | TPase | hAT-like transposase |

To rebuild, re-fetch each accession's `?annotation=hmm` (gzip), decompress, and concatenate in this order.
