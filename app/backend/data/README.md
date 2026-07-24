# te_domains.hmm — bundled TE protein-domain profiles

**File:** `te_domains.hmm` — 21 concatenated Pfam-A HMMER3 profiles used by `teagle_core.domains`
for TE protein-domain detection, covering the full retroviral/ERV **GAG–POL–ENV** architecture plus
chromodomain and the DNA-transposon transposases.

**Source:** InterPro (`https://www.ebi.ac.uk/interpro/api/entry/pfam/<ACC>?annotation=hmm`, gzip)
**Pfam release:** 37.0 (June 2024). Base set retrieved 2026-07-17; retroviral GAG/ENV set added 2026-07-23.
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
| PEG10_N-capsid | PF03732 | GAG | retrotransposon-gag / PEG10-type capsid (LTR-RT gag; NOT HERV-K capsid) |
| Gag_p10 | PF02337 | GAG | retroviral/ERV matrix (MA) |
| Gag_p24 | PF00607 | GAG | retroviral/ERV capsid (CA), N-terminal |
| Gag_p24_C | PF19317 | GAG | retroviral/ERV capsid (CA), C-terminal |
| zf-CCHC_5 | PF14787 | GAG | GAG-polyprotein nucleocapsid zinc-finger (NC) |
| HERV-K_env_2 | PF13804 | ENV | retroviral envelope glycoprotein (HERV-K) |
| GP41 | PF00517 | ENV | retroviral envelope, transmembrane (TM) |
| TLV_coat | PF00429 | ENV | retroviral/gypsy envelope, surface (SU) |
| Chromo | PF00385 | CHR | chromodomain (chromoviruses) |
| HTH_Tnp_Tc3_2 | PF01498 | TPase | Tc1/mariner transposase |
| DDE_1 | PF03184 | TPase | DDE transposase |
| DDE_3 | PF13358 | TPase | DDE transposase (Tc1/mariner) |
| Transposase_1 | PF01359 | TPase | mariner-type transposase |
| Dimer_Tnp_hAT | PF05699 | TPase | hAT transposase |
| hAT-like_RNase-H | PF14372 | TPase | hAT-like transposase |

Retroviral GAG/ENV set added 2026-07-23 (each fetched from the InterPro API, SHA-256 first-16 of the
downloaded HMM): PF00607 `7e7d5714bfce99b2`, PF19317 `e07954f49f4e9c59`, PF02337 `172800445c04c60f`,
PF14787 `0ec761b946995940`, PF13804 `6bd9941cb34e7e78`, PF00517 `7b0cf440ed9a999c`, PF00429 `0ad0c4497fe02ff1`.
The gag/env models close the HERV-K(HML-2) GAG/ENV gap (PF03732 does not annotate betaretroviral capsid;
the bundle previously carried no env model). To rebuild, re-fetch each accession's `?annotation=hmm` (gzip),
decompress, and concatenate in this order.
