# te_domains.hmm — bundled TE protein-domain profiles

**File:** `te_domains.hmm` — 30 concatenated Pfam-A HMMER3 profiles used by `teagle_core.domains`
for TE protein-domain detection, covering the retroviral/ERV **GAG–POL–ENV** architecture, LINE ORF1p and
endonuclease, tyrosine recombinase, Helitron helicase, chromodomain and the DNA-transposon transposases.

**Source:** InterPro (`https://www.ebi.ac.uk/interpro/api/entry/pfam/<ACC>?annotation=hmm`, gzip)
**Pfam release:** not recorded per profile. The `DATE` line of each model (last column below) is the build date Pfam stamped on it and runs from May to December 2025, so the set is later than Pfam 37.0 (June 2024); the bundled file, not a release number, is the reference. Base set retrieved 2026-07-17; retroviral GAG/ENV set added 2026-07-23; LINE, YR, Helitron, CACTA, MULE and IS4 set added 2026-07-28.
**License:** Pfam data is **CC0** (public domain) — redistributable.
**Cite:** Mistry J. et al. (2021) *Pfam: The protein families database in 2021.* NAR 49:D412–D419. doi:10.1093/nar/gkaa913

| Profile | Pfam (versioned) | Code | Used for | HMM built |
|---|---|---|---|---|
| RVT_1 | PF00078.33 | RT | reverse transcriptase | May 28 2025 |
| RVT_2 | PF07727.20 | RT | reverse transcriptase (copia-type) | May 28 2025 |
| RVT_3 | PF13456.13 | RNaseH | RNase H-like fold (Pfam: reverse transcriptase-like); coded RNaseH, not RT | May 28 2025 |
| rve | PF00665.33 | INT | integrase (Copia/Gypsy discriminator) | May 28 2025 |
| RNase_H | PF00075.31 | RNaseH | RNase H | May 28 2025 |
| RVP | PF00077.26 | PR | aspartic protease | May 28 2025 |
| PEG10_N-capsid | PF03732.25 | GAG | retrotransposon-gag / PEG10-type capsid (LTR-RT gag; not HERV-K capsid) | Dec 19 2025 |
| Chromo | PF00385.31 | CHR | chromodomain (chromoviruses) | May 28 2025 |
| HTH_Tnp_Tc3_2 | PF01498.25 | TPase | Tc1/mariner transposase | May 28 2025 |
| DDE_1 | PF03184.26 | TPase | DDE transposase | May 28 2025 |
| DDE_3 | PF13358.13 | TPase | DDE transposase (Tc1/mariner) | May 28 2025 |
| Transposase_1 | PF01359.25 | TPase | mariner-type transposase | May 28 2025 |
| Dimer_Tnp_hAT | PF05699.20 | TPase | hAT transposase | May 28 2025 |
| hAT-like_RNase-H | PF14372.13 | TPase | hAT-like transposase | Aug 7 2025 |
| Gag_p24 | PF00607.26 | GAG | retroviral/ERV capsid (CA), N-terminal | May 28 2025 |
| Gag_p24_C | PF19317.5 | GAG | retroviral/ERV capsid (CA), C-terminal | May 28 2025 |
| Gag_p10 | PF02337.23 | GAG | retroviral/ERV matrix (MA) | May 28 2025 |
| zf-CCHC_5 | PF14787.12 | GAG | GAG-polyprotein nucleocapsid zinc-finger (NC) | May 28 2025 |
| HERV-K_env_2 | PF13804.13 | ENV | retroviral envelope glycoprotein (HERV-K) | Sep 12 2025 |
| GP41 | PF00517.23 | ENV | retroviral envelope, transmembrane (TM) | May 28 2025 |
| TLV_coat | PF00429.25 | ENV | retroviral/gypsy envelope, surface (SU) | May 28 2025 |
| Transposase_22 | PF02994.21 | ORF1 | LINE ORF1p, RNA-binding domain | May 28 2025 |
| Tnp_22_trimer | PF17489.8 | ORF1 | LINE ORF1p, trimerisation domain | May 28 2025 |
| Tnp_22_dsRBD | PF17490.8 | ORF1 | LINE ORF1p, dsRBD-like domain | May 28 2025 |
| Exo_endo_phos | PF03372.30 | EN | APE-like endonuclease (ORF2p; broad family, counted only N-terminal to an RT) | May 29 2025 |
| Phage_integrase | PF00589.28 | YR | tyrosine recombinase (DIRS signal; broad family) | May 28 2025 |
| Helitron_like_N | PF14214.12 | HEL | Helitron helicase-like domain | May 29 2025 |
| Transposase_24 | PF03004.20 | TPase | CACTA/En-Spm transposase | May 28 2025 |
| MULE | PF10551.16 | TPase | MULE/Mutator transposase | May 28 2025 |
| DDE_Tnp_1_7 | PF13843.13 | TPase | IS4-like DDE transposase (piggyBac-like) | May 28 2025 |

Retroviral GAG/ENV set added 2026-07-23 (each fetched from the InterPro API, SHA-256 first-16 of the
downloaded HMM): PF00607 `7e7d5714bfce99b2`, PF19317 `e07954f49f4e9c59`, PF02337 `172800445c04c60f`,
PF14787 `0ec761b946995940`, PF13804 `6bd9941cb34e7e78`, PF00517 `7b0cf440ed9a999c`, PF00429 `0ad0c4497fe02ff1`.
The gag/env models close the HERV-K(HML-2) GAG/ENV gap (PF03732 does not annotate betaretroviral capsid;
the bundle previously carried no env model).

RVT_3 (PF13456, named "Reverse transcriptase-like" in Pfam) belongs to the RNase H clan CL0219 and is integrated
in InterPro under IPR002156 (Ribonuclease H domain), so TEagle codes it RNaseH rather than RT.

To rebuild, re-fetch each accession's `?annotation=hmm` (gzip), decompress, and concatenate the 30 profiles in
the table order. InterPro serves the current version of each model, so a fresh fetch can differ from the
versioned accessions above; the file shipped here is the reference, and a rebuilt file reproduces a sealed
result only if it matches the `HMM_SHA256` recorded in that run's manifest.
