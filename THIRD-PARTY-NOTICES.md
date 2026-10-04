# Third-party notices

TEagle is licensed under the GNU Affero General Public License, version 3 or later (see [LICENSE](LICENSE)).
It incorporates or depends on the components below, each under its own licence. Versions are the pinned ones
in `app/backend/requirements.txt`; the provenance manifest of every run records the versions actually used.

## Bundled in the released application

| Component | Version | Licence | Role |
|---|---|---|---|
| [primer3-py](https://github.com/libnano/primer3-py) | 2.3.0 | GPL-2.0-**or-later** | Primer design and nearest-neighbour thermodynamics |
| [PySide6](https://doc.qt.io/qtforpython/) | 6.11.1 | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only | Qt user interface |
| [shiboken6](https://doc.qt.io/qtforpython/shiboken6/) | 6.11.1 | LGPL-3.0-only OR GPL-2.0-only OR GPL-3.0-only | PySide6 binding runtime |
| [pyhmmer](https://github.com/althonos/pyhmmer) | 0.12.1 | MIT AND BSD-3-Clause AND BSD-2-Clause | HMMER3 profile search, in-process; the wheel statically links HMMER and Easel (BSD-3-Clause, see below) |
| [openpyxl](https://openpyxl.readthedocs.io/) | 3.1.5 | MIT | XLSX table export |
| [psutil](https://github.com/giampaolo/psutil) | 7.2.2 | BSD-3-Clause | Process and memory inspection for the resource gates |
| et_xmlfile | 2.0.0 (openpyxl dependency) | MIT | XLSX writing |
| [certifi](https://github.com/certifi/python-certifi) | 2026.7.22 | MPL-2.0 | CA bundle for NCBI/Ensembl requests |
| [pypdf](https://github.com/py-pdf/pypdf) | 6.15.0 | BSD-3-Clause | PDF metadata minimisation while preserving document navigation |
| Pfam-A TE domain profiles (`data/te_domains.hmm`) | 30 profiles, listed below | CC0-1.0 | Protein-domain panel |
| CPython runtime | 3.12 | PSF-2.0 | Python interpreter and standard library frozen into the executable |
| OpenSSL libraries linked by CPython | as shipped with CPython 3.12 | Apache-2.0 | TLS for NCBI/Ensembl requests |
| [Cascadia Code / Cascadia Mono](https://github.com/microsoft/cascadia-code) | — | SIL OFL 1.1 | Bundled UI typeface |
| [Roboto](https://github.com/googlefonts/roboto) | — | Apache-2.0 | Bundled UI typeface |

primer3-py is distributed under "version 2 of the License, or (at your option) any later version", which is
what allows it to be combined with AGPL-3.0 code. PySide6 and shiboken6 are used under their LGPL-3.0 option.
MPL-2.0 and Apache-2.0 are both compatible with GPL-3.0 / AGPL-3.0.

### Licence texts installed with the application

The installer places `LICENSE` (AGPL-3.0-or-later), this file, `LGPL-3.0.txt`, `GPL-3.0.txt`, `HMMER-BSD-3-Clause.txt` and `Easel-BSD-2-Clause.txt` in the `licenses` folder of the installation directory. The LGPL-3.0 text is required because Qt and PySide6 are redistributed under that option; the wheels carry only the commercial-licence notice. The Qt and PySide6 libraries are used unmodified and are loaded as separate shared libraries (DLLs) in the installation folder, so they can be replaced by a user-built version. Their source code is available from <https://download.qt.io/official_releases/QtForPython/> and <https://code.qt.io/cgit/pyside/pyside-setup.git/>. primer3-py's GPL-2.0 text, psutil's, certifi's, pypdf's, openpyxl's and pyhmmer's licence files travel inside the bundled `*.dist-info` directories.

### HMMER and Easel

pyhmmer statically links HMMER 3 and its Easel library. HMMER is distributed under the BSD three-clause licence and Easel under the BSD two-clause licence, which accounts for the licence expression declared by the pyhmmer wheel. The licence files are installed verbatim as `HMMER-BSD-3-Clause.txt` and `Easel-BSD-2-Clause.txt` in the `licenses` folder. They were taken from the bioconda hmmer 3.3.2 package, not from the HMMER release that pyhmmer 0.12.1 compiles in, so their copyright years and contributor lists may be older than that release's own files. Replace them with the files of the bundled release once those are verified.

### Pfam profiles in the bundled panel

The bundled file holds 30 HMMER3 profiles from Pfam-A, retrieved through InterPro (initial set 2026-07-17, retroviral GAG/ENV set 2026-07-23, the remaining profiles later). The Pfam release number is not asserted here; each profile is identified by the accession and version stored in the file itself: RVT_1 PF00078.33, RVT_2 PF07727.20, RVT_3 PF13456.13, rve PF00665.33, RNase_H PF00075.31, RVP PF00077.26, PEG10_N-capsid PF03732.25, Chromo PF00385.31, HTH_Tnp_Tc3_2 PF01498.25, DDE_1 PF03184.26, DDE_3 PF13358.13, Transposase_1 PF01359.25, Dimer_Tnp_hAT PF05699.20, hAT-like_RNase-H PF14372.13, Gag_p24 PF00607.26, Gag_p24_C PF19317.5, Gag_p10 PF02337.23, zf-CCHC_5 PF14787.12, HERV-K_env_2 PF13804.13, GP41 PF00517.23, TLV_coat PF00429.25, Transposase_22 PF02994.21, Tnp_22_trimer PF17489.8, Tnp_22_dsRBD PF17490.8, Exo_endo_phos PF03372.30, Phage_integrase PF00589.28, Helitron_like_N PF14214.12, Transposase_24 PF03004.20, MULE PF10551.16, DDE_Tnp_1_7 PF13843.13.

## Deliberately NOT bundled

**[ViennaRNA](https://www.tbi.univie.ac.at/RNA/)** (2.7.2) is an optional, user-installed dependency and is
**not** shipped in the installer. Its licence permits research, educational and commercial use but forbids
redistribution "for any fee, other than media costs" and asks that inclusion in a commercial product be
cleared with the authors. A copyleft licence cannot carry those added restrictions, so ViennaRNA cannot be
distributed inside an AGPL work.

TEagle's primer secondary-structure QC therefore runs Primer3's thermodynamics alone unless the user has
installed ViennaRNA themselves, in which case the independent second engine and the cross-engine agreement
flag become available. The interface reports which engines actually ran, and the provenance manifest records
it. Install it with `pip install ViennaRNA` or `conda install -c bioconda viennarna`, or let the in-app
backend installer add it as the one-click **ViennaRNA (primer QC)** component — this triggers a bioconda
download into a separate WSL environment (`teagle-vrna`) on the user's machine and redistributes nothing.

## Installed by the user into the optional WSL2 backend

These are fetched by the in-app installer and run as separate processes. TEagle does not redistribute
them; the installer pins their versions separately from the native Python requirements, and they are
listed because their terms affect what you may do with the results.

| Tool | Licence | Note |
|---|---|---|
| [RepeatMasker](https://www.repeatmasker.org/) | Open Software License 2.1 | Family naming |
| [Dfam](https://dfam.org/) | CC0-1.0 (curated partitions) | Family library |
| [minimap2](https://github.com/lh3/minimap2) | MIT | Splice alignment |
| [isPcr](https://genome.ucsc.edu/) | UCSC — free for academic, non-profit and personal use | **Commercial use requires a licence from UCSC.** Relevant if TEagle is ever hosted as a paid service. |
| [faToTwoBit](https://genome.ucsc.edu/) | UCSC — free for academic, non-profit and personal use | Optional compact-genome-cache conversion. **Commercial use requires a licence from UCSC.** |
| [ncbi-datasets-cli](https://www.ncbi.nlm.nih.gov/datasets/) | US Government public domain | Assembly download |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | Own licence (see *Deliberately NOT bundled*) | Optional second primer-QC engine (`teagle-vrna` env). **Redistribution for a fee is forbidden** — relevant if TEagle is ever hosted as a paid service. |
| [micromamba](https://github.com/mamba-org/mamba) | BSD-3-Clause | Environment manager |

RepBase is **not** used, bundled, or downloaded by TEagle at any point; it is not redistributable.

## Reference data retrieved at run time

Sequences and assemblies fetched from NCBI, ENA and Ensembl are retrieved by the user's own installation on
demand and are not redistributed with TEagle. NCBI data is generally not subject to copyright within the
United States; users remain responsible for the terms attaching to any record they retrieve.
