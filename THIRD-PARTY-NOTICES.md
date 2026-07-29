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
| [pyhmmer](https://github.com/althonos/pyhmmer) | 0.12.1 | MIT | HMMER3 profile search, in-process |
| [openpyxl](https://openpyxl.readthedocs.io/) | 3.1.5 | MIT | XLSX table export |
| et_xmlfile | (openpyxl dependency) | MIT | XLSX writing |
| [certifi](https://github.com/certifi/python-certifi) | 2026.6.17 | MPL-2.0 | CA bundle for NCBI/Ensembl requests |
| Pfam-A TE domain profiles (`data/te_domains.hmm`) | Pfam 37.x | CC0-1.0 | Protein-domain panel |
| [Cascadia Code / Cascadia Mono](https://github.com/microsoft/cascadia-code) | — | SIL OFL 1.1 | Bundled UI typeface |
| [Roboto](https://github.com/googlefonts/roboto) | — | Apache-2.0 | Bundled UI typeface |

primer3-py is distributed under "version 2 of the License, or (at your option) any later version", which is
what allows it to be combined with AGPL-3.0 code. PySide6 and shiboken6 are used under their LGPL-3.0 option.
MPL-2.0 and Apache-2.0 are both compatible with GPL-3.0 / AGPL-3.0.

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

These are fetched from bioconda by the in-app installer and run as separate processes. TEagle does not
redistribute them; they are listed because their terms affect what you may do with the results.

| Tool | Licence | Note |
|---|---|---|
| [RepeatMasker](https://www.repeatmasker.org/) | Open Software License 2.1 | Family naming |
| [Dfam](https://dfam.org/) | CC0-1.0 (curated partitions) | Family library |
| [minimap2](https://github.com/lh3/minimap2) | MIT | Splice alignment |
| [miniprot](https://github.com/lh3/miniprot) | MIT | Homology alignment (currently dormant) |
| [isPcr](https://genome.ucsc.edu/) | UCSC — free for academic, non-profit and personal use | **Commercial use requires a licence from UCSC.** Relevant if TEagle is ever hosted as a paid service. |
| [ncbi-datasets-cli](https://www.ncbi.nlm.nih.gov/datasets/) | US Government public domain | Assembly download |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | Own licence (see *Deliberately NOT bundled*) | Optional second primer-QC engine (`teagle-vrna` env). **Redistribution for a fee is forbidden** — relevant if TEagle is ever hosted as a paid service. |
| [micromamba](https://github.com/mamba-org/mamba) | BSD-3-Clause | Environment manager |

RepBase is **not** used, bundled, or downloaded by TEagle at any point; it is not redistributable.

## Reference data retrieved at run time

Sequences and assemblies fetched from NCBI, ENA and Ensembl are retrieved by the user's own installation on
demand and are not redistributed with TEagle. NCBI data is generally not subject to copyright within the
United States; users remain responsible for the terms attaching to any record they retrieve.
