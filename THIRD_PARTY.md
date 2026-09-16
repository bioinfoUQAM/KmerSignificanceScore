# Third-party material

`LICENSE` covers the code, documentation and data files this project authored. It does **not**
cover the third-party material listed here, and nothing in this repository purports to relicense
any of it. Each entry below says what the material is, where it comes from, whether this
repository redistributes it, and under what terms it may be used.

Three categories are distinguished throughout, because they carry different obligations:

- **Redistributed** — the file itself is committed here.
- **Downloaded and verified** — the file is not committed; a build script fetches it from its
  upstream and checks it against a pinned SHA-256 before use. The copies under
  `scripts/added_matrix_sources/` are ignored local caches.
- **Derived values only** — no upstream file is committed. A matrix in
  `src/substitution_matrices/` carries values that match a published matrix, or is computed by this
  project from a published model.

---

## Redistributed in this repository

### Manuscript typesetting files (3 files)

`LICENSE` states that this document lists every item of third-party material, and these three were
missing from it. They live under `manuscript/`, which the public repository does not carry, so the
omission had no practical effect there, but the claim was false of the repository as it stands.

| File | Origin | Terms |
|---|---|---|
| `manuscript/plos2025.bst` | the Vancouver BibTeX style of Folkert van der Beek, <https://gitlab.com/fvdbeek/vancouver.bst> | LaTeX Project Public License 1.3 or later. Redistribution is permitted and the notice travels inside the file |
| `manuscript/plos_latex_template.tex` | the PLOS submission template, version 3.7 | supplied by PLOS for preparing submissions to their journals |
| `manuscript/plos_bibtex_sample.bib` | the bibliography of this manuscript, built on the PLOS sample file | the entries are bibliographic records; the sample scaffolding comes from the same PLOS template |

### UniProtKB entries (11 files)

`src/uniprot/` contains eleven UniProtKB/Swiss-Prot entries, one per analysed gene, retrieved from
the UniProt REST API (release 2025_02) and totalling 4,691,689 bytes:

| Accession | Entry | Gene |
|---|---|---|
| `P04578` | ENV_HV1H2 | HIV-1 *env* |
| `P04585` | POL_HV1H2 | HIV-1 *pol* |
| `P04591` | GAG_HV1H2 | HIV-1 *gag* |
| `P06473` | GB_HCMVA | HCMV *UL55* |
| `P0DTC2` | SPIKE_SARS2 | SARS-CoV-2 *S* |
| `P0DTC4` | VEMP_SARS2 | SARS-CoV-2 *E* |
| `P0DTC5` | VME1_SARS2 | SARS-CoV-2 *M* |
| `P0DTC9` | NCAP_SARS2 | SARS-CoV-2 *N* |
| `P0DTD1` | R1AB_SARS2 | SARS-CoV-2 *ORF1ab* |
| `P16795` | GN_HCMVA | HCMV *UL73* |
| `P69332` | US28_HCMVA | HCMV *US28* |

**Terms.** UniProt data are distributed under the Creative Commons Attribution 4.0 International
licence (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/>. **These eleven files are
outside the scope of `LICENSE`** and remain under CC BY 4.0.

**Attribution.** The UniProt Consortium, UniProtKB/Swiss-Prot, <https://www.uniprot.org/>.

**Why they are committed rather than fetched.** The protein component is computed from these
records. Fetching them at run time made the published scores depend on the service answering, and
on 29 July 2026 an offline run scored every SARS-CoV-2 protein at zero. The service has since moved
on: it now serves later entry versions than the ones the published scores were computed from, so it
no longer reproduces them. The committed copies are the records the results were computed from.

**Verification.** Each of the eleven is byte-for-byte identical to the copy deposited at
<https://doi.org/10.5281/zenodo.21937072> and matches the SHA-256 recorded for it in
`results/protein_score_validation/protein_score_provenance.csv`.

> **`src/uniprot/index.json` is not covered by this attribution.** It is a 286-byte mapping from
> `taxid|gene` to accession, written by `src/protein_score.py`, and contains no UniProt content. It
> is part of this project and falls under `LICENSE`.

---

## Downloaded and verified, not redistributed

None of the files below is committed, but they are not all obtained the same way.

The five matrix sources are fetched automatically: `scripts/build_added_matrices.py` and
`scripts/build_viral_matrices.py` download each one from a pinned upstream revision and verify a
pinned SHA-256 before reading. Run both builders once after a fresh clone; no manual download is
needed.

The two external datasets are **obtained manually**. `scripts/external_inputs.py` does not download
them: it holds their expected SHA-256 digests and refuses to read a file that does not match,
reporting where to obtain it. Follow its instructions once, then the analyses that read them will
run.

| File | Upstream | Terms |
|---|---|---|
| `iqtree_modelprotein.cpp` | IQ-TREE 2, <https://github.com/iqtree/iqtree2> | GPL-2.0-or-later |
| `raxml_models.c` | Standard RAxML, <https://github.com/stamatak/standard-RAxML> | GPL-2.0 |
| `VTML200_spaln.txt` | spaln, <https://github.com/ogotoh/spaln>, `table/vtml200` | spaln is GPL-2.0-or-later; the table itself was produced by K. Kneutgen and T. Müller and carries no licence of its own |
| `FLAVI_authors.PAML` | <https://github.com/thulekm/flavi> | **No licence is published**, which grants no redistribution right by default |
| `wag.dat` | N. Goldman's server at the EBI, <https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat> | The authors' own file; it carries a citation request and no licence statement |
| COV2Var annotation layers (7 files, **obtained manually**) | <https://biomedbdc.wchscu.cn/COV2Var/download/> | **No explicit licence identified.** The download page carries only "© COV2Var. All Rights Reserved."; not redistributed by this repository. The snapshot the published numbers were computed from is described, file by file, in `results/functional_benchmark/cov2var_input_manifest.csv`: name, size, SHA-256, the members of the selection archive, and what each layer is used for. Archiving the bytes themselves would require the authors' written permission, which has not been sought |
| Amino acid fitness estimates (**obtained manually**) | `jbloomlab/SARS2-mut-fitness`, `results/aa_fitness/aa_fitness.csv`, at commit `067fce16` | **MIT License** at that commit (`LICENSE.md`, Copyright (c) 2023 Jesse D Bloom and Richard A Neher) |

Citations: IQ-TREE (Minh et al.), RAxML (Stamatakis), spaln (Gotoh), VTML (Müller and Vingron
2000; Müller, Spang and Vingron 2002), FLAVI (Le and Vinh 2020), WAG (Whelan and Goldman 2001),
COV2Var (Feng et al. 2024), amino acid fitness (Bloom and Neher 2023). Full references are in the
manuscript bibliography.

---

## Derived values only, no upstream file redistributed

The directory `src/substitution_matrices/` contains 38 matrices: **37 candidate matrices based on
published sources**, documented below, and **`MIYATA_EVO`, an original output generated by this
project** through genetic-algorithm optimization initialized from MIYATA. `MIYATA_EVO` is covered by
the project's MIT License; Miyata et al. (1979) are cited as the source of the starting matrix.

The 37 others are JSON files written by this project, but they do not all relate to their source in
the same way. Thirty reproduce published scoring or distance matrices. `VTML200` reproduces its
upstream score table. `WAG` and the five viral rate models (`HIVb`, `HIVw`, `rtREV`, `FLU`, `FLAVI`)
are **converted by this project into score matrices** from documented upstream rate parameters, by
the construction described in `scripts/build_added_matrices.py`; their numerical values appear in no
publication.

The sources are credited below. Whether a derived numerical matrix is itself covered by the source's
terms is not settled by this document; what follows is attribution, not a licence claim.

### Biopython

Twenty-three matrices agree numerically with those distributed with Biopython
(`Bio.Align.substitution_matrices`), to within 1e-9 across all 400 ordered amino-acid pairs:
BENNER22, BENNER6, BENNER74, BLASTP, BLOSUM45, BLOSUM50, BLOSUM80, BLOSUM90, DAYHOFF, FENG, GENETIC,
GONNET1992, JOHNSON, JONES, LEVIN, MCLACHLAN, MDM78, PAM250, PAM30, PAM70, RAO, RISLER, STR.

The agreement is numerical, not a file identity: ours are JSON and Biopython's are tab-delimited
text, so the files differ in bytes and in size. Biopython's own data files cite the primary
publication for each matrix rather than a database.

No Biopython file is redistributed here. Its permission notice is nevertheless reproduced in full
below, since the Biopython License Agreement asks that it accompany supporting documentation.

### AAindex

Sixteen matrices carry values matching entries of the AAindex2 database: four for which AAindex2 is
the only source found (`BLOSUM62` → HENS920102, `GRANTHAM` → GRAR740104, `MIYATA` → MIYT790101,
`PAM120` → ALTS910101), two also present in the spaln distribution (`VTML160` → MUET020101,
`VTML250` → MUET020102), and ten also present in Biopython (BENNER22, BLOSUM45, BLOSUM80, FENG,
GENETIC, GONNET1992, LEVIN, RAO, RISLER, STR). Values alone cannot say which source a transcription
came from, and no acquisition script survives, so the overlaps are stated rather than resolved.

**No explicit licence or redistribution terms were found**, neither in the `aaindex2` file nor on
the AAindex site; the database asks to be cited. Citation: Kawashima et al., AAindex, and the
primary publication of each entry.

### Sneath (1966)

`SNEATH` carries the dissimilarity index *D* of Table 2, pages 162-163, of Sneath P.H.A. (1966),
"Relations between chemical structure and biological activity in peptides", *Journal of Theoretical
Biology* 12(2), 157-195. The values are the published integers; no conversion is applied.

---

## Sequence data

**The analysis datasets are not redistributed.** The 278,738 SARS-CoV-2, 12,223 HIV-1 and 399 to 646
HCMV sequences the results are computed from are public records of GenBank and of the Los Alamos
National Laboratory HIV Sequence Database, listed by accession under `data/accessions/` so each
dataset can be rebuilt.

**Two kinds of sequence file are committed**, and both are public NCBI records:

| What | Where | Count |
|---|---|---|
| GenBank reference records, one per analysed gene directory | `data/*/*/` | 12 files, 4 distinct accessions: `NC_045512`, `AF033819`, `FJ527563`, `FJ616285` |
| The quick-start toy dataset, 16 HCMV *US28* sequences | `data/toy/US28/sequences.fasta` | 1 file |

NCBI places no restrictions on the use or distribution of GenBank data, while noting that some
submitters may claim patent, copyright or other intellectual property rights in their submissions.
See GenBank Data Usage, <https://www.ncbi.nlm.nih.gov/genbank/#genbank-data-usage>, and the NCBI
policies statement, <https://www.ncbi.nlm.nih.gov/home/about/policies/>.

---

## Biopython License Agreement

Reproduced in full, as its terms request. Biopython is not redistributed by this project.

```
Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
```
