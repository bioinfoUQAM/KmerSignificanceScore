# KmerSignificance Score (KSS)

[![Python 3.10-3.14](https://img.shields.io/badge/python-3.10--3.14-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A discriminative and biologically informed framework for viral k-mer prioritization. KSS integrates discriminative power, mutational impact, and protein-level characterization depth into a standardized [0,1] score reported on a common bounded numerical axis.

## Overview

KSS combines three complementary components:

- **Discriminative score**: Information-theoretic measure of strain-distinguishing capacity with adaptive complexity scaling, producing bounded [0,1] scores on a common numerical axis, with an explicit complexity term for the number of classes; the limits of transferring a raw value between datasets are quantified in the manuscript
- **Mutational score**: Biophysical impact assessment using MIYATA_EVO, an optimized amino acid substitution matrix (+28.4% improvement over the original MIYATA matrix)
- **Protein score**: Characterization depth quantified from UniProt annotations (Gene Ontology, protein existence, structural features, interactions, literature)

## Installation

```bash
git clone https://github.com/bioinfoUQAM/KmerSignificanceScore.git
cd KmerSignificanceScore

python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

pip install -e .                 # scoring pipeline
pip install -e ".[analysis]"     # plus the packages the reported analyses need
```

Requires Python >=3.10,<3.15. Every dependency is bounded at both ends in `pyproject.toml`.
To recreate the environment the reported numbers were produced in, install
`requirements-lock.txt`, which gives the exact resolved versions (Python 3.14.6).

## Quick start

The repository ships with a toy dataset: 16 HCMV `US28` sequences across four genotypes,
small enough to score end to end in a few seconds. Run it to check your install and to
see what KSS produces:

```bash
python main.py data/toy/config.yaml
```

It prints its progress and finishes with:

```
✓ Completed: toy
```

and writes the scores to `data/toy/US28/results/US28_compiled_results.json`. Each genomic
position gets a final `kss` score in [0, 1] and the three components it combines. The
highest-scoring positions in the toy are:

| Position | KSS | Discriminative | Mutational | Protein |
|---:|---:|---:|---:|---:|
| 55  | 0.717 | 0.923 | 0.342 | 0.885 |
| 622 | 0.632 | 0.149 | 0.861 | 0.885 |
| 73  | 0.609 | 0.801 | 0.141 | 0.885 |

`kss` is the weighted mean of the three components (equal weights by default): a position
scores high when its k-mers discriminate the genotypes well (**discriminative**), the amino
acid changes they carry are biophysically important (**mutational**), and they fall in a
well-characterized protein (**protein**, constant per gene). To score your own data, copy
`data/toy/` and replace the config, the `sequences.fasta` (headers `>accession|dataset|gene|class`)
and the GenBank reference.

Add `--report N` to also write a PDF summarizing the top `N` positions. The report draws with
matplotlib, so it needs the analysis extra: `pip install -e ".[analysis]"`.

```bash
python main.py data/toy/config.yaml --report 12
```

The report (`data/toy/toy_report.pdf`) has a summary table of the top positions and, per gene,
the discriminative, mutational and KSS scores along the gene with the protein's UniProt link.

To redraw it from results already on disk, without rescoring:

```bash
python main.py data/toy/config.yaml --report-only
```

## Usage

```bash
# Analyze all configured datasets
python main.py

# Analyze a specific dataset
python main.py data/Human_betaherpesvirus_5/config.yaml

# Analyze multiple datasets
python main.py data/*/config.yaml
```

### Recalculate scores without realignment

To adjust KSS parameters (weights, thresholds) without rerunning the costly alignment step:

```bash
python recalculate_scores.py
```

### Configuration

Each dataset has a `config.yaml`:

```yaml
dataset:
  name: "Human_betaherpesvirus_5"
  genes: ["UL55", "UL73", "US28"]

parameters:
  k: 9
  threshold: 25
  weights:
    discriminative: 1
    mutational: 1
    protein: 1
  scoring:
    mutational_matrix: "MIYATA_EVO"
    substitution_matrix: "BLOSUM62"
  alignment:
    open_gap_score: -10
    extend_gap_score: -1
```

## Datasets

The repository includes reference annotations, configuration files, and the exact accession lists used in the manuscript for three families of viruses. **Raw FASTA files are not redistributed** because of their size; they can be reproduced from the accession lists.

| Virus | Sequences | Classes | Source |
|-------|-----------|---------|--------|
| SARS-CoV-2 | 278,738 | 19 variants | [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) |
| HIV-1 | 12,223 | 15 subtypes | [Los Alamos National Laboratory HIV Sequence Database](https://www.hiv.lanl.gov/) |
| HCMV | 399-646 | 4-8 genotypes | [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) |

Accession lists live in `data/accessions/`. To regenerate the FASTA files from public databases, see [`data/README.md`](data/README.md):

```bash
python scripts/fetch_sequences.py --email you@example.com --virus hcmv
```

FASTA files use the header format: `>accession|Virus|Gene|class`

### Partial sequences

Partial records are processed rather than rejected. Regions absent relative to the reference are
represented as alignment gaps, kept as a variant, and scored as indels, which take the configured
indel score, 1.0 by default. KSS does not distinguish a biological deletion from incomplete
sequence coverage, so supply sequences covering the whole target CDS unless that reading is
intended, or narrow the reference CDS to the region every sequence covers.

Internal in-frame insertions and deletions relative to the reference are supported, and are what
the mutational component is meant to score.

## Project Structure

```
KmerSignificanceScore/
├── main.py                          # Main analysis pipeline
├── recalculate_scores.py            # Re-score without realignment
├── pyproject.toml                   # Packaging, bounded dependencies, `kss` entry point
├── requirements.txt                 # The same ranges, for `pip install -r`
├── requirements-lock.txt            # Exact versions of the reported environment
├── CONTRIBUTING.md                  # How the modules fit together, and what is public
├── LICENSE
│
├── src/
│   ├── kanalyzer.py                 # K-mer analysis engine
│   ├── kss.py                       # KSS score computation
│   ├── discriminative_score.py      # Discriminative scoring component
│   ├── mutation_score.py            # Mutational impact scoring
│   ├── protein_score.py             # Protein-level scoring
│   ├── pipeline.py                  # Configuration loading and orchestration
│   ├── report.py                    # Optional PDF report, --report N
│   ├── utils.py                     # Utility functions
│   └── substitution_matrices/       # 38 matrices as JSON, MIYATA_EVO among them
│
├── data/                            # References, configs, and accession lists
│   ├── README.md                    # How to reproduce the FASTA files
│   ├── accessions/                  # Per-virus accession TSVs
│   ├── Severe_acute_respiratory_syndrome_coronavirus_2/
│   ├── Human_immunodeficiency_virus_1/
│   └── Human_betaherpesvirus_5/
│
├── scripts/                         # Everything the revision added, one file per analysis
│   ├── extract_accessions.py        # Build accession lists from FASTA files
│   ├── fetch_sequences.py           # Download sequences from NCBI Entrez
│   ├── kss_selection.py             # The ranking and the selection Fig 3 draws, kept apart
│   ├── statistical_comparison.py    # Paired comparison against the baselines
│   ├── threshold_sensitivity.py     # Prevalence threshold sweep, with its two sweep drivers
│   ├── reference_matrix.py          # COV2Var reference layers and error rates
│   ├── floor_sensitivity.py         # What the 0.1 floor of Eq. (9) binds for and changes
│   ├── build_s7_table.py …          # One builder per supplementary workbook
│   └── …                            # 30 files in all; every one is cited where it is used
│
├── results/                         # One directory per analysis, regenerated by its script
├── supplementary/                   # S-numbered workbooks as submitted
│
└── notebooks/                       # Validation and evaluation
    ├── discriminative_score_validation.ipynb
    ├── matrix_evaluation.ipynb
    ├── matrix_optimization.ipynb
    ├── generate_kss_figure.ipynb
    ├── *_results/                   # Pre-computed validation results
    └── publication_figures/          # Figures for the manuscript
```

To change the code rather than run it, [`CONTRIBUTING.md`](CONTRIBUTING.md) gives the flow through
the modules, what each one owns, the public API, and how to add a substitution matrix.

## Validation Results

KSS was validated on all three families of viruses:

- **Discriminative component**: Mean F1 = 0.881 across all datasets, similar to the strongest of seven established feature selection methods with no conclusive difference either way (chi-squared, odds ratio, NMI, MI, ANOVA, Cramer's V, normalized symmetric chi-squared), and significantly above three of the remaining six
- **MIYATA_EVO matrix**: Composite agreement score of 4.578 across eight amino acid characteristics, vs 3.566 for the original MIYATA matrix (+28.4%), optimized via genetic algorithm over 625 generations
- **Protein characterization score**: Spearman rho = 0.900, Kendall tau = 0.775 against UniProt annotation levels on 17,471 viral proteins, which calibrates the component rather than validating it, both being computed from the same records
- **Functional validation**: Top-ranked positions across all three families of viruses correspond to established variant-defining mutations, drug resistance sites, immune escape loci, and genotype markers documented in independent studies

Detailed results are available in `notebooks/*_results/` directories.

## Citation

If you use KSS in your research, please cite:

> Lebatteux D, Corso F, Soudeyns H, Boucoiran I, Gantt S, Diallo AB. KmerSignificance Score: A discriminative and biologically informed framework for viral k-mer prioritization. *Submitted to PLOS Computational Biology*.

## License

MIT License - see [LICENSE](LICENSE) for details.

The MIT licence covers the code. It does not cover everything the repository carries:
[THIRD_PARTY.md](THIRD_PARTY.md) lists each item of third-party material, its origin, whether any
upstream file is redistributed, and the terms that govern it. Eleven UniProtKB entries ship under
CC BY 4.0, and thirty-seven of the thirty-eight substitution matrices are files written by this
project whose values reproduce published sources, so the MIT licence is not claimed over them.
Read it before reusing anything but the code.

## Contact

For questions or issues, please open an issue on GitHub.
