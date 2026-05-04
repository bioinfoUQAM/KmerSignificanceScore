# KmerSignificance Score (KSS)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A discriminative and biologically-informed framework for viral k-mer prioritization. KSS integrates discriminative power, mutational impact, and protein-level functional importance into a standardized [0,1] score enabling direct cross-dataset comparison.

## Overview

KSS combines three complementary components:

- **Discriminative score**: Information-theoretic measure of strain-distinguishing capacity with adaptive complexity scaling, producing bounded [0,1] scores comparable across datasets regardless of the number of classes
- **Mutational score**: Biophysical impact assessment using MIYATA_EVO, an optimized amino acid substitution matrix (+28.4% improvement over the original MIYATA matrix)
- **Protein score**: Functional importance quantification from UniProt annotations (Gene Ontology, protein existence, structural features, interactions, literature)

## Installation

```bash
git clone https://github.com/bioinfoUQAM/KmerSignificanceScore.git
cd KmerSignificanceScore

python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

pip install -r requirements.txt
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

The repository includes reference annotations, configuration files, and the exact accession lists used in the manuscript for three classes of viruses. **Raw FASTA files are not redistributed** because of their size; they can be reproduced from the accession lists.

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

## Project Structure

```
KmerSignificanceScore/
├── main.py                          # Main analysis pipeline
├── recalculate_scores.py            # Re-score without realignment
├── requirements.txt
├── LICENSE
│
├── src/
│   ├── kanalyzer.py                 # K-mer analysis engine
│   ├── kss.py                       # KSS score computation
│   ├── discriminative_score.py      # Discriminative scoring component
│   ├── mutation_score.py            # Mutational impact scoring
│   ├── protein_score.py             # Protein-level scoring
│   ├── pipeline.py                  # Pipeline orchestration
│   ├── utils.py                     # Utility functions
│   └── substitution_matrices/
│       └── MIYATA_EVO.pkl           # Optimized substitution matrix
│
├── data/                            # References, configs, and accession lists
│   ├── README.md                    # How to reproduce the FASTA files
│   ├── accessions/                  # Per-virus accession TSVs
│   ├── Severe_acute_respiratory_syndrome_coronavirus_2/
│   ├── Human_immunodeficiency_virus_1/
│   └── Human_betaherpesvirus_5/
│
├── scripts/
│   ├── extract_accessions.py        # Build accession lists from FASTA files
│   └── fetch_sequences.py           # Download sequences from NCBI Entrez
│
└── notebooks/                       # Validation and evaluation
    ├── discriminative_score_validation.ipynb
    ├── matrix_evaluation.ipynb
    ├── matrix_optimization.ipynb
    ├── protein_score_validation.ipynb
    ├── generate_kss_figure.ipynb
    ├── *_results/                   # Pre-computed validation results
    └── publication_figures/          # Figures for the manuscript
```

## Validation Results

KSS was validated on all three classes of viruses:

- **Discriminative component**: Mean F1 = 0.880 across all datasets, comparable to or above six established feature selection methods (chi-squared, odds ratio, NMI, MI, ANOVA, Cramer's V)
- **MIYATA_EVO matrix**: Composite biophysical correlation score of 4.578 vs 3.566 for the original MIYATA matrix (+28.4%), optimized via genetic algorithm over 625 generations
- **Protein score**: Spearman rho = 0.900, Kendall tau = 0.777 against UniProt annotation levels on 17,470 viral proteins
- **Functional validation**: Top-ranked positions across all three classes of viruses correspond to established variant-defining mutations, drug resistance sites, immune escape loci, and genotype markers documented in independent studies

Detailed results are available in `notebooks/*_results/` directories.

## Citation

If you use KSS in your research, please cite:

> Lebatteux D, Corso F, Soudeyns H, Boucoiran I, Gantt S, Diallo AB. KmerSignificance Score: A discriminative and biologically-informed framework for viral k-mer prioritization. *Submitted to PLOS Computational Biology*.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contact

For questions or issues, please open an issue on GitHub.
