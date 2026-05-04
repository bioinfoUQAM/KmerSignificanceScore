# Data — sequence accessions and reproduction

Raw nucleotide sequences are not redistributed in this repository (they exceed
GitHub size limits and, for SARS-CoV-2, are subject to GISAID terms of use).
Instead, the exact accession lists used in the manuscript are provided in
`accessions/` so the analysis can be reproduced from public databases.

## Accession files

| File | Rows | Source | Class label |
|---|---:|---|---|
| `accessions/sars_cov2.tsv` | 278,738 | NCBI GenBank | variant |
| `accessions/hiv1.tsv`      | 12,223  | Los Alamos National Laboratory HIV Sequence Database (also indexed in GenBank) | subtype |
| `accessions/hcmv_UL55.tsv` | 399     | NCBI GenBank | genotype |
| `accessions/hcmv_UL73.tsv` | 646     | NCBI GenBank | genotype |
| `accessions/hcmv_US28.tsv` | 443     | NCBI GenBank | genotype |

For SARS-CoV-2 and HIV-1 the accession lists are shared across all analysed
genes (one TSV per virus). For HCMV the genotype assignment depends on the
gene, so each gene has its own TSV.

## Reproducing the FASTA files

A helper script downloads the genomes and extracts each gene region using the
GenBank annotation already shipped under `data/<Virus>/<Gene>/`:

```bash
pip install -r requirements.txt           # includes biopython
python scripts/fetch_sequences.py --email you@example.com --virus hcmv
python scripts/fetch_sequences.py --email you@example.com --virus hiv1
python scripts/fetch_sequences.py --email you@example.com --virus sars_cov2
```

Outputs are written to `data/<Virus>/<Gene>/sequences.fasta` with the header
format expected by the pipeline:

```
>accession|Virus|Gene|class
```

NCBI rate limits apply (3 requests/sec without an API key, 10/sec with). For
the SARS-CoV-2 dataset (≈278k sequences) consider supplying an API key via
`--api-key` or, for production-scale runs, the NCBI `datasets` CLI or the
EBI ENA bulk download.

## Reference annotations

Each gene folder contains the GenBank annotation file used as the alignment
target:

| Virus | Reference accession | File |
|---|---|---|
| SARS-CoV-2 (Wuhan-Hu-1) | NC_045512 | `data/Severe_acute_respiratory_syndrome_coronavirus_2/<Gene>/NC_045512.gb` |
| HIV-1 (HXB2) | AF033819 | `data/Human_immunodeficiency_virus_1/<Gene>/AF033819.gb` |
| HCMV (Merlin) | FJ616285 | `data/Human_betaherpesvirus_5/<Gene>/FJ616285.gb` |

## Extraction date

Accession lists were extracted on 2026-02-27. Public databases evolve; some
accessions may have been updated, suppressed, or merged since this date.
