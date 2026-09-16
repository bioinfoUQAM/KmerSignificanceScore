"""
Assemble S3 Table, the inventory of the classes used in each dataset.

Answers the first part of Reviewer #2 comment 1, which asks how each set of classes was
originally established. The table lists every class, its number of sequences, and the
nomenclature it belongs to.

Usage
-----
    python notebooks/build_s3_table.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# workbook_format lives in scripts/, which is not on the path for a notebook-side builder.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from workbook_format import format_sheets  # noqa: E402

HERE = Path(__file__).resolve().parent
ACCESSIONS = HERE.parent / "data" / "accessions"
OUTPUT = HERE.parent / "supplementary" / "S3_Table.xlsx"

DATASETS = [
    ("SARS-CoV-2", "genome-wide (5 genes analysed)", "sars_cov2.tsv", "variant",
     "Pango lineage labels"),
    ("HIV-1", "genome-wide (3 genes analysed)", "hiv1.tsv", "subtype",
     # The summary named subtypes and recombinant forms only, while group O is one of the 15
     # labels and sits outside group M. Sheet 2 always said so; sheet 1 now does too.
     "Los Alamos HIV Sequence Database subtype field: group M subtypes, circulating recombinant "
     "forms and group O"),
    ("HCMV", "UL55", "hcmv_UL55.tsv", "genotype", "gB (UL55) genotyping scheme"),
    ("HCMV", "UL73", "hcmv_UL73.tsv", "genotype", "gN (UL73) genotyping scheme"),
    ("HCMV", "US28", "hcmv_US28.tsv", "genotype", "US28 genotyping scheme"),
]

# Classes that are not simple monophyletic units, flagged so that the table documents the
# heterogeneity of the nomenclatures rather than hiding it.
RECOMBINANT_PREFIXES = ("01_", "02_", "07_", "08_", "20_", "63_", "85_")


def annotate(dataset: str, label: str) -> str:
    text = str(label)
    if dataset == "HIV-1":
        if text.startswith(RECOMBINANT_PREFIXES):
            return "Circulating recombinant form"
        if text == "O":
            return "Group O (outside group M)"
        return "Group M subtype"
    if dataset == "SARS-CoV-2":
        if text.startswith("X"):
            return "Recombinant lineage"
        if "_" in text:
            return "Two Pango lineages grouped under one label"
        return "Pango lineage"
    return "Locus-specific genotype"


def main() -> None:
    rows = []
    for dataset, scope, filename, column, scheme in DATASETS:
        table = pd.read_csv(ACCESSIONS / filename, sep="\t")
        counts = table[column].value_counts()
        for label, n in counts.items():
            rows.append({
                "Dataset": dataset,
                "Scope of the classification": scope,
                "Class label": label,
                "Sequences": int(n),
                "Nomenclature": scheme,
                "Nature of the class": annotate(dataset, label),
            })

    inventory = pd.DataFrame(rows)

    summary = (
        inventory.groupby(["Dataset", "Scope of the classification", "Nomenclature"])
        .agg(Classes=("Class label", "count"), Sequences=("Sequences", "sum"))
        .reset_index()
    )

    print(summary.to_string(index=False))
    print("\nnature of the classes:")
    print(inventory["Nature of the class"].value_counts().to_string())

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="1_Summary", index=False)
        inventory.to_excel(writer, sheet_name="2_Classes", index=False)
        format_sheets(writer)

    print(f"\nWrote {OUTPUT}")


if __name__ == "__main__":
    main()
