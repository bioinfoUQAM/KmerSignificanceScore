"""
Assemble S7 Table, the supporting workbook on selection overlap between methods.

Run `topk_overlap.py` first, which produces the CSV files consumed here.

Usage
-----
    python scripts/build_s7_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from workbook_format import format_sheets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results" / "statistical_comparison"
OUTPUT = ROOT / "supplementary" / "S7_Table.xlsx"


# Source column -> workbook header, one map per sheet. The harness imports these to compare
# each sheet with its CSV, so a header renamed here cannot go unnoticed there.
SUMMARY_COLUMNS = {"top_k": "Selection size (m)"}

BY_DATASET_COLUMNS = {
    "dataset_gene": "Evaluation unit",
    "top_k": "Selection size (m)",
    "baseline": "Baseline method",
    "shared_positions": "Positions shared with the component",
    "jaccard": "Jaccard index",
}

CORRELATION_COLUMNS = {
    "dataset_gene": "Evaluation unit",
    "baseline": "Baseline method",
    "n_positions": "Positions analysed",
    "spearman_rho": "Spearman rho with the component",
}

SHEETS = {
    "1_Jaccard summary": ("topk_jaccard_summary.csv", SUMMARY_COLUMNS),
    "2_Jaccard by dataset": ("topk_jaccard_by_dataset.csv", BY_DATASET_COLUMNS),
    "3_Rank correlations": ("score_rank_correlations.csv", CORRELATION_COLUMNS),
}


def main() -> None:
    summary = pd.read_csv(RESULTS / "topk_jaccard_summary.csv")
    by_dataset = pd.read_csv(RESULTS / "topk_jaccard_by_dataset.csv")
    correlations = pd.read_csv(RESULTS / "score_rank_correlations.csv")

    summary = summary.rename(columns=SUMMARY_COLUMNS)
    by_dataset = by_dataset.rename(columns=BY_DATASET_COLUMNS)
    correlations = correlations.rename(columns=CORRELATION_COLUMNS)

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="1_Jaccard summary", index=False)
        by_dataset.to_excel(writer, sheet_name="2_Jaccard by dataset", index=False)
        correlations.to_excel(writer, sheet_name="3_Rank correlations", index=False)
        format_sheets(writer)

    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
