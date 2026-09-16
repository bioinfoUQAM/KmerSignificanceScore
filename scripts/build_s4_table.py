"""
Assemble S4 Table, the prevalence-threshold sensitivity analysis (Reviewer #1, comment 2.2).

The reviewer asks for a sensitivity analysis over t in [0.1, 0.5], reporting score variance and
method stability, because the default t = 0.25 was raised to 0.33 for HIV-1. This table reports,
for every dataset and every threshold, the retained positions and variant k-mers, the spread of
the KSS distribution, and two stability measures against the threshold used in the manuscript:
the overlap of the top-ranked positions (top-12 and top-25) and the Spearman rank correlation of
the complete score vector.

Reads results/threshold_sensitivity/threshold_sensitivity.csv, produced by
threshold_sensitivity.py from the sweep directories. Run that first.

Usage
-----
    python scripts/build_s4_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / "results" / "threshold_sensitivity" / "threshold_sensitivity.csv"
OUTPUT = ROOT / "supplementary" / "S4_Table.xlsx"

COLUMNS = {
    "dataset": "Dataset",
    "threshold": "Threshold t",
    "is_published_value": "Used in manuscript",
    "positions": "Positions",
    "variants": "k-mers retained",
    "alternatives": "Of which alternative to the reference",
    "kss_mean": "KSS mean",
    "kss_sd": "KSS standard deviation",
    "kss_variance": "KSS variance",
    "selection_rule": "Published selection rule",
    "selection_size": "Positions in the published selection",
    "shared_selection": "Published selection shared with manuscript t",
    "jaccard_selection": "Jaccard, published selection",
    "shared_top12": "Top-12 shared with manuscript t",
    "jaccard_top12": "Jaccard, top-12",
    "shared_top25": "Top-25 shared with manuscript t",
    "jaccard_top25": "Jaccard, top-25",
    "spearman_full_vs_published": "Spearman of full ranking vs manuscript t",
}
ORDER = ["HIV-1", "HCMV", "SARS-CoV-2"]


def _format_sheet(sheet, table: pd.DataFrame) -> None:
    """Make the workbook readable: headers that wrap, columns that fit, a frozen header row.

    The column names carry the meaning here, and several are long enough that a viewer clips
    them to something ambiguous, "Top-12 shared with manuscript..." reading as if it were the
    count for every threshold at once.
    """
    header = Font(bold=True)
    wrapped = Alignment(wrap_text=True, vertical="bottom", horizontal="center")
    for cell in sheet[1]:
        cell.font = header
        cell.alignment = wrapped
    sheet.row_dimensions[1].height = 46

    for index, name in enumerate(table.columns, start=1):
        longest_value = max((len(str(v)) for v in table[name]), default=0)
        longest_word = max(len(word) for word in str(name).split())
        width = min(28, max(11, longest_value + 2, longest_word + 2))
        sheet.column_dimensions[get_column_letter(index)].width = width

    sheet.freeze_panes = "A2"
    sheet.auto_filter.ref = sheet.dimensions


def main() -> None:
    table = pd.read_csv(SOURCE)
    table["threshold"] = table["threshold"] / 100.0
    # The retained count includes each position's reference k-mer, which essentially always
    # meets the threshold. Reporting the alternatives separately is what shows the filter at
    # work, and it explains why the number of positions does not move with t.
    table["alternatives"] = table["variants"] - table["positions"]
    # yes/no rather than yes/blank: an empty cell reads as missing data in a viewer, and as
    # NaN through pandas, where the answer is simply "no".
    table["is_published_value"] = table["is_published_value"].map({True: "yes", False: "no"})
    table["dataset"] = pd.Categorical(table["dataset"], categories=ORDER, ordered=True)
    table = table.sort_values(["dataset", "threshold"])
    table = table[list(COLUMNS)].rename(columns=COLUMNS).round(4)

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name="Threshold sensitivity", index=False)
        _format_sheet(writer.sheets["Threshold sensitivity"], table)

    for dataset in ORDER:
        g = table[table["Dataset"] == dataset]
        print(f"{dataset}: Spearman {g['Spearman of full ranking vs manuscript t'].min():.4f}"
              f"-{g['Spearman of full ranking vs manuscript t'].max():.4f}, "
              f"KSS sd {g['KSS standard deviation'].min():.4f}-{g['KSS standard deviation'].max():.4f}")
    print(f"\nWrote {OUTPUT} ({len(table)} rows)")


if __name__ == "__main__":
    main()
