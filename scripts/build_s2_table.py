#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Assemble S2 Table, the scoring rubric of the protein characterization component.

S2 was the last supplementary workbook with no producer, and the only one carrying an Office
identity rather than openpyxl's: it had been edited by hand. That showed. Its category labels sat
in ten vertically merged cells in column A, and merged text is anchored at the top of its range,
so scrolling a few rows down left a reader looking at "3 terms | 15" with no way to tell which
category it belonged to. The blocks are written out row by row here instead. Nothing is merged,
the sheet carries a filter, and both the category and the field it is read from repeat on every
line, which is also what makes the file usable by a script rather than only by an eye.

The rubric is a literal in this file, deliberately not imported from `src/protein_score.py`.
The cross-check compares the two, and a builder that read the scorer would turn that comparison
into the scorer agreeing with itself. Published values belong on one side of a check and the code
on the other.

Usage
-----
    python scripts/build_s2_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl.styles import Alignment, Font

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = ROOT / "supplementary" / "S2_Table.xlsx"
SHEET = "Protein characterization score"

# category -> (UniProt field or cross-reference, [(level, points), ...]).
# The first level of each category is its ceiling, and the ten ceilings sum to 100, which is the
# denominator of Eq (10). The cross-check holds all three of those to the scorer and the equation.
RUBRIC: dict[str, tuple[str, list[tuple[str, float]]]] = {
    "Protein existence": ("proteinExistence", [
        ("Experimental evidence at protein level", 15),
        ("Experimental evidence at transcript level", 10),
        ("Inferred from homology", 5),
        ("Predicted", 0),
        ("Uncertain", 0),
    ]),
    "Molecular function": ("GO cross-references, F: terms", [
        ("4+ terms", 20), ("3 terms", 15), ("2 terms", 10), ("1 term", 5), ("0 terms", 0),
    ]),
    "Biological processes": ("GO cross-references, P: terms; Reactome cross-references", [
        ("4+ terms", 20), ("3 terms", 15), ("2 terms", 10), ("1 term", 5), ("0 terms", 0),
    ]),
    "Cellular component": ("GO cross-references, C: terms", [
        ("3+ terms", 5), ("2 terms", 3.5), ("1 term", 2), ("0 terms", 0),
    ]),
    "Structural/functional annotations": ("features, 14 types listed in the note below", [
        ("7+ terms", 10), ("5-6 terms", 7.5), ("3-4 terms", 5), ("1-2 terms", 2.5),
        ("0 terms", 0),
    ]),
    "Protein-protein interactions": (
        "IntAct, BioGRID, STRING, ComplexPortal cross-references",
        [("Present", 5), ("Absent", 0)]),
    "Post-translational modifications": ("keywords, category PTM only",
                                        [("Present", 5), ("Absent", 0)]),
    "3D structures": ("PDB cross-references", [("Present", 5), ("Absent", 0)]),
    "Drug target interactions": ("DrugBank, ChEMBL, BindingDB cross-references",
                                 [("Present", 2.5), ("Absent", 0)]),
    "Literature references": ("references", [
        ("8+ references", 12.5), ("6-7 references", 10), ("4-5 references", 7.5),
        ("2-3 references", 5), ("0-1 references", 0),
    ]),
}

NOTE = (
    "Note: structural/functional annotations count UniProt features of these 14 types only: "
    "topological domain, transmembrane, intramembrane, domain, repeat, zinc finger, dna binding, "
    "region, coiled coil, motif, compositional bias, active site, binding site, site."
)

COLUMNS = ["Evaluation Criteria", "Level/Category", "Points", "UniProt field or cross-reference"]
WIDTHS = {"A": 31, "B": 41, "C": 13, "D": 58}


def build() -> pd.DataFrame:
    """One row per level, with the category and its source field repeated on each."""
    rows = []
    for category, (field, levels) in RUBRIC.items():
        for level, points in levels:
            rows.append({COLUMNS[0]: category, COLUMNS[1]: level,
                         COLUMNS[2]: points, COLUMNS[3]: field})
    return pd.DataFrame(rows, columns=COLUMNS)


def main() -> int:
    """Write the workbook and report the totals a reader can check against the manuscript."""
    table = build()
    ceilings = {category: levels[0][1] for category, (_, levels) in RUBRIC.items()}
    total = sum(ceilings.values())

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name=SHEET, index=False)
        sheet = writer.book[SHEET]

        note_row = len(table) + 3
        sheet.cell(note_row, 1, NOTE).alignment = Alignment(vertical="top", wrap_text=True)

        for letter, width in WIDTHS.items():
            sheet.column_dimensions[letter].width = width
        for cell in sheet[1]:
            cell.font = Font(bold=True)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
        for row in sheet.iter_rows(min_row=2, max_row=len(table) + 1):
            for cell in row:
                cell.alignment = Alignment(vertical="top", wrap_text=True)
        # Column A stays put as well as the header, so the category a level belongs to is on
        # screen whichever way the reader scrolls. The filter is what S4 already offers.
        sheet.freeze_panes = "B2"
        sheet.auto_filter.ref = f"A1:D{len(table) + 1}"

    print(f"{len(table)} rows, {len(RUBRIC)} categories, ceilings summing to {total:g}")
    for category, ceiling in ceilings.items():
        print(f"  {category:36} {ceiling:6g}")
    print(f"\nwritten to {OUTPUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
