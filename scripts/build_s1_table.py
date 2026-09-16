"""
Assemble S1 Table, the complete list of substitution matrices evaluated for the mutational
component.

Answers Reviewer #2 comment 6, which asks which 30 AAindex matrices were tested. The ranking
file also contains MIYATA_EVO and an internal working copy of it, neither of which belongs in
the list of candidate matrices, so both are separated out here.

Usage
-----
    python scripts/build_s1_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from workbook_format import format_sheets

HERE = Path(__file__).resolve().parent
RANKING = (HERE.parent / "results" / "matrix_evaluation"
           / "substitution_matrices_ranking_complete.csv")
OUTPUT = HERE.parent / "supplementary" / "S1_Table.xlsx"

# Produced by the genetic algorithm, not an AAindex candidate.
OPTIMIZED = "MIYATA_EVO"
# Internal working copy kept by the optimization notebook, not a distinct matrix.
INTERNAL = "MIYATA_EVO_columnwise_backup"

FAMILIES = {
    "Physicochemical": {"MIYATA", "GRANTHAM", "SNEATH", "RISLER", "MCLACHLAN", "RAO",
                        "LEVIN", "FENG", "JOHNSON", "STR"},
    "Evolutionary substitution or replacement rate": {
        "DAYHOFF", "PAM30", "PAM70", "PAM120", "PAM250", "JONES", "WAG", "GONNET1992",
        "BENNER6", "BENNER22", "BENNER74", "VTML160", "VTML200", "VTML250", "MDM78"},
    "Alignment scoring": {"BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90",
                          "BLASTP"},
    "Genetic code": {"GENETIC"},
    # Kept apart from the general evolutionary models, which is the distinction reviewer 2.10
    # is drawing: these were estimated from viral alignments rather than from a broad set of
    # protein families. rtREV, HIVb and HIVw come from reverse-transcribing RNA viruses.
    "Virus-derived replacement rate": {"HIVb", "HIVw", "rtREV", "FLU", "FLAVI"},
}


def family_of(name: str) -> str:
    for family, members in FAMILIES.items():
        if name in members:
            return family
    return "Unclassified"


# Where each matrix's values were matched, and the source they were verified against.
#
# The column says "matched to", not "read from". No acquisition script survives for the thirty
# matrices of the base panel: they were in the repository before the JSON files were tracked, so
# the route by which their numbers arrived is unrecoverable. What is measurable is where the same
# numbers are found today, and that is what this table reports. Agreement means numerical equality
# to 1e-9 over the 400 ordered amino-acid pairs, not identical files: ours are JSON, the sources are
# not. Ten of them agree with both Biopython and AAindex2, which no comparison of values can
# disambiguate; those say so.
#
# Until 28 August 2026 this dict held only the seven matrices added for comment 2.10 and every
# other matrix fell through `get(name, ("AAindex", ""))` to a silent "AAindex". That was wrong for
# twenty-three of them: they agree numerically with Biopython's distribution, whose data files cite
# the primary literature and not AAindex. The default is gone rather than corrected, so a matrix
# added without provenance now stops the build instead of being attributed to whatever the default
# happened to say.
SOURCES = {
    # Values found in Biopython's distribution and not in AAindex2.
    "BENNER6": ("Biopython", ""),
    "BENNER74": ("Biopython", ""),
    "BLASTP": ("Biopython", ""),
    "BLOSUM50": ("Biopython", ""),
    "BLOSUM90": ("Biopython", ""),
    "DAYHOFF": ("Biopython", ""),
    "JOHNSON": ("Biopython", ""),
    "JONES": ("Biopython", ""),
    "MCLACHLAN": ("Biopython", ""),
    "MDM78": ("Biopython", ""),
    "PAM250": ("Biopython", ""),
    "PAM30": ("Biopython", ""),
    "PAM70": ("Biopython", ""),
    # Values found in AAindex2 and in neither Biopython nor the spaln tables.
    "BLOSUM62": ("AAindex2 (HENS920102)", ""),
    "GRANTHAM": ("AAindex2 (GRAR740104)", ""),
    "MIYATA": ("AAindex2 (MIYT790101)", ""),
    "PAM120": ("AAindex2 (ALTS910101)", ""),
    # Values found in AAindex2 and in the spaln tables; the route cannot be told apart.
    "VTML160": ("AAindex2 (MUET020101) or Spaln distribution", ""),
    "VTML250": ("AAindex2 (MUET020102) or Spaln distribution", ""),
    # Values identical in Biopython and AAindex2; the route cannot be told apart.
    "BENNER22": ("Biopython or AAindex2", ""),
    "BLOSUM45": ("Biopython or AAindex2", ""),
    "BLOSUM80": ("Biopython or AAindex2", ""),
    "FENG": ("Biopython or AAindex2", ""),
    "GENETIC": ("Biopython or AAindex2", ""),
    "GONNET1992": ("Biopython or AAindex2", ""),
    "LEVIN": ("Biopython or AAindex2", ""),
    "RAO": ("Biopython or AAindex2", ""),
    "RISLER": ("Biopython or AAindex2", ""),
    "STR": ("Biopython or AAindex2", ""),
    # In neither; the values are those published by Sneath.
    "SNEATH": ("Sneath (1966), Table 2", ""),
    # The seven added for comment 2.10, where the file actually read is known and pinned.
    "VTML200": ("Spaln distribution", "SeqAn"),
    "WAG": ("Authors' file, EBI server", "PAML"),
    "HIVb": ("IQ-TREE", "RAxML"),
    "HIVw": ("IQ-TREE", "RAxML"),
    "rtREV": ("IQ-TREE", "RAxML"),
    "FLU": ("IQ-TREE", "RAxML"),
    "FLAVI": ("IQ-TREE", "Authors' distribution"),
}


def source_of(name: str) -> tuple[str, str]:
    """Provenance of one matrix. Unknown names stop the build rather than defaulting."""
    try:
        return SOURCES[name]
    except KeyError:
        raise SystemExit(
            f"{name}: no provenance recorded in SOURCES.\n"
            "  Add where its values were matched before regenerating S1. There is deliberately no\n"
            "  default: a silent one attributed twenty-three Biopython matrices to AAindex."
        ) from None


def main() -> None:
    ranking = pd.read_csv(RANKING)
    divergence = pd.read_csv(RANKING.parent / "divergence_sensitivity.csv")

    internal = ranking[ranking["Matrix"] == INTERNAL]
    candidates = ranking[~ranking["Matrix"].isin({OPTIMIZED, INTERNAL})].copy()
    optimized = ranking[ranking["Matrix"] == OPTIMIZED]

    candidates["Family"] = candidates["Matrix"].map(family_of)
    candidates["Values matched to"] = candidates["Matrix"].map(lambda name: source_of(name)[0])
    candidates["Verification source"] = candidates["Matrix"].map(lambda name: source_of(name)[1])
    candidates = candidates.sort_values("Composite", ascending=False)
    candidates.insert(0, "Rank among candidates", range(1, len(candidates) + 1))
    candidates = candidates.drop(columns=["Rank"])

    unclassified = candidates.loc[candidates["Family"] == "Unclassified", "Matrix"].tolist()
    if unclassified:
        print(f"warning: no family assigned to {unclassified}")

    # The manuscript describes a base panel of thirty plus the seven added for comment 2.10. The
    # split is checked, but no longer by counting an "AAindex" label: that label was the silent
    # default and counting it only confirmed the default had been applied thirty times.
    added = {"VTML200", "WAG", "HIVb", "HIVw", "rtREV", "FLU", "FLAVI"}
    base = candidates[~candidates["Matrix"].isin(added)]
    if len(base) != 30 or len(candidates) - len(base) != 7:
        raise SystemExit(
            f"{len(base)} base matrices and {len(candidates) - len(base)} added, "
            "the text says 30 and 7")
    print(f"candidate matrices evaluated : {len(candidates)} "
          f"({len(base)} base panel, {len(candidates) - len(base)} added for comment 2.10)")
    print("\nby provenance:")
    print(candidates["Values matched to"].value_counts().to_string())
    print(f"optimized matrix rows        : {len(optimized)}")
    print(f"internal working copies      : {len(internal)} (excluded)")
    print("\nby family:")
    print(candidates["Family"].value_counts().to_string())

    # A rate model has no score matrix until a divergence is chosen, so the sheet that
    # shows what that choice costs belongs beside the ranking it moves.
    divergence = divergence.rename(columns={
        "model": "Rate model",
        "divergence_delta": "Divergence (expected substitutions per site)",
        "is_reference_value": "Reference value",
        "composite": "Composite",
        "rank_overall": "Rank overall",
        "rank_among_candidates": "Rank among candidates",
        "distinct_values": "Distinct off-diagonal values",
        "spearman_vs_reference_delta": "Spearman against the reference divergence",
    })

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        candidates.to_excel(writer, sheet_name="1_Matrices evaluated", index=False)
        optimized.to_excel(writer, sheet_name="2_Optimized matrix", index=False)
        divergence.to_excel(writer, sheet_name="3_Divergence sensitivity", index=False)
        format_sheets(writer)

    print(f"\nWrote {OUTPUT}")


if __name__ == "__main__":
    main()
