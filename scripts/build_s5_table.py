"""
Assemble S5 Table, the supporting workbook for the statistical comparison of the
discriminative component.

Run `statistical_comparison.py` first, which produces the CSV files consumed here.

Usage
-----
    python scripts/build_s5_table.py
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from workbook_format import format_sheets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results" / "statistical_comparison"
OUTPUT = ROOT / "supplementary" / "S5_Table.xlsx"

METHOD_LABELS = {
    "KSS_discriminative": "Discriminative component",
    "Chi2": "Chi-squared",
    "MI": "Mutual information",
    "NMI": "NMI",
    "OddsRatio": "Odds ratio",
    "CramerV": "Cramer's V",
    "SymChi2": "Normalized symmetric chi-squared divergence",
    "ANOVA": "ANOVA",
    "Random": "Random baseline",
}

CONTRAST_COLUMNS = {
    "baseline": "Baseline method",
    "kss_mean": "Discriminative component mean F1",
    "baseline_mean": "Baseline mean F1",
    "mean_difference": "Mean paired difference",
    "se_naive": "SE (uncorrected)",
    "se_corrected": "SE (corrected)",
    "ci95_low": "95% CI lower",
    "ci95_high": "95% CI upper",
    "t_corrected": "t (corrected)",
    "p_corrected_t": "p (corrected t)",
    "p_corrected_t_holm": "p (Holm-adjusted)",
    "p_wilcoxon": "p (Wilcoxon, uncorrected)",
    "cohen_dz": "Cohen dz",
    "wins": "Iterations won",
    "ties": "Iterations tied",
    "losses": "Iterations lost",
    "n_pairs": "Paired observations",
    "equivalent_margin_0.01": "Equivalent (margin 0.01)",
    "equivalent_margin_0.02": "Equivalent (margin 0.02)",
}


# Source column -> workbook header, one map per sheet. The harness imports these to compare
# each sheet with the CSV it was written from, so a header renamed here cannot go unnoticed
# there. Sheets 1 and 2 share CONTRAST_COLUMNS above.
TOPK_COLUMNS = {
    "scope": "Scope",
    "top_k": "Selection size (m)",
    "kss_mean": "Discriminative component mean F1",
    "baseline_mean": "Chi-squared mean F1",
    "mean_difference": "Mean paired difference",
    "se_corrected": "SE (corrected)",
    "ci95_low": "95% CI lower",
    "ci95_high": "95% CI upper",
    "p_corrected_t": "p (corrected t)",
    "p_holm": "p (Holm-adjusted)",
}

RANKING_COLUMNS = {
    "method": "Method",
    "mean_rank": "Mean rank",
    "mean_score": "Mean F1",
    "rank_gap_to_best": "Rank gap to best method",
}

CONVERGENCE_COLUMNS = {
    "n_iterations": "Iterations used",
    "mean_difference": "Mean paired difference (discriminative component - chi-squared)",
    "se_naive": "SE (uncorrected)",
    "se_corrected": "SE (corrected)",
    "p_naive": "p (uncorrected)",
    "p_corrected": "p (corrected)",
}


def format_contrasts(table: pd.DataFrame) -> pd.DataFrame:
    """Select and rename the reporting columns of a contrast table."""
    columns = [c for c in CONTRAST_COLUMNS if c in table.columns]
    formatted = table[columns].rename(columns=CONTRAST_COLUMNS)
    formatted["Baseline method"] = formatted["Baseline method"].map(
        lambda m: METHOD_LABELS.get(m, m)
    )
    return formatted


def main() -> None:
    contrasts = pd.read_csv(RESULTS / "paired_contrasts_f1.csv")
    topk = pd.read_csv(RESULTS / "contrast_by_topk_f1.csv")
    ranking = pd.read_csv(RESULTS / "descriptive_mean_ranks_f1.csv")
    convergence = pd.read_csv(RESULTS / "iteration_convergence_f1.csv")
    with open(RESULTS / "descriptive_mean_ranks_summary.json", encoding="utf-8") as handle:
        ranks_info = json.load(handle)

    is_global = contrasts["scope"].str.startswith("Global")
    global_contrasts = format_contrasts(contrasts[is_global])

    per_dataset = contrasts[~is_global].copy()
    per_dataset = format_contrasts(per_dataset).assign(
        **{"Evaluation unit": contrasts.loc[~is_global, "scope"].to_numpy()}
    )
    per_dataset = per_dataset[
        ["Evaluation unit"] + [c for c in per_dataset.columns if c != "Evaluation unit"]
    ].sort_values(["Evaluation unit", "Baseline method"])

    # The averaged stratification and the same contrast inside each unit belong on one sheet:
    # reading the first alone is what let the manuscript conclude that no selection size carried
    # the US28 difference, when averaging the five units is what hides it.
    topk_by_unit = pd.read_csv(RESULTS / "contrast_by_topk_unit_f1.csv")
    topk = pd.concat([
        topk.assign(scope="All units, averaged"),
        topk_by_unit.rename(columns={"dataset_gene": "scope"}),
    ], ignore_index=True)

    topk_out = topk.rename(columns=TOPK_COLUMNS).drop(
        columns=["baseline", "se_naive", "t_corrected"], errors="ignore")
    topk_out = topk_out[["Scope"] + [c for c in topk_out.columns if c != "Scope"]]

    ranking_out = ranking.rename(columns=RANKING_COLUMNS)
    ranking_out["Method"] = ranking_out["Method"].map(lambda m: METHOD_LABELS.get(m, m))
    ranking_out.loc[len(ranking_out)] = [""] * ranking_out.shape[1]
    ranking_out.loc[len(ranking_out)] = [
        f"Mean ranks over {ranks_info['n_blocks']} blocks and "
        f"{ranks_info['n_methods']} methods. The blocks are nested, the nine values of m within "
        f"an evaluation unit selecting from the same sequences, so this ranking is descriptive "
        f"and nothing inferential is computed from it.",
        "", "", "",
    ]

    convergence_out = convergence.rename(columns=CONVERGENCE_COLUMNS)

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        global_contrasts.to_excel(writer, sheet_name="1_Global contrasts", index=False)
        per_dataset.to_excel(writer, sheet_name="2_Per dataset", index=False)
        topk_out.to_excel(writer, sheet_name="3_By top-m", index=False)
        ranking_out.to_excel(writer, sheet_name="4_Descriptive ranks", index=False)
        convergence_out.to_excel(writer, sheet_name="5_Iteration convergence", index=False)
        format_sheets(writer)

    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
