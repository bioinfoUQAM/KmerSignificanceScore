"""
Assemble S6 Table, the COV2Var reference-layer matrix (Reviewer #1, comment 2.11).

The reviewer asks for false positive and false negative rates against a large-scale empirical
database. No layer of COV2Var answers that on its own, so the manuscript reports every
admissible layer at every selection size and this table carries the detail the main text
summarizes. Sheet 1 is the matrix. Sheet 2 records which layers were admitted and which were
excluded, with the reason, so that the choice of reference is auditable rather than implied.

Reads results/functional_benchmark/reference_matrix.csv, produced by
reference_matrix.py from the COV2Var download. Run that first.

Usage
-----
    python scripts/build_s6_table.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from workbook_format import format_sheets

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "results" / "functional_benchmark" / "reference_matrix.csv"
OUTPUT = HERE.parent / "supplementary" / "S6_Table.xlsx"

# The workbook headers say m, as the manuscript does: k is the k-mer length everywhere in this
# work, and the number of positions a method may select needed its own symbol. The CSV keeps
# the shorter column ids, where the code makes the distinction by name (K_NT for the length).
COLUMNS = {
    "reference": "Reference layer",
    "nature": "What the layer is",
    "admitted": "Admitted as evidence",
    "universe": "Evaluable k-mers",
    "prevalence": "Prevalence of positives",
    "expected_flat": "Positives expected at m, random draw",
    "expected_gene_matched": "Positives expected at m, gene-matched draw",
    "k": "Selection size m",
    "TP": "True positives",
    "FP": "False positives",
    "FN": "False negatives",
    "TN": "True negatives",
    "precision_at_k": "Precision at m",
    "fdr_at_k": "False discovery rate at m",
    "fpr_at_k": "False positive rate at m",
    "recall_at_k": "Recall at m",
    "fnr_at_k": "False negative rate at m",
    "enrichment": "Enrichment over prevalence",
    # Both probabilities are nominal: the hypergeometric null treats the evaluable k-mers as
    # exchangeable, and the selection is spatially aggregated, so the column heading says under
    # what assumption the number holds rather than leaving a reader to assume it is calibrated.
    "p_hypergeometric": "P(at least this many positives by chance), nominal under an "
                        "exchangeable-k-mer null",
    "p_perfect_by_chance": "P(all m positive by chance), nominal under an exchangeable-k-mer null",
    "perfect_prefix": "Longest prefix with no false positive",
    "auroc": "AUROC over the full ranking",
    "auroc_p": "AUROC p-value",
}

# Kept in the table itself rather than only in the legend, because the exclusions are part of
# the result: a benchmark against a prediction would establish nothing, and a layer that never
# looked at a gene must not score that gene as negative.
# Six decimals keeps the small exact probabilities readable without inventing precision: the
# p-values span 1e-1 to 1e-300 and are therefore left in full. The harness reads this list, so
# a column that stops being rounded here stops being compared as rounded there.
ROUNDED = ["Prevalence of positives", "Positives expected at m, random draw",
           "Positives expected at m, gene-matched draw", "Precision at m",
           "False discovery rate at m", "False positive rate at m", "Recall at m",
           "False negative rate at m", "Enrichment over prevalence",
           "AUROC over the full ranking"]


ADMISSIBILITY = [
    ("COV2Var common-mutation membership", "admitted",
     "Frequency-defined membership of the 9,832 common mutations, of which 40.7% are "
     "synonymous. Reported as coverage: its prevalence is 91.5% and a random selection of "
     "twelve matches it in 34.3% of draws."),
    ("COV2Var selection, MEME", "admitted",
     "Phylogenetic test for episodic diversifying selection on observed sequences. Restricted "
     "to the k-mers the selection analysis reports on, since COV2Var distributes only the "
     "significant rows of this test and the sites it examined without calling are unknown."),
    ("COV2Var selection, FEL", "admitted",
     "Phylogenetic test for pervasive selection, two-sided, so it is restricted to the sites "
     "where the non-synonymous rate exceeds the synonymous one: 74.9% of its significant calls "
     "are purifying and belong to the opposite question. Same universe restriction as MEME."),
    ("COV2Var selection, FUBAR", "admitted",
     "Phylogenetic test for pervasive selection, a different inference procedure from FEL, "
     "taken on its diversifying posterior Prob[alpha<beta]. 48.5% of its calls are also FEL "
     "calls, so the two agree on less than half of what they find."),
    ("COV2Var selection, union and intersection", "admitted",
     "Sensitivity analyses over a sensitive and a strict definition of selection."),
    ("COV2Var DMS antibody escape", "admitted, restricted",
     "Experimental measurement, receptor-binding domain only. Restricted to k-mers carrying a "
     "measured residue and to substitutions, which is what the assay varies: 68 of the 189 "
     "residues between spike 334 and 522 were assayed, and unmeasured or non-substituted "
     "k-mers leave the denominator rather than counting as negatives."),
    ("COV2Var DMS ACE2 binding", "admitted, restricted",
     "Experimental measurement, same restriction. Each assay measured 68 residues but not the "
     "same 68: their union is 69 and their intersection 67, which is why the two evaluable "
     "universes differ by one k-mer."),
    ("Union of the admitted layers", "admitted",
     "Every admitted layer is keyed on the same 9,832 catalogued mutations, so this union "
     "reproduces the membership layer. The union excluding membership is reported beside it, "
     "and is the one that says how much of a selection this resource speaks to beyond the "
     "fact that a mutation is frequent."),
    ("COV2Var protein function (MutPred2)", "reported, not admitted",
     "Software prediction of functional impact, not a measurement, and our own mutational "
     "term is a prediction too, so agreement would establish that two predictors agree. "
     "Computed and reported because it is one of the few genome-wide layers that is "
     "selective rather than saturated, restricted to substitutions and to the residues it "
     "scored."),
    ("COV2Var protein stability (iMutant 2.0)", "reported, not admitted",
     "Software prediction of stability change, same treatment as MutPred2."),
    ("COV2Var antigenicity and immunogenicity (VaxiJen, IEDB)", "reported, not admitted",
     "Software prediction, and spike only: all 671 scored mutations are in S, so the "
     "evaluable universe is 231 k-mers rather than the genome."),
    ("COV2Var physicochemical parameters (ProtParam)", "excluded",
     "Calculated descriptors of the sequence, not measurements and not effect estimates."),
    ("COV2Var mutation and patient metadata (age, sex, status)", "not applicable",
     "Associates a mutation with host metadata, confounded by the epidemiology of the "
     "lineages carrying it and by sampling. It says nothing about whether a position "
     "deserves prioritization, so it is neither an admissible reference nor a rate we quote."),
    ("COV2Var co-mutation correlation", "not applicable",
     "Measures linkage between pairs of mutations within a lineage, which is the lineage "
     "haplotype structure the discriminative component already exploits, on the same corpus."),
    ("COV2Var animal host records", "not applicable",
     "Records the host a sequence came from, not a property of the position."),
    ("COV2Var literature curation (PubMed)", "not obtainable",
     "The 615 mutations reported in 2,587 papers are the closest thing the resource has to a "
     "relevance list, and it is the one category of the fifteen that the bulk download does "
     "not distribute: it exists per mutation on the web interface only."),
    ("Bloom and Neher fitness effects", "exploratory characterization; not used for validation",
     "Estimates evolutionary constraint from observed against expected counts, which is a "
     "different quantity from the one ranked here. Ranking its sites by effect magnitude "
     "coincides with ranking them by constraint, so the reference cannot separate whether a "
     "site has a measurable effect from whether it is constrained. Examined while candidate "
     "references were being surveyed and not used to validate the integrated score; the "
     "characterization it yields is written to "
     "results/functional_benchmark/fitness_characterisation.csv."),
]


def main() -> None:
    matrix = pd.read_csv(SOURCE)
    matrix = matrix[list(COLUMNS)].rename(columns=COLUMNS)
    for column in ROUNDED:
        matrix[column] = matrix[column].round(6)

    admissibility = pd.DataFrame(ADMISSIBILITY,
                                 columns=["Layer", "Status", "Reason"])

    with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
        matrix.to_excel(writer, sheet_name="Reference matrix", index=False)
        admissibility.to_excel(writer, sheet_name="Layer admissibility", index=False)
        format_sheets(writer)

    layers = matrix["Reference layer"].nunique()
    sizes = sorted(matrix["Selection size m"].unique())
    print(f"{layers} layers, selection sizes {sizes}")
    print(f"admissibility sheet: {len(admissibility)} layers, "
          f"{(admissibility['Status'].str.startswith('admitted')).sum()} admitted, "
          f"{(~admissibility['Status'].str.startswith('admitted')).sum()} not")
    print(f"\nWrote {OUTPUT} ({len(matrix)} rows)")


if __name__ == "__main__":
    main()
