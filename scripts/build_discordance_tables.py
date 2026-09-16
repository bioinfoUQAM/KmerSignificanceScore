"""S8 and S9 Tables, the extreme discordance cases (Reviewer #1, comment 2.5).

The two tables list the proteins where the protein component and the normalized UniProt
annotation score disagree most, and the manuscript reads their annotation content to say what
the component sees that the annotation score does not.

They had no generator. The November 2025 predecessor files of the current S8 and S9 still held
the pre-correction set of twenty positive cases, produced by a float comparison that admitted
proteins sitting exactly on the threshold; the strict decimal rule of
regenerate_protein_scores.py leaves eighteen. The published share of positive cases with a
resolved structure, 70%, is 14/20 on that superseded set and 12/18 on the correct one.

Selection is therefore stated once here and nowhere else:
  - positive cases: every protein with delta strictly greater than +0.30, which is all 18 of
    them, not a top-20 that cannot exist;
  - negative cases: the 20 most negative of the 1,703 with delta strictly below -0.30.

Annotation content is read from the same UniProt snapshot the scores were computed from, so
the tables cannot drift from the scores.

Usage
-----
    python scripts/build_discordance_tables.py
"""

from __future__ import annotations

import json
from decimal import Decimal
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
JSON_DIR = ROOT / "notebooks" / "uniprot_viral_proteins"
RESULTS = ROOT / "results" / "protein_score_validation"
# The two tables also ship under their publication names, beside the other
# supplementary files. Files carrying those names existed with no generator and
# went stale: they still held the twenty pre-correction positives in July 2026.
SUPPLEMENTARY = ROOT / "supplementary"
SCORES = RESULTS / "complete_protein_scores.csv"
THRESHOLD = Decimal("0.30")
NEGATIVE_CASES = 20

DELTA = "Delta_Protein_Component_minus_UniProt_norm"

# Copied from get_structural_features_score in src/protein_score.py, so the table counts what
# the score counts.
INTERACTION_DATABASES = {"IntAct", "BioGRID", "STRING", "ComplexPortal"}
STRUCTURAL_FEATURE_TYPES = {
    "topological_domain", "transmembrane", "intramembrane", "domain", "repeat", "zinc_finger",
    "dna_binding", "region", "coiled_coil", "motif", "compositional_bias", "active_site",
    "binding_site", "site",
}


def annotation_content(accession: str) -> dict:
    """The fields the manuscript comments on, read from the entry the score was computed from."""
    data = json.loads((JSON_DIR / f"{accession}.json").read_text(encoding="utf-8"))
    cross = data.get("uniProtKBCrossReferences", [])
    go = [x for x in cross if x.get("database") == "GO"]

    def aspect(letter: str) -> int:
        # Distinct lowercased terms, as the scorer counts them: an entry can cite the same GO
        # term through several sources, and counting the citations would inflate the column.
        return len({p["value"].lower() for x in go for p in x.get("properties", [])
                    if p.get("key") == "GoTerm" and p.get("value", "").startswith(f"{letter}:")})

    # The two structural columns follow src/protein_score.py rather than a definition of their
    # own: the same feature types it counts, and PTM read from the keyword category it reads.
    features = data.get("features", [])
    structural = [f for f in features
                  if f.get("type", "").lower().replace(" ", "_") in STRUCTURAL_FEATURE_TYPES]
    ptm = any(k.get("category", "").strip() == "PTM" for k in data.get("keywords", []))
    return {
        "Protein_Name": (data.get("proteinDescription", {}).get("recommendedName", {})
                         .get("fullName", {}).get("value", "")),
        "Organism": data.get("organism", {}).get("scientificName", ""),
        # UniProt writes "1: Evidence at protein level"; the level is the leading digit.
        "Protein_Existence": int(str(data.get("proteinExistence", "0"))[0]),
        "Num_References": len(data.get("references", [])),
        "GO_MF": aspect("F"),
        "GO_BP": aspect("P"),
        "GO_CC": aspect("C"),
        "Has_3D_Structure": any(x.get("database") == "PDB" for x in cross),
        "Has_PTM": ptm,
        "Has_Interactions": any(x.get("database") in INTERACTION_DATABASES for x in cross),
        "Num_Structural_Functional_Features": len(structural),
    }


def main() -> None:
    scores = pd.read_csv(SCORES)
    delta = scores[DELTA].map(lambda v: Decimal(str(v)))

    positive = scores[delta > THRESHOLD].sort_values(DELTA, ascending=False)
    negative = scores[delta < -THRESHOLD].sort_values(DELTA).head(NEGATIVE_CASES)

    SUPPLEMENTARY.mkdir(parents=True, exist_ok=True)
    for name, published, table in (("S2_positive_discordance.csv", "S8_Table.csv", positive),
                                   ("S3_negative_discordance.csv", "S9_Table.csv", negative)):
        enriched = table.join(pd.DataFrame(
            [annotation_content(a) for a in table["Accession"]], index=table.index))
        enriched.to_csv(RESULTS / name, index=False, encoding="utf-8")
        enriched.to_csv(SUPPLEMENTARY / published, index=False, encoding="utf-8")
        structures = int(enriched["Has_3D_Structure"].sum())
        print(f"{name}: {len(enriched)} proteins, {structures} with a resolved structure "
              f"({100 * structures / len(enriched):.1f}%), "
              f"protein existence 1 for {int((enriched['Protein_Existence'] == 1).sum())}")

    print(f"\nwritten to {RESULTS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
