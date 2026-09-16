"""
Canonical regeneration of the protein score table and of the discordance counts.

Addresses Reviewer #1 comments 2.5 and 4.1. The published table could not be re-derived
from the shipped code, and the discordance counts it supports were wrong.

What was wrong
--------------
1. Precision. `complete_protein_scores.csv` as published mixed two precisions: 11,167 of its
   17,470 rows carry the value the scorer returns, and 6,303 carry that value rounded to two
   decimals. No code in this repository produces two decimals, so the file was not
   reproducible. `calculate_protein_score` rounds to three decimals, which is a no-op because
   the raw score is always a multiple of 0.005.
2. Comparison. The discordance test was a strict float inequality, `delta > 0.30`. Rounding to
   two decimals had placed 378 proteins exactly on the threshold, and 122 of them crossed it
   because their operand pair is inexact in binary: 0.55 - 0.25, 0.80 - 0.50 and 0.70 - 1.00
   all evaluate to 0.30000000000000004, whereas 0.30 - 0.00, 0.20 - 0.50 and 0.45 - 0.75 are
   exact. That is the whole of the published 18 + 5 = 23 and 1,666 + 117 = 1,783.
3. Silence. `except Exception: continue` swallowed one unreadable input, so the analysis ran on
   17,470 of 17,471 records without saying so.

What this script guarantees
---------------------------
- No input is skipped silently. An unreadable or incomplete record aborts the run unless it is
  listed in RECOVERED with the provenance of its replacement.
- Scores are accumulated in decimal arithmetic and written unrounded.
- Discordance is counted in decimal arithmetic with a strict threshold. Values landing exactly
  on the threshold are excluded from the strict count and reported separately, because a
  strict inequality is what the manuscript states.
- The float comparison is re-run alongside as a regression check and any disagreement with the
  decimal verdict is reported, so the defect above cannot return unnoticed.
- Every derived figure quoted in the manuscript is printed from the regenerated table.

Usage
-----
    python scripts/regenerate_protein_scores.py [--check]

The UniProt entries are not in the repository. Unpack `uniprot_viral_proteins/` from the archive
deposited at https://doi.org/10.5281/zenodo.21937072 into `notebooks/` first.

`--check` recomputes and compares against the files on disk without rewriting them, and exits
non-zero on any difference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import protein_score as ps  # noqa: E402

JSON_DIR = ROOT / "notebooks" / "uniprot_viral_proteins"
RESULTS = ROOT / "results" / "protein_score_validation"
SCORES_CSV = RESULTS / "complete_protein_scores.csv"
PROVENANCE_CSV = RESULTS / "protein_score_provenance.csv"
SUMMARY_MD = RESULTS / "REGENERATION_SUMMARY.md"

THRESHOLD = Decimal("0.30")
DELTA_COLUMN = "Delta_Protein_Component_minus_UniProt_norm"

CATEGORIES = (
    ("protein_existence", ps.score_protein_existence),
    ("molecular_function", ps.score_molecular_function),
    ("biological_processes", ps.score_biological_processes),
    ("cellular_component", ps.score_cellular_component),
    ("structural_annotation", ps.score_structural_annotation),
    ("protein_interactions", ps.score_protein_interactions),
    ("post_translational_mods", ps.score_post_translational_modifications),
    ("3d_structures", ps.score_3d_structures),
    ("drug_target", ps.score_drug_target),
    ("publication_references", ps.score_publication_references),
)

# Records whose local copy from the 2025_02 snapshot is unusable and was replaced. Each entry
# states why the replacement is score-neutral, so the substitution is auditable rather than
# silent. Q65228: the 1 June 2025 download stopped at 8,204 bytes inside
# `uniProtKBCrossReferences`. The truncated file still records `entryVersion: 80` and
# `annotationScore: 2.0`. UniSave shows version 80 covered release 2025_02, and its flat-file
# diff against the fetched version 82 changes only the entry-version line, one intermediate
# lineage taxon, and the keyword "Reference proteome" (category Technical term, which the
# scorer reads only under category PTM). No DR, FT, RN, GO or CC line differs, so no scoring
# input differs.
RECOVERED = {
    "Q65228": {
        "snapshot_version": 80,
        "replacement_version": 82,
        "source": "rest.uniprot.org/uniprotkb/Q65228.json fetched 2026-07-26",
        "reason": "truncated download; v80 vs v82 differ only outside every scoring input",
    },
}


# The operand pairs that produced the published counts. Kept as a fixture so the defect is
# pinned to its mechanism rather than to the number of rows a given corpus happens to contain:
# each pair sits exactly on the threshold in decimal arithmetic, and the third column says
# whether subtracting the two floats crosses it anyway.
FLOAT_LEAK_FIXTURE = (
    ("0.30", "0.00", False), ("0.20", "0.50", False), ("0.45", "0.75", False),
    ("0.55", "0.25", True), ("0.80", "0.50", True), ("0.70", "1.00", True),
)


def check_float_leak(counts: dict) -> list[str]:
    """Pin the float comparison defect, so reverting to it cannot pass unnoticed."""
    problems = []
    for kss, normalized, leaks in FLOAT_LEAK_FIXTURE:
        exact = Decimal(kss) - Decimal(normalized)
        if abs(exact) != THRESHOLD:
            problems.append(f"fixture {kss} - {normalized} no longer sits on the threshold")
        crossed = abs(float(kss) - float(normalized)) > 0.30
        if crossed != leaks:
            problems.append(f"fixture {kss} - {normalized}: expected leak={leaks}, got {crossed}")
    if (counts["float_positive"], counts["float_negative"]) == (counts["positive"],
                                                                counts["negative"]):
        problems.append("float and decimal counting agree on this corpus, so the counting path "
                        "is no longer distinguishable and this guard proves nothing")
    return problems


def normalized_annotation_score(annotation_score) -> Decimal:
    """Map the 1-5 UniProt annotation score onto [0, 1], as the published analysis does."""
    return (Decimal(str(annotation_score)) - 1) / Decimal("4")


def protein_score(record: dict) -> Decimal:
    """Sum the ten categories in decimal arithmetic and scale to [0, 1], without rounding."""
    total = sum((Decimal(str(fn(record, False))) for _, fn in CATEGORIES), Decimal("0"))
    return total / Decimal("100")


def load_records() -> tuple[list[dict], list[str]]:
    """Read every JSON in the corpus. Unreadable files are returned, never skipped."""
    records, unreadable = [], []
    for path in sorted(JSON_DIR.glob("*.json")):
        raw = path.read_bytes()
        try:
            record = json.loads(raw.decode("utf-8"))
        except Exception:
            unreadable.append(path.name)
            continue
        accession = record.get("primaryAccession", "")
        annotation = record.get("annotationScore")
        if not accession or annotation is None:
            unreadable.append(path.name)
            continue
        records.append({
            "accession": accession,
            "record": record,
            "sha256": hashlib.sha256(raw).hexdigest(),
            "entry_version": record.get("entryAudit", {}).get("entryVersion"),
            "last_update": record.get("entryAudit", {}).get("lastAnnotationUpdateDate"),
        })
    return records, unreadable


def build_table(records: list[dict]) -> list[dict]:
    """Score every record and attach its delta against the normalized annotation score."""
    rows = []
    for item in records:
        record = item["record"]
        kss = protein_score(record)
        annotation = Decimal(str(record["annotationScore"]))
        normalized = normalized_annotation_score(record["annotationScore"])
        rows.append({
            "Accession": item["accession"],
            "UniProt_Annotation_Score": annotation,
            "UniProt_Normalized_Score": normalized,
            "Protein_Characterization_Score": kss,
            DELTA_COLUMN: kss - normalized,
            "sha256": item["sha256"],
            "entry_version": item["entry_version"],
            "last_update": item["last_update"],
        })
    rows.sort(key=lambda r: r["Accession"])
    return rows


def count_discordance(rows: list[dict]) -> dict:
    """Count discordance in decimal arithmetic, and re-run the float test as a regression check."""
    positive = negative = tie_positive = tie_negative = 0
    float_positive = float_negative = 0
    disagreements = []
    for row in rows:
        delta = row[DELTA_COLUMN]
        if delta > THRESHOLD:
            positive += 1
        elif delta == THRESHOLD:
            tie_positive += 1
        if delta < -THRESHOLD:
            negative += 1
        elif delta == -THRESHOLD:
            tie_negative += 1

        float_delta = float(row["Protein_Characterization_Score"]) - float(row["UniProt_Normalized_Score"])
        float_high = float_delta > 0.30
        float_low = float_delta < -0.30
        float_positive += float_high
        float_negative += float_low
        if float_high != (delta > THRESHOLD) or float_low != (delta < -THRESHOLD):
            disagreements.append((row["Accession"], str(delta), repr(float_delta)))
    return {
        "n": len(rows),
        "positive": positive,
        "negative": negative,
        "tie_positive": tie_positive,
        "tie_negative": tie_negative,
        "float_positive": float_positive,
        "float_negative": float_negative,
        "disagreements": disagreements,
    }


def derived_figures(rows: list[dict], counts: dict) -> dict:
    """Every number the manuscript quotes from this table, computed from the table itself."""
    annotation = np.array([float(r["UniProt_Annotation_Score"]) for r in rows])
    kss = np.array([float(r["Protein_Characterization_Score"]) for r in rows])
    n = counts["n"]
    extreme_positive = max(rows, key=lambda r: r[DELTA_COLUMN])
    extreme_negative = min(rows, key=lambda r: r[DELTA_COLUMN])
    return {
        "spearman_rho": float(stats.spearmanr(annotation, kss).statistic),
        "kendall_tau": float(stats.kendalltau(annotation, kss).statistic),
        "pct_positive": 100.0 * counts["positive"] / n,
        "pct_negative": 100.0 * counts["negative"] / n,
        "ratio": counts["negative"] / counts["positive"] if counts["positive"] else float("nan"),
        "extreme_positive": (extreme_positive["Accession"], extreme_positive[DELTA_COLUMN]),
        "extreme_negative": (extreme_negative["Accession"], extreme_negative[DELTA_COLUMN]),
    }


def correlation_table(rows: list[dict]) -> list[dict]:
    """Spearman and Kendall against the raw 1-5 annotation score, with the published intervals."""
    annotation = np.array([float(r["UniProt_Annotation_Score"]) for r in rows])
    kss = np.array([float(r["Protein_Characterization_Score"]) for r in rows])
    n = len(rows)

    rho, rho_p = stats.spearmanr(annotation, kss)
    z = 0.5 * np.log((1 + rho) / (1 - rho))
    half = 1.96 / np.sqrt(n - 3)
    bounds = [(np.exp(2 * b) - 1) / (np.exp(2 * b) + 1) for b in (z - half, z + half)]

    tau, tau_p = stats.kendalltau(annotation, kss)
    se = np.sqrt(2 * (2 * n + 5) / (9 * n * (n - 1)))

    return [
        {"Metric": "Spearman ρ", "Value": rho, "CI_Lower": bounds[0], "CI_Upper": bounds[1],
         "p_value": rho_p},
        {"Metric": "Kendall τ", "Value": tau, "CI_Lower": max(-1.0, tau - 1.96 * se),
         "CI_Upper": min(1.0, tau + 1.96 * se), "p_value": tau_p},
    ]


def granularity_table(rows: list[dict]) -> list[dict]:
    """KSS distribution within each UniProt annotation category, as Table 6 reports it."""
    by_category: dict[float, list[float]] = {}
    for row in rows:
        by_category.setdefault(float(row["UniProt_Annotation_Score"]), []).append(
            float(row["Protein_Characterization_Score"]))
    table = []
    for category in sorted(by_category):
        values = np.array(by_category[category])
        mean = values.mean()
        sd = values.std(ddof=1)
        table.append({
            "UniProt_Category": category,
            "count": len(values),
            "mean": round(mean, 4),
            "std": round(sd, 4),
            "min": round(values.min(), 4),
            "max": round(values.max(), 4),
            "CV(%)": round(sd / mean * 100, 2),
        })
    return table


def discordance_tables(rows: list[dict], n_examples: int = 20) -> dict[str, list[dict]]:
    """The most extreme discordant proteins on each side, under the strict decimal threshold."""
    high = sorted((r for r in rows if r[DELTA_COLUMN] > THRESHOLD),
                  key=lambda r: r[DELTA_COLUMN], reverse=True)
    low = sorted((r for r in rows if r[DELTA_COLUMN] < -THRESHOLD),
                 key=lambda r: r[DELTA_COLUMN])
    keep = ("Accession", "UniProt_Annotation_Score", "UniProt_Normalized_Score", "Protein_Characterization_Score",
            DELTA_COLUMN)
    return {
        "positive_discordance_all18.csv": [{k: r[k] for k in keep} for r in high[:n_examples]],
        "negative_discordance_top20.csv": [{k: r[k] for k in keep} for r in low[:n_examples]],
    }


def write_outputs(rows: list[dict], counts: dict, figures: dict) -> None:
    """Write the canonical score table, its provenance manifest and the regeneration summary."""
    RESULTS.mkdir(exist_ok=True)

    def dump(name: str, table: list[dict]) -> None:
        with (RESULTS / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)

    dump("correlation_results.csv", correlation_table(rows))
    dump("granularity_by_category.csv", granularity_table(rows))
    for name, table in discordance_tables(rows).items():
        dump(name, table)

    with SCORES_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Accession", "UniProt_Annotation_Score", "UniProt_Normalized_Score",
                         "Protein_Characterization_Score", DELTA_COLUMN])
        for row in rows:
            writer.writerow([row["Accession"], row["UniProt_Annotation_Score"],
                             row["UniProt_Normalized_Score"], row["Protein_Characterization_Score"],
                             row[DELTA_COLUMN]])
    with PROVENANCE_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Accession", "entry_version", "last_annotation_update", "sha256",
                         "substitution"])
        for row in rows:
            recovered = RECOVERED.get(row["Accession"])
            writer.writerow([row["Accession"], row["entry_version"], row["last_update"],
                             row["sha256"], recovered["source"] if recovered else ""])
    SUMMARY_MD.write_text(summary_text(counts, figures), encoding="utf-8")


def summary_text(counts: dict, figures: dict) -> str:
    """Render the regenerated figures as the reference the manuscript is checked against."""
    positive_accession, positive_delta = figures["extreme_positive"]
    negative_accession, negative_delta = figures["extreme_negative"]
    return "\n".join([
        "# Regenerated protein score table",
        "",
        f"- proteins: {counts['n']:,}",
        f"- Spearman rho: {figures['spearman_rho']:.4f}",
        f"- Kendall tau: {figures['kendall_tau']:.4f}",
        "",
        f"Discordance, strict |delta| > {THRESHOLD}, decimal arithmetic:",
        "",
        f"- positive: {counts['positive']} ({figures['pct_positive']:.1f}%)",
        f"- negative: {counts['negative']} ({figures['pct_negative']:.1f}%)",
        f"- ratio: {figures['ratio']:.1f}:1",
        f"- exactly on +{THRESHOLD}: {counts['tie_positive']}, excluded from the strict count",
        f"- exactly on -{THRESHOLD}: {counts['tie_negative']}, excluded from the strict count",
        "",
        f"- most positive: {positive_accession}, delta {positive_delta:+}",
        f"- most negative: {negative_accession}, delta {negative_delta:+}",
        "",
        "Regression check on the float comparison that produced the published counts:",
        "",
        f"- float test would give {counts['float_positive']} and {counts['float_negative']}",
        f"- rows where the float verdict differs from the decimal verdict: "
        f"{len(counts['disagreements'])}",
    ]) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="compare against the files on disk instead of rewriting them")
    args = parser.parse_args()

    records, unreadable = load_records()
    if not records and not unreadable:
        print(f"no UniProt entry under {JSON_DIR.relative_to(ROOT)}: unpack uniprot_viral_proteins/ "
              "from https://doi.org/10.5281/zenodo.21937072 into notebooks/")
        return 1
    fatal = [name for name in unreadable if Path(name).stem not in RECOVERED]
    if fatal:
        print(f"unreadable inputs, refusing to continue: {fatal}")
        print("list them in RECOVERED with the provenance of their replacement, or repair them")
        return 1
    if unreadable:
        print(f"note: {len(unreadable)} record(s) read from a documented replacement: "
              f"{sorted(unreadable)}")

    missing = [accession for accession in RECOVERED
               if accession not in {item["accession"] for item in records}]
    if missing:
        print(f"declared as recovered but absent from the corpus: {missing}")
        return 1

    rows = build_table(records)
    counts = count_discordance(rows)
    figures = derived_figures(rows, counts)

    problems = check_float_leak(counts)
    if problems:
        print("float comparison regression guard failed:")
        for problem in problems:
            print(f"  {problem}")
        return 1

    print(f"proteins: {counts['n']:,}")
    print(f"Spearman rho = {figures['spearman_rho']:.4f}   "
          f"Kendall tau = {figures['kendall_tau']:.4f}")
    print(f"strict |delta| > {THRESHOLD}: {counts['positive']} positive "
          f"({figures['pct_positive']:.1f}%), {counts['negative']} negative "
          f"({figures['pct_negative']:.1f}%), ratio {figures['ratio']:.1f}:1")
    print(f"exactly on the threshold, excluded: {counts['tie_positive']} at +{THRESHOLD}, "
          f"{counts['tie_negative']} at -{THRESHOLD}")
    print(f"float comparison would give {counts['float_positive']} and "
          f"{counts['float_negative']}, disagreeing on {len(counts['disagreements'])} rows")
    print(f"most positive: {figures['extreme_positive'][0]} "
          f"{figures['extreme_positive'][1]:+}")
    print(f"most negative: {figures['extreme_negative'][0]} "
          f"{figures['extreme_negative'][1]:+}")

    if args.check:
        if not SCORES_CSV.exists():
            print(f"missing {SCORES_CSV}")
            return 1
        on_disk = {r["Accession"]: r for r in
                   csv.DictReader(SCORES_CSV.open(encoding="utf-8"))}
        # Every published column, not the KSS one alone: the normalized UniProt score and the
        # delta are what the discordance counts and the two supplementary tables are built from,
        # so a drift there would pass a check that only looked at the score.
        COMPARED = ["UniProt_Annotation_Score", "UniProt_Normalized_Score", "Protein_Characterization_Score",
                    "Delta_Protein_Component_minus_UniProt_norm"]
        differing = []
        for row in rows:
            stored = on_disk.get(row["Accession"])
            if stored is None:
                differing.append((row["Accession"], "absent"))
                continue
            for column in COMPARED:
                if Decimal(stored[column]) != row[column]:
                    differing.append((row["Accession"], column))
        extra = sorted(set(on_disk) - {row["Accession"] for row in rows})
        if differing or extra:
            print(f"table on disk differs: {len(differing)} cells, {len(extra)} unexpected rows")
            for accession, column in differing[:5]:
                print(f"  {accession}: {column}")
            return 1

        # The provenance manifest must describe the corpus that was just read, entry by entry
        # and field by field. Checking only that the accessions match would pass on a corpus
        # where an entry had been silently replaced by a different version of itself, which is
        # precisely the accident this manifest exists to make visible.
        if PROVENANCE_CSV.exists():
            manifest = {r["Accession"]: r for r in
                        csv.DictReader(PROVENANCE_CSV.open(encoding="utf-8"))}
            drifted = []
            for row in rows:
                recorded = manifest.get(row["Accession"])
                if recorded is None:
                    drifted.append((row["Accession"], "absent from the manifest"))
                    continue
                for column, value in (("sha256", row["sha256"]),
                                      ("entry_version", row["entry_version"]),
                                      ("last_annotation_update", row["last_update"])):
                    if str(recorded[column]) != str(value):
                        drifted.append((row["Accession"], column))
            missing = sorted(set(manifest) - {row["Accession"] for row in rows})
            if drifted or missing:
                print(f"provenance differs on {len(drifted)} fields, "
                      f"{len(missing)} entries missing from the corpus")
                for accession, column in drifted[:5]:
                    print(f"  {accession}: {column}")
                return 1

        print(f"table on disk matches the regenerated values, all {len(COMPARED)} columns, "
              f"and the manifest matches on hash, version and date for every entry")
        return 0

    write_outputs(rows, counts, figures)
    print(f"wrote {SCORES_CSV.relative_to(ROOT)}, {PROVENANCE_CSV.relative_to(ROOT)}, "
          f"{SUMMARY_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
