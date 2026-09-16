"""
Computational cost of the KSS scoring components (Reviewer #1, comment 2.7).

    "The manuscript lacks a comparative analysis of computational complexity (O-notation) and
    execution time. Given that the chi-squared method performed equally well in several
    scenarios, the authors must address the trade-off between the computational cost of the KSS
    pipeline and simpler, traditional feature selection algorithms."

Two facts frame the answer. First, the preprocessing that turns raw sequences into the
position-by-position k-mer matrix (translation, alignment, k-mer extraction) is paid once and is
shared by every feature selection method that operates on that matrix, chi-squared included, so
it is not a cost that distinguishes KSS. Second, the differential cost is the scoring step, and
this script measures it on the real per-position matrices: the KSS discriminative score against
chi-squared, mutual information, and the mutational lookup, over identical inputs.

Each per-position score is a reduction of the same k-mer-by-class contingency table, so as an
algorithm the discriminative component is the same order as chi-squared: one pass to build the
table and O(K C) to reduce it, for N sequences, K distinct k-mers and C classes. The reference
implementation in src/discriminative_score.py does not reach that bound, and the manuscript says
so: it deduplicates the k-mer configurations with a row sort, then takes one pass over the labels
per distinct k-mer and per class, which costs O(NV log N + N(V+K)) and dominates the O(V K)
reduction at the small K these data carry. The mutational component is a
bounded number of dictionary lookups in the scaled substitution matrix, and the protein component
is a per-gene constant computed once from UniProt. This script reports the measured constants
behind those statements.

The claim that preprocessing dominates is measured rather than asserted (--stages), by timing the
two shipped entry points on a real gene: kanalyzer.analyze_records for translation, alignment and
k-mer extraction, then kss.compile_results and kss.compute_kss_scores for the whole scoring layer.
Nothing is written back to data/: the saving step of main.py is not called.

Usage
-----
    python scripts/computational_cost.py            # per-position scoring cost, one row per gene
    python scripts/computational_cost.py --stages   # preprocessing against scoring (needs the
                                                      # raw sequences.fasta of the gene measured)
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
import time
from pathlib import Path

# Same guard as main.py: the pipeline prints check marks, which a Windows console or a redirected
# stdout encodes as cp1252 and refuses.
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

import numpy as np
import pandas as pd
from sklearn.feature_selection import chi2, mutual_info_classif

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

# Package form, as main.py uses it: kanalyzer imports its siblings relatively.
from src import kanalyzer, kss, pipeline, protein_score  # noqa: E402
from src.discriminative_score import get_discriminative_score  # noqa: E402
from src.mutation_score import get_mutational_scores  # noqa: E402

BASE_PATH = ROOT / "data"
RESULTS = ROOT / "results" / "computational_cost"

# Every analysed gene. Until 25 August 2026 the four large SARS-CoV-2 genes were skipped to keep
# the run short, which left the highest sequence count measured only at the lowest k-mer diversity
# of the panel: E carries 278,738 sequences and 1.16 distinct k-mers per position. Cost grows with
# the product of the two, so that corner said nothing about the regime it was meant to cover.
#
# The raw per-sequence results are 2.3 GB for ORF1ab and 2.6 GB for S. The genes are measured one
# per process (--gene), so the peak is one gene rather than their sum, and each partial CSV is
# merged afterwards (--merge).
#
# One row per gene, not per dataset. Pooling positions across genes attributes a single sequence
# count to a mean taken over genes that do not share one: HCMV UL55, UL73 and US28 carry 399, 646
# and 443 sequences and 4, 8 and 6 classes. The three HIV-1 genes do share their 12,223 sequences,
# so their rows are directly comparable to each other.
MEASURE = {
    "Human_betaherpesvirus_5": ["UL55", "UL73", "US28"],
    "Human_immunodeficiency_virus_1": ["gag", "pol", "env"],
    "Severe_acute_respiratory_syndrome_coronavirus_2": ["E", "M", "N", "S", "ORF1ab"],
}
LABELS = {
    "Human_betaherpesvirus_5": "HCMV",
    "Human_immunodeficiency_virus_1": "HIV-1",
    "Severe_acute_respiratory_syndrome_coronavirus_2": "SARS-CoV-2",
}
REPEATS = 3
MUTATIONAL_MATRIX = "MIYATA_EVO"

# Genes used for the preprocessing-against-scoring split, at the two extremes of dataset size.
# They need their raw sequences.fasta, which the repository does not ship.
STAGES = [
    ("Human_betaherpesvirus_5", "US28"),
    ("Severe_acute_respiratory_syndrome_coronavirus_2", "E"),
]


def load_gene(dataset: str, gene: str) -> dict:
    """Per-position feature matrices and amino acid changes for one gene (real data)."""
    gene_path = BASE_PATH / dataset / gene / "results"
    compiled = json.loads((gene_path / f"{gene}_compiled_results.json").read_text())[gene]
    results = json.loads((gene_path / f"{gene}_results.json").read_text())
    gene_raw = results["genes"][gene]

    X, kmers_variants, labels, _ = kss.build_feature_matrix(gene_raw["sequences"], compiled)
    y = np.array(labels)

    positions = []
    for pos_str, record in compiled.items():
        if pos_str not in kmers_variants:
            continue
        cols = kmers_variants[pos_str]["col_indices"]
        if not len(cols):
            continue
        aa_changes = {}
        for var in record.get("alts", {}).values():
            for mutation in var.get("amino_acid_changes", {}):
                aa_changes[mutation] = 0
        positions.append((X[:, cols], y, aa_changes))
    return {"positions": positions, "n_sequences": len(y), "n_classes": len(np.unique(y))}


def timed(fn, positions) -> tuple[float, float, float]:
    """Median, minimum and maximum wall-clock seconds to score every position once.

    One uncounted warm-up pass first, so the first run does not carry the import, allocation and
    cache costs the later ones do not pay. The manuscript reports observed execution time and a
    practical trade-off rather than a noise floor, so the median leads and the range is published
    beside it; the best of three, used until 25 August 2026, gave neither.
    """
    for args in positions:
        fn(*args)
    elapsed = []
    for _ in range(REPEATS):
        start = time.perf_counter()
        for args in positions:
            fn(*args)
        elapsed.append(time.perf_counter() - start)
    return float(np.median(elapsed)), min(elapsed), max(elapsed)


def measure_stages(dataset: str, gene: str) -> dict:
    """Preprocessing against scoring for one gene, through the shipped entry points.

    The protein component resolves its UniProt accession through the REST API on every run, only
    the entry JSON being cached, so leaving it inside the scoring stage would time a network round
    trip as if it were computation and would make the split depend on the connection. It is timed
    on its own and held out of the scoring stage, which is then local computation throughout.
    """
    dataset_name, info, parameters = pipeline.load_config(
        str(BASE_PATH / dataset / "config.yaml"))
    info = dict(info, cds_selection=gene)

    start = time.perf_counter()
    results = kanalyzer.analyze_records(info, parameters)
    preprocessing = time.perf_counter() - start

    gene_dir = BASE_PATH / dataset / gene
    gb_path = next(p for p in sorted(gene_dir.iterdir()) if p.suffix.lower() == ".gb")
    taxon_id = kss.get_taxon_id(str(gb_path))
    start = time.perf_counter()
    score = protein_score.get_protein_score(taxon_id, gene, verbose=False)
    uniprot = time.perf_counter() - start

    original = protein_score.get_protein_score
    protein_score.get_protein_score = lambda *args, **kwargs: score
    try:
        start = time.perf_counter()
        compiled = kss.compile_results(results, parameters, target_gene=gene)
        compiled, _ = kss.compute_kss_scores(compiled, results, info, parameters,
                                             target_gene=gene)
        scoring = time.perf_counter() - start
    finally:
        protein_score.get_protein_score = original

    n_sequences = len(results["genes"][gene]["sequences"])
    return {
        "Dataset": LABELS[dataset],
        "Gene": gene,
        "Sequences": n_sequences,
        "Positions": len(compiled[gene]),
        "Preprocessing (s)": round(preprocessing, 2),
        "Scoring (s)": round(scoring, 2),
        "UniProt retrieval (s)": round(uniprot, 2),
        "Preprocessing share": round(preprocessing / (preprocessing + scoring), 4),
    }


def run_stages(workers: int) -> None:
    """Time the two stages on each gene whose raw sequences are available locally.

    The worker count is fixed for the run and written into the table. It does not make the
    timings deterministic, machine load still moves them, but it fixes the resources so the two
    stages are compared under the same conditions and so a later run can be set up the same way.
    The published split was measured with a single worker, before the alignment step gained
    multiprocessing, and reruns need that value to be comparable with it.
    """
    RESULTS.mkdir(parents=True, exist_ok=True)
    previous = os.environ.get("KSS_WORKERS")
    os.environ["KSS_WORKERS"] = str(workers)
    try:
        _run_stages(workers)
    finally:
        if previous is None:
            os.environ.pop("KSS_WORKERS", None)
        else:
            os.environ["KSS_WORKERS"] = previous


def _run_stages(workers: int) -> None:
    rows = []
    for dataset, gene in STAGES:
        fasta = BASE_PATH / dataset / gene / "sequences.fasta"
        if not fasta.exists():
            print(f"skipped {LABELS[dataset]} {gene}: {fasta.relative_to(ROOT)} not present; "
                  f"regenerate it from data/accessions/ to reproduce this measurement")
            continue
        row = measure_stages(dataset, gene)
        row["Workers"] = workers
        rows.append(row)
        print(f"{row['Dataset']} {row['Gene']}: N={row['Sequences']}, "
              f"preprocessing {row['Preprocessing (s)']} s, scoring {row['Scoring (s)']} s, "
              f"preprocessing share {row['Preprocessing share']:.1%} "
              f"(UniProt round trip {row['UniProt retrieval (s)']} s, held out)")
        # Written after each gene: the large one takes long enough that a partial table is worth
        # keeping.
        pd.DataFrame(rows).to_csv(RESULTS / "pipeline_stage_cost.csv", index=False)

    if rows:
        print(f"written to {(RESULTS / 'pipeline_stage_cost.csv').relative_to(ROOT)}")


def main(only_gene: str | None = None) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)

    scorers = {
        "KSS discriminative": lambda X, y, _: get_discriminative_score(X, y),
        "Chi-squared": lambda X, y, _: np.sum(chi2(X, y)[0]),
        "Mutual information": lambda X, y, _: np.sum(
            mutual_info_classif(X, y, discrete_features=True)),
        "Mutational lookup": lambda X, y, aa: get_mutational_scores(aa, MUTATIONAL_MATRIX),
    }

    rows = []
    for dataset, genes in MEASURE.items():
        for gene in genes:
            if only_gene and gene != only_gene:
                continue
            data = load_gene(dataset, gene)
            positions = data["positions"]
            n_seq, n_cls = data["n_sequences"], data["n_classes"]
            kmers = np.mean([Xy[0].shape[1] for Xy in positions])

            timings = {}
            for name, fn in scorers.items():
                # Mutational needs the aa dict; the others ignore it.
                timings[name] = timed(fn, positions)

            median = {name: t[0] for name, t in timings.items()}
            per_pos_us = {name: t / len(positions) * 1e6 for name, t in median.items()}
            spread_us = {name: (timings[name][1] / len(positions) * 1e6,
                                timings[name][2] / len(positions) * 1e6) for name in scorers}
            ratio = median["KSS discriminative"] / median["Chi-squared"]
            rows.append({
                "Dataset": LABELS[dataset],
                "Gene": gene,
                "Positions": len(positions),
                "Sequences": n_seq,
                "Classes": n_cls,
                "Mean k-mers per position": round(kmers, 2),
                **{f"{name} (us/position)": round(per_pos_us[name], 2) for name in scorers},
                **{f"{name} (us/position, min)": round(spread_us[name][0], 2) for name in scorers},
                **{f"{name} (us/position, max)": round(spread_us[name][1], 2) for name in scorers},
                "KSS/chi-squared ratio": round(ratio, 2),
                "Timed repeats": REPEATS,
            })
            print(f"{LABELS[dataset]} {gene}: {len(positions)} positions, N={n_seq}, "
                  f"C={n_cls}, K~{kmers:.1f}")
            for name in scorers:
                print(f"    {name:22s} {per_pos_us[name]:8.2f} us/position "
                      f"[{spread_us[name][0]:.2f}, {spread_us[name][1]:.2f}]")
            print(f"    KSS discriminative / chi-squared = {ratio:.2f}x\n")

    table = pd.DataFrame(rows)
    if only_gene:
        target = RESULTS / f"computational_cost_{only_gene}.csv"
    else:
        target = RESULTS / "computational_cost.csv"
    table.to_csv(target, index=False)
    print(f"written to {target.relative_to(ROOT)}")


def merge_partials() -> None:
    """Assemble the per-gene CSVs into the published table, in the order MEASURE declares."""
    order = [(LABELS[d], g) for d, genes in MEASURE.items() for g in genes]
    frames = []
    missing = []
    for _, gene in order:
        partial = RESULTS / f"computational_cost_{gene}.csv"
        if partial.exists():
            frames.append(pd.read_csv(partial))
        else:
            missing.append(gene)
    if missing:
        raise SystemExit(f"no partial CSV for: {', '.join(missing)}")
    table = pd.concat(frames, ignore_index=True)
    table.to_csv(RESULTS / "computational_cost.csv", index=False)
    print(f"merged {len(frames)} genes into computational_cost.csv")


if __name__ == "__main__":
    # A real parser rather than a membership test on sys.argv: the test ran the whole
    # measurement, minutes of it, for any argument including --help, which is how a check that
    # a moved script still starts ended up rewriting the published timings.
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stages", action="store_true",
                        help="time preprocessing against scoring instead of the per-position "
                             "cost; needs the raw sequences.fasta of the genes measured")
    parser.add_argument("--workers", type=int, default=1,
                        help="worker processes for the alignment step of --stages, written into "
                             "the table; 1 reproduces the published measurement, which predates "
                             "the multiprocessing of that step")
    parser.add_argument("--gene", default=None,
                        help="measure one gene and write computational_cost_<gene>.csv; the raw "
                             "results of the large SARS-CoV-2 genes are gigabytes, so one process "
                             "per gene keeps the peak at one of them rather than their sum")
    parser.add_argument("--merge", action="store_true",
                        help="assemble the per-gene CSVs into computational_cost.csv")
    arguments = parser.parse_args()
    if arguments.merge:
        merge_partials()
    elif arguments.stages:
        run_stages(arguments.workers)
    else:
        main(arguments.gene)
