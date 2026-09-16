"""
The threshold sweep for SARS-CoV-2, loading each gene once instead of once per threshold.

Shelling out to `recalculate_scores.py` for every combination, as `threshold_sweep_fast.py`
does, reads the raw results afresh each time. That is fine for HIV-1 and HCMV, whose files are
small, but SARS-CoV-2 keeps 2.3 GB for ORF1ab and 2.6 GB for S, and re-reading them nine times
dominates everything else: the first threshold alone had not finished a single gene after forty
minutes.

The threshold is only consumed at the compilation step, well after loading, so the loop can be
turned inside out. Each gene is read once and then compiled at every threshold, which divides
the reading cost by the number of thresholds.

Results go to <output-dir>/t<threshold>/<dataset>/<gene>/results/, the same layout the other
sweep scripts produce, so `threshold_sensitivity.py` reads them without changes. Nothing is
written to the input tree.

Usage
-----
    python scripts/threshold_sweep_sars.py --output-dir <directory outside the repository>
"""

from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))

import kss  # noqa: E402
import utils  # noqa: E402
from pipeline import load_config  # noqa: E402

THRESHOLDS = [10, 15, 20, 25, 30, 33, 40, 45, 50]
CONFIG = ROOT / "data" / "Severe_acute_respiratory_syndrome_coronavirus_2" / "config.yaml"

# Smallest first, so a failure surfaces before the expensive genes are attempted.
GENE_ORDER = ["E", "M", "N", "ORF1ab", "S"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--thresholds", type=int, nargs="*", default=THRESHOLDS)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dataset_name, dataset_info, parameters = load_config(str(CONFIG))
    base_folder = dataset_info["input_folder"]
    started = time.time()

    for gene in GENE_ORDER:
        raw_path = Path(base_folder) / gene / "results" / f"{gene}_results.json"
        if not raw_path.exists():
            print(f"[{gene}] raw results missing, skipping", flush=True)
            continue

        pending = [t for t in args.thresholds
                   if not (args.output_dir / f"t{t:02d}" / dataset_name / gene /
                           "results" / f"{gene}_compiled_results.json").exists()]
        if not pending:
            print(f"[{gene}] all thresholds already present, skipping", flush=True)
            continue

        size_gb = raw_path.stat().st_size / 1024 ** 3
        print(f"[{gene}] loading {size_gb:.1f} GB "
              f"({(time.time() - started) / 60:.0f} min elapsed)", flush=True)
        results = utils.load_data_from_json(str(raw_path))
        print(f"[{gene}] loaded, compiling at {len(pending)} thresholds", flush=True)

        for threshold in pending:
            target = args.output_dir / f"t{threshold:02d}" / dataset_name
            gene_dir = target / gene / "results"
            gene_dir.mkdir(parents=True, exist_ok=True)
            local = dict(parameters)
            local["threshold"] = threshold
            infos = {"input_folder": base_folder, "output_folder": str(target)}
            try:
                compiled = kss.compile_results(results, local, target_gene=gene)
                compiled, _ = kss.compute_kss_scores(compiled, results, infos, local,
                                                     target_gene=gene)
                payload = {gene: compiled[gene]} if gene in compiled else {}
                utils.save_data_as_json(
                    payload, str(gene_dir / f"{gene}_compiled_results.json"))
                print(f"[{gene}] t={threshold} done "
                      f"({(time.time() - started) / 60:.0f} min)", flush=True)
            except Exception as error:
                print(f"[{gene}] t={threshold} FAILED: {error}", flush=True)

        del results
        gc.collect()

    # A threshold directory is complete once every gene wrote its file there.
    for threshold in args.thresholds:
        target = args.output_dir / f"t{threshold:02d}"
        produced = list(target.glob(f"*/*/results/*_compiled_results.json")) \
            if target.exists() else []
        if len(produced) == len(GENE_ORDER):
            (target / "COMPLETE").write_text(f"threshold={threshold}\n", encoding="utf-8")

    print(f"\nfinished in {(time.time() - started) / 60:.0f} minutes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
