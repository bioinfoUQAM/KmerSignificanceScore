"""
Sensitivity of the prioritization to the prevalence threshold (Reviewer #1, comment 2.2).

The reviewer objects that the default threshold t = 0.25 was raised to 0.33 for HIV-1 without
justification, and asks for a sensitivity analysis over t in [0.1, 0.5].

No realignment is needed. The raw per-sequence k-mer results are already stored, and
`recalculate_scores.py` replays the scoring from them, so a whole grid costs a few hours of
recomputation rather than a full rerun of the pipeline.

Every run writes to its own directory under --output-dir. The published results are never
touched, which is enforced by the command-line guard in `recalculate_scores.py`: passing
--threshold without an output directory is rejected.

Usage
-----
    python scripts/threshold_sweep.py --output-dir <directory>
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# The manuscript uses 25 for SARS-CoV-2 and HCMV and 33 for HIV-1; both are in the grid so
# the published configuration is reproduced inside the sweep rather than assumed.
THRESHOLDS = [10, 15, 20, 25, 30, 33, 40, 45, 50]

CONFIGS = [
    ROOT / "data" / "Human_betaherpesvirus_5" / "config.yaml",
    ROOT / "data" / "Human_immunodeficiency_virus_1" / "config.yaml",
    ROOT / "data" / "Severe_acute_respiratory_syndrome_coronavirus_2" / "config.yaml",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    for threshold in THRESHOLDS:
        target = args.output_dir / f"t{threshold:02d}"
        marker = target / "COMPLETE"
        if marker.exists():
            print(f"[t={threshold}] already done, skipping", flush=True)
            continue
        target.mkdir(parents=True, exist_ok=True)
        complete = True
        for config in CONFIGS:
            if not config.exists():
                print(f"[t={threshold}] missing config {config}", flush=True)
                complete = False
                continue
            label = config.parent.name
            print(f"[t={threshold}] {label} starting "
                  f"({(time.time() - started) / 60:.0f} min elapsed)", flush=True)
            outcome = subprocess.run(
                [sys.executable, str(ROOT / "recalculate_scores.py"), str(config),
                 "--threshold", str(threshold), "-o", str(target)],
                cwd=str(ROOT), capture_output=True, text=True)
            if outcome.returncode != 0:
                print(f"[t={threshold}] {label} FAILED", flush=True)
                print(outcome.stdout[-2000:], flush=True)
                print(outcome.stderr[-2000:], flush=True)
                complete = False
            else:
                print(f"[t={threshold}] {label} done", flush=True)

        # Only mark a threshold done when every dataset produced results. Marking it
        # unconditionally would let a rerun skip a threshold that is missing a dataset, and the
        # gap would then be invisible: the ingestion step would simply average what it found.
        if complete:
            marker.write_text(f"threshold={threshold}\n", encoding="utf-8")
        else:
            print(f"[t={threshold}] incomplete, not marked as done", flush=True)

    print(f"\nsweep finished in {(time.time() - started) / 60:.0f} minutes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
