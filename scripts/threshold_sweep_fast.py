"""
The part of the threshold sweep that answers the comment, without the part that costs hours.

Reviewer #1 objects to the default t = 0.25 having been "arbitrarily increased to t = 0.33 for
HIV-1". The objection is about the HIV-1 exception, and HIV-1 recomputes in about a minute
while SARS-CoV-2 takes an hour per threshold because two of its genes hold multi-gigabyte raw
results. This runs the two datasets that recompute cheaply and produced the published values;
`threshold_sweep_sars.py` runs the third, and `threshold_sweep.py` runs all three in one pass.

Usage
-----
    python scripts/threshold_sweep_fast.py --output-dir <directory outside the repository>
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

THRESHOLDS = [10, 15, 20, 25, 30, 33, 40, 45, 50]
CONFIGS = [
    ROOT / "data" / "Human_immunodeficiency_virus_1" / "config.yaml",
    ROOT / "data" / "Human_betaherpesvirus_5" / "config.yaml",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    for threshold in THRESHOLDS:
        target = args.output_dir / f"t{threshold:02d}"
        if (target / "COMPLETE").exists():
            print(f"[t={threshold}] already done", flush=True)
            continue
        target.mkdir(parents=True, exist_ok=True)
        ok = True
        for config in CONFIGS:
            outcome = subprocess.run(
                [sys.executable, str(ROOT / "recalculate_scores.py"), str(config),
                 "--threshold", str(threshold), "-o", str(target)],
                cwd=str(ROOT), capture_output=True, text=True)
            if outcome.returncode != 0:
                ok = False
                print(f"[t={threshold}] {config.parent.name} FAILED", flush=True)
                print(outcome.stderr[-1500:], flush=True)
        if ok:
            (target / "COMPLETE").write_text(f"threshold={threshold}\n", encoding="utf-8")
        print(f"[t={threshold}] done ({(time.time() - started) / 60:.1f} min elapsed)",
              flush=True)

    print(f"\nfinished in {(time.time() - started) / 60:.1f} minutes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
