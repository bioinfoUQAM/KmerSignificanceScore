# Prevalence-threshold sensitivity

Written by `scripts/threshold_sweep_fast.py` (HIV-1, HCMV) and `scripts/threshold_sweep_sars.py`
(SARS-CoV-2), whose sweeps are then assembled by `scripts/threshold_sensitivity.py`.

| File | What it holds |
|---|---|
| `threshold_sensitivity.csv` | one row per dataset and threshold: positions, variants, mean, standard deviation and variance of the score, standard deviation of the discriminative component, overlap of the reported selection sizes with the published threshold, and the rank correlation of the complete score vector |

The number of analyzed positions does not change with the threshold. Only the number of variants
retained inside each position does, so the threshold filters rare variants within a fixed search
space rather than redefining it, which is what makes the score vectors comparable across
thresholds.

No realignment is involved. The raw per-sequence k-mers are already stored and
`recalculate_scores.py` replays the scoring from them, which is why a whole grid is affordable.
With `--threshold`, `recalculate_scores.py` requires `--output-dir`, so a sweep cannot overwrite the
published results.

## Reproduce

```
python scripts/threshold_sweep_fast.py --output-dir <a directory outside the repository>
python scripts/threshold_sweep_sars.py --output-dir <another directory outside the repository>
python scripts/threshold_sensitivity.py --sweep-dir <first directory> --sweep-dir <second>
```

Both sweeps are idempotent, marking each finished threshold, so a relaunch resumes where it
stopped. HIV-1 and HCMV take minutes; SARS-CoV-2 takes hours, its two largest genes exceeding a
gigabyte, which is why it runs separately with a script that loads each gene once instead of once
per threshold. `scripts/threshold_sweep.py` sweeps all three datasets in one pass, rereading the
raw results for every combination; it was used for an independent verification run rather than
for the published values.
