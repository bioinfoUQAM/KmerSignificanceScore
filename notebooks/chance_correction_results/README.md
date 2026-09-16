# Chance-corrected comparability of the score across datasets

Written by `notebooks/chance_correction.py`, which is deterministic under a fixed seed.

The discriminative score is a pure function of the k-mer by class contingency table, which the
script verifies by recomputing the published value from those tables. Permuting the class labels
therefore amounts to resampling a table at fixed margins, which it does by multivariate
hypergeometric draws. No sequence is reread and no classifier is retrained. Sequences whose k-mer
did not pass the prevalence threshold appear as an all-zero row, a configuration in its own right
that weighs in the conditional entropy but not in the purity term; that residual row is part of
the table and is resampled with the rest.

| File | What it holds |
|---|---|
| `per_position.csv` | one row per analyzed position: observed score, mean and standard deviation of the null, its 95th percentile and maximum, the empirical p value, and the z score |
| `null_by_dataset.csv` | the same distribution aggregated per dataset |
| `chance_ceiling_sensitivity.csv` | the null threshold recomputed at reduced effective sample sizes, as a stress test over an order of magnitude rather than an estimate of phylogenetic structure |

The chance ceilings do not coincide across the three datasets, so a raw score low in the range does
not mean the same thing everywhere. What the manuscript states from these files is bounded
accordingly: above the highest ceiling observed anywhere in these data, a score separates from
chance in every dataset, and the positions the framework prioritizes score well above that
ceiling. `per_position.csv` carries the observed score and the ceiling for every position, so
the margin can be read there rather than quoted here. Separation from chance is not comparability of meaning, and the manuscript concludes
that neither a universal threshold nor a statistical or biological interpretation transfers from
one dataset to another.

The lowest ceiling belongs to the largest dataset, which is a sample-size effect rather than a
class-count effect: the complexity factor corrects for the number of classes, not for the number
of sequences.

## Reproduce

```
python notebooks/chance_correction.py --permutations 1000
```

About an hour. A dry run at fifty permutations agrees to the third decimal and takes minutes.
