# What the protein weight controls

Written by `scripts/protein_weight_sensitivity.py`. No rerun of the pipeline: the compiled
results already store the three components per position and the score is their weighted mean.
The script reports, for each dataset, how many positions disagree with the published score at
the default weight, and whether the protein score is constant within every gene.

The protein component measures how thoroughly a protein has been characterised, not how important
it is, so its weight is a deliberate user setting rather than a biological correction. One
structural fact drives the outcome: the protein score is constant within a gene, so the weight
cannot reorder positions inside a gene. It only changes the balance between genes, and therefore
has no effect at all on a dataset ranked gene by gene.

| File | What it holds |
|---|---|
| `topk_overlap_by_weight.csv` | per dataset, unit, selection size and weight: positions shared with the default, Jaccard index, and the number changed |
| `topk_composition_by_weight.csv` | how the selected positions distribute over genes, with each gene's protein score |

Both files report two distinct objects, because they behave differently: the raw ranking by score,
and the selection actually drawn, which applies the cap on indel positions described in the
Methods.

## Reproduce

```
python scripts/protein_weight_sensitivity.py
```
