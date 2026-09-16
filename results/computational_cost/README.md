# Computational cost

Written by `scripts/computational_cost.py`. Two measurements that answer different questions and
do not have the same scope.

| File | What it holds |
|---|---|
| `computational_cost.csv` | per-position scoring cost, one row per analyzed gene, with the ratio to chi-squared and to mutual information |
| `pipeline_stage_cost.csv` | the split between preprocessing and scoring over a complete run, for two genes, with the worker count in a column |

The stage split is measured with a single worker. That is the condition the published timings were
taken in, since the alignment step was parallelised afterwards. Fixing the worker count does not
make timings deterministic, machine load still moves them; it makes the two stages comparable
under the same resources. The split is measured on SARS-CoV-2 E and HCMV US28 only. On a longer
gene the number of scored positions grows with length while pairwise alignment grows with its
square, so the preprocessing share is not expected to fall; that is a direction given by the
complexities, not a measurement.

## Provenance of the whole-run claim

The Conclusion states, as the submitted version already did, that the complete analysis of more
than 290,000 sequences over eleven genes ran in under 24 hours using a single worker, on a laptop
with an Intel Core i9-14900HX and 64 GB of RAM.

That note is retrospective. It records the conditions of the run reported, and the raw execution
log was not kept, so the duration cannot be regenerated independently from the versioned
artefacts. It is the only published timing in that situation, which is why the Conclusion presents
it as an observation rather than as a demonstration of efficiency.

- It is the historical run of the submitted version, before the alignment step was parallelised.
- The one large-scale pipeline measurement that survives covers ORF1ab: 0.25 s per sequence,
  19.5 h on one core, and 136 min of alignment on twelve workers. It was recorded in the body of
  commit `39e0e68` and nowhere else, and is repeated here so that it outlives that history.
- ORF1ab dominates the total, alignment cost growing with the square of the gene length. That
  consistency makes the claim plausible; it does not demonstrate it and must not be presented as
  its proof.

Any re-measurement would give a different number, the code having changed twice since submission,
the window fix and the parallelisation.

## Reproduce

```
python scripts/computational_cost.py --gene <gene>   # one gene, writes a partial table
python scripts/computational_cost.py --merge         # assembles computational_cost.csv
python scripts/computational_cost.py --stages --workers 1
```

One process per gene keeps the peak memory at one gene rather than their sum: the raw results of
the large SARS-CoV-2 genes are gigabytes. The per-position run reads each gene's compiled results,
which ship under `data/`, and its raw per-sequence results, which a pipeline run writes and the
repository does not carry. `--stages` needs the raw `sequences.fasta` of the two genes it times.
