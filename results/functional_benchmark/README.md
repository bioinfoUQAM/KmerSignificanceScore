# External reference layers

Written by `scripts/reference_matrix.py`, which builds the published matrix, and by
`scripts/functional_benchmark.py`, which characterises the selection against amino acid fitness
estimates. The protocol was fixed before any of it was computed, and every departure from it is
recorded with its reason in the pre-specification kept with the manuscript's revision record.

| File | What it holds |
|---|---|
| `reference_matrix.csv` | one row per annotation layer and selection size: the evaluable universe, the layer's prevalence, the four confusion counts, precision, false discovery rate, false positive rate, recall, false negative rate, enrichment, exact hypergeometric probabilities, and the area under the ROC curve over the complete ranking. Source of S6 Table |
| `fitness_characterisation.csv` | the sign composition and the effect magnitude of the retained positions against the whole set |
| `cov2var_input_manifest.csv` | the COV2Var snapshot the published rates were computed from, file by file: name, size, SHA-256 and the members of the selection archive |

Every rate is relative to a pair of layer and selection size. A selected k-mer the layer does not
annotate counts as a false positive, and an annotated k-mer the selection does not reach counts as
a false negative, so the false negative rate mainly measures the size of the selection against the
number of positives a genome-wide layer carries. Coverage never becomes a negative in either
direction: a layer that assayed a set of residues is evaluated on the k-mers carrying one, and a
layer that scores substitutions only is evaluated on the ranking without its indel positions,
which are not evaluable against it and not wrong.

What the analysis concludes is that no available resource is a complete, genome-wide and
non-saturated reference for the quantity this prioritization ranks. Each layer fails for a
different measured reason, which is why the layers are reported side by side, each with its own
prevalence next to its rates. The prediction layers are reported under that label and are not read
as evidence, the mutational component being itself a prediction.

## Reproduce

The third-party inputs are not carried here. `scripts/external_inputs.py` names the source of each
and pins its SHA-256, and the analyses stop rather than run on a different download.

```
python scripts/reference_matrix.py --cov2var <directory of COV2Var downloads>
python scripts/build_s6_table.py
python scripts/functional_benchmark.py --fitness <aa_fitness.csv>
python scripts/build_cov2var_manifest.py
```

Cite Bloom and Neher 2023 (Virus Evolution 9(2):vead055) for the fitness estimates and COV2Var
(Nucleic Acids Research 2024;52(D1):D701-D713) for the annotation layers.
