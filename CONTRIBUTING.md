# Contributing to KSS

This document is the map a reader needs before changing anything: how a run flows through the
modules, what each one owns, and what is public. For installing and running the tool, see the
[README](README.md).

## How a run flows

```
main.py                        parses the command line, loops over datasets and genes
  ├─ pipeline.py               finds and validates config.yaml, resolves output paths
  ├─ kanalyzer.py              aligns each sequence to the reference, extracts the windows
  │    └─ mutation_score.py    impact of the amino acid changes it identifies
  ├─ kss.py                    combines the components, compiles the results
  │    ├─ discriminative_score.py   class separation at a window
  │    ├─ mutation_score.py         impact of the amino acid changes
  │    └─ protein_score.py          characterization depth of the protein
  ├─ utils.py                  reads and writes the JSON results
  └─ report.py                 optional PDF of the highest-scoring positions, --report N
```

`main.py` calls these itself: `pipeline` hands it the configuration, and the rest run in the
order above. `pipeline` also supports rescoring, which `recalculate_scores.py` drives: it loads
saved results through `utils` and calls `kss` without rerunning `kanalyzer`.

The three scoring modules do not call one another, and `kss` is where the components are
combined. `mutation_score` is the one called from two places, by `kanalyzer` as it identifies
the changes and by `kss` as it scores them.

## What each module owns

| Module | Owns |
|---|---|
| `pipeline.py` | Configuration: finding `config.yaml`, validating parameters, output paths. Also rescoring from results already on disk |
| `kanalyzer.py` | Alignment, window extraction, mutation identification |
| `discriminative_score.py` | The discriminative component and its complexity scaling |
| `mutation_score.py` | Substitution matrices, their loading, inversion and rescaling |
| `protein_score.py` | UniProt retrieval, caching, the ten annotation categories |
| `kss.py` | Weighted combination, per-position records, compiled output |
| `report.py` | Reads compiled results and writes the optional PDF, without altering the scores |
| `utils.py` | JSON read and write |

## Public API

What `import src` exposes, and what a downstream user is expected to call:

```python
from src import (
    compute_kss_scores,      # compiled results -> per-position scores
    compile_results,         # per-gene analysis -> compiled results
    get_taxon_id,            # taxon of a GenBank reference
    analyze_records,         # sequences + reference -> per-gene analysis
    load_data_from_json,
    save_data_as_json,
)
```

Anything else, including every name starting with an underscore, is internal and may change
without notice. `main.py` is the command-line surface; `pipeline.py` is an internal configuration
and recalculation helper, not a library API.

## Development install and tests

```bash
pip install -e ".[analysis]"             # core plus what the analyses and the report need
python tests/test_matrix_formats.py      # substitution matrices survive a JSON round trip
python tests/test_docstring_coverage.py  # every definition of src/ and main.py documents itself
python tests/test_weight_bounds.py       # admissible weights keep KSS inside [0, 1]
python tests/test_offline_and_failure.py # the quick start runs offline, a missing input fails loudly
python main.py data/toy/config.yaml      # the toy dataset, deterministic, about a second
```

The toy run is a fast regression check for the code paths it exercises, and the README gives the
scores it is expected to print.

## Adding a substitution matrix

Matrices live in `src/substitution_matrices/` as JSON of amino acid pairs and values, one file
per matrix, named after the matrix. Write them through the project's own writer, which does not
round:

```python
from src.mutation_score import save_substitution_matrix, get_substitution_matrix

save_substitution_matrix(matrix, "src/substitution_matrices/MY_MATRIX")
get_substitution_matrix("MY_MATRIX")   # loads it back by name
```

`tests/test_matrix_formats.py` then covers it automatically: it reads every matrix in the
directory and checks that a round trip changes no value. Rounding on write is what made an
earlier matrix differ from its source in the fifteenth decimal, so the test is the guard against
reintroducing it.

## Changing a component

Each component returns a value in `[0, 1]`, and `kss.py` combines them under weights that default
to equal thirds. Two things are worth knowing before editing one:

- Each component is normalized to `[0, 1]`, an invariant their weighted combination requires. A
  component that can leave that interval breaks it.
- `--report-only` redraws an existing compiled result. It does not incorporate scoring changes
  that have not been rerun, so it shows what a change did only after the dataset is rescored.

## Docstrings

Every function, method and class in `src/` and `main.py` carries one, and
`tests/test_docstring_coverage.py` fails if a definition is added without it. Keep them to a
summary, the arguments and the return: the contract, not a narrative.
