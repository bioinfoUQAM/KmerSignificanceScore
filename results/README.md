# Result tables

Each directory holds the committed output of one analysis. The numbers reported in the article
are read from these files, so a table can be regenerated and compared against the copy here.
Scripts are deterministic: rerunning one reproduces its directory unchanged.

| Directory | Produced by | What the article reports from it |
|---|---|---|
| `statistical_comparison/` | `scripts/statistical_comparison.py`, `scripts/topk_overlap.py` | Paired comparison of the discriminative component against the eight baselines, and the overlap between their selections (Table 4, S5 and S7 Tables) |
| `matrix_evaluation/` | `scripts/build_matrix_ranking.py`, `scripts/divergence_sensitivity.py` | Ranking of the 37 candidate substitution matrices over the eight amino acid characteristics, and its sensitivity to the conversion divergence (Table 5, S1 Table) |
| `matrix_optimization/` | `notebooks/matrix_optimization.ipynb` | The genetic algorithm run that produced MIYATA\_EVO, with the generation it stopped at and its composite value |
| `protein_score_validation/` | `scripts/regenerate_protein_scores.py`, `scripts/build_discordance_tables.py` | The protein characterization score of 17,471 UniProt entries, its calibration against the UniProt annotation score, and the discordant entries (Table 6, S8 and S9 Tables) |
| `protein_weight_sensitivity/` | `scripts/protein_weight_sensitivity.py` | Effect of the protein weight on the reported selections |
| `threshold_sensitivity/` | `scripts/threshold_sweep_fast.py`, `scripts/threshold_sweep_sars.py`, `scripts/threshold_sensitivity.py` | Sweep of the prevalence threshold from 0.10 to 0.50 on the three datasets (S4 Table) |
| `floor_sensitivity/` | `scripts/floor_sensitivity.py` | Which amino acid pairs the mutational floor binds for, and over what range the reported selections hold |
| `functional_benchmark/` | `scripts/reference_matrix.py`, `scripts/functional_benchmark.py` | Error rates and enrichment of the prioritization against the COV2Var annotation layers (S6 Table) |
| `computational_cost/` | `scripts/computational_cost.py` | Per-position scoring cost per gene, and the split between preprocessing and scoring |

The twelve supplementary tables the article ships are collected under `supplementary/`, and each
names its producer. Eight are written from these directories by `scripts/build_s1_table.py`,
`scripts/build_s4_table.py`, `scripts/build_s5_table.py`, `scripts/build_s6_table.py`,
`scripts/build_s7_table.py` and `scripts/build_discordance_tables.py`, the last writing both S8
and S9. `notebooks/build_s3_table.py` writes S3 from `data/accessions/`, and
`notebooks/generate_kss_figure.ipynb` writes the three validation workbooks S10, S11 and S12 from
the compiled results under `data/`. `scripts/build_s2_table.py` writes
`supplementary/S2_Table.xlsx`, the ten annotation categories and their weights. It was the last
workbook without one, and it showed: the categories sat in ten vertically merged cells, so
scrolling a few rows down left a reader unable to tell which category a level belonged to. The
rubric is a literal in that script and is deliberately not imported from `src/protein_score.py`,
since the cross-check compares the two and a builder reading the scorer would turn that
comparison into the scorer agreeing with itself.

None of this appears in the supplementary legends: they describe what each file holds, and the
commands that rebuild it live here.

Two further result directories live beside the code that writes them:
`notebooks/discriminative_score_validation_results/`, from
`notebooks/discriminative_score_validation.ipynb`, and `notebooks/chance_correction_results/`,
from `notebooks/chance_correction.py`.

The external files two of these analyses read are not redistributed here. They are pinned by
checksum in `scripts/external_inputs.py`, which names the source of each, and the analyses stop
rather than run on a different download. The UniProt entries the protein characterization
score is computed from are deposited at <https://doi.org/10.5281/zenodo.21937072>: unpack
`uniprot_viral_proteins/` from that archive into `notebooks/` before running the two scripts of
`protein_score_validation/`.
