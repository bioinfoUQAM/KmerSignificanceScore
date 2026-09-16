# Source data for the seven matrices added in reviewer response 2.10

Reviewer 2.10 asked us to compare against more recent or contextual substitution
matrices, naming three things: **VTML200**, **WAG**, and models derived for RNA
virus evolution. This directory is the local source cache for the seven matrices
that answer those three, so every derived file in
`../../src/substitution_matrices/` is auditable back to a documented upstream source.

| Matrix | Source read | Checked against | Built by |
|---|---|---|---|
| VTML200 | `VTML200_spaln.txt`, from the Spaln distribution | SeqAn's `Vtml200` | `../build_added_matrices.py` |
| WAG | `wag.dat`, the authors' file on the EBI server | PAML's copy | `../build_added_matrices.py` |
| HIVb, HIVw | `iqtree_modelprotein.cpp` | `raxml_models.c` | `../build_viral_matrices.py` |
| rtREV, FLU | `iqtree_modelprotein.cpp` | `raxml_models.c` | `../build_viral_matrices.py` |
| FLAVI | `iqtree_modelprotein.cpp` | `FLAVI_authors.PAML`, from its authors | `../build_viral_matrices.py` |

**None of the five source files is redistributed here.** Each one belongs to
someone else, and each for a different reason:

| File | Why it is not ours to redistribute |
|---|---|
| `iqtree_modelprotein.cpp` | complete IQ-TREE source, carries a GPL-2.0-or-later notice |
| `raxml_models.c` | complete RAxML source, carries a GPL-2.0 notice |
| `VTML200_spaln.txt` | travels inside the GPL-licensed `spaln` distribution; the table itself is Kneutgen and Mueller's, with no licence of its own |
| `FLAVI_authors.PAML` | `thulekm/flavi` publishes **no licence at all**, which grants no redistribution right by default |
| `wag.dat` | the authors' own file on the EBI server, carrying its citation request but no licence statement |

They are downloaded on first use and verified against pinned SHA-256 digests,
and the copies in this directory are ignored local caches. A digest mismatch is
fatal in both directions: a bad download is never cached, and a tampered cache is
never read. Where the upstream is a git repository the URL is pinned to an
immutable commit; the EBI serves a static file, so there the digest is the pin.

After a fresh clone, run both builders once; no manual download is required:

```
python scripts/build_added_matrices.py     # VTML200, WAG
python scripts/build_viral_matrices.py     # HIVb, HIVw, rtREV, FLU, FLAVI
```

Removing the raw files from the repository does not settle every question about
the derived matrices in `../../src/substitution_matrices/`, which remain under
the repository licence. It settles the plain one: we no longer ship copies of
other people's files.

Every model is required to agree with its second source before it is written. For
the six rate models the check covers both the exchangeabilities and the
equilibrium frequencies, the latter within the precision each source prints:
RAxML reports four decimals, so a disagreement up to 5e-5 is expected there and
nowhere else. RAxML predates FLAVI, published in 2020, which is why that model is
checked against its authors' own distribution instead, a better source in any
case.

Sensitivity to the divergence `t`, the one free parameter introduced when a rate
model is converted into a score matrix, is computed for all six rate models by
`../divergence_sensitivity.py`. The conversion itself is checked by independent
routes in `../verify_viral_matrices.py`.

## VTML200 (`VTML200_spaln.txt`)

Integer score matrix of the VTML family, at the 200 evolutionary-time point
(analogous to PAM250 within the PAM series). Units: third-bits.

- The VTML200 numerical scoring matrix was obtained from the `spaln`
  distribution and independently cross-checked against the VTML200
  implementation distributed with SeqAn. The spaln copy was taken from the
  GeneMark-ETP repository (`gatech-genemark/GeneMark-ETP`, path
  `bin/gmes/ProtHint/dependencies/spaln_table/vtml200`); its header states it
  was produced by Kai Kneutgen and Tobias Mueller's own VTML scripts (May 2002).
- **Cross-checked** entry by entry against the independent `Vtml200` matrix in
  the SeqAn C++ library (`seqan/seqan`,
  `include/seqan/score/score_matrix_data.h`). The two agree exactly on all 190
  off-diagonal amino acid pairs and the 20 diagonal entries (e.g. W:W = 15,
  C:C = 12, A:R = -2, I:L = 3, D:E = 3).
- We keep only the 20 standard amino acids, dropping the ambiguity columns
  (B, Z, X, *).
- Primary reference: Mueller T, Vingron M (2000), *Modeling amino acid
  replacement*, J Comput Biol 7(6):761-776. Method companion: Mueller T,
  Spang R, Vingron M (2002), *Estimating amino acid substitution models*,
  Mol Biol Evol 19(1):8-13.

## WAG (`wag.dat`)

An empirical model of protein evolution, distributed as a symmetric
exchangeability matrix plus equilibrium amino acid frequencies (the `wag.dat`
format), rather than as a score matrix.

- **Taken verbatim** from the authors' own distribution on Nick Goldman's
  server at the EBI: `https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat`. We
  verified that PAML's redistribution (`abacus-gene/paml`, `dat/wag.dat`) is
  numerically identical to it (maximum absolute difference 0.0 on both the
  exchangeabilities and the equilibrium frequencies).
- Reference: Whelan S, Goldman N (2001), *A general empirical model of protein
  evolution derived from multiple protein families using a maximum-likelihood
  approach*, Mol Biol Evol 18(5):691-699.
- Because WAG is a rate model and the benchmark scores similarity matrices, we
  convert it to a log-odds score matrix by the standard construction (see
  `build_added_matrices.py`): the instantaneous rate is `Q = R diag(pi)` with
  `Q_ii = -sum_j Q_ij`, normalised to one expected substitution per site;
  `P(t) = exp(Q t)`; and `S_ij = log(P(t)_ij / pi_j)`, which is symmetric by
  reversibility. Scores are expressed in third-bit units (multiply the natural
  log-odds by `3/ln 2`) and rounded to integers, matching the VTML200 unit. We
  use `t = 2.5`, a divergence comparable to the 250-PAM matrices already in the
  benchmark (PAM250, VTML250, GONNET).

  The choice of `t` has no material effect on the result, for two reasons the
  code checks explicitly. First, the evaluation ranks matrices by Spearman rank
  correlation over the 190 off-diagonal pairs, so it is invariant to the scale
  of the entries; the off-diagonal rank order of the log-odds is stable above
  rho = 0.99 for t in [2, 4] (asserted in `build_added_matrices.py`). That
  figure describes the continuous log-odds, before the rounding to integer
  third-bits that produces the committed matrix; rounding collapses the 190
  pairs onto 11 distinct values, so between the matrices actually evaluated the
  correlation is 0.94 to 0.97. Second, and this is the argument that carries the
  point, recomputing the full biophysical composite with WAG derived at
  t = 2.0, 2.5, 3.0 and 4.0 gives 3.259, 3.291, 3.313 and 3.133 respectively, a
  spread of 0.18, and WAG stays mid-ranked among the candidate matrices in every
  case, so no reasonable choice of t changes the conclusion.

  The same analysis now covers all six rate models, not WAG alone, and is
  computed by `../divergence_sensitivity.py`, which writes
  `results/matrix_evaluation/divergence_sensitivity.csv` and reports for each
  model and each divergence the composite, both ranks, the number of distinct
  entry values after rounding, and the rank correlation against the published t.
  Before it existed the four WAG composites appeared only as prose, in this file
  and in the response letter, with nothing behind them a reader could recompute.
