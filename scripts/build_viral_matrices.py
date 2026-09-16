"""Build the virus-derived replacement models named by reviewer 2.10.

Point 2.10 asks for three things: VTML200, WAG, and "modern models specifically derived for
RNA virus evolution". The first two are handled by `build_added_matrices.py`; this covers the
third, with the five empirical models that answer the description:

  HIVb    HIV-1, between-host    Nickle, Heath, Jensen, Gilbert, Mullins & Kosakovsky Pond,
                                 PLoS ONE 2007;2(6):e503
  HIVw    HIV-1, within-host     same reference
  rtREV   retroviral reverse transcriptase   Dimmic, Rest, Mindell & Goldstein,
                                 J Mol Evol 2002;55(1):65-73
  FLU     influenza, predominantly influenza A   Dang, Le, Le & Le Si Quang,
                                 BMC Evol Biol 2010;10:99
  FLAVI   flaviviruses           Le & Vinh, J Mol Evol 2020;88(5):445-452

HIVb and HIVw matter most here, since HIV-1 is one of the three datasets analysed. FLAVI is the
most recent of the family, and its own paper positions it against the other four.

All five are rate models, like WAG, so they are converted to log-odds score matrices by the same
construction and at the same divergence, and `build_added_matrices.wag_logodds` does the work.
Nothing in the conversion is specific to WAG.

Provenance. Every source is read from a local cache under `added_matrix_sources/`, and each model
is cross-checked against a second, independent source, exactly as VTML200 was cross-checked
between spaln and SeqAn. The complete GPL-licensed IQ-TREE and RAxML sources are deliberately not
redistributed in this MIT-licensed repository: they are downloaded from immutable upstream
commits and verified by SHA-256. HIVb, HIVw, rtREV and FLU are read from IQ-TREE and checked
against RAxML; FLAVI is read from IQ-TREE and checked against the matrix its authors distribute,
RAxML predating it. Both the exchangeabilities and the equilibrium frequencies are compared, the
latter within the precision each source prints. Agreement is required, not assumed: the script
refuses to write a model whose two sources disagree.

Usage
-----
    python scripts/build_viral_matrices.py
    python scripts/build_viral_matrices.py --refresh   # re-download the upstream sources
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from build_added_matrices import (  # noqa: E402
    AA20, THIRD_BITS, ensure_source, wag_logodds, write_json,
)

SOURCES = HERE / "added_matrix_sources"
DIVERGENCE = 2.5          # the divergence WAG was converted at, kept identical

UPSTREAM = {
    "iqtree_modelprotein.cpp":
        "https://raw.githubusercontent.com/iqtree/iqtree2/"
        "a00094e03d1ae984e1497e16738f91514df8c366/model/modelprotein.cpp",
    "raxml_models.c":
        "https://raw.githubusercontent.com/stamatak/standard-RAxML/"
        "36ec36110631c34692abcd4f24ca7b3e2fea742a/models.c",
    "FLAVI_authors.PAML":
        "https://raw.githubusercontent.com/thulekm/flavi/"
        "8a974c9755106e13b4532eeb3568e1d26ba7d308/FLAVI.PAML",
}

# None of these three upstream files is ours. The two complete sources carry GPL notices, and
# the FLAVI rate file comes from a repository that publishes no licence at all, which grants no
# redistribution right by default. They are inputs, not part of KSS, so they are cached locally
# rather than redistributed. Pinning both the commit and the bytes makes removing them from git
# neutral for reproducibility.
PINNED_SHA256 = {
    "iqtree_modelprotein.cpp":
        "2d568ee0d6cec80f364cdaa34db64ffef345221c62f77bdc3843d3a373659901",
    "raxml_models.c":
        "b2fda59ef65c2b5a58d309ec5eee47c325f780b92ee96d0cb5f11addd8d8ab35",
    "FLAVI_authors.PAML":
        "e688785520dfd80fe89aa6baf92119cdbd2f062257aadf02f392172fa42faff8",
}

# IQ-TREE spells them upper case; the name on the left is what we publish. RAxML predates
# FLAVI, published in 2020, so that model is checked against its authors' own distribution
# instead, which is the better source of the two anyway.
MODELS = {"HIVB": "HIVb", "HIVW": "HIVw", "RTREV": "rtREV", "FLU": "FLU", "FLAVI": "FLAVI"}
SECOND_SOURCE = {name: "raxml_models.c" for name in MODELS}
SECOND_SOURCE["FLAVI"] = "FLAVI_authors.PAML"

# How far the frequencies of a second source may sit from IQ-TREE's before the model is
# refused. RAxML prints them to four decimals, which caps the disagreement at 5e-5; the
# authors' own file carries full precision and must match. A single loose bound would have
# accepted ten times the discrepancy any of these sources can legitimately produce.
FREQUENCY_TOLERANCE = {"raxml_models.c": 6e-5, "FLAVI_authors.PAML": 1e-9}


def parse_iqtree(text: str, model: str) -> tuple[np.ndarray, np.ndarray]:
    """Exchangeabilities and equilibrium frequencies of one IQ-TREE model block.

    The block is PAML's layout: nineteen rows of a strictly lower triangle, then a line of
    twenty frequencies closed by a semicolon.
    """
    match = re.search(rf"^model {model}\s*=\s*$(.*?);", text, re.MULTILINE | re.DOTALL)
    if match is None:
        raise LookupError(f"{model} not found in the IQ-TREE source")
    numbers = [float(value) for value in match.group(1).split()]

    expected = 19 * 20 // 2 + 20
    if len(numbers) != expected:
        raise ValueError(f"{model}: expected {expected} numbers, found {len(numbers)}")

    rates = np.zeros((20, 20))
    cursor = 0
    for i in range(1, 20):
        for j in range(i):
            rates[i, j] = rates[j, i] = numbers[cursor]
            cursor += 1
    frequencies = np.array(numbers[cursor:])
    return rates, frequencies / frequencies.sum()


def parse_paml(text: str) -> tuple[np.ndarray, np.ndarray]:
    """A bare PAML rate file: nineteen lower-triangle rows, then twenty frequencies."""
    numbers = [float(value) for value in text.split()]
    expected = 19 * 20 // 2 + 20
    if len(numbers) != expected:
        raise ValueError(f"expected {expected} numbers in the PAML file, found {len(numbers)}")

    rates = np.zeros((20, 20))
    cursor = 0
    for i in range(1, 20):
        for j in range(i):
            rates[i, j] = rates[j, i] = numbers[cursor]
            cursor += 1
    frequencies = np.array(numbers[cursor:])
    return rates, frequencies / frequencies.sum()


def parse_raxml(text: str, model: str) -> tuple[np.ndarray, np.ndarray]:
    """Exchangeabilities and frequencies of one RAxML model, from its `case` block."""
    match = re.search(rf"case {model}:(.*?)break;", text, re.DOTALL)
    if match is None:
        raise LookupError(f"{model} not found in the RAxML source")

    # The exponent has to be part of the capture: FLU writes entries such as 1.41e-02, and
    # stopping the number at the "e" silently multiplies them by a power of ten.
    rates = np.zeros((20, 20))
    for i, j, value in re.findall(
            r"daa\[\s*(\d+)\s*\*\s*20\s*\+\s*(\d+)\s*\]\s*=\s*([\d.]+(?:[eE][+-]?\d+)?)",
            match.group(1)):
        rates[int(i), int(j)] = rates[int(j), int(i)] = float(value)
    if not rates.any():
        raise ValueError(f"{model}: no exchangeabilities parsed from the RAxML source")

    frequencies = np.zeros(20)
    for index, value in re.findall(r"f\[\s*(\d+)\s*\]\s*=\s*([\d.]+(?:[eE][+-]?\d+)?)",
                                   match.group(1)):
        frequencies[int(index)] = float(value)
    return rates, frequencies / frequencies.sum()


def agree(a: np.ndarray, b: np.ndarray) -> tuple[bool, float]:
    """Whether two exchangeability matrices agree up to a common scale factor.

    A rate model is only defined up to scale, since the diagonal is fixed by the rows and the
    whole matrix is renormalised to one substitution per site. Proportionality is therefore the
    right notion of agreement, and it is what the two distributions are expected to share.
    """
    upper = np.triu_indices(20, 1)
    x, y = a[upper], b[upper]
    if not y.any():
        return False, float("nan")
    scale = float(np.sum(x * y) / np.sum(y * y))
    residual = float(np.max(np.abs(x - scale * y)) / max(np.max(x), 1e-12))
    return residual < 1e-6, residual


def fetch(refresh: bool) -> dict[str, str]:
    """Read cached sources, downloading missing ones and verifying pinned inputs."""
    return {
        name: ensure_source(SOURCES / name, url, PINNED_SHA256.get(name), refresh)
              .decode("utf-8", errors="replace")
        for name, url in UPSTREAM.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="re-download the upstream sources before building")
    args = parser.parse_args()

    SOURCES.mkdir(exist_ok=True)
    contents = fetch(args.refresh)

    print()
    for iqtree_name, published in MODELS.items():
        rates, frequencies = parse_iqtree(contents["iqtree_modelprotein.cpp"], iqtree_name)

        source_name = SECOND_SOURCE[iqtree_name]
        if source_name.endswith(".PAML"):
            other_rates, other_frequencies = parse_paml(contents[source_name])
        else:
            other_rates, other_frequencies = parse_raxml(contents[source_name], iqtree_name)

        ok, residual = agree(rates, other_rates)
        if not ok:
            raise SystemExit(f"{published}: the two sources disagree on the exchangeabilities "
                             f"(relative residual {residual:.2e}); refusing to write it")

        # Frequencies decide the log-odds as much as the rates do, so they are compared too,
        # against what the second source can legitimately differ by and no more.
        frequency_gap = float(np.max(np.abs(frequencies - other_frequencies)))
        if frequency_gap > FREQUENCY_TOLERANCE[source_name]:
            raise SystemExit(f"{published}: the two sources disagree on the equilibrium "
                             f"frequencies (max gap {frequency_gap:.2e}); refusing to write it")

        scores = wag_logodds(rates, frequencies, DIVERGENCE)
        if np.abs(scores - scores.T).max() > 1e-9:
            raise SystemExit(f"{published}: log-odds not symmetric")
        rounded = np.round(scores * THIRD_BITS)
        matrix = {(AA20[i], AA20[j]): int(rounded[i, j])
                  for i in range(20) for j in range(20)}
        write_json(matrix, published)
        print(f"{published:<7} agrees with {source_name} on rates (residual {residual:.1e}) "
              f"and frequencies (gap {frequency_gap:.1e}), written at t = {DIVERGENCE}")


if __name__ == "__main__":
    main()
