"""Build the two substitution matrices named by reviewer 2.10: VTML200 and WAG.

Both are amino acid replacement models of the same log-odds family already
benchmarked (JTT/JONES, DAYHOFF, GONNET, VTML160, VTML250). They are written
into src/substitution_matrices/ in the project's JSON format so the matrix
evaluation notebook picks them up automatically.

Sources (see added_matrix_sources/README.md for URLs and citations):

  VTML200  Canonical integer score matrix produced by T. Mueller's own VTML
           scripts (Mueller, Spang & Vingron 2002), read from the spaln
           distribution at a pinned revision. Units: third-bits. Cross-checked
           entry by entry against SeqAn's Vtml200 (independent source, agrees
           exactly).

  WAG      An empirical rate model (Whelan & Goldman 2001), read from the
           authors' own file on Nick Goldman's EBI server as exchangeabilities
           plus equilibrium frequencies, and checked against PAML's copy.

Neither source file is redistributed here. Both are downloaded on first use and
verified against a pinned SHA-256; the local copies under added_matrix_sources/
are ignored caches. Run this script once after a fresh clone.
           Unlike VTML200 it is a rate model rather than a score matrix, so we
           convert it to a log-odds score matrix by the standard construction:
           Q = R diag(pi) with the usual diagonal, normalised to one expected
           substitution per site, P(t) = exp(Q t), and S_ij = log(P(t)_ij/pi_j),
           which is symmetric by reversibility. t is set to 2.5 (a divergence
           comparable to the 250-PAM matrices already in the benchmark). The
           evaluation ranks matrices by Spearman rank correlation on the 190
           off-diagonal amino acid pairs, so it is invariant to the scale and
           to strictly monotone transforms of the entries; the off-diagonal
           rank order of the log-odds is stable above rho = 0.99 for t in [2, 4]
           (asserted below), so the choice of t has no material effect.

           That assertion concerns the continuous log-odds, before the rounding
           to integer third-bits that produces the committed matrix. Rounding
           collapses the 190 pairs onto 11 distinct values, so between the
           matrices actually evaluated the rank correlation is 0.94 to 0.97
           rather than above 0.99. The claim that t does not matter therefore
           rests on the composite itself, which spans 0.18 over t in [2, 4] and
           leaves WAG mid-ranked throughout. `divergence_sensitivity.py`
           computes that, and reproduces the committed ranking at t = 2.5.
"""

import argparse
import hashlib
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import spearmanr

sys_dir = Path(__file__).resolve().parent
SOURCES = sys_dir / "added_matrix_sources"
OUT = sys_dir.parent / "src" / "substitution_matrices"

# PAML / spaln amino acid order.
AA20 = list("ARNDCQEGHILKMFPSTWYV")
THIRD_BITS = 3.0 / np.log(2)  # natural-log-odds -> third-bit units

# Neither source file is ours. VTML200 travels inside the GPL-licensed spaln distribution, and
# wag.dat is the authors' own file served by the EBI with no licence statement of its own. They
# are inputs, not part of KSS, so they are downloaded and verified rather than redistributed.
# Pinning the bytes is what makes removing them from git neutral for reproducibility.
UPSTREAM = {
    "VTML200_spaln.txt":
        "https://raw.githubusercontent.com/ogotoh/spaln/"
        "35ae966a7da21ad8c1a37909508801386c39c0a9/table/vtml200",
    # The EBI serves a static file, so there is no revision to pin; the digest is the pin.
    "wag.dat":
        "https://www.ebi.ac.uk/goldman-srv/WAG/wag.dat",
}

PINNED_SHA256 = {
    "VTML200_spaln.txt":
        "81b48e787645f600e13dc8f63ebc6f4052b2b111b89488b8cd8837f90589fb96",
    "wag.dat":
        "616c2f0c8c9a87e00590a9886a6789d95ee12706f586a56477a24209bee4b38e",
}


def ensure_source(path: Path, url: str, expected: str | None, refresh: bool) -> bytes:
    """Return a cached upstream source, downloading it if needed and verifying its digest.

    Shared with build_viral_matrices.py so that one implementation guards every external
    input. A digest mismatch is fatal in both directions: a bad download is never cached,
    and a tampered cache is never read.
    """
    if refresh or not path.exists():
        import urllib.request
        print(f"downloading {url}")
        with urllib.request.urlopen(url, timeout=60) as response:
            downloaded = response.read()
        digest = hashlib.sha256(downloaded).hexdigest()
        if expected is not None and digest != expected:
            raise SystemExit(
                f"{path.name}: downloaded source failed SHA-256 verification\n"
                f"  expected {expected}\n"
                f"  found    {digest}\n"
                "refusing to cache or use it"
            )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(downloaded)

    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if expected is not None and digest != expected:
        raise SystemExit(
            f"{path.name}: cached source failed SHA-256 verification\n"
            f"  expected {expected}\n"
            f"  found    {digest}\n"
            "delete the cache and rerun, or use --refresh"
        )
    print(f"  {path.name}: {len(raw):,} bytes, sha256 {digest[:16]}...")
    return raw


def fetch(refresh: bool = False) -> None:
    """Make sure both source files are present and verified before anything reads them."""
    for name, url in UPSTREAM.items():
        ensure_source(SOURCES / name, url, PINNED_SHA256.get(name), refresh)


def write_json(matrix: dict, name: str) -> None:
    """Write a matrix dict to OUT/name.json in the project's "(A, R)" format."""
    import json
    string_matrix = {f"({a}, {b})": float(v) for (a, b), v in matrix.items()}
    with open(OUT / f"{name}.json", "w") as f:
        json.dump(string_matrix, f, indent=4, sort_keys=True)


def build_vtml200() -> dict:
    rows, header = {}, None
    for line in (SOURCES / "VTML200_spaln.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        toks = line.split()
        if header is None:
            header = toks
            continue
        rows[toks[0]] = [int(x) for x in toks[1:]]
    col = {aa: i for i, aa in enumerate(header)}
    matrix = {(a, b): rows[a][col[b]] for a in AA20 for b in AA20}
    asym = sum(matrix[(a, b)] != matrix[(b, a)] for a in AA20 for b in AA20)
    assert asym == 0, f"VTML200 not symmetric: {asym} entries"
    return matrix


def load_wag_rate():
    lines = [l for l in (SOURCES / "wag.dat").read_text().splitlines()
             if l.strip() and not l.strip().startswith("//")]
    R = np.zeros((20, 20))
    for i in range(1, 20):
        vals = [float(x) for x in lines[i - 1].split()]
        assert len(vals) == i
        for j in range(i):
            R[i, j] = R[j, i] = vals[j]
    pi = np.array([float(x) for x in lines[19].split()])
    assert len(pi) == 20
    return R, pi / pi.sum()


def wag_logodds(R, pi, t):
    Q = R * pi[None, :]
    for i in range(20):
        Q[i, i] = -Q[i].sum()
    Q /= -np.sum(pi * np.diag(Q))          # one expected substitution per site
    P = expm(Q * t)
    return np.log(P / pi[None, :])          # symmetric by reversibility


def build_wag() -> dict:
    R, pi = load_wag_rate()
    S = wag_logodds(R, pi, 2.5)
    assert np.abs(S - S.T).max() < 1e-9, "WAG log-odds not symmetric"

    # Rank order of the 190 off-diagonal pairs is what the evaluation uses; confirm it barely
    # moves with the divergence t before committing to 2.5. This is measured on the continuous
    # log-odds; after the integer rounding below it drops to 0.94-0.97, which is why the
    # argument that t is immaterial rests on the composite, not on this number. See
    # divergence_sensitivity.py.
    idx = list(combinations(range(20), 2))
    base = np.array([S[i, j] for i, j in idx])
    for t in (2.0, 3.0, 4.0):
        alt = np.array([wag_logodds(R, pi, t)[i, j] for i, j in idx])
        rho = spearmanr(base, alt).statistic
        assert rho > 0.99, f"WAG rank order unstable at t={t}: rho={rho:.3f}"

    Sint = np.round(S * THIRD_BITS)         # third-bit log-odds, like VTML200
    return {(AA20[i], AA20[j]): int(Sint[i, j]) for i in range(20) for j in range(20)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build VTML200 and WAG from their upstream sources.")
    parser.add_argument("--refresh", action="store_true",
                        help="re-download both sources even if the local cache is present")
    args = parser.parse_args()

    fetch(args.refresh)
    vtml = build_vtml200()
    wag = build_wag()
    write_json(vtml, "VTML200")
    write_json(wag, "WAG")
    print(f"VTML200 -> {OUT / 'VTML200.json'}  (diag mean {np.mean([vtml[(a, a)] for a in AA20]):.2f})")
    print(f"WAG     -> {OUT / 'WAG.json'}      (diag mean {np.mean([wag[(a, a)] for a in AA20]):.2f})")
    print("Both symmetric; WAG off-diagonal rank order stable for t in [2, 4] (rho > 0.99).")


if __name__ == "__main__":
    main()
