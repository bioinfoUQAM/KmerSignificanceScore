"""Every admissible COV2Var reference layer, at every selection size, in one table.

Comment 2.11 asks for false positive and false negative rates against an empirical database.
No single layer of COV2Var answers that on its own: the membership layer is saturated, the
selection layers are conditioned on the same variation the discriminative component exploits,
and the deep mutational scanning layers cover only the receptor binding domain. Reporting one
of them and calling it the reference would be choosing a number after seeing it.

So this reports all of them, across the whole range of k, with each layer's own prevalence
next to its rates. The reader sees what the rate is a rate of.

Admissibility follows the rule fixed in the pre-specification: a layer qualifies when it is
an experimental measurement or a statistical test on observed sequences. Layers produced by
software that predicts an effect from a sequence (MutPred2, iMutant, VaxiJen and IEDB,
ProtParam) are named and excluded, since benchmarking a prediction against a prediction
establishes nothing.

Coverage is never converted into a negative. A k-mer no layer evaluates is dropped from that
layer's universe and counted separately, which is what makes the receptor binding domain
layers usable at all: they measure 68 residues between spike 334 and 522, not that whole
interval, and the prioritization mostly selects the N-terminal domain, so scoring the rest as
negatives would report an assay's coverage as a method's error. Taking the interval for the
coverage is not hypothetical: it turned 17 unmeasured k-mers into negatives here once.

Usage
-----
    python scripts/reference_matrix.py --cov2var data/external/COV2Var [--fitness aa_fitness.csv]
"""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from external_inputs import verify_cov2var, verify_fitness  # noqa: E402
from kss_selection import INDEL_MUTATIONAL_SCORE  # noqa: E402
from functional_benchmark import (  # noqa: E402
    COV2VAR_GENES, K_NT, MIN_EXPECTED_COUNT, confusion, label, label_residues, load_rankings,
    ranked,
)

RESULTS = HERE.parent / "results" / "functional_benchmark"
K_SWEEP = [10, 12, 25, 50, 75, 100]
SELECTION_ZIP = "results of natural selection analysis.zip"

# Layers excluded by the admissibility rule, kept here so the exclusion is reported rather
# than silent, with the software that produces each one.
PREDICTED_LAYERS = {
    "protein function": "MutPred2",
    "protein stability": "iMutant 2.0",
    "antigenicity and immunogenicity": "VaxiJen and IEDB",
    "physicochemical parameters": "ProtParam",
}


# --------------------------------------------------------------------------------------
# Reference layers
# --------------------------------------------------------------------------------------

def mutation_id_sites(directory: Path) -> dict[str, set[tuple[str, int]]]:
    """Mutation ID to the (gene, residue) pairs it covers, the mapping every layer reuses.

    Joining the other files on this identifier avoids reconstructing gene coordinates from
    genome positions, which fails on ORF1ab: the ribosomal frameshift breaks the affine
    relation between nucleotide position and residue number. A deletion covers a range of
    residues rather than one, so the value is a set; see `label_residues`.
    """
    table = pd.read_csv(directory / "mutations_summary.txt", sep="\t", dtype=str)
    table["kss_gene"] = table["Gene name"].map(COV2VAR_GENES)
    usable = table.dropna(subset=["kss_gene"])
    mapping = {}
    for identifier, gene, label in zip(usable["Mutation ID"], usable["kss_gene"],
                                       usable["Protein mutation-1 letter"]):
        sites = {(gene, residue) for residue in label_residues(label)}
        if sites:
            mapping[identifier] = sites
    return mapping


def selection_layers(directory: Path, id_to_site: dict) -> tuple[dict[str, set], set]:
    """The three selection tests, in one direction, with the universe they report on.

    Direction first. FUBAR gives two posteriors, and `Prob[alpha>beta]` is the one for
    *purifying* selection: filtering on it, as this did until 28 July 2026, keeps 20,586 rows
    of which every single one has beta < alpha. The diversifying posterior is
    `Prob[alpha<beta]`, and it is what belongs beside MEME. FEL reports a two-sided p-value and
    74.9% of its significant sites are purifying, so it is restricted by sign as well. All
    three layers are then the same question asked three ways, which is what the manuscript
    claims they are.

    Universe second. COV2Var ships only the significant rows of MEME and FEL, every p-value in
    both files being below 0.05, so the sites those tests examined and did not call cannot be
    recovered. FUBAR ships its non-significant rows. The universe returned here is the set of
    residues any of the three files reports on at all, which is a lower bound on what was
    tested and is still far better than treating every scored k-mer as evaluable.
    """
    archive = directory / SELECTION_ZIP
    if not archive.exists():
        return {}, set()

    def sites_of(ids) -> set:
        return {site for i in set(ids) for site in id_to_site.get(i, ())}

    layers: dict[str, set] = {}
    reported: set = set()
    with zipfile.ZipFile(archive) as bundle:
        names = {Path(n).name: n for n in bundle.namelist() if not n.endswith("/")}
        tables = {}
        for tag in ("meme", "fel", "fubar"):
            filename = f"selection_{tag}.txt"
            if filename not in names:
                return {}, set()
            with bundle.open(names[filename]) as handle:
                tables[tag] = pd.read_csv(
                    io.TextIOWrapper(handle, encoding="utf-8", errors="replace"),
                    sep="\t", dtype=str)
            reported |= sites_of(tables[tag]["Mutation ID"])

    meme = tables["meme"]
    layers["MEME (episodic)"] = sites_of(
        meme[pd.to_numeric(meme["P value"], errors="coerce") < 0.05]["Mutation ID"])

    fel = tables["fel"]
    fel_diversifying = fel[(pd.to_numeric(fel["P value"], errors="coerce") < 0.05)
                           & (pd.to_numeric(fel["Beta"], errors="coerce")
                              > pd.to_numeric(fel["Alpha"], errors="coerce"))]
    layers["FEL (pervasive)"] = sites_of(fel_diversifying["Mutation ID"])

    fubar = tables["fubar"]
    layers["FUBAR"] = sites_of(
        fubar[pd.to_numeric(fubar["Prob[alpha<beta]"], errors="coerce") > 0.95]["Mutation ID"])

    found = [layers[name] for name in ("MEME (episodic)", "FEL (pervasive)", "FUBAR")]
    layers["any of the three"] = set().union(*found)
    layers["all three"] = found[0] & found[1] & found[2]
    return layers, reported


def rbd_layers(directory: Path) -> dict[str, tuple[set, set]]:
    """The two measured layers, each with the set of residues its assay actually measured.

    The assayed set, not the interval it spans. Deep mutational scanning of the receptor
    binding domain covers 68 of the 189 residues between 334 and 522, so taking the endpoints
    as a proxy for coverage admits k-mers carrying no measurement at all and silently turns
    them into negatives, which is the one thing the evaluation universe must never do.
    """
    layers = {}

    # COV2Var ships these two with stray non-UTF-8 bytes in free-text columns.
    escape = pd.read_csv(directory / "rbd_immune_escape_score.txt", sep="\t", dtype=str,
                         encoding="utf-8", encoding_errors="replace")
    escape["site"] = pd.to_numeric(escape["Mutation site of protein"], errors="coerce")
    escape["value"] = pd.to_numeric(escape["Average mutation escape"], errors="coerce")
    escape = escape.dropna(subset=["site", "value"])
    assayed = {("S", int(s)) for s in escape["site"]}
    hits = escape[escape["value"] > 0.1]
    layers["DMS antibody escape"] = ({("S", int(s)) for s in hits["site"]}, assayed)

    binding = pd.read_csv(directory / "RBD binding affinity.txt", sep="\t", dtype=str,
                          encoding="utf-8", encoding_errors="replace")
    binding["site"] = pd.to_numeric(binding["Protein position of mutation"], errors="coerce")
    binding["p"] = pd.to_numeric(binding["P-value"], errors="coerce")
    binding = binding.dropna(subset=["site", "p"])
    assayed_b = {("S", int(s)) for s in binding["site"]}
    hits_b = binding[binding["p"] < 0.05]
    layers["DMS ACE2 binding"] = ({("S", int(s)) for s in hits_b["site"]}, assayed_b)
    return layers


def fitness_layer(path: Path, kmers: pd.DataFrame, quantile: float) -> set:
    """Bloom and Neher, top `quantile` of sites by median absolute effect.

    Kept because the prespecification kept it, and read as constraint rather than as
    functional relevance: the top decile by magnitude is entirely deleterious.
    """
    fitness = pd.read_csv(path)
    confident = fitness[fitness["expected_count"] >= MIN_EXPECTED_COUNT].copy()
    confident["magnitude"] = confident["fitness"].abs()
    sites = confident.groupby(["gene", "aa_site"])["magnitude"].median().reset_index()
    sites = sites[sites["gene"].isin(set(kmers["gene"]))]
    cutoff = sites["magnitude"].quantile(1 - quantile)
    hits = sites[sites["magnitude"] >= cutoff]
    return set(zip(hits["gene"], hits["aa_site"]))


# --------------------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------------------

def annotation_layers(directory: Path, id_to_site: dict) -> dict[str, tuple[set, set, str]]:
    """The per-mutation annotation files: the sites each one flags, the sites it scored at all,
    and the software that produces it.

    These are predictions, so they are reported and never admitted as evidence: our own
    mutational term is a prediction too, and two predictors agreeing says they agree. They are
    computed all the same, because they are the only genome-wide layers COV2Var offers that
    are selective rather than saturated, and a reader is entitled to see that.

    The scored set matters as much as the flagged one. MutPred2 and iMutant cover the missense
    mutations of every gene, but the antigenicity file is spike and nothing else, 671 rows all
    of them S, so evaluating it over the whole genome would score four genes as negatives for
    never having been looked at.
    """
    def sites_of(ids) -> set:
        return {site for i in set(ids) for site in id_to_site.get(i, ())}

    layers: dict[str, tuple[set, set, str]] = {}

    function = pd.read_csv(directory / "mutation induced protein function changes.txt",
                           sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
    score = pd.to_numeric(function["Score"], errors="coerce")
    layers["MutPred2 functional impact, score >= 0.5"] = (
        sites_of(function[score >= 0.5]["Mutation ID"]),
        sites_of(function["Mutation ID"]), "MutPred2")

    stability = pd.read_csv(directory / "mutation induced protein stability changes.txt",
                            sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
    column = next(c for c in stability.columns if "ddg" in c.lower())
    ddg = pd.to_numeric(stability[column], errors="coerce").abs()
    layers["iMutant stability change, |ddG| >= 1"] = (
        sites_of(stability[ddg >= 1]["Mutation ID"]),
        sites_of(stability["Mutation ID"]), "iMutant 2.0")

    antigen = pd.read_csv(directory / "mutation effects on antigenicity and immunogenicity.txt",
                          sep="\t", dtype=str, encoding="utf-8", encoding_errors="replace")
    flag = antigen["Significant changes in antigenicity"].astype(str).str.strip().str.lower()
    layers["VaxiJen antigenicity change, significant"] = (
        sites_of(antigen[flag.str.startswith(("y", "t", "sig"))]["Mutation ID"]),
        sites_of(antigen["Mutation ID"]), "VaxiJen")
    return layers


def substitutions_only(kmers: pd.DataFrame) -> pd.DataFrame:
    """The ranking with its insertion and deletion positions removed.

    Some layers annotate substitutions and nothing else, so against them an indel position is
    unevaluable, in the same sense as an unassayed residue. Dropping those k-mers from the
    universe and re-taking the top k is what keeps the comparison fair in both directions: it
    neither counts them as errors nor lets them occupy the selection without being judged.
    """
    return kmers[kmers["mutational_score"] < INDEL_MUTATIONAL_SCORE]


def restrict(kmers: pd.DataFrame, assayed: set) -> pd.DataFrame:
    """The k-mers an assay evaluates: at least one of the three codons carries a measurement.

    The prespecified rule, applied to coverage rather than to the interval: a k-mer with no
    measured residue is unevaluable and leaves the denominator, it does not become a negative.
    """
    evaluable = [any((row.gene, row.aa_first + offset) in assayed
                     for offset in range(K_NT // 3))
                 for row in kmers.itertuples()]
    return kmers[pd.Series(evaluable, index=kmers.index)]


def evaluate(kmers: pd.DataFrame, sites: set, name: str, method: str = "KSS") -> pd.DataFrame:
    positive = label(kmers, sites)
    total_positive = int(positive.sum())
    universe = len(kmers)
    if total_positive in (0, universe):
        return pd.DataFrame()

    order = ranked(kmers, method)
    auroc_test = mannwhitneyu(kmers[method][positive], kmers[method][~positive],
                              alternative="greater")
    auroc = auroc_test.statistic / (total_positive * (universe - total_positive))

    # A flat expectation assumes the selection is drawn from the universe at random, and it is
    # not: the universe is mostly ORF1ab while the top of the ranking is mostly spike, and
    # prevalence differs between genes. The gene-matched expectation draws the same number of
    # k-mers per gene as the selection contains, which is the null the confounder allows.
    by_gene = positive.groupby(kmers["gene"]).mean()

    rows = []
    # Half the universe, not all of it: past that a cut selects most of what exists and the
    # rates stop describing a prioritization. The receptor binding domain layers, with 63
    # evaluable k-mers, are the ones this protects.
    for k in [k for k in K_SWEEP if k <= universe // 2]:
        row = confusion(order, positive, k)
        # Exact, not Monte Carlo: the chance of drawing at least this many positives in k
        # draws without replacement from the universe.
        row["p_hypergeometric"] = float(
            hypergeom.sf(row["TP"] - 1, universe, total_positive, k))
        row["p_perfect_by_chance"] = float(hypergeom.pmf(k, universe, total_positive, k))
        selected = kmers.loc[order[:k]]
        row["expected_flat"] = round(total_positive / universe * k, 4)
        row["expected_gene_matched"] = round(
            float(sum(by_gene[gene] for gene in selected["gene"])), 4)
        row.update({"reference": name, "universe": universe, "method": method,
                    "auroc": auroc, "auroc_p": float(auroc_test.pvalue)})
        rows.append(row)

    hits = positive.loc[order].to_numpy()
    perfect = int((~hits).argmax()) if (~hits).any() else universe
    for row in rows:
        row["perfect_prefix"] = perfect
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cov2var", type=Path, required=True)
    parser.add_argument("--fitness", type=Path, default=None)
    parser.add_argument("--method", default="KSS")
    args = parser.parse_args()
    # Identity, not presence: the COV2Var site serves no version tag, so the checksums are
    # the only thing that says which download these numbers came from.
    problem = verify_cov2var(args.cov2var)
    if problem:
        print(problem)
        return 1
    # The fitness file is optional, but naming one and getting the layer silently dropped is
    # worse than not naming it: the run would still write a matrix, one layer short, with
    # nothing in the output saying so. A path that was given has to exist and has to match its
    # pin, which is what the letter says of every pinned input.
    if args.fitness is not None:
        if not args.fitness.exists():
            print(f"--fitness {args.fitness} does not exist")
            return 1
        problem = verify_fitness(args.fitness)
        if problem:
            print(problem)
            return 1

    RESULTS.mkdir(exist_ok=True)
    kmers = load_rankings()
    id_to_site = mutation_id_sites(args.cov2var)

    print("=" * 96)
    print("Reference layers, admissibility, and coverage")
    print("=" * 96)
    print(f"  universe: {len(kmers)} scored k-mers over {kmers['gene'].nunique()} genes")
    for layer, software in PREDICTED_LAYERS.items():
        print(f"  excluded: {layer:<34} predicted by {software}")

    frames = []
    substitution_universe = substitutions_only(kmers)
    print(f"  of these, {len(kmers) - len(substitution_universe)} carry an insertion or a "
          f"deletion, so {len(substitution_universe)} are evaluable against a layer that "
          f"annotates substitutions only")

    def add(sites, name, nature, admitted, universe=kmers):
        frame = evaluate(universe, sites, name, args.method)
        if not frame.empty:
            frame["nature"] = nature
            frame["admitted"] = admitted
        frames.append(frame)

    membership = set().union(*id_to_site.values())
    add(membership, "COV2Var membership", "frequency-defined catalogue", True)

    beyond_membership = set()
    selection, reported_on = selection_layers(args.cov2var, id_to_site)
    tested = restrict(kmers, reported_on)
    print(f"  selection analysis: reports on {len(reported_on)} residues, "
          f"{len(tested)} of {len(kmers)} k-mers evaluable")
    for name, sites in selection.items():
        add(sites, f"selection: {name}", "test on observed sequences", True, universe=tested)
        beyond_membership |= sites

    for name, (sites, assayed) in rbd_layers(args.cov2var).items():
        # Two restrictions, and both are coverage rather than judgement: the assay measured 68
        # residues, and it measured substitutions at them.
        covered = restrict(substitutions_only(kmers), assayed)
        residues = sorted(site for _, site in assayed)
        print(f"  {name}: {len(assayed)} residues measured between spike {residues[0]} and "
              f"{residues[-1]}, {len(covered)} of {len(kmers)} k-mers evaluable")
        add(sites, f"{name} (spike RBD, measured residues)", "measurement", True,
            universe=covered)
        beyond_membership |= sites

    for name, (sites, scored, software) in annotation_layers(args.cov2var, id_to_site).items():
        covered = restrict(substitution_universe, scored)
        print(f"  {software}: {len(covered)} of {len(kmers)} k-mers carry a scored residue")
        add(sites, name, f"software prediction ({software})", False, universe=covered)

    if args.fitness is not None:
        for quantile, tag in [(0.05, "top 5%"), (0.10, "top 10%"), (0.20, "top 20%")]:
            sites = fitness_layer(args.fitness, kmers, quantile)
            add(sites, f"fitness magnitude, {tag}", "constraint, reported separately", False)
    else:
        print("  fitness layer skipped: pass --fitness aa_fitness.csv to include it")

    # Two unions, and the second is the informative one. Every other layer is keyed on the
    # 9,832 catalogued mutations, so the first union is the membership layer again, to the
    # site. The second drops membership and asks what is left: how much of a selection this
    # resource says anything about beyond the fact that the mutation is frequent.
    add(membership | beyond_membership, "any admitted layer", "union of the layers above",
        True)
    add(beyond_membership, "any layer beyond membership",
        "union, frequency catalogue excluded", True)

    table = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    table.to_csv(RESULTS / "reference_matrix.csv", index=False)

    union = table[table["reference"] == "any admitted layer"]
    print("\n" + "=" * 96)
    print("Coverage of the selection itself: how many of the top k no layer annotates")
    print("=" * 96)
    print(f"  {'k':>5} {'annotated':>10} {'by nothing':>11}")
    for row in union.itertuples():
        print(f"  {row.k:>5} {row.TP:>10} {row.k - row.TP:>11}")

    print("\n" + "=" * 96)
    print(f"Rates for {args.method} at every selection size, per reference layer")
    print("=" * 96)
    for name, block in table.groupby("reference", sort=False):
        head = block.iloc[0]
        print(f"\n{name}")
        print(f"  universe {int(head['universe'])}, prevalence {head['prevalence']:.1%}, "
              f"AUROC {head['auroc']:.3f} (p {head['auroc_p']:.1e}), "
              f"perfect to k={int(head['perfect_prefix'])}")
        print(f"  {'k':>5} {'TP':>5} {'FP':>4} {'precision':>10} {'FDR':>7} {'FPR':>7} "
              f"{'recall':>7} {'enrich':>7} {'p(>=TP)':>10} {'p(k/k)':>9}")
        for row in block.itertuples():
            print(f"  {row.k:>5} {row.TP:>5} {row.FP:>4} {row.precision_at_k:>10.3f} "
                  f"{row.fdr_at_k:>7.3f} {row.fpr_at_k:>7.3f} {row.recall_at_k:>7.3f} "
                  f"{row.enrichment:>7.2f} {row.p_hypergeometric:>10.2e} "
                  f"{row.p_perfect_by_chance:>9.2e}")

    print(f"\nwritten: {RESULTS / 'reference_matrix.csv'}")
    print("\nFDR is 1 - precision@k. FPR is FP / (FP + TN) over the layer's own negatives.")
    print("They answer different questions and must not be quoted for one another.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
