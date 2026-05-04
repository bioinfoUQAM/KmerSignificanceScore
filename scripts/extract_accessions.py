"""Extract accession + class metadata from sequences.fasta files into accession TSV lists.

Output:
- data/accessions/sars_cov2.tsv      (one row per accession, variant)
- data/accessions/hiv1.tsv           (one row per accession, subtype)
- data/accessions/hcmv_<GENE>.tsv    (per-gene because HCMV genotypes vary by gene)
"""
from collections import OrderedDict
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_DIR = DATA_DIR / "accessions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

VIRUS_DIRS = {
    "sars_cov2": ("Severe_acute_respiratory_syndrome_coronavirus_2", "variant"),
    "hiv1":       ("Human_immunodeficiency_virus_1", "subtype"),
    "hcmv":       ("Human_betaherpesvirus_5", "genotype"),
}


def parse_headers(fasta_path: Path):
    """Yield (accession, class_label) from FASTA headers (>acc|virus|gene|class)."""
    with fasta_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\n").split("|")
            if len(parts) < 4:
                continue
            yield parts[0], parts[3]


def consolidate_per_virus(short: str, virus_folder: str, label: str):
    """For a virus where all gene FASTAs share the same accession+class set,
    produce a single TSV listing each accession once."""
    virus_dir = DATA_DIR / virus_folder
    if not virus_dir.exists():
        return None
    seen = OrderedDict()
    for gene_dir in sorted(virus_dir.iterdir()):
        if not gene_dir.is_dir():
            continue
        fasta = gene_dir / "sequences.fasta"
        if not fasta.exists():
            continue
        for acc, cls in parse_headers(fasta):
            if acc not in seen:
                seen[acc] = cls
    out = OUT_DIR / f"{short}.tsv"
    with out.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"accession\t{label}\n")
        for acc, cls in seen.items():
            fh.write(f"{acc}\t{cls}\n")
    return out, len(seen)


def per_gene(short: str, virus_folder: str, label: str):
    """For HCMV, produce one TSV per gene (genotypes differ by gene)."""
    virus_dir = DATA_DIR / virus_folder
    results = []
    for gene_dir in sorted(virus_dir.iterdir()):
        if not gene_dir.is_dir():
            continue
        fasta = gene_dir / "sequences.fasta"
        if not fasta.exists():
            continue
        out = OUT_DIR / f"{short}_{gene_dir.name}.tsv"
        n = 0
        with out.open("w", encoding="utf-8", newline="\n") as fh:
            fh.write(f"accession\t{label}\n")
            for acc, cls in parse_headers(fasta):
                fh.write(f"{acc}\t{cls}\n")
                n += 1
        results.append((out, gene_dir.name, n))
    return results


def main():
    # Wipe previous TSVs to remove the old per-gene files for SARS/HIV
    for old in OUT_DIR.glob("*.tsv"):
        old.unlink()

    # SARS-CoV-2 and HIV-1: one TSV per virus (accessions identical across genes)
    for short in ("sars_cov2", "hiv1"):
        folder, label = VIRUS_DIRS[short]
        out, n = consolidate_per_virus(short, folder, label)
        print(f"{out.name}: {n:,} accessions")

    # HCMV: one TSV per gene
    folder, label = VIRUS_DIRS["hcmv"]
    for out, gene, n in per_gene("hcmv", folder, label):
        print(f"{out.name}: {n:,} accessions ({gene})")


if __name__ == "__main__":
    main()
