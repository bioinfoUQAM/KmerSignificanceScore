"""Download nucleotide sequences from NCBI Entrez for the accession lists in
data/accessions/ and write them as data/<Virus>/<Gene>/sequences.fasta with the
header format expected by KSS: ``>accession|Virus|Gene|class``.

Usage
-----
Set your NCBI contact email (required by Entrez) and optionally an API key, then
fetch one virus at a time::

    python scripts/fetch_sequences.py --email you@example.com --virus hcmv
    python scripts/fetch_sequences.py --email you@example.com --virus hiv1
    python scripts/fetch_sequences.py --email you@example.com --virus sars_cov2

Notes
-----
Each accession resolves to a complete genome record. Gene regions are extracted
from the matching reference annotation provided in
``data/<Virus>/<Gene>/<reference>.gb`` via pairwise alignment.

For SARS-CoV-2 (~278k sequences) this is intentionally slow. NCBI rate limits
apply (3 requests/sec without an API key, 10/sec with). Consider using the NCBI
``datasets`` CLI or an EBI ENA bulk download for production-scale runs; the
script is provided for transparent, reproducible small-to-medium downloads.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parent.parent
ACCESSIONS_DIR = ROOT / "data" / "accessions"
DATA_DIR = ROOT / "data"

VIRUSES = {
    "sars_cov2": {
        "folder": "Severe_acute_respiratory_syndrome_coronavirus_2",
        "label_in_header": "Severe_acute_respiratory_syndrome_coronavirus_2",
        "tsv": "sars_cov2.tsv",
        "genes": ["E", "M", "N", "ORF1ab", "S"],
        "per_gene": False,
    },
    "hiv1": {
        "folder": "Human_immunodeficiency_virus_1",
        "label_in_header": "HIV-1",
        "tsv": "hiv1.tsv",
        "genes": ["env", "gag", "pol"],
        "per_gene": False,
    },
    "hcmv": {
        "folder": "Human_betaherpesvirus_5",
        "label_in_header": "Human_betaherpesvirus_5",
        "tsv_pattern": "hcmv_{gene}.tsv",
        "genes": ["UL55", "UL73", "US28"],
        "per_gene": True,
    },
}


def read_accession_tsv(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            acc, cls = line.split("\t", 1)
            rows.append((acc, cls))
    return rows


def chunked(seq: list, n: int) -> Iterable[list]:
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def fetch_genomes(accessions: list[str], email: str, api_key: str | None,
                  batch: int = 200) -> dict[str, str]:
    """Return ``{accession: full-length nucleotide sequence}`` from NCBI."""
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        sys.exit("Biopython is required: pip install biopython")

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    delay = 0.11 if api_key else 0.34  # respect rate limits
    out: dict[str, str] = {}
    total = len(accessions)
    for i, group in enumerate(chunked(accessions, batch)):
        print(f"  fetching {i*batch + 1}-{min((i+1)*batch, total)} / {total}", flush=True)
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(group), rettype="fasta", retmode="text"
        )
        for rec in SeqIO.parse(handle, "fasta"):
            key = rec.id.split(".")[0]  # match accessions with or without version
            out[key] = str(rec.seq)
        handle.close()
        time.sleep(delay)
    return out


def extract_gene(genome: str, reference_gb: Path, gene: str) -> str | None:
    """Locate ``gene`` in ``genome`` by aligning to the gene region of the
    reference annotation. Returns the matching subsequence or None."""
    try:
        from Bio import SeqIO
        from Bio.Align import PairwiseAligner
    except ImportError:
        sys.exit("Biopython is required: pip install biopython")

    ref = SeqIO.read(reference_gb, "genbank")
    target = None
    for feat in ref.features:
        if feat.type != "gene":
            continue
        names = feat.qualifiers.get("gene", []) + feat.qualifiers.get("locus_tag", [])
        if any(n.lower() == gene.lower() for n in names):
            target = str(feat.extract(ref.seq))
            break
    if target is None:
        return None

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    try:
        aln = aligner.align(genome.upper(), target.upper())[0]
    except Exception:
        return None
    start, end = aln.aligned[0][0][0], aln.aligned[0][-1][1]
    return genome[start:end]


def write_fasta(path: Path, virus_label: str, gene: str,
                rows: list[tuple[str, str]], sequences: dict[str, str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for acc, cls in rows:
            key = acc.split(".")[0]
            seq = sequences.get(key)
            if not seq:
                continue
            fh.write(f">{acc}|{virus_label}|{gene}|{cls}\n{seq}\n")
            n_written += 1
    return n_written


def fetch_virus(short: str, email: str, api_key: str | None) -> None:
    cfg = VIRUSES[short]
    virus_dir = DATA_DIR / cfg["folder"]
    print(f"\n=== {short} ===")
    if cfg["per_gene"]:
        for gene in cfg["genes"]:
            tsv = ACCESSIONS_DIR / cfg["tsv_pattern"].format(gene=gene)
            rows = read_accession_tsv(tsv)
            print(f"[{gene}] {len(rows):,} accessions")
            seqs = fetch_genomes([a for a, _ in rows], email, api_key)
            ref_dir = virus_dir / gene
            ref_gb = next(ref_dir.glob("*.gb"), None)
            extracted: dict[str, str] = {}
            for acc, _ in rows:
                key = acc.split(".")[0]
                if key not in seqs:
                    continue
                if ref_gb is not None:
                    sub = extract_gene(seqs[key], ref_gb, gene)
                    if sub:
                        extracted[key] = sub
                else:
                    extracted[key] = seqs[key]
            out = ref_dir / "sequences.fasta"
            n = write_fasta(out, cfg["label_in_header"], gene, rows, extracted)
            print(f"  wrote {n:,} sequences -> {out}")
    else:
        tsv = ACCESSIONS_DIR / cfg["tsv"]
        rows = read_accession_tsv(tsv)
        print(f"{len(rows):,} accessions, fetching once and reusing across {len(cfg['genes'])} genes")
        seqs = fetch_genomes([a for a, _ in rows], email, api_key)
        for gene in cfg["genes"]:
            ref_dir = virus_dir / gene
            ref_gb = next(ref_dir.glob("*.gb"), None)
            extracted: dict[str, str] = {}
            for acc, _ in rows:
                key = acc.split(".")[0]
                if key not in seqs:
                    continue
                if ref_gb is not None:
                    sub = extract_gene(seqs[key], ref_gb, gene)
                    if sub:
                        extracted[key] = sub
                else:
                    extracted[key] = seqs[key]
            out = ref_dir / "sequences.fasta"
            n = write_fasta(out, cfg["label_in_header"], gene, rows, extracted)
            print(f"  [{gene}] wrote {n:,} sequences -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--email", required=True, help="NCBI contact email")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    ap.add_argument(
        "--virus", choices=list(VIRUSES.keys()) + ["all"], default="all",
        help="Which virus to fetch (default: all)",
    )
    args = ap.parse_args()
    targets = list(VIRUSES) if args.virus == "all" else [args.virus]
    for short in targets:
        fetch_virus(short, args.email, args.api_key)


if __name__ == "__main__":
    main()
