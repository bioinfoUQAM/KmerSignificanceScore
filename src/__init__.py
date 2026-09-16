"""
K-mer Significance Score (KSS) Package

A comprehensive toolkit for computing k-mer significance scores in genomic sequences,
combining mutational impact, discriminative power, and protein characterization depth.

Main modules:
    - pipeline: Configuration loading and per-gene orchestration
    - kanalyzer: Sequence alignment and mutation identification
    - discriminative_score: Class discrimination metrics
    - mutation_score: Amino acid substitution scoring
    - protein_score: Protein characterization depth from UniProt
    - kss: Component combination, compilation and final scores
    - report: Optional PDF report of the highest-scoring positions
    - utils: Data I/O utilities

See CONTRIBUTING.md for how these fit together.
"""

# Kept equal to the version pyproject.toml declares, which the cross-check enforces.
__version__ = "1.1.0"
__author__ = "KSS Development Team"

# Import main functions for convenient top-level access
from .kss import compute_kss_scores, compile_results, get_taxon_id
from .kanalyzer import analyze_records
from .utils import load_data_from_json, save_data_as_json

__all__ = [
    "compute_kss_scores",
    "compile_results",
    "get_taxon_id",
    "analyze_records",
    "load_data_from_json",
    "save_data_as_json",
]
