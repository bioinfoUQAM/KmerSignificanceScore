#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Amino Acid Utilities Module

This module provides utility functions for amino acid analysis, including
codon distance calculations and property-based comparisons.

NOTE: This module is used exclusively by the evaluation notebooks
(matrix_evaluation.ipynb and matrix_optimization.ipynb). It is NOT part of
the main KSS application (src/).
"""

import numpy as np
from itertools import combinations

# Standard genetic code (NCBI translation table 1)
CODON_TABLE = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT']
}

# Standard amino acids (sorted for consistency)
AMINO_ACIDS = sorted(CODON_TABLE.keys())

# Mapping between 3-letter and 1-letter codes
AA_3_TO_1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
}

AA_1_TO_3 = {v: k for k, v in AA_3_TO_1.items()}


def compute_codon_distance(codon1, codon2):
    """
    Calculate Hamming distance between two codons.

    Args:
        codon1 (str): First codon (3 nucleotides)
        codon2 (str): Second codon (3 nucleotides)

    Returns:
        int: Number of differing positions (0-3)

    Examples:
        >>> compute_codon_distance('AAA', 'AAG')
        1
        >>> compute_codon_distance('AAA', 'GGG')
        3
    """
    if len(codon1) != 3 or len(codon2) != 3:
        raise ValueError("Codons must be exactly 3 nucleotides long")
    return sum(b1 != b2 for b1, b2 in zip(codon1, codon2))


def get_min_codon_distance(aa1, aa2, codon_table=None):
    """
    Get minimum codon distance between two amino acids.

    Args:
        aa1 (str): First amino acid (1-letter code)
        aa2 (str): Second amino acid (1-letter code)
        codon_table (dict, optional): Custom codon table. Defaults to standard genetic code.

    Returns:
        int: Minimum Hamming distance between any codon pair (1-3)

    Examples:
        >>> get_min_codon_distance('K', 'R')  # AAA->AGA
        1
        >>> get_min_codon_distance('M', 'W')  # ATG->TGG
        2
    """
    if codon_table is None:
        codon_table = CODON_TABLE

    if aa1 not in codon_table or aa2 not in codon_table:
        raise ValueError(f"Amino acids must be valid 1-letter codes: {aa1}, {aa2}")

    codons1 = codon_table[aa1]
    codons2 = codon_table[aa2]

    return min(
        compute_codon_distance(c1, c2)
        for c1 in codons1
        for c2 in codons2
    )


def calculate_property_distances(property_values, pairs):
    """
    Calculate absolute differences in property values for amino acid pairs.

    Args:
        property_values (dict): Mapping of amino acid (1-letter) to property value
        pairs (list): List of amino acid pairs as strings (e.g., ['AC', 'AD', ...])

    Returns:
        numpy.ndarray: Array of absolute property differences

    Examples:
        >>> props = {'A': 1.0, 'C': 2.5, 'D': 3.0}
        >>> pairs = ['AC', 'AD', 'CD']
        >>> calculate_property_distances(props, pairs)
        array([1.5, 2.0, 0.5])
    """
    if not property_values:
        raise ValueError("property_values cannot be empty")

    distances = []
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError(f"Each pair must be exactly 2 amino acids: {pair}")

        aa1, aa2 = pair[0], pair[1]

        if aa1 not in property_values or aa2 not in property_values:
            raise KeyError(f"Amino acid not found in properties: {pair}")

        distance = abs(property_values[aa1] - property_values[aa2])
        distances.append(distance)

    return np.array(distances)


def generate_amino_acid_pairs(amino_acids=None):
    """
    Generate all unique pairs of amino acids.

    Args:
        amino_acids (list, optional): List of amino acids. Defaults to standard 20.

    Returns:
        list: List of amino acid pairs as strings (e.g., ['AC', 'AD', ...])

    Examples:
        >>> pairs = generate_amino_acid_pairs(['A', 'C', 'D'])
        >>> len(pairs)
        3
        >>> 'AC' in pairs
        True
    """
    if amino_acids is None:
        amino_acids = AMINO_ACIDS

    return [''.join(pair) for pair in combinations(sorted(amino_acids), 2)]


def precompute_codon_distances(amino_acids=None, codon_table=None):
    """
    Precompute minimum codon distances for all amino acid pairs.

    Args:
        amino_acids (list, optional): List of amino acids. Defaults to standard 20.
        codon_table (dict, optional): Custom codon table. Defaults to standard genetic code.

    Returns:
        dict: Mapping of amino acid pairs to minimum codon distances

    Examples:
        >>> distances = precompute_codon_distances(['A', 'C', 'D'])
        >>> len(distances)
        3
        >>> distances['AC']  # Example distance
        2
    """
    if amino_acids is None:
        amino_acids = AMINO_ACIDS

    if codon_table is None:
        codon_table = CODON_TABLE

    pairs = generate_amino_acid_pairs(amino_acids)

    return {
        pair: get_min_codon_distance(pair[0], pair[1], codon_table)
        for pair in pairs
    }
