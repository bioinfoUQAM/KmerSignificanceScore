#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Discriminative Scoring Module - KSS Method

This module calculates discriminative scores for genomic positions using
k-mer based classification. The score quantifies how well k-mer patterns
at a position discriminate between different classes (e.g., viral strains).

Scoring Formula:
    KSS = tanh(sqrt(NMI × weighted_purity²) × ln(n_classes + 1))

    where:
    - NMI = Normalized Mutual Information between k-mers and classes
    - weighted_purity² = Σ P(kmer_j) × [max_c P(class_c | kmer_j)]²
    - sqrt() = Geometric mean to reduce over-penalization
    - ln(n_classes + 1) = Class complexity adjustment

The final score is normalized to [0, 1] range, where:
    - 1.0 = Perfect discrimination (k-mers uniquely identify classes)
    - 0.0 = No discriminative power (k-mers independent of classes)

Main function:
    get_discriminative_score: Calculate discriminative score for a genomic position
"""

import numpy as np
from typing import Dict, Any


def get_discriminative_score(
    X: np.ndarray,
    y: np.ndarray
) -> Dict[str, Any]:
    """Discriminative score for one position: how well its k-mers separate the classes.

    Combines normalized mutual information with class-weighted k-mer purity, scaled by
    ln(n_classes + 1) and mapped into [0, 1] by tanh (see the module docstring for the
    formula). Needs at least two classes; returns 0 otherwise.

    Args:
        X: binary k-mer matrix (n_sequences, n_kmers), at most one 1 per row. A sequence whose
            k-mer did not reach the prevalence threshold has no column and contributes a zero row.
        y: class label per sequence (n_sequences,).

    Returns:
        dict with 'raw_score', 'normalized_score' (the [0, 1] score), 'n_classes',
        'n_sequences' and 'n_kmers'.
    """
    X = np.asarray(X)
    y = np.asarray(y)

    n_sequences, n_kmers = X.shape
    unique_classes = np.unique(y)
    n_classes = len(unique_classes)

    # Edge cases
    if n_kmers == 0 or n_classes < 2:
        return {
            'raw_score': 0.0,
            'normalized_score': 0.0,
            'n_classes': n_classes,
            'n_sequences': n_sequences,
            'n_kmers': n_kmers
        }

    # Class complexity factor: ln(n_classes + 1)
    class_complexity = np.log(n_classes + 1)

    # ========================================================================
    # STEP 1: Calculate H(Y) - Entropy of classes
    # ========================================================================
    class_probs = np.array([np.mean(y == c) for c in unique_classes])
    h_class = -np.sum(class_probs * np.log(class_probs + 1e-10))

    if h_class == 0:
        return {
            'raw_score': 0.0,
            'normalized_score': 0.0,
            'n_classes': n_classes,
            'n_sequences': n_sequences,
            'n_kmers': n_kmers
        }

    # ========================================================================
    # STEP 2: Calculate H(Y|X) - Conditional entropy by k-mer configurations
    # ========================================================================
    unique_rows, inverse = np.unique(X, axis=0, return_inverse=True)
    h_conditional = 0.0

    for i in range(len(unique_rows)):
        mask = (inverse == i)
        p_config = np.mean(mask)
        if p_config > 0:
            y_subset = y[mask]
            subset_probs = np.array([np.mean(y_subset == c) for c in unique_classes])
            h_subset = -np.sum(subset_probs * np.log(subset_probs + 1e-10))
            h_conditional += p_config * h_subset

    # Normalized Mutual Information
    nmi = float(np.clip((h_class - h_conditional) / h_class, 0, 1))

    # ========================================================================
    # STEP 3: Calculate boosted weighted purity
    # ========================================================================
    weighted_purity_boosted = 0.0

    for j in range(n_kmers):
        mask_kmer = (X[:, j] == 1)
        freq_j = np.mean(mask_kmer)

        if freq_j > 0:
            y_kmer = y[mask_kmer]
            probs_j = np.array([np.mean(y_kmer == c) for c in unique_classes])
            purity_j = np.max(probs_j)
            weighted_purity_boosted += freq_j * (purity_j ** 2)

    # ========================================================================
    # STEP 4: Calculate final score
    # ========================================================================
    raw_score = np.sqrt(nmi * weighted_purity_boosted)
    x = raw_score * class_complexity
    normalized_score = float(np.clip(np.tanh(x), 0, 1))

    return {
        'raw_score': round(float(raw_score), 4),
        'normalized_score': round(normalized_score, 3),
        'n_classes': n_classes,
        'n_sequences': n_sequences,
        'n_kmers': n_kmers
    }
