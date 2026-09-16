#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mutational impact scoring from amino acid substitution matrices.

Scores an amino acid change by the biochemical dissimilarity of the two residues, read from a
substitution matrix (MIYATA_EVO, BLOSUM, GRANTHAM, PAM, ...) scaled to [0, 1]. Matrices load
from and save to JSON. Main entry point: get_mutational_scores.
"""

import os
import json
import numpy as np
from typing import Dict, Any, List, Tuple, Optional, Union


# Cache for substitution matrices to avoid repeated file I/O
_MATRIX_CACHE = {}


def load_substitution_matrix(filename: str) -> Dict[Tuple[str, str], float]:
    """
    Load a substitution matrix from a JSON file.

    Args:
        filename: Path to the matrix file

    Returns:
        Substitution matrix keyed by amino acid pair. Example: {('A', 'G'): 0.45, ...}

    Raises:
        FileNotFoundError: If the file does not exist
        RuntimeError: If the file cannot be parsed

    Notes:
        JSON is the only format the matrices ship in. Keys are converted from strings to
        tuples, and both '(A, G)' and 'AG' key formats are accepted.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    return _load_matrix_json(filename)


def _load_matrix_json(filename: str) -> Dict[Tuple[str, str], float]:
    """
    Load substitution matrix from JSON file.

    Args:
        filename: Path to JSON (.json) file

    Returns:
        Substitution matrix dictionary with tuple keys

    Raises:
        RuntimeError: If JSON loading or parsing fails

    Notes:
        - Converts string keys to tuple format
        - Handles both '(A, G)' and 'AG' key formats
    """
    try:
        with open(filename, 'r') as f:
            matrix_data = json.load(f)
        
        matrix = {}
        for key, value in matrix_data.items():
            if ',' in key:
                aa_pair = tuple(aa.strip("'\" ()") for aa in key.strip().split(','))
            else:
                aa_pair = tuple(key.strip())
            
            matrix[aa_pair] = float(value)
        
        return matrix
        
    except (FileNotFoundError, json.JSONDecodeError) as e:
        raise RuntimeError(f"Error loading JSON file {filename}: {e}")


def save_substitution_matrix(matrix: Dict[Tuple[str, str], float],
                            base_filename: str) -> str:
    """
    Save a substitution matrix as JSON.

    JSON is the only format written: it is open, readable and language-agnostic, and it
    stores float64 without loss.

    Args:
        matrix: Substitution matrix dictionary with tuple keys
        base_filename: Base filename without extension (e.g., 'matrices/BLOSUM62')

    Returns:
        Path of the file created. Parent directories are created if needed.
    """
    os.makedirs(os.path.dirname(base_filename), exist_ok=True)

    json_path = f"{base_filename}.json"
    _save_matrix_json(matrix, json_path)

    return json_path


def _save_matrix_json(matrix: Dict[Tuple[str, str], float], filepath: str) -> None:
    """
    Save matrix in JSON format with string keys.

    Args:
        matrix: Substitution matrix dictionary
        filepath: Complete path including .json extension

    Notes:
        - Converts tuple keys to string format '(A, G)'
        - Values are written unrounded. json.dump emits the shortest decimal string that
          reads back as the same float64, so the round trip is exact. An earlier version
          rounded to 15 decimal places, which silently altered 216 of the 400 entries of
          MIYATA_EVO by up to 4e-16. See tests/test_matrix_formats.py.
    """
    string_matrix = {}
    for (aa1, aa2), value in matrix.items():
        key = f"({aa1}, {aa2})"
        string_matrix[key] = float(value)

    with open(filepath, "w") as outfile:
        json.dump(string_matrix, outfile, indent=4, sort_keys=True)


def adjust_and_scale_substitution_matrix(matrix: Dict[Tuple[str, str], float],
                                        matrix_type: str) -> Dict[Tuple[str, str], float]:
    """
    Adjust and scale a substitution matrix to [0, 1] range.

    Distance-based matrices (MIYATA, GRANTHAM, etc.) are inverted so that
    similar amino acids receive higher scores. All matrices are then scaled
    using min-max normalization to ensure values fall in [0, 1] range.

    Args:
        matrix: Original substitution matrix with tuple keys
        matrix_type: Name/type of the matrix (e.g., "MIYATA", "BLOSUM62")
                    Used to determine if inversion is needed

    Returns:
        Scaled substitution matrix with values in [0, 1] range.
        Higher values indicate more conservative substitutions (or smaller
        distances for inverted matrices).

    Notes:
        - Distance matrices are inverted: new_value = max_value - old_value
        - Ensures symmetry: if (A,G) exists, (G,A) will have same value
        - Uses global min-max normalization to preserve symmetry
        - Distance matrix types: MIYATA, GRANTHAM, SNEATH, and variants
    """
    # Handle distance-based matrices by inverting values
    distance_matrices = {
        "MIYATA", "GRANTHAM", "SNEATH", "MIYATA_EVO", 
        "MIYATA_BEST_GLOBAL", "MIYATA_OPTIMIZED"
    }
    
    if matrix_type.upper() in distance_matrices:
        max_value = max(matrix.values())
        matrix = {key: max_value - value for key, value in matrix.items()}
    
    # Extract all amino acids from the matrix
    amino_acids = sorted(set(aa for pair in matrix for aa in pair))
    
    # Build a 2D numpy array for scaling
    matrix_array = np.zeros((len(amino_acids), len(amino_acids)))
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Fill the array with substitution values
    for (aa1, aa2), value in matrix.items():
        i, j = aa_to_index[aa1], aa_to_index[aa2]
        matrix_array[i, j] = value
        # Ensure symmetry if not already present
        matrix_array[j, i] = value
    
    # Scale values to [0, 1] range using global min-max normalization
    # Global scaling preserves symmetry: score(A→G) == score(G→A)
    min_val = matrix_array.min()
    max_val = matrix_array.max()
    scaled_array = (matrix_array - min_val) / (max_val - min_val)
    
    # Convert back to dictionary
    scaled_matrix = {}
    for aa1 in amino_acids:
        for aa2 in amino_acids:
            i, j = aa_to_index[aa1], aa_to_index[aa2]
            scaled_matrix[(aa1, aa2)] = scaled_array[i, j]
    
    return scaled_matrix


def get_substitution_matrix(matrix_type: str) -> Dict[Tuple[str, str], float]:
    """Load a substitution matrix from src/substitution_matrices/, scaled to [0, 1] and cached.

    Reads the JSON file, then applies adjust_and_scale_substitution_matrix.

    Args:
        matrix_type: matrix name, e.g. 'MIYATA_EVO' or 'BLOSUM62'.

    Returns:
        Substitution matrix with values scaled to [0, 1].

    Raises:
        FileNotFoundError: if no .json exists for matrix_type.
    """
    # OPTIMIZATION: Check cache first to avoid repeated file loading
    if matrix_type in _MATRIX_CACHE:
        return _MATRIX_CACHE[matrix_type]

    current_dir = os.path.dirname(os.path.abspath(__file__))

    json_path = os.path.join(current_dir, f'substitution_matrices/{matrix_type}.json')
    if not os.path.exists(json_path):
        raise FileNotFoundError(f"No file found for {matrix_type} (.json)")

    raw_matrix = load_substitution_matrix(json_path)

    processed_matrix = adjust_and_scale_substitution_matrix(raw_matrix, matrix_type)

    # Cache the processed matrix for future use
    _MATRIX_CACHE[matrix_type] = processed_matrix

    return processed_matrix


def get_mutational_scores(changes: Union[List[str], Dict[str, Any]],
                         substitution_matrix_type: str = "MIYATA_EVO",
                         custom_matrix: Optional[Dict[Tuple[str, str], float]] = None,
                         min_score: float = 0.1,
                         indel_score: float = 1.0) -> Dict[str, float]:
    """Score amino acid changes by biochemical dissimilarity (higher = more impactful).

    Each score is 1 - scaled_similarity from the substitution matrix, floored at min_score;
    insertions and deletions (a '-' in the mutation string) receive indel_score. Lookup is
    symmetric, so (A,G) and (G,A) score the same.

    Args:
        changes: mutation strings "{orig}{pos}{mut}" (e.g. "A123G", "D100-"), list or dict keys.
        substitution_matrix_type: matrix name (default "MIYATA_EVO"); ignored if custom_matrix set.
        custom_matrix: use this matrix (tuple keys) instead of loading one from file.
        min_score: floor applied to every substitution (default 0.1).
        indel_score: score for insertions and deletions (default 1.0).

    Returns:
        {mutation: score} with score in [min_score, indel_score], rounded to 3 decimals.
    """
    if not 0.0 <= indel_score <= 1.0:
        raise ValueError(
            f"indel_score must lie in [0, 1], got {indel_score}. "
            f"An indel takes this value directly, so anything outside the interval moves the "
            f"mutational score, and KSS with it, outside the [0, 1] they are defined on."
        )

    # Handle both list and dict inputs
    if isinstance(changes, dict):
        mutation_list = list(changes.keys())
    else:
        mutation_list = changes

    # Get processed substitution matrix
    if custom_matrix is not None:
        # Use the provided custom matrix and process it
        matrix = adjust_and_scale_substitution_matrix(custom_matrix, substitution_matrix_type)
    else:
        # Load and process matrix from its JSON file
        matrix = get_substitution_matrix(substitution_matrix_type)

    # Calculate impact scores for each mutation
    scores = {}
    for mutation in mutation_list:
        # Assign indel_score for deletions/insertions
        if "-" in mutation:
            score = indel_score
        else:
            # Extract original and mutated amino acids
            original_aa = mutation[0]
            mutated_aa = mutation[-1]
            
            # Get substitution score (symmetric lookup)
            pair = (original_aa, mutated_aa)
            reverse_pair = (mutated_aa, original_aa)
            
            substitution_score = matrix.get(pair, matrix.get(reverse_pair, 0))
            
            # Calculate mutation score (higher = more significant change)
            raw_score = 1.0 - substitution_score
            
            # Apply minimum score to avoid zero impact
            score = max(raw_score, min_score)

        scores[mutation] = round(score, 3)
    
    return scores