import os
import json
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def load_substitution_matrix(filename):
    """
    Load the substitution matrix from a JSON file and return it as a dictionary
    with tuple keys and float values.
    """
    with open(filename, 'r') as f:
        matrix = json.load(f)
    new_matrix = {}
    for key, value in matrix.items():
        # Adjust this parsing based on the actual format of your keys
        aa_pair = tuple(key.strip().split(','))
        aa_pair = tuple(aa.strip("'\" ()") for aa in aa_pair)
        new_matrix[aa_pair] = float(value)
    return new_matrix

def adjust_and_scale_substitution_matrix(matrix, matrix_type):
    """
    Adjust and scale the substitution matrix values based on the matrix type.
    """
    # Determine if the matrix is distance-based and needs adjustment
    if matrix_type.upper() in {"MIYATA", "GRANTHAM", "SNEATH"}:
        # Invert the values for distance-based matrices to convert them into similarity scores
        max_value = max(matrix.values())
        adjusted_matrix = {key: max_value - value for key, value in matrix.items()}
    else:
        adjusted_matrix = matrix.copy()

    # Get the set of amino acids involved in the matrix
    amino_acids = sorted(set(a for pair in adjusted_matrix for a in pair))

    # Construct a 2D array of substitution values
    matrix_array = np.array([
        [adjusted_matrix.get((aa1, aa2), adjusted_matrix.get((aa2, aa1), 0))
         for aa2 in amino_acids] for aa1 in amino_acids
    ])

    # Scale the matrix values between 0 and 1 using MinMaxScaler
    scaler = MinMaxScaler()
    scaled_values = scaler.fit_transform(matrix_array)

    # Rebuild the dictionary with scaled values
    scaled_matrix = {
        (aa1, aa2): scaled_values[i, j]
        for i, aa1 in enumerate(amino_acids)
        for j, aa2 in enumerate(amino_acids)
    }

    return scaled_matrix

def get_mutational_scores(changes, substitution_matrix_type="MIYATA", categorize=True):
    """
    Compute mutational scores for a list of amino acid changes based on a substitution matrix.
    """
    # Define the bin edges for categorization
    bins = [0.2, 0.4, 0.6, 0.8]

    # Load and adjust the substitution matrix
    current_dir = os.path.dirname(os.path.abspath(__file__))
    substitution_matrix_path = os.path.join(
        current_dir, f'substitution_matrices/{substitution_matrix_type}.json'
    )

    # Load the substitution matrix
    substitution_matrix = load_substitution_matrix(substitution_matrix_path)

    # Adjust and scale the substitution matrix
    scaled_matrix = adjust_and_scale_substitution_matrix(
        substitution_matrix, substitution_matrix_type
    )

    scores = {}
    for entry in changes:
        # Handle deletions represented by a '-' in the entry
        if "-" in entry:
            score = 1.0  # Assign the maximum score for deletions
        else:
            # Extract the original and mutated amino acids from the entry
            original_aa = entry[0]
            mutated_aa = entry[-1]
            # Retrieve the scaled substitution value
            substitution_score = scaled_matrix.get(
                (original_aa, mutated_aa),
                scaled_matrix.get((mutated_aa, original_aa), 0)
            )
            # Compute the mutational score
            score = 1.0 - substitution_score

        if categorize:
            # Categorize the score into bins and assign a category from 1 to 5
            category = np.digitize(score, bins) + 1
            scores[entry] = category
        else:
            scores[entry] = score

    return scores
