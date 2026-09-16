#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
K-mer Significance Score (KSS) Module

This module provides the core functionality for computing k-mer significance scores
in genomic sequences. It combines mutational, discriminative, and protein-level scores
to identify significant genetic variations that may correlate with phenotypic classes.

The KSS scoring system evaluates:
- Mutational impact based on amino acid substitution matrices
- Discriminative power for class separation
- Protein characterization depth from UniProt annotations

Main functions:
    compile_results: Aggregate sequence variations by position/k-mer
    compute_kss_scores: Calculate comprehensive KSS scores
    get_taxon_id: Extract taxonomic information from GenBank files
"""

import os
from typing import Dict, Any, List, Tuple, Optional, Union
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import math

# Import internal modules (relative imports for package, absolute for standalone)
try:
    from . import protein_score
    from . import mutation_score
    from . import discriminative_score
except ImportError:
    # Fallback for running as standalone script
    import protein_score
    import mutation_score
    import discriminative_score


# Cache for GenBank taxon IDs to avoid repeated file I/O
_TAXON_ID_CACHE = {}


def compile_results(results: Dict[str, Any],
                   parameters: Dict[str, Any],
                   target_gene: Optional[str] = None,
                   verbose: bool = False) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """Aggregate the per-sequence mutations by position and k-mer, then filter by prevalence.

    Groups the raw mutations by genomic position and variant k-mer, counts class occurrences,
    and keeps only variants (k-mers of ACGT and gap characters) reaching the prevalence
    threshold in at least one class. The reference k-mer passes the same filter as the others
    and is dropped when it does not reach the threshold; amino acid changes start at score 0
    (filled in later by compute_kss_scores).

    Args:
        results: the incremental structure returned by kanalyzer.analyze_records.
        parameters: needs 'threshold' (minimum prevalence in percent, 25 == t = 0.25) and 'k'.
        target_gene: if given, compile only this gene; otherwise all genes.
        verbose: print progress.

    Returns:
        {gene: {(position, ref_kmer): {'variations': {alt_kmer:
        {'amino_acid_changes': {mutation: 0}, 'class_counts': {class: count}}}}}}.
    """
    compiled_results = {}
    threshold = parameters["threshold"]

    def is_valid_nucleotide_sequence(sequence: str) -> bool:
        """Check if a sequence contains only valid ACGT nucleotides or gaps."""
        return all(nucleotide in 'ACGT-' for nucleotide in sequence.upper())

    # Extract genes from new structure
    genes_data = results.get("genes", {})

    # Filter to target gene if specified (OPTIMIZATION: process only 1 gene instead of all)
    if target_gene:
        if target_gene not in genes_data:
            return {}
        genes_to_process = {target_gene: genes_data[target_gene]}
    else:
        genes_to_process = genes_data

    for gene, gene_info in genes_to_process.items():
        compiled_results[gene] = {}
        sequences = gene_info.get("sequences", {})
        num_sequences = len(sequences)

        if verbose:
            print(f"  Compiling {gene}: {num_sequences} sequences...", end='', flush=True)

        # Count classes
        class_counts = {}
        for seq_id, seq_data in sequences.items():
            class_id = seq_data.get("class", "unknown")
            class_counts[class_id] = class_counts.get(class_id, 0) + 1

        # Build position-based structure: (position, ref_kmer) -> variations
        position_mutations = {}

        # Track all positions that exist (from gene metadata or mutations)
        all_positions = {}

        # First pass: collect all mutations
        for seq_id, seq_data in sequences.items():
            class_id = seq_data.get("class", "unknown")
            mutations = seq_data.get("mutations", [])

            for mutation in mutations:
                position = mutation["position"]
                ref_kmer = mutation["ref_kmer"]
                alt_kmer = mutation["alt_kmer"]
                aa_changes = mutation.get("aa_changes", [])

                # Validate nucleotide sequence
                if not is_valid_nucleotide_sequence(alt_kmer):
                    continue

                key = (position, ref_kmer)

                # Track position
                if key not in all_positions:
                    all_positions[key] = ref_kmer

                # Initialize position if not seen
                if key not in position_mutations:
                    position_mutations[key] = {}

                # Initialize alt_kmer if not seen
                if alt_kmer not in position_mutations[key]:
                    position_mutations[key][alt_kmer] = {
                        "class_counts": {cls: 0 for cls in class_counts},
                        "amino_acid_changes": {}
                    }

                # Increment class count
                position_mutations[key][alt_kmer]["class_counts"][class_id] += 1

                # Collect amino acid changes (initialize with score 0)
                for aa_change in aa_changes:
                    notation = aa_change["notation"]
                    if notation not in position_mutations[key][alt_kmer]["amino_acid_changes"]:
                        position_mutations[key][alt_kmer]["amino_acid_changes"][notation] = 0

        # Second pass: count reference k-mers (sequences without mutation at this position)
        # OPTIMIZATION: Build position→mutated_seqs mapping once (instead of per-position lookup)
        position_to_mutated_seqs = defaultdict(set)
        for seq_id, seq_data in sequences.items():
            mutations = seq_data.get("mutations", [])
            for mutation in mutations:
                position_to_mutated_seqs[mutation["position"]].add(seq_id)

        # Pre-compute sequence→class mapping for O(1) class lookup
        seq_to_class = {seq_id: seq_data.get("class", "unknown")
                       for seq_id, seq_data in sequences.items()}
        all_seq_ids = set(sequences.keys())

        # Now count reference k-mers efficiently using pre-computed mapping (O(1) lookup per position)
        for key, ref_kmer in all_positions.items():
            position, _ = key

            # Initialize reference if not present
            if key not in position_mutations:
                position_mutations[key] = {}

            if ref_kmer not in position_mutations[key]:
                position_mutations[key][ref_kmer] = {
                    "class_counts": {cls: 0 for cls in class_counts},
                    "amino_acid_changes": {}
                }

            # Get sequences that have a mutation at this position (O(1) lookup instead of O(n×m))
            mutated_at_position = position_to_mutated_seqs[position]

            # Sequences with reference k-mer = all sequences - mutated sequences (O(n) set operation)
            unmutated_seq_ids = all_seq_ids - mutated_at_position

            # Count by class using pre-computed lookup
            for seq_id in unmutated_seq_ids:
                class_id = seq_to_class[seq_id]
                position_mutations[key][ref_kmer]["class_counts"][class_id] += 1

        if verbose:
            print(" Done")

        # Build final compiled structure.
        #
        # The aggregation above is keyed by (position, reference k-mer) but the output below is
        # keyed by the position alone, so two reference k-mers at one position would silently
        # overwrite one another, and which one survived would depend on dictionary order. That
        # is how a defect stayed invisible: no count changed, no error was raised, one group of
        # observed variants simply never reached the feature matrix. Upstream, kanalyzer no
        # longer emits the truncated terminal windows that produced those collisions, so this
        # should now be unreachable. It is checked rather than assumed.
        seen_references = {}
        for position, ref_kmer in position_mutations:
            if seen_references.setdefault(position, ref_kmer) != ref_kmer:
                raise RuntimeError(
                    f"{gene} position {position} carries two reference k-mers, "
                    f"{seen_references[position]!r} and {ref_kmer!r}. Writing them under the "
                    f"position alone would discard one of them. Fix the window they come from "
                    f"rather than letting this pass.")

        for key, variations in position_mutations.items():
            position, ref_kmer = key

            # Apply threshold filter (minimum percentage in at least one class)
            if threshold > 0:
                filtered_vars = {}
                for var_kmer, details in variations.items():
                    # Check if variation meets threshold in at least one class
                    meets_threshold = False
                    for cls, count in details["class_counts"].items():
                        if class_counts.get(cls, 0) > 0:
                            percentage = (count / class_counts[cls]) * 100
                            if percentage >= threshold:
                                meets_threshold = True
                                break

                    if meets_threshold:
                        filtered_vars[var_kmer] = details

                # Only include position if it has variations passing threshold
                if filtered_vars:
                    compiled_results[gene][str(position)] = {
                        "ref": ref_kmer,
                        "alts": filtered_vars
                    }
            else:
                compiled_results[gene][str(position)] = {
                    "ref": ref_kmer,
                    "alts": variations
                }

    return compiled_results


def get_taxon_id(gb_path: str) -> Optional[str]:
    """
    Extract NCBI taxonomy ID from a GenBank file.

    Args:
        gb_path: Path to the GenBank (.gb) file

    Returns:
        NCBI taxonomy ID as a string, or None if not found

    Notes:
        - Searches for 'db_xref' qualifier in source features
        - Returns the first taxonomy ID found
        - Results are cached to avoid repeated file I/O
    """
    # OPTIMIZATION: Check cache first to avoid repeated GenBank parsing
    if gb_path in _TAXON_ID_CACHE:
        return _TAXON_ID_CACHE[gb_path]

    for record in SeqIO.parse(gb_path, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                db_xrefs = feature.qualifiers.get("db_xref", [])
                for db_xref in db_xrefs:
                    parts = db_xref.split(":")
                    if len(parts) >= 2 and parts[0] == "taxon":
                        taxon_id = parts[1]
                        _TAXON_ID_CACHE[gb_path] = taxon_id
                        return taxon_id

    _TAXON_ID_CACHE[gb_path] = None
    return None


def calculate_shannon_entropy_kmers(position_details: Dict[str, Any]) -> float:
    """
    Calculate Shannon entropy from k-mer distribution at a genomic position.

    **NOTE:** This function is currently not used in the main pipeline but is kept
    for potential future analysis or notebook use.

    Shannon entropy measures the variability/diversity of k-mers at a position:
    - H = 0: Only one k-mer variant (conserved position)
    - H high: Many k-mer variants with similar frequencies (variable position)

    This is distinct from KSS discriminative score:
    - Shannon measures variability (how many variants exist)
    - KSS measures discrimination (are variants class-specific)

    Args:
        position_details: Dictionary containing 'alts' with k-mer class_counts
                         from compiled_results[gene][position]

    Returns:
        Shannon entropy in bits (log2)

    Notes:
        - Uses k-mer frequency distribution from compiled_results
        - Calculated per position (same granularity as KSS)
        - Reference k-mer included if present in alts
    """
    alts = position_details.get('alts', {})

    if not alts:
        return 0.0  # No k-mer variants = no entropy (conserved)

    # Calculate total count across all k-mer variants
    total_count = 0
    for variant_info in alts.values():
        class_counts = variant_info.get('class_counts', {})
        total_count += sum(class_counts.values())

    if total_count == 0:
        return 0.0  # No observations

    # Calculate frequency for each k-mer variant
    frequencies = []
    for variant_info in alts.values():
        class_counts = variant_info.get('class_counts', {})
        variant_count = sum(class_counts.values())
        frequency = variant_count / total_count
        frequencies.append(frequency)

    # Calculate Shannon entropy: H = -Σ p(i) * log₂(p(i))
    entropy = 0.0
    for freq in frequencies:
        if freq > 0:  # log(0) undefined
            entropy -= freq * math.log2(freq)

    return round(entropy, 4)


def _get_all_positions(gene_metadata: Dict[str, Any], k: int) -> List[int]:
    """
    Generate all possible k-mer positions for a gene.

    Args:
        gene_metadata: Gene metadata containing reference sequence info
        k: K-mer size

    Returns:
        List of all positions (1-indexed, step=k)
    """
    # Get gene length from metadata
    aa_length = gene_metadata.get("reference_aa_length", 0)
    coding_length = aa_length * 3 if aa_length else gene_metadata.get("reference_length", 0)
    if coding_length == 0:
        return []

    # Generate positions: 1, 1+k, 1+2k, ..., until position+k-1 <= coding_length
    positions = []
    pos = 1
    while pos + k - 1 <= coding_length:
        positions.append(pos)
        pos += k

    return positions


def build_feature_matrix(gene_sequences: Dict[str, Any],
                        compiled_gene_data: Dict[str, Any]) -> Tuple[np.ndarray, Dict[str, Dict], List[str], List[str]]:
    """
    Build feature matrix for discriminative score calculation.

    This function constructs a binary presence matrix where each row represents
    a sequence and each column represents a k-mer variant at a specific position.
    Only k-mers that passed the threshold filter (in compiled_gene_data) are included.

    Args:
        gene_sequences: Dictionary of sequences with mutations and classes
            Format: {seq_id: {"class": str, "mutations": [...]}}
        compiled_gene_data: Compiled results for a single gene (after threshold filtering)
            Format: {position_str: {"ref": str, "alts": {...}}}

    Returns:
        Tuple of (X_gene_all, kmers_variants_X, y, seq_ids_ordered):
            - X_gene_all: Binary matrix (n_sequences, n_total_kmers)
            - kmers_variants_X: Metadata dict {position_str: {"var_keys": [...], "col_indices": [...]}}
            - y: List of class labels (aligned with X_gene_all rows)
            - seq_ids_ordered: List of sequence IDs (aligned with X_gene_all rows)

    Notes:
        - Sequences with filtered k-mers (not in var_keys) will have all-zero rows
        - Matrix is built by horizontal stacking of position-specific matrices
        - Row order is critical for y alignment
    """
    # Extract class labels (CRITICAL: maintain order for X_local alignment)
    seq_ids_ordered = []
    y = []
    for seq_id, seq_data in gene_sequences.items():
        seq_ids_ordered.append(seq_id)
        y.append(seq_data.get("class", "unknown"))

    # Pre-compute sequence mutation index (OPTIMIZATION)
    seq_mutation_index = {}
    for seq_id in seq_ids_ordered:
        seq_data = gene_sequences[seq_id]
        mutation_dict = {}
        for mutation in seq_data.get("mutations", []):
            mutation_dict[mutation["position"]] = mutation["alt_kmer"]
        seq_mutation_index[seq_id] = mutation_dict

    # Build feature matrix with threshold-filtered k-mers
    X_gene_all = []
    kmers_variants_X = {}
    col_idx = 0

    for position_str, details in compiled_gene_data.items():
        position = int(position_str)
        ref_kmer = details['ref']
        var_keys = sorted(list(details.get('alts', {}).keys()))

        if not var_keys:
            continue

        X_local = []

        # Build binary presence matrix (CRITICAL: use seq_ids_ordered for y alignment)
        for seq_id in seq_ids_ordered:
            alt_at_pos = seq_mutation_index[seq_id].get(position)
            var_presence = []

            for variation in var_keys:
                if alt_at_pos is None:
                    # No mutation → reference k-mer present
                    is_present = 1 if variation == ref_kmer else 0
                else:
                    # Mutation exists → check if matches this variation
                    is_present = 1 if alt_at_pos == variation else 0

                var_presence.append(is_present)

            X_local.append(var_presence)

        if not X_local:
            continue

        X_local_np = np.array(X_local)

        if len(X_gene_all) == 0:
            X_gene_all = X_local_np
        else:
            X_gene_all = np.hstack([X_gene_all, X_local_np])

        kmers_variants_X[position_str] = {
            'var_keys': var_keys,
            'col_indices': list(range(col_idx, col_idx + len(var_keys)))
        }
        col_idx += len(var_keys)

    return X_gene_all, kmers_variants_X, y, seq_ids_ordered


def _compute_position_scores(details: Dict[str, Any],
                            current_protein_score: float,
                            parameters: Dict[str, Any],
                            weight_sum: float) -> Dict[str, float]:
    """
    Compute mutational and protein scores for a single position.

    Helper function to avoid code duplication in compute_kss_scores.

    Args:
        details: Position details containing 'alts' with amino acid changes
        current_protein_score: Protein characterization score for this gene
        parameters: Configuration dictionary with mutational_matrix and
                   optionally indel_score (default: 1.0)
        weight_sum: Pre-computed sum of score weights

    Returns:
        Dictionary with 'protein_score' and 'mutational_score'
    """
    details['protein_score'] = current_protein_score

    # Collect all amino acid changes across all variations
    all_aa_changes = {}
    for variation, var_details in details["alts"].items():
        for mutation in var_details["amino_acid_changes"]:
            if mutation not in all_aa_changes:
                all_aa_changes[mutation] = 0

    # Calculate mutational scores
    indel_score = parameters.get("indel_score", 1.0)
    mutation_scores = mutation_score.get_mutational_scores(
        all_aa_changes, parameters["mutational_matrix"],
        indel_score=indel_score
    )
    details['mutational_score'] = max(mutation_scores.values()) if mutation_scores else 0

    # Update individual mutation scores in variations
    for variation, var_details in details["alts"].items():
        for mutation in var_details["amino_acid_changes"]:
            if mutation in mutation_scores:
                var_details["amino_acid_changes"][mutation] = mutation_scores[mutation]

    return details


def _compute_final_kss(details: Dict[str, Any],
                      parameters: Dict[str, Any],
                      weight_sum: float) -> float:
    """
    Compute final weighted KSS score.

    Args:
        details: Position details with discriminative_score, mutational_score, protein_score
        parameters: Configuration dictionary with weight parameters
        weight_sum: Pre-computed sum of weights

    Returns:
        Final KSS score (rounded to 3 decimals)
    """
    avg_score = (
        parameters["discriminative_weight"] * details['discriminative_score'] +
        parameters["mutational_weight"] * details['mutational_score'] +
        parameters["protein_weight"] * details['protein_score']
    ) / weight_sum
    return round(avg_score, 3)


def compute_kss_scores(compiled_results: Dict[str, Dict[str, Dict[str, Any]]],
                       results: Dict[str, Any],
                       infos: Dict[str, Any],
                       parameters: Dict[str, Any],
                       target_gene: Optional[str] = None) -> Tuple[Dict, Dict]:
    """Add the three component scores and the final weighted KSS to every position.

    Each position gets a discriminative, mutational and protein score in [0, 1], then their
    weighted mean as the final KSS. Reference-only positions get discriminative_score = 0.

    Args:
        compiled_results: output of compile_results().
        results: the incremental structure (used for class labels and sequences).
        infos: needs 'input_folder' (to locate the GenBank file per gene).
        parameters: needs 'mutational_matrix', the three component weights and 'k'.
        target_gene: if given, score only this gene; otherwise all genes.

    Returns:
        (compiled_results enriched with scores, {gene: {position: discriminative_score}}).
    """
    discriminative_scores = {}

    # Extract genes data from new structure
    genes_data = results.get("genes", {})

    # Filter to target gene if specified
    if target_gene:
        if target_gene not in compiled_results:
            return {}, {}
        genes_to_process = {target_gene: compiled_results[target_gene]}
    else:
        genes_to_process = compiled_results

    # Pre-compute weight_sum once (constant across all genes/positions)
    weight_sum = (
        parameters["discriminative_weight"] +
        parameters["mutational_weight"] +
        parameters["protein_weight"]
    )

    for gene, data in genes_to_process.items():
        print(f"Processing gene: {gene}")

        total_positions = len(data)
        processed_positions = 0

        discriminative_scores[gene] = {}

        # Load GenBank file for protein scoring
        gb_dir = os.path.join(infos["input_folder"], gene)
        gb_files = [f for f in os.listdir(gb_dir) if f.lower().endswith('.gb')]
        if not gb_files:
            print(f"  Warning: No GenBank file found for {gene}")
            continue

        gb_path = os.path.join(gb_dir, gb_files[0])
        # UniProt is queried only when the component is weighted. The lookup aborts the run when
        # it cannot reach UniProt, deliberately, because scoring a protein as uncharacterized on
        # a network failure is a statement about the network. At w_p = 0 nothing is scored from
        # it, so that guard would refuse a run whose result does not depend on the answer, which
        # is the setting for a protein that has no UniProt entry to begin with.
        if parameters["protein_weight"] == 0:
            current_protein_score = 0.0
        else:
            taxon_id = get_taxon_id(gb_path)
            current_protein_score = protein_score.get_protein_score(taxon_id, gene, verbose=False)

        # Build feature matrix using refactored function
        gene_sequences = genes_data.get(gene, {}).get("sequences", {})
        X_gene_all, kmers_variants_X, y, seq_ids_ordered = build_feature_matrix(
            gene_sequences, data
        )

        # Handle edge case: no valid features
        if isinstance(X_gene_all, list) or len(X_gene_all) == 0 or X_gene_all.shape[1] == 0:
            print(f"  ⚠ Warning: Gene {gene} has 0 variant features. Using default scores.")

            for position_str, details in data.items():
                # Compute protein and mutational scores
                _compute_position_scores(details, current_protein_score, parameters, weight_sum)

                # Discriminative score (default for no features)
                details['discriminative_score'] = 0.0
                discriminative_scores[gene][position_str] = 0.0

                # Final KSS
                details['kss'] = _compute_final_kss(details, parameters, weight_sum)

            continue

        # Normal processing for positions with features
        for position_str, details in data.items():
            processed_positions += 1
            if processed_positions % 50 == 0 or processed_positions == total_positions:
                pct = (processed_positions / total_positions) * 100
                print(f"\r  Scoring positions: {processed_positions}/{total_positions} ({pct:.1f}%)", end='', flush=True)

            # Compute protein and mutational scores
            _compute_position_scores(details, current_protein_score, parameters, weight_sum)

            # Discriminative score
            var_keys = kmers_variants_X[position_str]['var_keys']
            col_indices = kmers_variants_X[position_str]['col_indices']
            X_local = X_gene_all[:, col_indices]

            try:
                disc_result = discriminative_score.get_discriminative_score(
                    X_local, y
                )
                details['discriminative_score'] = disc_result['normalized_score']
            except Exception as e:
                print(f"\n  Warning: Failed discriminative score for {gene}:{position_str}: {e}")
                details['discriminative_score'] = 0.0

            # Store discriminative score
            discriminative_scores[gene][position_str] = details['discriminative_score']

            # Final KSS score
            details['kss'] = _compute_final_kss(details, parameters, weight_sum)

        # Complete progress line if positions were processed
        if total_positions > 0:
            print()  # Newline after progress bar

        # Calculate scores for reference-only positions (invariant)
        print(f"  Calculating scores for all positions...")
        gene_metadata = genes_data.get(gene, {}).get("metadata", {})
        all_positions = _get_all_positions(gene_metadata, parameters["k"])

        positions_to_process = [pos for pos in all_positions if str(pos) not in discriminative_scores[gene]]
        total_ref_positions = len(positions_to_process)
        processed_ref_positions = 0

        for pos in all_positions:
            pos_str = str(pos)

            if pos_str in discriminative_scores[gene]:
                continue

            processed_ref_positions += 1
            if processed_ref_positions % 100 == 0 or processed_ref_positions == total_ref_positions:
                pct = (processed_ref_positions / total_ref_positions) * 100 if total_ref_positions > 0 else 0
                print(f"\r  Reference positions: {processed_ref_positions}/{total_ref_positions} ({pct:.1f}%)", end='', flush=True)

            # Invariant positions → discriminative_score = 0
            discriminative_scores[gene][pos_str] = 0.0

        # Complete progress line if reference positions were processed
        if total_ref_positions > 0:
            print()  # Newline after progress bar

    return compiled_results, discriminative_scores