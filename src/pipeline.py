#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
KSS Pipeline Module

This module provides reusable workflow functions for the K-mer Significance Score
analysis pipeline. It centralizes common operations such as configuration loading,
score recalculation, and dataset processing to promote code reuse and maintainability.

Main functions:
    find_config_files: Find configuration files from arguments or defaults
    load_config: Load and parse configuration from YAML files
    recalculate_scores_for_gene: Recalculate scores for a single gene
    recalculate_scores_from_config: Recalculate scores for entire dataset
"""

import os
import sys
import glob
import yaml
from typing import Dict, Any, Tuple, Optional, List

# Import KSS modules
try:
    from . import utils
    from . import kss
except ImportError:
    # Fallback for running as standalone script
    import utils
    import kss


def find_config_files(args: List[str], input_root: Optional[str] = None) -> List[str]:
    """
    Find configuration files from command line arguments or default locations.

    Args:
        args: Command line arguments (e.g., sys.argv[1:])
        input_root: Directory to search for <root>/*/config.yaml when no explicit
                    path or pattern is given. Defaults to "data".

    Returns:
        List of configuration file paths

    Raises:
        SystemExit: If no configuration files are found

    Notes:
        - Supports glob patterns (e.g., "data/*/config.yaml")
        - If no arguments provided, searches for data/*/config.yaml
        - Prints usage information and exits if no files found

    Examples:
        >>> files = find_config_files(['data/CMV/config.yaml'])
        >>> files = find_config_files(['data/*/config.yaml'])
        >>> files = find_config_files([])  # Uses default pattern
    """
    config_files = []

    if args:
        # Use provided config files or patterns
        for arg in args:
            matched = glob.glob(arg)
            config_files.extend(matched)
    else:
        # Default: find all config.yaml in the subdirectories of the input root
        root = input_root or "data"
        config_files = glob.glob(os.path.join(root, "*", "config.yaml"))

    if not config_files:
        root = input_root or "data"
        print(f"Error: No configuration files found in '{root}'!")
        print("\nUsage:")
        print("  python main.py [config_files...] [-i INPUT_DIR] [-o OUTPUT_DIR]")
        print("\nExamples:")
        print("  python main.py data/Human_betaherpesvirus_5/config.yaml")
        print("  python main.py data/*/config.yaml")
        print("  python main.py                       # all config.yaml under data/")
        print("  python main.py -i /path/to/data -o /path/to/results")
        sys.exit(1)

    return config_files


def validate_component_weights(weights: Dict[str, float]) -> None:
    """Refuse the weights KSS is not defined on.

    Materials and methods states that the weights are non-negative and not all zero, and that
    KSS therefore stays in [0, 1]. Nothing enforced it: three zeros divide by zero when the
    weighted average is formed, and a negative weight carries the average outside the interval
    silently. A guarantee the code does not hold to is worse than none, since it is the one a
    reader of a published score relies on.

    Args:
        weights: The discriminative, mutational and protein weights, as read from the config.

    Raises:
        ValueError: If any weight is negative, or if all three are zero.
    """
    values = {name: weights[name] for name in ("discriminative", "mutational", "protein")}
    if any(value < 0 for value in values.values()):
        raise ValueError(
            f"component weights must be non-negative, got {values}. "
            f"A negative weight moves KSS outside the [0, 1] interval it is defined on."
        )
    if sum(values.values()) == 0:
        raise ValueError(
            "at least one component weight must be non-zero; all three are zero, "
            "which leaves the weighted average undefined."
        )


def validate_indel_score(indel_score: float) -> None:
    """Refuse an indel score that would carry the mutational component outside [0, 1].

    The mutational score of a window is the largest impact it carries, and an insertion or a
    deletion takes this configured value directly. Anything outside [0, 1] therefore leaves the
    bound the manuscript states for every component, and does so without any error: 2.0 gives a
    KSS of 1.333 at equal weights, and -1.0 gives -1.0.

    Args:
        indel_score: The score insertions and deletions receive, as read from the config.

    Raises:
        ValueError: If the score lies outside [0, 1].
    """
    if not 0.0 <= indel_score <= 1.0:
        raise ValueError(
            f"indel_score must lie in [0, 1], got {indel_score}. "
            f"The mutational score is the largest impact at a window, so a value outside the "
            f"interval moves KSS outside the [0, 1] it is defined on."
        )


def load_config(config_path: str,
                output_root: Optional[str] = None) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """
    Load configuration from YAML file.

    When output_root is given, results are written under
    <output_root>/<dataset_name>/ instead of next to the input data, which keeps the
    input directory read-only.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Tuple of (dataset_name, dataset_info, parameters)
            - dataset_name: Name of the dataset
            - dataset_info: Dictionary with input_folder, output_folder, cds_selection
            - parameters: Dictionary with all analysis parameters

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is malformed
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Extract dataset info
    dataset = config['dataset']
    dataset_name = dataset['name']

    # Build dataset info dictionary
    base_folder = os.path.dirname(config_path)

    # Results are written next to the input data unless an output root is supplied, in which
    # case each dataset gets its own subdirectory there and the input tree is left untouched.
    if output_root:
        output_folder = os.path.join(output_root, dataset['name'])
        os.makedirs(output_folder, exist_ok=True)
    else:
        output_folder = base_folder

    dataset_info = {
        "input_folder": base_folder,
        "output_folder": output_folder,
        "cds_selection": ",".join(dataset['genes'])
    }

    # Extract parameters
    params = config['parameters']

    # The k-mer window must cover whole codons: the nucleotide window [3j, 3j+k) is mapped to
    # the amino acids [j, j + k/3) by integer division. A k that is not a multiple of 3 would
    # silently truncate that mapping and misalign the mutational component, so it is rejected
    # rather than tolerated.
    if params['k'] % 3 != 0:
        raise ValueError(
            f"k must be a multiple of 3 so that each k-mer spans whole codons, got k={params['k']}. "
            f"The mutational component maps a k-nucleotide window to k/3 amino acids."
        )

    validate_component_weights(params['weights'])
    validate_indel_score(params['scoring'].get('indel_score', 1.0))

    parameters = {
        "k": params['k'],
        "threshold": params['threshold'],
        "discriminative_weight": params['weights']['discriminative'],
        "mutational_weight": params['weights']['mutational'],
        "protein_weight": params['weights']['protein'],
        "mutational_matrix": params['scoring']['mutational_matrix'],
        "substitution_matrix": params['scoring']['substitution_matrix'],
        "indel_score": params['scoring'].get('indel_score', 1.0),
        "open_gap_score": params['alignment']['open_gap_score'],
        "extend_gap_score": params['alignment']['extend_gap_score'],
        "save_raw": params['output'].get('save_raw', False)
    }

    return dataset_name, dataset_info, parameters


def recalculate_scores_for_gene(gene: str,
                                base_folder: str,
                                parameters: Dict[str, Any],
                                verbose: bool = True,
                                output_folder: Optional[str] = None) -> bool:
    """
    Recalculate KSS scores for a single gene from existing raw results.

    This function loads pre-computed raw results (mutations and alignments)
    and recalculates only the compiled results and KSS scores. It does not
    re-run sequence alignments or mutation detection.

    Args:
        gene: Gene name to process
        base_folder: Folder the raw results are read from
        parameters: Analysis parameters dictionary
        verbose: If True, prints progress information
        output_folder: Where to write the recalculated results. Defaults to base_folder,
            which overwrites in place. Pass a separate directory to sweep a parameter
            without touching the results that produced the published analysis.

    Returns:
        True if successful, False if raw results not found or error occurred

    Notes:
        - Requires existing {gene}_results.json file
        - Overwrites {gene}_compiled_results.json in output_folder
        - Uses same parameters as original analysis
        - Much faster than full analysis (no alignments)
    """
    if output_folder is None:
        output_folder = base_folder
    if verbose:
        print(f"\nProcessing gene: {gene}")

    # Load existing raw results
    results_path = os.path.join(base_folder, gene, "results", f"{gene}_results.json")

    if not os.path.exists(results_path):
        if verbose:
            print(f"  ✗ Raw results not found: {results_path}")
        return False

    if verbose:
        print(f"  Loading raw results: {results_path}")

    try:
        results = utils.load_data_from_json(results_path)
    except Exception as e:
        if verbose:
            print(f"  ✗ Error loading raw results: {e}")
        return False

    # Raw results and reference records are read from the input tree; everything produced
    # is written to the output tree, which may be the same directory.
    infos = {
        "input_folder": base_folder,
        "output_folder": output_folder
    }

    try:
        # Step 1: Compile results
        if verbose:
            print(f"  [1/2] Compiling results...")
        compiled_results = kss.compile_results(results, parameters, target_gene=gene)

        # Step 2: Compute scores
        if verbose:
            print(f"  [2/2] Computing KSS scores...")
        compiled_results, discriminative_scores = kss.compute_kss_scores(
            compiled_results,
            results,
            infos,
            parameters,
            target_gene=gene
        )

        # An empty compilation used to be written as {} and reported as a success, the defect
        # main.py carried on the analysis side: a file that looks like a result but holds none.
        # Nothing is written and the gene is counted as failed.
        if not compiled_results.get(gene):
            if verbose:
                print(f"  ✗ {gene}: the raw results carry no scored position, nothing written")
            return False

        # Save compiled results
        gene_output_dir = os.path.join(output_folder, gene, "results")
        os.makedirs(gene_output_dir, exist_ok=True)

        compiled_path = os.path.join(gene_output_dir, f"{gene}_compiled_results.json")

        utils.save_data_as_json({gene: compiled_results[gene]}, compiled_path)

        if verbose:
            print(f"  ✓ Saved: {compiled_path}")

        return True

    except Exception as e:
        if verbose:
            import traceback
            print(f"  ✗ Error recalculating scores: {e}")
            traceback.print_exc()
        return False


def recalculate_scores_from_config(config_path: str,
                                   target_gene: Optional[str] = None,
                                   verbose: bool = True,
                                   output_root: Optional[str] = None,
                                   overrides: Optional[Dict[str, Any]] = None) -> Dict[str, bool]:
    """Recalculate scores from stored raw results, without realigning, for one config file.

    Replays scoring for genes that already have raw results, so weights or thresholds can be
    changed cheaply. Only genes with existing raw results are processed.

    Args:
        config_path: path to the YAML configuration.
        target_gene: if given, process only this gene; otherwise all genes.
        verbose: print progress.
        output_root: write under <output_root>/<dataset_name>/ instead of next to the input.
            Without it the existing compiled results are overwritten (destructive in a sweep).
        overrides: parameter values applied after loading, e.g. {"threshold": 10}; recorded in
            the output so a sweep directory is self-describing.

    Returns:
        {gene: success (True/False)}.
    """
    # Load configuration
    if verbose:
        print(f"\nLoading configuration: {config_path}")

    try:
        dataset_name, dataset_info, parameters = load_config(config_path,
                                                             output_root=output_root)
    except Exception as e:
        if verbose:
            print(f"✗ Error loading configuration: {e}")
        return {}

    if overrides:
        parameters.update(overrides)
        if verbose:
            shown = ", ".join(f"{k}={v}" for k, v in sorted(overrides.items()))
            print(f"Parameter overrides: {shown}")

    base_folder = dataset_info['input_folder']
    output_folder = dataset_info['output_folder']

    # Determine genes to process
    all_genes = dataset_info['cds_selection'].split(',')
    genes_to_process = [target_gene] if target_gene else all_genes

    if verbose:
        print(f"\n{'='*70}")
        print(f"Recalculating KSS scores: {dataset_name}")
        print(f"{'='*70}")
        print(f"Genes to process: {', '.join(genes_to_process)}")
        print(f"{'='*70}")

    # Recalculate scores for each gene
    results = {}
    successful = 0
    failed = 0

    for gene in genes_to_process:
        success = recalculate_scores_for_gene(gene, base_folder, parameters,
                                              verbose=verbose,
                                              output_folder=output_folder)
        results[gene] = success

        if success:
            successful += 1
        else:
            failed += 1

    # Summary
    if verbose:
        print(f"\n{'='*70}")
        print(f"Summary")
        print(f"{'='*70}")
        print(f"  Successful: {successful}/{len(genes_to_process)}")
        if failed > 0:
            print(f"  Failed: {failed}/{len(genes_to_process)}")
        print(f"{'='*70}\n")

    return results
