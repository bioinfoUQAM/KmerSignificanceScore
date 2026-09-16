#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
K-mer Significance Score (KSS) - Score Recalculation Tool

This script recalculates KSS scores from existing raw results without re-running
sequence alignments. This is useful when you want to:
    - Test different scoring parameters (weights, thresholds, matrices)
    - Update protein annotations from UniProt
    - Regenerate compiled results after parameter changes

The script loads pre-computed mutation data and recalculates:
    - Compiled results (aggregated by position/k-mer)
    - Discriminative scores (class separation)
    - Mutational scores (amino acid substitution impact)
    - Protein scores (characterization depth)
    - Final weighted KSS scores

Usage:
    python recalculate_scores.py [config_files...] [--gene GENE]

    Examples:
        python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml
        python recalculate_scores.py data/*/config.yaml
        python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml --gene UL55

Requirements:
    - Existing *_results.json files (raw results from previous analysis)
    - Valid config.yaml configuration file
    - Same dependencies as main analysis pipeline

Performance:
    - Much faster than full analysis (no alignments)
    - Typical runtime: 1-2 minutes per gene (vs 10-15 minutes for full analysis)
"""

import os
import sys
import glob
import argparse
import io

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Import KSS modules
sys.path.insert(0, 'src')
from src import pipeline


def find_config_files(patterns: list) -> list:
    """
    Find configuration files from patterns or use defaults.

    Args:
        patterns: List of file patterns or paths

    Returns:
        List of configuration file paths

    Notes:
        - Supports glob patterns (e.g., "data/*/config.yaml")
        - If no patterns provided, searches data/*/config.yaml
        - Exits with error if no files found
    """
    config_files = []

    if patterns:
        # Use provided config files or patterns
        for pattern in patterns:
            matched = glob.glob(pattern)
            config_files.extend(matched)
    else:
        # Default: find all config.yaml in data/ subdirectories
        config_files = glob.glob("data/*/config.yaml")

    if not config_files:
        print("Error: No configuration files found!")
        print("\nUsage:")
        print("  python recalculate_scores.py [config_files...] [options]")
        print("\nExamples:")
        print("  python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml")
        print("  python recalculate_scores.py data/*/config.yaml")
        print("  python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml --gene UL55")
        sys.exit(1)

    return config_files


def main():
    """
    Main execution function for KSS score recalculation.

    Parses command line arguments and processes datasets to recalculate
    KSS scores from existing raw results.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Recalculate KSS scores from existing raw results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Recalculate all genes in dataset
  python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml

  # Recalculate specific gene
  python recalculate_scores.py data/Human_betaherpesvirus_5/config.yaml --gene UL55

  # Recalculate multiple datasets
  python recalculate_scores.py data/*/config.yaml

Notes:
  - Requires existing *_results.json files
  - Much faster than full analysis (no sequence alignments)
  - Overwrites existing compiled results
        """
    )
    parser.add_argument(
        'config_files',
        nargs='*',
        help='Configuration file(s) or glob pattern (default: data/*/config.yaml)'
    )
    parser.add_argument(
        '--gene',
        type=str,
        default=None,
        help='Specific gene to recalculate (default: all genes in config)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=str,
        default=None,
        metavar='DIR',
        help='Write results to <DIR>/<dataset_name>/ instead of overwriting them next to '
             'the input data. Use this when sweeping a parameter, so that the results '
             'behind an existing analysis are left intact.'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=None,
        metavar='PCT',
        help='Override the prevalence threshold from the configuration, in percent. '
             'Intended for sensitivity analyses; requires --output-dir so that the '
             'sweep cannot overwrite the configured results.'
    )

    args = parser.parse_args()

    if args.threshold is not None and not args.output_dir:
        parser.error("--threshold changes the results, so --output-dir is required to "
                     "keep the existing ones intact.")

    overrides = {} if args.threshold is None else {"threshold": args.threshold}

    # Header
    print("\n" + "="*70)
    print("K-mer Significance Score (KSS) - Score Recalculation")
    print("="*70)

    # Find configuration files
    config_files = find_config_files(args.config_files)

    print(f"\nFound {len(config_files)} configuration file(s):")
    for config_file in config_files:
        print(f"  • {config_file}")

    # Process each configuration file
    all_results = {}
    total_successful = 0
    total_failed = 0

    for idx, config_file in enumerate(config_files, 1):
        print(f"\n[Dataset {idx}/{len(config_files)}]")

        try:
            results = pipeline.recalculate_scores_from_config(
                config_file,
                target_gene=args.gene,
                verbose=True,
                output_root=args.output_dir,
                overrides=overrides or None
            )

            all_results[config_file] = results

            # Count successes and failures
            successful = sum(1 for success in results.values() if success)
            failed = sum(1 for success in results.values() if not success)

            total_successful += successful
            total_failed += failed

        except Exception as e:
            import traceback
            print(f"\n✗ Error processing {config_file}:")
            print(f"  {str(e)}")
            traceback.print_exc()
            total_failed += 1
            continue

    # Final summary
    print("\n" + "="*70)
    print("Overall Summary")
    print("="*70)
    print(f"  Datasets processed: {len(config_files)}")
    print(f"  Genes successful: {total_successful}")
    if total_failed > 0:
        print(f"  Genes failed: {total_failed}")
    print("="*70)

    if total_failed == 0:
        print("\n✓ All score recalculations completed successfully!\n")
    else:
        print(f"\n⚠ {total_failed} gene(s) failed. Check error messages above.\n")

    # Exit with appropriate code
    sys.exit(0 if total_failed == 0 else 1)


if __name__ == "__main__":
    main()
