#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import os
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
from Bio.Align import PairwiseAligner, substitution_matrices
from . import mutation_score

# Parallel alignment is an optimisation, not a requirement: joblib belongs to the analysis
# extra, and the scoring pipeline must keep running on a core install, which declares five
# runtime dependencies and not six. Without joblib the sequences are aligned one after the
# other, and the results are the same either way.
try:
    from joblib import Parallel, delayed
except ImportError:  # pragma: no cover - exercised only on a core install
    Parallel = delayed = None

# Suppress BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning)


def initialize_aligner(parameters: dict) -> PairwiseAligner:
    """
    Initialize a pairwise sequence aligner with specified parameters.

    Args:
        parameters: Dictionary containing alignment parameters:
            - substitution_matrix: Name of substitution matrix (e.g., 'BLOSUM62')
            - open_gap_score: Penalty for opening a gap
            - extend_gap_score: Penalty for extending a gap

    Returns:
        Configured PairwiseAligner object for global sequence alignment
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(parameters["substitution_matrix"])
    aligner.open_gap_score = parameters["open_gap_score"]
    aligner.extend_gap_score = parameters["extend_gap_score"]

    return aligner


def align_nucleotides(protein_a: str, protein_b: str, nucleotide_a: str, nucleotide_b: str) -> tuple:
    """Place nucleotide sequences on the coordinate system of their protein alignment.

    Each codon is copied verbatim from the query sequence and each amino acid gap becomes a
    three-nucleotide gap, preserving the reading frame. The genetic code is never inverted,
    so codon degeneracy introduces no ambiguity: the codons carried into the alignment are
    those actually observed, not codons inferred from an amino acid.

    Args:
        protein_a, protein_b: Aligned protein sequences, gaps written '-'
        nucleotide_a, nucleotide_b: The corresponding unaligned nucleotide sequences

    Returns:
        The two aligned nucleotide sequences.
    """
    nucleotide_a = Seq(nucleotide_a)
    nucleotide_b = Seq(nucleotide_b)

    aligned_nuc_a = []
    aligned_nuc_b = []
    index_a = index_b = 0

    for aa_a, aa_b in zip(protein_a, protein_b):
        if aa_a == '-':
            aligned_nuc_a.append('---')
        else:
            aligned_nuc_a.append(str(nucleotide_a[index_a:index_a+3]))
            index_a += 3

        if aa_b == '-':
            aligned_nuc_b.append('---')
        else:
            aligned_nuc_b.append(str(nucleotide_b[index_b:index_b+3]))
            index_b += 3

    return ''.join(aligned_nuc_a), ''.join(aligned_nuc_b)


def replace_first_three_gaps(sequence: str, replacement: str = '') -> str:
    """
    Replace the first three gap characters '-' in a sequence with a replacement string.

    Used for handling insertions in nucleotide alignment (one codon = 3 nucleotides).

    Args:
        sequence: Input sequence containing gap characters
        replacement: String to replace gaps with (default: empty string)

    Returns:
        Modified sequence with first three gaps replaced
    """
    count = 0
    modified_sequence = []
    for char in sequence:
        if char == '-' and count < 3:
            modified_sequence.append(replacement)
            count += 1
        else:
            modified_sequence.append(char)
    return ''.join(modified_sequence)


def identify_mutations(infos: dict, parameters: dict, ref_nuc_seq: str,
                      ref_aa_seq: str, query_nuc_seq: str, query_aa_seq: str, aligner=None) -> list:
    """
    Identify nucleotide and amino acid mutations between reference and query sequences.

    This function performs k-mer based mutation analysis by:
    1. Aligning protein sequences
    2. Placing the observed codons on that alignment, never inverting the genetic code
    3. Extracting k-mers and identifying variations
    4. Annotating amino acid changes

    Only positions with actual mutations are returned (incremental storage).

    Args:
        infos: Dataset information dictionary
        parameters: Analysis parameters including k-mer size ('k')
        ref_nuc_seq: Reference nucleotide sequence
        ref_aa_seq: Reference amino acid sequence
        query_nuc_seq: Query nucleotide sequence to compare
        query_aa_seq: Query amino acid sequence to compare

    Returns:
        List of mutation objects, each containing:
        - position: Nucleotide position (1-indexed)
        - ref_kmer: Reference k-mer sequence
        - alt_kmer: Alternative k-mer sequence
        - aa_changes: List of amino acid changes with impact scores
    """
    mutations = []
    skipped = []
    k = parameters["k"]

    # Perform sequence alignments
    if aligner is None:
        aligner = initialize_aligner(parameters)
    alignments = aligner.align(ref_aa_seq, query_aa_seq)
    aligned_ref_aa, aligned_query_aa = alignments[0]
    aligned_ref_nuc, aligned_query_nuc = align_nucleotides(aligned_ref_aa, aligned_query_aa, ref_nuc_seq, query_nuc_seq)

    n_insertion = 0  # Track insertions

    # Analyze mutations - only store if different from reference
    for i in range(0, len(aligned_ref_nuc) - k + 1, k):
        adjusted_k = k
        current_pos = i + (n_insertion * 3)

        # Get current k-mers
        ref_kmer = aligned_ref_nuc[current_pos:current_pos + adjusted_k]
        query_kmer = aligned_query_nuc[current_pos:current_pos + adjusted_k]

        # Handle gaps in reference sequence
        while '-' in ref_kmer and current_pos + adjusted_k < len(aligned_ref_nuc):
            next_block = current_pos + adjusted_k
            ref_kmer = replace_first_three_gaps(ref_kmer) + aligned_ref_nuc[next_block:next_block + 3]
            adjusted_k += 3
            query_kmer = aligned_query_nuc[current_pos:current_pos + adjusted_k]

        # Update insertion tracking
        current_insertion = (adjusted_k - k) // 3
        n_insertion += current_insertion

        # A window is a k-mer position only if the reference actually supplies k nucleotides
        # for it. At the end of a gene, once the accumulated insertions have pushed the window
        # past the last complete codon, the extension loop above runs out of sequence and hands
        # back a truncated reference: fewer than k non-gap nucleotides, padded with gaps. Such a
        # window is not a position, and what would be recorded against it is the stop codon and
        # the tail of the alignment. It is therefore not emitted.
        #
        # The insertion accounting above is deliberately left untouched, so that skipping this
        # window changes nothing for any window that follows it.
        #
        # It is reported rather than dropped in silence. A silent skip here would be the same
        # fault this fix exists to remove: something observed disappearing without trace.
        reference_length = len(ref_kmer.replace('-', ''))
        if reference_length != k:
            skipped.append({"position": current_pos + 1,
                            "reference_nucleotides": reference_length})
            continue

        # Skip if k-mers are identical (incremental storage)
        if ref_kmer == query_kmer:
            continue

        # Calculate amino acid positions
        aa_start = current_pos // 3
        aa_end = min(aa_start + (adjusted_k // 3), len(aligned_ref_aa))

        ref_aa_kmer = aligned_ref_aa[aa_start:aa_end]
        query_aa_kmer = aligned_query_aa[aa_start:aa_end]

        # Process amino acid changes
        aa_changes = []
        aa_position = aa_start - n_insertion
        processed_insertions = 0

        for idx, (ref_aa, query_aa) in enumerate(zip(ref_aa_kmer, query_aa_kmer)):
            if ref_aa != query_aa:
                # Calculate mutation position
                mut_pos = max(0, aa_position + idx + 1 + current_insertion - processed_insertions)

                # Create mutation annotation
                notation = f"{ref_aa}{mut_pos}{query_aa}"
                if ref_aa == '-':
                    processed_insertions += 1

                # Handle duplicate positions with suffix
                original_notation = notation
                suffix = ""
                while any(change["notation"] == notation for change in aa_changes):
                    suffix += "*"
                    notation = original_notation + suffix

                aa_changes.append({
                    "notation": notation,
                    "mut_score": 0  # Placeholder, will be updated below
                })

        # Store mutation only if there are changes
        if aa_changes or ref_kmer != query_kmer:
            mutations.append({
                "position": i + 1,
                "ref_kmer": ref_kmer,
                "alt_kmer": query_kmer,
                "aa_changes": aa_changes
            })

    # Calculate mutation scores for all amino acid changes using substitution matrix
    if mutations:
        # Collect all unique mutation notations
        all_notations = {}
        for mutation in mutations:
            for aa_change in mutation["aa_changes"]:
                notation = aa_change["notation"]
                # Remove suffix for score calculation (e.g., "A10G*" -> "A10G")
                clean_notation = notation.rstrip("*")
                if clean_notation not in all_notations:
                    all_notations[clean_notation] = 0

        # Get mutation scores from substitution matrix
        mutation_scores = mutation_score.get_mutational_scores(
            all_notations, parameters["mutational_matrix"]
        )

        # Update mutation scores in aa_changes
        for mutation in mutations:
            for aa_change in mutation["aa_changes"]:
                notation = aa_change["notation"]
                clean_notation = notation.rstrip("*")
                if clean_notation in mutation_scores:
                    aa_change["mut_score"] = mutation_scores[clean_notation]

    return mutations, skipped


SEQUENCES_PER_BATCH = 500


MAX_DEFAULT_WORKERS = 16


def _worker_count() -> int:
    """Workers to use, honouring KSS_WORKERS, and 1 when joblib is not installed.

    The default is capped well below the core count. Each worker holds its own alignment,
    which on the largest reference here is a 7,097 by 7,097 dynamic programming problem
    costing on the order of 140 MB, while the parent accumulates every result before writing
    them. Taking every core would trade a modest gain for the risk of exhausting memory
    partway through a run of several hours.
    """
    if Parallel is None:
        return 1
    requested = os.environ.get("KSS_WORKERS")
    if requested:
        return max(1, int(requested))
    return max(1, min((os.cpu_count() or 1) - 2, MAX_DEFAULT_WORKERS))


def _batched(records, size: int):
    """Read the FASTA in batches, so the whole file is never held in memory at once."""
    batch = []
    for record in records:
        batch.append((record.id, str(record.seq)))
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _process_batch(batch: list, infos: dict, parameters: dict,
                   ref_nuc: str, ref_aa: str) -> list:
    """One batch of sequences, aligned and scanned exactly as the serial loop did.

    Defined at module level, and taking only picklable arguments, because a worker process on
    Windows is started fresh rather than forked: it re-imports this module and rebuilds its own
    aligner from `parameters`, so nothing about the parent's state can leak into the result.
    """
    aligner = initialize_aligner(parameters)
    processed = []

    for query_id, query_nuc in batch:
        # Biopython's PairwiseAligner does not support 'J' (ambiguous Leu/Ile), which is
        # normalised to 'L', the closest standard amino acid. Rare, and only in sequences with
        # a non-standard translation.
        query_aa = str(Seq(query_nuc).translate()).replace('J', 'L')

        # Class label from the sequence header, format >ID|Virus|Gene|CLASS
        class_label = query_id.split('|')[-1] if '|' in query_id else "unknown"

        mutations, skipped = identify_mutations(
            infos, parameters, ref_nuc, ref_aa, query_nuc, query_aa, aligner)
        processed.append((query_id, class_label, mutations, skipped))

    return processed


def calculate_similarity(seq_a: str, seq_b: str) -> float:
    """
    Calculate the percentage similarity between two sequences.

    Args:
        seq_a: First sequence
        seq_b: Second sequence

    Returns:
        Percentage similarity (0-100)
    """
    matches = sum(a == b for a, b in zip(seq_a, seq_b))
    total = len(seq_a)
    return (100 * matches / total) if total else 0


def analyze_records(infos: dict, parameters: dict) -> dict:
    """Align query sequences to the GenBank reference and record their mutations.

    Main entry point for sequence analysis: loads the reference and the FASTA queries for
    each selected gene and identifies mutations, storing only positions that carry one.

    Args:
        infos: needs 'input_folder' and 'cds_selection' (comma-separated gene list).
        parameters: analysis parameters (k-mer size, alignment penalties, ...).

    Returns:
        {"_metadata": {...}, "genes": {gene: {"metadata": {...},
        "sequences": {seq_id: {"class": label, "mutations": [...]}}}}}.
    """
    # Initialize results structure with metadata
    results = {
        "_metadata": {
            "version": "2.0.0",
            "tool": "K-mer Significance Score (KSS)",
            "k": parameters["k"]
        },
        "genes": {}
    }

    cds_list = infos["cds_selection"].split(",")
    total_sequences = 0

    for cds in cds_list:
        # Get reference data
        gb_dir = os.path.join(infos["input_folder"], cds)
        gb_files = [f for f in os.listdir(gb_dir) if f.lower().endswith('.gb')]
        if not gb_files:
            print(f"  ⚠ No GenBank file found for {cds}")
            continue

        gb_path = os.path.join(gb_dir, gb_files[0])
        gb_record = list(SeqIO.parse(gb_path, "genbank"))[0]

        # Extract reference sequences. A gene qualifier does not always identify one CDS:
        # NC_045512 annotates both the ORF1ab and the ORF1a polyprotein with gene="ORF1ab",
        # 7,096 and 4,405 residues, the second being the product of the same reading frame
        # without the ribosomal frameshift. Taking the first match made the choice depend on
        # the order of the features in the record. The longest translation is taken instead,
        # which is the whole gene and is what the reported results were computed from.
        matches = [feature for feature in gb_record.features
                   if feature.type == "CDS"
                   and feature.qualifiers.get("gene", [""])[0] == cds
                   and feature.qualifiers.get("translation")]
        ref_nuc = ""
        ref_aa = ""
        if matches:
            feature = max(matches, key=lambda f: len(f.qualifiers["translation"][0]))
            ref_nuc = str(feature.location.extract(gb_record.seq))
            ref_aa = feature.qualifiers["translation"][0]

        if not ref_nuc or not ref_aa:
            print(f"  ⚠ CDS {cds} not found in GenBank file")
            continue

        # Get query data
        fasta_dir = os.path.join(infos["input_folder"], cds)
        fasta_files = [f for f in os.listdir(fasta_dir) if f.lower().endswith('.fasta')]
        if not fasta_files:
            print(f"  ⚠ No FASTA file found for {cds}")
            continue

        fasta_path = os.path.join(fasta_dir, fasta_files[0])

        # Count sequences first without loading them all
        print(f"  Counting {cds} sequences...", end='', flush=True)
        num_sequences = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
        print(f" {num_sequences} sequences found")

        # Initialize gene structure
        k = parameters["k"]
        # From the translated reference, as the scoring grid does: len(ref_nuc) counted the stop codon.
        num_positions = (len(ref_aa) * 3) // k

        results["genes"][cds] = {
            "metadata": {
                "reference_length": len(ref_nuc),
                "reference_aa_length": len(ref_aa),
                "genbank_id": gb_record.id,
                "total_positions": num_positions,
                "num_sequences": num_sequences
            },
            "sequences": {}
        }

        # position -> how many reference nucleotides that window could actually be given
        incomplete_windows: dict[int, int] = {}

        # Pre-initialize aligner once (reuse across all sequences)
        aligner = initialize_aligner(parameters)

        # Process each query sequence - stream instead of loading all at once
        print(f"  Processing {cds} sequences: ", end='', flush=True)

        # Sequences are independent: each one is aligned against the same reference and shares
        # nothing with the others, so the work splits cleanly. It is split in batches rather
        # than one task per sequence, which is what makes the difference here: a single
        # alignment takes a few milliseconds, while starting a worker on Windows and shipping
        # its arguments costs far more than that, so per-sequence tasks spend their time being
        # distributed. A batch amortises that cost over thousands of alignments.
        #
        # The work each batch does is bit for bit the work the serial loop did, and the batches
        # are reassembled in the order they were read, so the results do not depend on how many
        # workers ran or on the order they finished in.
        workers = _worker_count()
        batches = _batched(SeqIO.parse(fasta_path, "fasta"), SEQUENCES_PER_BATCH)

        if workers > 1:
            print(f"\r  Processing {cds} sequences: {num_sequences:,} in batches of "
                  f"{SEQUENCES_PER_BATCH} across {workers} workers", end='', flush=True)
            # The batches stay a generator and `pre_dispatch` bounds how many are held at once,
            # so the file is still streamed rather than loaded whole. Materialising them would
            # cost, on the largest gene here, 278,738 sequences of some 21 kb held together
            # with every result they produce. Order is preserved by joblib regardless of which
            # worker finishes first.
            # `verbose` matters on the largest gene: without it a stage that runs for over an
            # hour prints nothing at all, which is indistinguishable from a hung process.
            processed = Parallel(n_jobs=workers, backend="loky",
                                 pre_dispatch="2 * n_jobs", verbose=5)(
                delayed(_process_batch)(batch, infos, parameters, ref_nuc, ref_aa)
                for batch in batches)
        else:
            processed = (_process_batch(batch, infos, parameters, ref_nuc, ref_aa)
                         for batch in batches)

        for batch_results in processed:
            for query_id, class_label, mutations, skipped in batch_results:
                for window in skipped:
                    incomplete_windows[window["position"]] = window["reference_nucleotides"]
                results["genes"][cds]["sequences"][query_id] = {
                    "class": class_label,
                    "mutations": mutations
                }
                total_sequences += 1

        # Print summary for this gene
        print(f"\n  ✓ {cds}: {num_sequences} sequences processed")

        # Windows the reference could not fill are recorded, not merely skipped: a reader of
        # these results should be able to see that a position was declined and why, without
        # having to notice its absence.
        if incomplete_windows:
            positions = sorted(incomplete_windows)
            results["genes"][cds]["metadata"]["windows_without_full_reference"] = {
                str(position): incomplete_windows[position] for position in positions
            }
            where = (str(positions[0]) if len(positions) == 1
                     else f"{positions[0]}-{positions[-1]}")
            print(f"    {len(positions)} window(s) declined at {where}: "
                  f"reference incomplete (see metadata)")

    # Update global metadata
    results["_metadata"]["sequences_analyzed"] = total_sequences

    return results