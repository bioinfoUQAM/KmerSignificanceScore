import re
import json
import warnings
from Bio.Seq import Seq
from Bio import BiopythonWarning
from Bio.Align import PairwiseAligner, substitution_matrices

# Suppress BiopythonWarning from Biopython
warnings.filterwarnings("ignore", category=BiopythonWarning)

def load_data_from_json(filename):
    """
    Load data from a JSON file.

    Parameters:
        filename (str): Path to the JSON file.

    Returns:
        dict: Data loaded from the JSON file.

    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file contains invalid JSON.
    """
    try:
        with open(filename, 'r') as file:
            return json.load(file)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from {filename}: {e}")
        raise
    except FileNotFoundError as e:
        print(f"File not found: {filename}")
        raise

def find_orfs(dna_sequence):
    """
    Find all open reading frames (ORFs) in the given DNA sequence.

    Parameters:
        dna_sequence (str): DNA sequence.

    Returns:
        list: List of translated amino acid sequences from all ORFs.
    """
    seq = Seq(dna_sequence)

    def translate_frames(sequence):
        return [str(sequence[frame:].translate(to_stop=False)) for frame in range(3)]

    orfs = translate_frames(seq) + translate_frames(seq.reverse_complement())
    return orfs

def initialize_aligner(matrix, open_gap_score, extend_gap_score, end_open_gap_score=None, end_extend_gap_score=None):
    """
    Initialize a pairwise aligner with specified parameters.

    Parameters:
        mode (str): Alignment mode.
        matrix (str): Substitution matrix name.
        open_gap_score (float): Gap opening penalty.
        extend_gap_score (float): Gap extension penalty.
        end_open_gap_score (float, optional): Gap opening penalty at sequence ends.
        end_extend_gap_score (float, optional): Gap extension penalty at sequence ends.

    Returns:
        PairwiseAligner: Configured aligner.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    aligner.end_open_gap_score = end_open_gap_score if end_open_gap_score is not None else open_gap_score
    aligner.end_extend_gap_score = end_extend_gap_score if end_extend_gap_score is not None else extend_gap_score
    return aligner

def calculate_similarity(seqA, seqB):
    """
    Calculate the percentage similarity between two sequences.

    Parameters:
        seqA (str): First sequence.
        seqB (str): Second sequence.

    Returns:
        float: Percentage of matching characters between the sequences.
    """
    matches = sum(a == b for a, b in zip(seqA, seqB))
    total = len(seqA)
    return (100 * matches / total) if total else 0

def align_nucleotides(protein_a, protein_b, nucleotide_a, nucleotide_b):
    """
    Align nucleotide sequences based on aligned protein sequences.

    Parameters:
        protein_a (str): Aligned protein sequence A.
        protein_b (str): Aligned protein sequence B.
        nucleotide_a (str): Original nucleotide sequence A.
        nucleotide_b (str): Original nucleotide sequence B.

    Returns:
        tuple: Aligned nucleotide sequences (aligned_nuc_a, aligned_nuc_b).
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

def replace_first_three_gaps(sequence, replacement=''):
    """
    Replace the first three gap characters '-' in a sequence with a replacement string.

    Parameters:
        sequence (str): The sequence containing gaps.
        replacement (str): The string to replace gaps with.

    Returns:
        str: Modified sequence with first three gaps replaced.
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

def adjust_mutation(mutation, increment):
    """
    Adjust the numerical position in a mutation string by a specified increment.

    Parameters:
        mutation (str): Mutation string in the format 'A123B'.
        increment (int): Value to add to the position number.

    Returns:
        str: Adjusted mutation string.

    Raises:
        ValueError: If the mutation string is not in the expected format.
    """
    parts = re.split(r'(\d+)', mutation)
    if len(parts) < 3:
        raise ValueError(f"Mutation string '{mutation}' is not in expected format.")
    new_number = int(parts[1]) + increment
    return f'{parts[0]}{new_number}{parts[2]}'


def add_unique_mutation(mutation, mutations_list):
    original_mutation = mutation
    suffix = ""
    while mutation in mutations_list:
        suffix += "*"
        mutation = original_mutation + suffix
    mutations_list.append(mutation)

def identify_mutations(reference_nucleotide_sequence, reference_amino_acid_sequence,
                       query_amino_acid_sequence, query_nucleotide_sequence, k, infos, gene):
    """
    Identify mutations between reference and query sequences.

    Parameters:
        reference_nucleotide_sequence (str): Reference nucleotide sequence.
        reference_amino_acid_sequence (str): Reference amino acid sequence.
        query_amino_acid_sequence (str): Query amino acid sequence.
        query_nucleotide_sequence (str): Query nucleotide sequence.
        k (int): K-mer size.
        infos (dict): Alignment parameters.
        gene (str): Gene name.

    Returns:
        dict: Dictionary of identified mutations.
    """
    if True:
        # Extract alignment parameters
        matrix = infos["aligner_matrix"]
        open_gap_score = infos["aligner_open_gap_score"]
        extend_gap_score = infos["aligner_extend_gap_score"]

        # Initialize the aligner
        aligner = initialize_aligner(matrix, open_gap_score, extend_gap_score)

        # Prepare sequences for alignment
        ref_aa_seq_clean = reference_amino_acid_sequence.replace('*', '')
        query_aa_seq_clean = query_amino_acid_sequence.replace('*', '').replace('J', 'L')

        # Perform amino acid alignment
        amino_acid_alignments = aligner.align(ref_aa_seq_clean, query_aa_seq_clean)
        aligned_ref_aa_seq, aligned_query_aa_seq = amino_acid_alignments[0]

        # Align nucleotide sequences based on amino acid alignment
        aligned_ref_nuc_seq, aligned_query_nuc_seq = align_nucleotides(
            aligned_ref_aa_seq,
            aligned_query_aa_seq,
            reference_nucleotide_sequence,
            query_nucleotide_sequence
        )

        # Initialize results dictionary with appropriate offsets
        if gene == "gag-pol":
            offset = 1636
            mutation_offset = 434
        else:
            offset = 0
            mutation_offset = 0

        results = {}
        for i in range(0, len(reference_nucleotide_sequence) - k + 1, k):
            kmer = reference_nucleotide_sequence[i:i + k]
            key = (i + 1 + offset, kmer)
            results[key] = {"variations": "", "amino_acid_changes": []}

        n_insertion = 0  # Total number of insertions encountered

        # Loop over positions in the aligned nucleotide sequences
        for i in range(0, len(aligned_ref_nuc_seq) - k + 1, k):
            adjusted_k = k
            idx = i + (n_insertion * 3)

            # Extract k-mer from aligned sequences
            kmer_ref_nuc = aligned_ref_nuc_seq[idx:idx + adjusted_k]
            kmer_query_nuc = aligned_query_nuc_seq[idx:idx + adjusted_k]
            
            # Adjust k-mer size if gaps are present in the reference nucleotide sequence
            while '-' in kmer_ref_nuc:
                temp_kmer_ref_nuc = kmer_ref_nuc
                next_indices = idx + adjusted_k + 3
                next_nucleotides = aligned_ref_nuc_seq[idx + adjusted_k:next_indices]
                kmer_ref_nuc = replace_first_three_gaps(temp_kmer_ref_nuc) + next_nucleotides
                adjusted_k += 3
                kmer_query_nuc = aligned_query_nuc_seq[idx:idx + adjusted_k]
                if idx + adjusted_k > len(aligned_ref_nuc_seq):
                    break
            # Update insertion counts
            current_insertion = (adjusted_k - k) // 3
            n_insertion += current_insertion
            
            # Determine amino acid positions in the alignment
            start_aa = (idx // 3) #- current_insertion
            end_aa = start_aa + (adjusted_k // 3)
            #print("current_insertion", current_insertion, "n_insertion", n_insertion, "start_aa", start_aa, "end_aa", end_aa)
            if end_aa > len(aligned_ref_aa_seq):
                end_aa = len(aligned_ref_aa_seq)

            kmer_ref_aa = aligned_ref_aa_seq[start_aa:end_aa]
            kmer_query_aa = aligned_query_aa_seq[start_aa:end_aa]

            # Identify mutations
            mutations = []
            position = start_aa - n_insertion
            considered_insertion = 0

            for idx_aa, (aa_ref, aa_query) in enumerate(zip(kmer_ref_aa, kmer_query_aa)):
                if aa_ref != aa_query:
                    # Calculate mutation position
                    if position + idx_aa + 1 < 0:
                        mutation_pos = 0
                    else:
                        mutation_pos = position + idx_aa + 1 + current_insertion - considered_insertion
                    mutation = f"{aa_ref}{mutation_pos}{aa_query}"
                    if aa_ref == '-':
                        considered_insertion += 1
                    original_mutation = mutation
                    suffix = ""
                    while mutation in mutations:
                        suffix += "*"
                        mutation = original_mutation + suffix
                    mutations.append(mutation)

            # Update the results
            key = (i + 1 + offset, kmer_ref_nuc)
            if key in results:
                if gene == "gag-pol":
                    adjusted_mutations = [adjust_mutation(mutation, mutation_offset) for mutation in mutations]
                    results[key] = {
                        "variations": kmer_query_nuc,
                        "amino_acid_changes": adjusted_mutations
                    }
                else:
                    results[key] = {
                        "variations": kmer_query_nuc,
                        "amino_acid_changes": mutations
                    }
        
        return results

def process_sequence(start_motifs, ref_nuc_seq, ref_aa_seq, query_nuc_seq, orfs, k, infos, gene, seq_id):
    """
    Process a single sequence to identify mutations.

    Parameters:
        start_motifs (list): List of start motifs to look for.
        ref_nuc_seq (str): Reference nucleotide sequence.
        ref_aa_seq (str): Reference amino acid sequence.
        query_nuc_seq (str): Query nucleotide sequence.
        orfs (list): List of open reading frames.
        k (int): K-mer size.
        infos (dict): Alignment and scoring parameters.
        gene (str): Gene name.
        seq_id (str): Sequence identifier.

    Returns:
        dict: Identified mutations or None if not processed.
    """
    for i, orf in enumerate(orfs):
        is_reverse = i >= 3
        frame_offset = i % 3

        for motif in start_motifs:
            match_start_index = orf.find(motif)
            if match_start_index != -1:
                dna_start_index = 3 * match_start_index + frame_offset
                if is_reverse:
                    dna_start_index = len(query_nuc_seq) - (dna_start_index + 3 * len(orf[match_start_index:]))

                # Find the stop codon in the ORF
                stop_codon_index = orf.find('*', match_start_index)

                # Extract the sequence up to the stop codon
                if stop_codon_index != -1:
                    query_aa_seq = orf[match_start_index:stop_codon_index + 1]
                    dna_length = (stop_codon_index + 1 - match_start_index) * 3
                    query_nuc_seq = query_nuc_seq[dna_start_index:dna_start_index + dna_length]
                else:
                    query_aa_seq = orf[match_start_index:]
                    dna_length = len(query_aa_seq) * 3
                    query_nuc_seq = query_nuc_seq[dna_start_index:dna_start_index + dna_length]

                output = identify_mutations(
                    ref_nuc_seq,
                    ref_aa_seq,
                    query_aa_seq,
                    query_nuc_seq,
                    k,
                    infos,
                    gene
                )
                return output
    print(f"{gene} {seq_id} not processed")
    return None

def analyze_records(records, infos):
    """
    Analyze multiple sequences to identify mutations across records.

    Parameters:
        records (list): List of sequence records.
        infos (dict): Information and parameters for analysis.

    Returns:
        dict: Results containing identified mutations for each gene and sequence.
    """
    k = infos['k']
    infos_file = load_data_from_json(infos['infos_file'])
    organism = next(iter(infos_file))
    results = {}

    for record in records:
        query_nuc_seq = str(record.seq.upper())
        seq_id = record.description
        orfs = find_orfs(query_nuc_seq)
        if True:
            for gene, gene_info in infos_file[organism].items():
                start_motifs = gene_info["start_motifs"]
                ref_nuc_seq = gene_info["CDS"]
                ref_aa_seq = gene_info["translation"]

                if gene not in results:
                    results[gene] = {}
                results[gene][seq_id] = process_sequence(
                    start_motifs,
                    ref_nuc_seq,
                    ref_aa_seq,
                    query_nuc_seq,
                    orfs,
                    k,
                    infos,
                    gene,
                    seq_id
                )
    return results
