import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import itertools

# mostly just to ensure that the input is the same
def set_sequences(seq1, seq2):
    """
        Clean and normalize two input sequences.

        Removes whitespace, newlines, and tabs from each sequence,
        and converts them to uppercase. Intended to ensure consistent formatting
        before alignment.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.

        Returns:
            tuple: A tuple containing the cleaned versions of (seq1, seq2).
    """

    # Remove spaces, tabs, and newlines
    seq1 = seq1.replace(" ", "").replace("\n", "").replace("\t", "").upper()
    seq2 = seq2.replace(" ", "").replace("\n", "").replace("\t", "").upper()

    return seq1, seq2


def matrix_building(seq1, seq2, gap_value):
    """
        Initialize the scoring matrix with gap penalties for the Needleman-Wunsch algorithm.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            gap_value (float): Penalty value for introducing a gap.

        Returns:
            np.ndarray: Initialized scoring matrix with gap penalties.
    """

    m = len(seq1)
    n = len(seq2)
    zero_matrix = np.zeros((m+1,n+1))

    # Fill out first column
    for i in range(0, m + 1):
        zero_matrix[i][0] = gap_value * i

    # Fill out first row
    for j in range(0, n + 1):
        zero_matrix[0][j] = gap_value * j

    return zero_matrix


def algorithm(seq1, seq2, array, gap_value, match_value, mismatch_value):
    """
        Fill in the scoring matrix using the Needleman-Wunsch algorithm.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            array (np.ndarray): Initialized scoring matrix.
            gap_value (float): Penalty for a gap.
            match_value (float): Reward for a match.
            mismatch_value (float): Penalty for a mismatch.

        Returns:
            np.ndarray: Completed scoring matrix.
    """

    m = len(seq1)
    n = len(seq2)

    # recursion algorithm that is found on https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # for the row - will be penalty
            if seq1[i-1] == seq2[j-1]:
                score = array[i-1][j-1] + match_value
            else:
                score =  array[i-1][j-1] + mismatch_value

            score2 = array[i, j - 1] + gap_value
            score3 = array[i - 1, j] + gap_value
            # getting the max score
            array[i][j] = max(score, score2, score3)

    return array

def reconstruct_alignment(seq1, seq2, path):
    """
        Reconstruct two aligned sequences from a traceback path.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            path (list): Traceback path from the alignment matrix.

        Returns:
            tuple: (aligned_seq1, aligned_seq2), both containing possible gaps.
    """

    aligned_seq1 = ""
    aligned_seq2 = ""

    for idx in range(1, len(path)):
        i_curr, j_curr = path[idx - 1]
        i_next, j_next = path[idx]

        if i_next == i_curr + 1 and j_next == j_curr + 1:
            # diagonal: match/mismatch
            aligned_seq1 += seq1[i_curr]
            aligned_seq2 += seq2[j_curr]
        elif i_next == i_curr + 1 and j_next == j_curr:
            # up: gap in seq2
            aligned_seq1 += seq1[i_curr]
            aligned_seq2 += "-"
        elif i_next == i_curr and j_next == j_curr + 1:
            # left: gap in seq1
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j_curr]

    return aligned_seq1, aligned_seq2

def merge_alignment(prev_center, new_center, new_seq, aligned_others):
    """
        Merge a new sequence into an existing MSA using a new center sequence.

        Args:
            prev_center (str): Previously aligned central sequence.
            new_center (str): Newly aligned central sequence (with gaps).
            new_seq (str): Newly aligned sequence corresponding to new_center.
            aligned_others (list): List of other sequences already aligned to prev_center.

        Returns:
            tuple: (new_center, realigned_new_seq, updated_others)
    """

    # Project both the previous center and new_seq onto new_center
    updated_prev = project_onto_master(new_center, prev_center)
    updated_new = project_onto_master(new_center, new_seq)

    # Project all previously aligned sequences (aligned to prev_center) onto new_center
    updated_others = [project_onto_master(new_center, s) for s in aligned_others]

    # Insert the realigned previous center at the beginning of others (it’s part of the alignment)
    updated_others.append(updated_new)


    return new_center, updated_new, updated_others


def traceback(score_matrix, seq1, seq2, gap, match, mismatch):
    """
        Trace back through the scoring matrix to recover the optimal alignment path.

        Args:
            score_matrix (np.ndarray): Completed scoring matrix.
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            gap (float): Gap penalty.
            match (float): Match reward.
            mismatch (float): Mismatch penalty.

        Returns:
            list: List of (i, j) tuple indices representing the optimal alignment path.
    """

    i, j = len(seq1), len(seq2)
    path = [(i, j)]

    while i > 0 or j > 0:
        current_score = score_matrix[i, j]

        if i > 0 and j > 0:
            diag = score_matrix[i-1, j-1]
            score_diag = match if seq1[i-1] == seq2[j-1] else mismatch
            if current_score == diag + score_diag:
                i -= 1
                j -= 1
                path.append((i, j))
                continue

        if i > 0 and current_score == score_matrix[i-1, j] + gap:
            i -= 1
            path.append((i, j))
            continue

        if j > 0 and current_score == score_matrix[i, j-1] + gap:
            j -= 1
            path.append((i, j))
            continue

    path.reverse()
    return path

def project_onto_master(master: str, seq):
    """
        Project an unaligned sequence onto a gapped master sequence.

        Inserts gaps into `seq` where `master` has them to maintain column consistency.

        Args:
            master (str): Gapped reference sequence.
            seq (str): Sequence to align with the master.

        Returns:
            str: Gapped version of `seq` aligned to `master`.
    """

    out = []
    p = 0
    for c in master:
        if c == "-":

            if len(master) > len(seq):
                out.append("-")
            else:
                continue
        else:
            if p < len(seq):
                out.append(seq[p])
                p += 1

    while p < len(seq):
        out.append(seq[p])
        p += 1

    return "".join(out)


def calculate_identity_percentage(msa_sequences):
    """
        Calculate average pairwise identity percentage for an MSA.

        Args:
            msa_sequences (list): List of aligned sequences.

        Returns:
            float: Average identity percentage (0–100), rounded to 2 decimal places.
    """

    num_sequences = len(msa_sequences)
    alignment_length = len(msa_sequences[0])
    total_identity = 0
    pair_count = 0

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            matches = 0
            valid_positions = 0
            for a, b in zip(msa_sequences[i], msa_sequences[j]):
                if a != '-' and b != '-':
                    valid_positions += 1
                    if a == b:
                        matches += 1
            if valid_positions > 0:
                identity = (matches / valid_positions) * 100
                total_identity += identity
                pair_count += 1

    average_identity = total_identity / pair_count if pair_count > 0 else 0
    return round(average_identity, 2)


# Compute MSA score
def calculate_msa_score(msa, match, mismatch, gap):
    """
        Calculate the total MSA score based on pairwise scoring.

        Args:
            msa (list): List of aligned sequences.
            match (float): Score for a match.
            mismatch (float): Penalty for mismatch.
            gap (float): Penalty for a gap.

        Returns:
            float: Total score for the MSA.
    """

    total_score = 0
    num_seqs = len(msa)
    seq_length = len(msa[0])

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = msa[i]
            seq2 = msa[j]
            for k in range(seq_length):
                a = seq1[k]
                b = seq2[k]
                if a == '-' or b == '-':
                    total_score += gap
                elif a == b:
                    total_score += match
                else:
                    total_score += mismatch
    return total_score

def count_msa_statistics(msa_sequences):
    """
        Count matches, mismatches, and gaps across all pairs in the MSA.

        Args:
            msa_sequences (list): List of aligned sequences.

        Returns:
            tuple: (matches, mismatches, gaps)
    """

    alignment_length = len(msa_sequences[0])
    num_sequences = len(msa_sequences)

    matches = 0
    mismatches = 0
    gaps = 0

    for i in range(alignment_length):
        column = [seq[i] for seq in msa_sequences]
        for a, b in itertools.combinations(column, 2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1

    return matches, mismatches, gaps
