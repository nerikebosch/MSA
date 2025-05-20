import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import itertools

# mostly just to ensure that the input is the same
def set_sequences(seq1, seq2):
    """
        Cleans and validates input sequences.
        Removes whitespace and ensures only valid letters (no digits or special characters).
        Raises Streamlit error if validation fails.
        Convert both input sequences to uppercase.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.

        Returns:
            tuple: Uppercased first and second sequences.
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
    Reconstructs aligned sequences (with gaps), based on an alignment path.

    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        path (list): List of (i, j) positions along the alignment path.

    Returns:
        tuple: (aligned_seq1, aligned_seq2) with gaps inserted.
    """

    #seq1 = seq1.replace("-","")
    aligned_seq1 = ""
    aligned_seq2 = ""
    st.write("aligned sequence 1" + seq1 + " " + seq2)

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
    Simplified star‐merge:
      • new_center is now the one and only master.
      • For each old seq in aligned_others, and for new_seq, rebuild them
        by inserting a '-' wherever new_center has one, otherwise
        consuming the next character from the old seq.
    """
    st.write(">>> merge_alignment()")
    st.write(" prev_center:", prev_center)
    st.write(" new_center :", new_center)
    st.write(" new_seq    :", new_seq)
    st.write(" #others    :", len(aligned_others))
    for k, s in enumerate(aligned_others):
        st.write(f"  other[{k}]: {s}")

    # Project both the previous center and new_seq onto new_center
    updated_prev = project_onto_master(new_center, prev_center)
    updated_new = project_onto_master(new_center, new_seq)

    # Project all previously aligned sequences (aligned to prev_center) onto new_center
    updated_others = [project_onto_master(new_center, s) for s in aligned_others]

    # Insert the realigned previous center at the beginning of others (it’s part of the alignment)
    updated_others.append(updated_new)

    # 3) The new_center itself becomes the next prev_center
    return new_center, updated_new, updated_others


def traceback(score_matrix, seq1, seq2, gap, match, mismatch):
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
    Given a gapped master template and a raw (ungapped) sequence,
    return a new string of the same length as `master` that
    places each character of raw_seq under the non-gap columns,
    and inserts '-' wherever master has a gap.
    """
    st.write(" -- project_onto_master --")
    st.write("  master:", master)
    st.write("  input seq:", seq)

    out = []
    p = 0
    for c in master:
        if c == "-":
            if len(master) > len(seq):
                out.append("-")
                st.write(f"   col 1: master gap → output '-'")
            else:
                # Skip adding a gap if seq is already long enough
                st.write(f"   col 1: master gap ignored — seq already aligned")
        else:
            if p < len(seq):
                out.append(seq[p])
                st.write(f"   col 2: master '{c}' → take seq[{p}]='{seq[p]}'")
                p += 1
    #         else:
    #             out.append("-")
    #             st.write(f"   col 3: master '{c}' but seq exhausted → '-'")
    #
    while p < len(seq):
        out.append(seq[p])
    #     master += '-'
        p += 1

    return "".join(out)


def calculate_identity_percentage(msa_sequences):
    num_sequences = len(msa_sequences)
    alignment_length = len(msa_sequences[0])
    total_identity = 0
    pair_count = 0

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            matches = sum(
                1 for a, b in zip(msa_sequences[i], msa_sequences[j]) if a == b
            )
            identity = (matches / alignment_length) * 100
            total_identity += identity
            pair_count += 1

    average_identity = total_identity / pair_count if pair_count > 0 else 0
    return round(average_identity, 2)
