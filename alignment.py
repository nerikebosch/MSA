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

    # Check for invalid characters (anything that's not a letter)
    if not seq1.isalpha():
        st.error("First sequence contains invalid characters. Only letters (A-Z) are allowed.")
        st.stop()

    if not seq2.isalpha():
        st.error("Second sequence contains invalid characters. Only letters (A-Z) are allowed.")
        st.stop()

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

