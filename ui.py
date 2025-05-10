from alignment import *
from file_utils import *
import streamlit as st
import numpy as np

def app_creation():
    st.title("MSA Star Algorithm")

    # Initialize session state
    if "sequences" not in st.session_state:
        st.session_state.sequences = [""]

    # Function to add a new sequence
    def add_sequence():
        st.session_state.sequences.append("")

    # Function to remove a sequence by index
    def remove_sequence(index):
        if len(st.session_state.sequences) > 1:
            st.session_state.sequences.pop(index)

    # UI rendering
    new_sequences = []

    st.markdown("### Input Sequences")

    for i, seq in enumerate(st.session_state.sequences):
        cols = st.columns([6, 1])
        # Text input
        new_value = cols[0].text_input(f"Sequence {i + 1}", value=seq, key=f"seq_{i}")
        new_sequences.append(new_value)

        # Show remove button only if more than one input remains
        if len(st.session_state.sequences) > 1:
            if cols[1].button("❌", key=f"remove_{i}"):
                remove_sequence(i)
                st.rerun()

    # Update the sequence values
    st.session_state.sequences = new_sequences

    # Add button
    if st.button("➕ Add a new sequence"):
        add_sequence()
        st.rerun()

    # Display current sequences
    st.markdown("### Current Sequences")
    st.write(st.session_state.sequences)

    st.divider()

    col1, col2, col3 = st.columns(3)
    with col1:
        gap_value = st.number_input("Gap penalty", value=-2.0, step=1.0, format="%.2f")

    with col2:
        match_value = st.number_input("Match reward", value=1.0, step=1.0, format="%.2f")

    with col3:
        mismatch_value = st.number_input("Mismatch penalty", value=-1.0, step=1.0, format="%.2f")

    # Run button
    if st.button("Run MSA Star"):
        st.success("Running MSA Star with these sequences:")
        st.write(st.session_state.sequences)

        final_scores_list = []
        final_score_matrix = np.zeros((len(st.session_state.sequences), len(st.session_state.sequences)))
        i = 0
        j = 1
        for seq1, seq2 in itertools.combinations(new_sequences, 2):
            try:
                seq1_clean, seq2_clean = set_sequences(seq1, seq2)
                matrix = matrix_building(seq1_clean, seq2_clean, gap_value)
                score_matrix = algorithm(seq1_clean, seq2_clean, matrix, gap_value, match_value, mismatch_value)

                final_score = score_matrix[-1, -1]

                length = len(new_sequences)

                if i < length and j < length:
                    final_score_matrix[i, j] = final_score
                    final_score_matrix[j, i] = final_score
                    j += 1
                    if j == length:
                        i += 1
                        j = i + 1

                final_scores_list.append(final_score)

            except st.exception:
                st.warning(f"Skipping pair ({seq1}, {seq2}) due to input error.")


        st.write(final_scores_list)
        st.write(final_score_matrix)


