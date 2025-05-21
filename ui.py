import datetime

from alignment import *
from file_utils import *
from graph import *
import streamlit as st
import numpy as np
import itertools

def app_creation():
    """
        Initialize and run the Streamlit-based MSA (Multiple Sequence Alignment) Star Algorithm application.

        This function sets up the user interface using Streamlit, including:
        - Dynamic sequence input (manual or via FASTA file)
        - User-defined scoring parameters (match, mismatch, gap)
        - Execution of the Star MSA algorithm
        - Calculation and display of alignment statistics
        - Saving the final results to a text file

        The Star alignment method identifies a central sequence with the highest pairwise score sum,
        and aligns all other sequences to it.

        Uses:
            - load_fasta_sequences: To load sequences from a file
            - set_sequences, matrix_building, algorithm, traceback, reconstruct_alignment,
              merge_alignment, calculate_identity_percentage, calculate_msa_score,
              count_msa_statistics: For MSA processing
            - get_text, save_to_text_file: For result formatting and saving

        Returns:
            None
    """

    st.set_page_config(layout="wide")
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

    new_sequences = []

    st.markdown("### Input Sequences")
    col1, col2 = st.columns(2)
    with col1:
        for i, seq in enumerate(st.session_state.sequences):
            cols = st.columns([6, 1])
            new_value = cols[0].text_input(
                f"Sequence {i + 1}", value=seq, key=f"seq_{i}"
            )
            new_sequences.append(new_value)
            if len(st.session_state.sequences) > 1:
                if cols[1].button(
                    "❌", key=f"remove_{i}"
                ):
                    remove_sequence(i)
                    st.rerun()
        st.session_state.sequences = new_sequences
        if st.button("➕ Add a new sequence"):
            add_sequence()
            st.rerun()

    with col2:
        seq_file = st.file_uploader(
            "Choose a file for all the sequences", type="fasta", key="seq_file"
        )
        if seq_file is not None:
            st.session_state.sequences = load_fasta_sequences(seq_file)

    st.markdown("### Current Sequences")
    st.write(st.session_state.sequences)
    st.divider()

    col1, col2, col3 = st.columns(3)
    with col1:
        gap_value = st.number_input(
            "Gap penalty", value=-2.0, step=1.0, format="%.2f"
        )
    with col2:
        match_value = st.number_input(
            "Match reward", value=1.0, step=1.0, format="%.2f"
        )
    with col3:
        mismatch_value = st.number_input(
            "Mismatch penalty", value=-1.0, step=1.0, format="%.2f"
        )

    if st.button("Run MSA Star"):
        # Check for empty or invalid sequences
        valid_chars = set("ACGTURYKMSWBDHVN-")  # Adjust depending on DNA/RNA/protein
        invalid_sequences = []
        empty_sequences = []

        for i, seq in enumerate(new_sequences):
            seq_upper = seq
            if not seq_upper:
                empty_sequences.append(i + 1)

        if empty_sequences:
            st.error(f"Sequence(s) {empty_sequences} are empty. Please provide valid input.")
        elif invalid_sequences:
            bad_info = "\n".join([f"Sequence {i}: {s}" for i, s in invalid_sequences])
            st.error(f"The following sequence(s) contain invalid characters:\n{bad_info}")
        else:
            st.success("Running MSA Star with these sequences:")
            st.write(st.session_state.sequences)

            # Build pairwise score matrix
            length = len(new_sequences)
            final_score_matrix = np.zeros((length, length))
            idx_i = 0
            for i, j in itertools.combinations(range(length), 2):
                seq1_clean, seq2_clean = set_sequences(
                    new_sequences[i], new_sequences[j]
                )
                mat = matrix_building(seq1_clean, seq2_clean, gap_value)
                score_mat = algorithm(
                    seq1_clean,
                    seq2_clean,
                    mat,
                    gap_value,
                    match_value,
                    mismatch_value,
                )
                score = score_mat[-1, -1]
                final_score_matrix[i, j] = score
                final_score_matrix[j, i] = score

            row_sums = final_score_matrix.sum(axis=1)
            max_row_index = int(np.argmax(row_sums))

            # Star alignment
            center_seq = new_sequences[max_row_index]
            aligned_center = center_seq
            aligned_others = []
            order = [i for i in range(length) if i != max_row_index]

            for idx in order:
                seq = new_sequences[idx]

                seq1, seq2 = set_sequences(aligned_center, seq)
                seq2 = project_onto_master(seq1,seq2)
                mat = matrix_building(seq1, seq2, gap_value)
                score_mat = algorithm(
                    seq1,
                    seq2,
                    mat,
                    gap_value,
                    match_value,
                    mismatch_value,
                )
                path = traceback(
                    score_mat,
                    seq1,
                    seq2,
                    gap_value,
                    match_value,
                    mismatch_value,
                )

                new_c, new_s = reconstruct_alignment(seq1, seq2, path)
                aligned_center, realigned_s, aligned_others = merge_alignment(
                    aligned_center, new_c, new_s, aligned_others
                )

            # Final MSA display
            msa = [aligned_center] + aligned_others

            identity_percent = calculate_identity_percentage(msa)
            scoring = calculate_msa_score(msa, match_value, mismatch_value, gap_value)
            st.markdown(f"### MSA Score")
            st.write(f"Total Alignment Score: **{scoring}**")

            st.markdown(f"### Identity Percentage")
            st.write(f"Average Identity: **{identity_percent}%**")

            matches, mismatches, gaps = count_msa_statistics(msa)
            st.markdown("### Alignment Statistics")
            st.write(f"Matches: **{matches}**")
            st.write(f"Mismatches: **{mismatches}**")
            st.write(f"Gaps: **{gaps}**")

            full_text = ""
            full_text += get_text(
                scoring,
                identity_percent,
                matches,
                gaps,
                mismatches,
                msa,
                match_value,
                mismatch_value,
                gap_value
            )

            date = datetime.datetime.now().strftime("%Y-%m-%d")
            filename = f"alignment_{date}.txt"

            new_text = save_to_text_file(filename, full_text)

            # Create a display for all sequences with proper styling
            msa_display = ""
            for idx, row in enumerate(msa, start=1):
                # Join the characters with spaces for better readability
                spaced = " ".join(row)
                msa_display += f"s{idx}: {spaced}\n"

            st.markdown("### Final Multiple Sequence Alignment")
            st.code(msa_display, language=None)

            st.markdown("""
            <style>
                .stCodeBlock {
                    max-width: 100% !important;
                    overflow-x: auto;
                    white-space: pre;
                    font-family: monospace;
                }
            </style>
            """, unsafe_allow_html=True)

            # Show heatmap
            st.markdown("### MSA Heatmap")
            fig = plot_msa_heatmap(msa)
            st.pyplot(fig)

            st.divider()

            col1, col2, col3 = st.columns(3)
            with col2:
                st.download_button(
                    label='Download Results',
                    data=new_text,
                    file_name=filename,
                )