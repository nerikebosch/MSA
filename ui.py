from alignment import *
from file_utils import *
import streamlit as st
import numpy as np
import itertools

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
        # center_ungapped = center_seq
        # st.write("center_ungapped:", center_ungapped)
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
            st.write("sequence 1 before reconstruct", seq1)
            st.write("sequence 2 before reconstruct", seq2)


            new_c, new_s = reconstruct_alignment(seq1, seq2, path)
            st.write("sequence 1 after reconstruct", new_c)
            st.write("sequence 2 after reconstruct", new_s)
            aligned_center, realigned_s, aligned_others = merge_alignment(
                aligned_center, new_c, new_s, aligned_others
            )

            #aligned_others.append(realigned_s)
            st.write("── Step merge idx", idx, "──")
            st.write(" merged_center after    : ", aligned_center)
            st.write(" realigned_s            : ", realigned_s)
            st.write(" aligned_others so far  :")
            for j, o in enumerate(aligned_others):
                st.write(f"   seq[{order[j]}]: {o}")

        # Final MSA display
        msa = [aligned_center] + aligned_others

        st.write("=== FINAL ALIGNMENT ===")
        st.write("aligned_center:", aligned_center)
        for i, raw in enumerate(new_sequences):
            st.write(f"s{i + 1} raw      :", raw)
            st.write(f"s{i + 1} aligned  :", msa[i])

        identity_percent = calculate_identity_percentage(msa)

        st.markdown(f"### Identity Percentage")
        st.write(f"Average Identity: **{identity_percent}%**")

        st.markdown("### Final Multiple Sequence Alignment")
        width = len(aligned_center)  # all rows are now this length
        for idx, row in enumerate(msa, start=1):
            spaced = " ".join(row)  # already exactly width columns
            st.text(f"s{idx}: {spaced}")

        st.write("=== FINAL ALIGNMENT ===")
        st.write("aligned_center:", aligned_center)
        for i, raw in enumerate(new_sequences):
            st.write(f"s{i + 1} raw      :", raw)
            st.write(f"s{i + 1} aligned  :", msa[i])