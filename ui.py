from alignment import *
from file_utils import *
import streamlit as st

def app_creation():
    """
        Main function to create and manage the Streamlit web application interface
        for the Needleman-Wunsch sequence alignment algorithm.

        Loads sequences from user input or file uploads, gathers gap/match/mismatch
        parameters, triggers calculation, and displays the alignment results including
        heatmaps and download functionality.

        Uses session state to persist results across Streamlit reruns.
    """


    if 'result_text' not in st.session_state:
        st.session_state.result_text = None
    if 'seq1_text' not in st.session_state:
        st.session_state.seq1_text = ""
    if 'seq2_text' not in st.session_state:
        st.session_state.seq2_text = ""

    st.title('Needleman-Wunsch Algorithm')

    col1, col2 = st.columns(2)

    with col1:
        seq1_file = st.file_uploader("Choose a file for the first sequence", type="fasta", key="seq1_file")
        if seq1_file is not None:
            st.session_state.seq1_text = load_fasta_sequences(seq1_file)

        seq1 = st.text_input("Enter the first sequence", value=st.session_state.seq1_text, key="seq1_input")

    with col2:
        seq2_file = st.file_uploader("Choose a file for the second sequence", type="fasta", key="seq2_file")
        if seq2_file is not None:
            st.session_state.seq2_text = load_fasta_sequences(seq2_file)

        seq2 = st.text_input("Enter the second sequence", value=st.session_state.seq2_text, key="seq2_input")

    seq1, seq2 = set_sequences(seq1, seq2)

    if not seq1 or not seq2:
        st.error("Please provide at least one sequence.")
        st.stop()

    st.divider()

    col1,col2,col3 = st.columns(3)
    with col1:
        gap_value = st.number_input("Gap penalty", value=-2.0, step=1.0, format="%.2f")

    with col2:
        match_value = st.number_input("Match reward", value=1.0, step=1.0, format="%.2f")

    with col3:
        mismatch_value = st.number_input("Mismatch penalty", value=-1.0, step=1.0, format="%.2f")

    if st.button("Calculate", use_container_width=True):
        zero_matrix = matrix_building(seq1, seq2, gap_value)
        # create the matrix
        completed_matrix = algorithm(seq1, seq2, zero_matrix, gap_value, match_value, mismatch_value)
        scoring = completed_matrix[-1,-1]
        # create all the paths for these specific sequences
        all_paths = find_all_paths(completed_matrix, seq1, seq2, gap_value, match_value, mismatch_value)

        st.session_state.completed_matrix = completed_matrix
        st.session_state.all_paths = all_paths
        st.session_state.scoring = scoring
        st.session_state.gap_penalty = gap_value
        st.session_state.match_reward = match_value
        st.session_state.mismatch_reward = mismatch_value
        st.session_state.seq1 = seq1
        st.session_state.seq2 = seq2

    if 'all_paths' in st.session_state and st.session_state.all_paths:
        st.subheader(f"Found {len(st.session_state.all_paths)} Optimal Paths!")

        # list all the found paths
        selected_idx = st.selectbox(
            "Select which path to display",
            list(range(1, len(st.session_state.all_paths) + 1)),
            format_func=lambda x: f"Path {x}"
        )

        selected_path = st.session_state.all_paths[selected_idx - 1]

        aligned_seq1, aligned_seq2, matching_seq = reconstruct_alignment(st.session_state.seq1, st.session_state.seq2, selected_path)

        # display the results and related traceback
        st.text(aligned_seq1)
        st.text(matching_seq)
        st.text(aligned_seq2)

        # draw the heatmap for this selected path
        st.caption(f"Graphical Interface for Path {selected_idx}")
        graph_path = draw_graphical_interface(st.session_state.completed_matrix, selected_path, st.session_state.seq1, st.session_state.seq2)
        st.pyplot(graph_path)

        # create the text for the file
        full_text = f"Sequence 1: {st.session_state.seq1}\nSequence 2: {st.session_state.seq2}\n\n"
        full_text += f"Gap penalty: {st.session_state.gap_penalty}\n\n"
        full_text += f"Match reward: {st.session_state.match_reward}\n\n"
        full_text += f"Mismatch penalty: {st.session_state.mismatch_reward}\n\n"
        full_text += f"Scoring: {st.session_state.scoring}\n\n"

        # calculate the statistics to save in the text file for all the paths
        for idx, path in enumerate(st.session_state.all_paths, start=1):
            aligned_seq1, aligned_seq2, matching_seq = reconstruct_alignment(seq1, seq2, path)

            matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
            gaps = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == '-' or b == '-')
            total = len(aligned_seq1)

            match_percentage = (matches / total) * 100
            gap_percentage = (gaps / total) * 100

            full_text += get_text(
                matching_length=total,
                match_percentage=match_percentage,
                gap_percentage=gap_percentage,
                aligned_seq1=aligned_seq1,
                aligned_seq2=aligned_seq2,
                matching_seq=matching_seq,
                idx=idx
            )

        seq1_short = st.session_state.seq1[:5]  # First 5 characters
        seq2_short = st.session_state.seq2[:5]  # First 5 characters

        filename = f"alignment_{seq1_short}_vs_{seq2_short}.txt"

        # save the file to the program to show how to save file
        save_to_text_file(filename, full_text)
        st.session_state.result_text = full_text


    st.divider()

    # also add a manually saving the file to your computer
    if st.session_state.result_text is not None:
        seq1_short = st.session_state.seq1[:5]  # First 5 characters
        seq2_short = st.session_state.seq2[:5]  # First 5 characters

        filename = f"alignment_{seq1_short}_vs_{seq2_short}.txt"

        col1, col2, col3 = st.columns(3)
        with col2:
            st.download_button(
                label='Download Results',
                data=st.session_state.result_text,
                file_name=filename,
            )