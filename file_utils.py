from io import StringIO
from Bio import SeqIO

# load sequence from fasta file
def load_fasta_sequences(input_file):
    """
        Load and parse a single sequence from a FASTA file.

        Args:
            input_file (UploadedFile): Uploaded file object from Streamlit's file uploader.

        Returns:
            str: The loaded sequence as a string.
    """

    stringio = StringIO(input_file.getvalue().decode('utf-8'))
    fasta_sequences = SeqIO.read(stringio, "fasta")
    fasta_sequences = str(fasta_sequences.seq)

    return fasta_sequences


def get_text(matching_length, match_percentage, gap_percentage, aligned_seq1, aligned_seq2, matching_seq, idx=None):
    """
        Generate a formatted text block summarizing the alignment results.

        Args:
            matching_length (int): Total length of the alignment.
            match_percentage (float): Percentage of identical positions.
            gap_percentage (float): Percentage of gaps in the alignment.
            aligned_seq1 (str): Aligned version of the first sequence.
            aligned_seq2 (str): Aligned version of the second sequence.
            matching_seq (str): String representing match/mismatch/gap positions.
            idx (int, optional): Path index number for labeling multiple paths. Defaults to None.

        Returns:
            str: Formatted text block summarizing alignment results.
    """

    header = f"--- Path {idx} ---\n" if idx is not None else ""
    text = f"""{header}
    {aligned_seq1}
    {matching_seq}
    {aligned_seq2}
    Matching Length: {matching_length}
    Percentage of Identical Positions: {match_percentage:.2f}%
    Percentage of Gaps: {gap_percentage:.2f}%

    """
    return text

def save_to_text_file(filename, text):
    """
        Save the provided text to a local file.

        Args:
            filename (str): Local file name.
            text (str): Text content to save.

        Returns:
            str: The same text that was saved.

    """

    with open(filename, 'w') as output_file:
        output_file.write(text)

    return text