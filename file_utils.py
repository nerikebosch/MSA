from io import StringIO
from Bio import SeqIO

# load sequence from fasta file
def load_fasta_sequences(input_file):
    """
        Load and parse sequences from a FASTA file.

        Args:
            input_file (UploadedFile): Uploaded file object from Streamlit's file uploader.

        Returns:
            list[str]: List of sequences as strings.
    """

    fasta_sequences = []

    stringio = StringIO(input_file.getvalue().decode('utf-8'))
    for sequence in SeqIO.parse(stringio, "fasta"):
        fasta_sequences.append(str(sequence.seq))

    return fasta_sequences


def get_text(scoring, identity_percentage, match_amount, gap_amount, mismatch_amount, sequences, match_value, mismatch_value, gap_value):
    """
        Generate a formatted string summarizing sequence alignment results.

        Args:
            scoring (str): Description or method of the scoring system used.
            identity_percentage (float): Percentage of identical positions in the alignment.
            match_amount (int): Number of matches found in the alignment.
            gap_amount (int): Number of gaps in the alignment.
            mismatch_amount (int): Number of mismatches in the alignment.
            sequences (list[str]): List of sequences involved in the alignment.
            match_value (float): Score assigned to a match.
            mismatch_value (float): Score assigned to a mismatch.
            gap_value (float): Score assigned to a gap.

        Returns:
            str: A formatted string containing alignment summary and scoring details.
    """


    sequences_text = "\n".join([f"s{i + 1}: {seq}" for i, seq in enumerate(sequences)])

    text = f"""
    {sequences_text}
    
    Scoring: {scoring}
    Identity Percentage: {identity_percentage:.2f}%
    
    Match value: {match_value:.2f}
    Gap value: {gap_value:.2f}
    Mismatch value: {mismatch_value:.2f}
    
    Number of Matches: {match_amount}
    Number of Gaps: {gap_amount}
    Number of Mismatches: {mismatch_amount}

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