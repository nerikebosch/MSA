import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter

def plot_msa_heatmap(msa):
    """
    Plot a heatmap of MSA matches/mismatches/gaps.
    - Match (to majority): 0
    - Mismatch (from majority): 1
    - Gap: 2
    """
    num_seqs = len(msa)
    length = len(msa[0])
    matrix = np.zeros((num_seqs, length))

    for j in range(length):
        column = [msa[i][j] for i in range(num_seqs)]
        # Exclude gaps to find the most common base
        non_gaps = [c for c in column if c != "-"]
        if non_gaps:
            majority_base, _ = Counter(non_gaps).most_common(1)[0]
        else:
            majority_base = "-"  # fallback

        for i in range(num_seqs):
            char = msa[i][j]
            if char == "-":
                matrix[i][j] = 2  # gap
            elif char == majority_base:
                matrix[i][j] = 0  # match
            else:
                matrix[i][j] = 1  # mismatch

    fig, ax = plt.subplots(figsize=(length * 0.3, num_seqs * 0.6))
    sns.heatmap(
        matrix,
        cmap=sns.color_palette(["#aaffaa", "#ffaaaa", "#ddddff"]),  # match, mismatch, gap
        cbar_kws={
            "ticks": [0.5, 1.5, 2.5],
        },
        linewidths=0.5,
        linecolor="gray"
    )
    ax.set_yticks(np.arange(num_seqs) + 0.5)
    ax.set_yticklabels([f"s{i+1}" for i in range(num_seqs)], rotation=0)
    ax.set_xticks([])  # hide x-axis for clarity
    return fig
