#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Prepare TCR clone data for IMNET by merging with motif and cluster information."
    )

    # Input arguments
    parser.add_argument(
        "--input-clones",
        required=True,
        help="CSV file containing TCR clones for a sample (expanded, random, or other)."
    )
    parser.add_argument(
        "--motifs-csv",
        required=True,
        help="CSV file containing TCR motifs (constant across samples)."
    )
    parser.add_argument(
        "--clusters-csv",
        required=True,
        help="CSV file containing TCR clusters (constant across samples)."
    )

    # Output argument
    parser.add_argument(
        "--out-file",
        required=True,
        help="Path to write the final exploded CSV file for IMNET input."
    )

    args = parser.parse_args()

    input_path = Path(args.input_clones)
    motifs_path = Path(args.motifs_csv)
    clusters_path = Path(args.clusters_csv)
    out_path = Path(args.out_file)

    # -------------------------------
    # Load input data
    # -------------------------------
    clones_df = pd.read_csv(input_path)

    motifs_df = pd.read_csv(motifs_path, sep=",", index_col=0)
    motifs_df['cluster'] = motifs_df.index

    clusters_df = pd.read_csv(clusters_path, sep=",")

    # -------------------------------
    # Merge clones with motifs
    # -------------------------------
    clones_with_motifs = clones_df.merge(
        motifs_df[['motif', 'cluster']],
        how='left',
        left_on='amino_acid',
        right_on='motif'
    ).drop(columns=['motif'])

    # -------------------------------
    # Merge with clusters and prepare final exploded table
    # -------------------------------
    exploded_df = clones_with_motifs.merge(
        clusters_df,
        on='cluster',
        how='left'
    ).drop(columns=['amino_acid']).rename(columns={'junction_aa': 'amino_acid'})

    exploded_df = exploded_df.reset_index(drop=True)

    # -------------------------------
    # Write output
    # -------------------------------
    exploded_df.to_csv(out_path, index=False)

    print(f"Exploded TCR table written to {out_path}")

if __name__ == "__main__":
    main()
