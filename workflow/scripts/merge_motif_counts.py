#!/usr/bin/env python3
import pandas as pd
import argparse
import os

def main(count_file, annotation_file, motif_file, output_file):
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Load files
    count_df = pd.read_csv(count_file, index_col=0)           # Rows: sequences, Columns: samples
    annotation_df = pd.read_csv(annotation_file)              # Columns: 'junction_aa', 'cluster'
    motif_df = pd.read_csv(motif_file, index_col=0)           # Index: cluster, Columns: 'size', 'motif'

    # Merge counts with annotation
    # Step 1: Add cluster to each sequence in count_df
    count_with_cluster = count_df.merge(annotation_df, left_index=True, right_on='junction_aa')

    # Step 2: Add motif to each sequence via cluster
    count_with_cluster_motif = count_with_cluster.merge(motif_df[['motif']], left_on='cluster', right_index=True)

    # Step 3: Group by motif and sum counts across samples
    motif_count_df = count_with_cluster_motif.groupby('motif')[count_df.columns].sum()

    # Optional: Sort motifs by total count across all samples
    motif_count_df['total'] = motif_count_df.sum(axis=1)
    motif_count_df = motif_count_df.sort_values('total', ascending=False).drop(columns='total')

    # Save output
    motif_count_df.to_csv(output_file)
    print(f"Saved motif counts to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge TCR counts with cluster annotations and motifs, aggregate by motif.")
    parser.add_argument("-c", "--counts", required=True, help="Input counts CSV file")
    parser.add_argument("-a", "--annotation", required=True, help="TCR clusters CSV file")
    parser.add_argument("-m", "--motifs", required=True, help="TCR motifs CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for motif counts")
    args = parser.parse_args()

    main(args.counts, args.annotation, args.motifs, args.output)
