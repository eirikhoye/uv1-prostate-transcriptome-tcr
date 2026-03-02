import pandas as pd

# Load files
count_df = pd.read_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/data/merged_counts.csv", index_col=0)  # Rows: sequences, Columns: samples
annotation_df = pd.read_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/TCR_clusters.csv")      # Columns: 'junction_aa', 'cluster'
motif_df = pd.read_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/TCR_motifs.csv", index_col=0)  # Index: cluster, Columns: 'size', 'motif'

# Merge counts with annotation
# Step 1: Add cluster to each sequence in count_df
count_with_cluster = count_df.merge(annotation_df, left_index=True, right_on='junction_aa')

# Step 2: Add motif to each sequence via cluster
count_with_cluster_motif = count_with_cluster.merge(motif_df[['motif']], left_on='cluster', right_index=True)

# Step 3: Group by motif and sum the counts across samples
motif_count_df = count_with_cluster_motif.groupby('motif')[count_df.columns].sum()

# Optional: Sort motifs by total count across all samples
motif_count_df['total'] = motif_count_df.sum(axis=1)
motif_count_df = motif_count_df.sort_values('total', ascending=False).drop(columns='total')

# Save output
motif_count_df.to_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts.csv")
