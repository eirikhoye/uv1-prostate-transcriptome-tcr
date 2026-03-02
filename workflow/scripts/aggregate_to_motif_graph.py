#!/usr/bin/env python3

import argparse
import pickle
import networkx as nx
import pandas as pd
import collections
import numpy as np

# ------------------------------------------------------------
# Function to build motif-level network from sequence graph
# ------------------------------------------------------------
def build_motif_graph(graph_file, cluster_file, motif_file, out_file):
    # ------------------------------------------------------------
    # Load sequence-level graph (full TCR similarity network)
    # ------------------------------------------------------------
    with open(graph_file, 'rb') as f:
        G = pickle.load(f)

    # ------------------------------------------------------------
    # Load cluster and motif annotations
    # ------------------------------------------------------------
    df_clusters = pd.read_csv(cluster_file)
    df_motifs = (
        pd.read_csv(motif_file)
        .rename(columns={"Unnamed: 0": "cluster"})
    )

    # Merge cluster → motif information
    cluster_map = df_clusters.merge(df_motifs, on='cluster')

    # Build dict: sequence (junction_aa) → motif
    seq_to_motif = dict(zip(cluster_map.junction_aa, cluster_map.motif))

    # ------------------------------------------------------------
    # Count edges between motifs based on edges between sequences
    # ------------------------------------------------------------
    motif_edge_counts = collections.Counter()

    for u, v in G.edges():
        m1 = seq_to_motif.get(u)
        m2 = seq_to_motif.get(v)

        # Skip sequences without motif annotation    
        if m1 is None or m2 is None:
            continue

        # Sort motifs so edges are undirected and counted once
        key = tuple(sorted((m1, m2)))
        motif_edge_counts[key] += 1

    # ------------------------------------------------------------
    # Build normalized motif–motif matrix
    # ------------------------------------------------------------
    motif_list = df_motifs.motif.tolist()
    motif_index = {m: i for i, m in enumerate(motif_list)}
    sizes = dict(zip(df_motifs['motif'], df_motifs['size']))

    N = len(motif_list)
    norm_matrix = np.zeros((N, N))

    for (m1, m2), count in motif_edge_counts.items():
        i = motif_index[m1]
        j = motif_index[m2]

        # Normalize by motif sizes to adjust for cluster size imbalance    
        norm = count / (sizes[m1] * sizes[m2])

        norm_matrix[i, j] = norm
        norm_matrix[j, i] = norm  # symmetry

    # ------------------------------------------------------------
    # Build motif-level graph with raw and normalized edge weights
    # ------------------------------------------------------------
    H = nx.Graph()

    for (m1, m2), count in motif_edge_counts.items():
        H.add_edge(
            m1, 
            m2, 
            weight=count,                                # raw edge count
            normalized=count / (sizes[m1] * sizes[m2])  # size-normalized connectivity
        )

    # ------------------------------------------------------------
    # Add node attributes to motif graph
    # ------------------------------------------------------------
    for n in H.nodes:
        # 'size': number of sequences in the motif cluster
        H.nodes[n]['size'] = H.nodes[n].get('size', 1)

        # 'degree': number of connections to other motifs
        # Useful for initial color/size mapping in Gephi
        H.nodes[n]['degree'] = H.degree(n)

    # ------------------------------------------------------------
    # Add or reinforce edge attributes
    # ------------------------------------------------------------
    for u, v in H.edges:
        # 'weight': raw count of edges between motifs    
        H.edges[u, v]['weight'] = H.edges[u, v].get('weight', 1)

        # 'normalized': edge count normalized by motif sizes    
        H.edges[u, v]['normalized'] = H.edges[u, v].get('normalized', 0.01)

    # ------------------------------------------------------------
    # Export motif graph for visualization
    # ------------------------------------------------------------
    nx.write_gexf(H, out_file)
    print(f"Motif network exported to {out_file}")

# TODO Possible confounder, how "degenerate" are the "motifs" i.e. a very 

# TODO Possible simplification is to uncollapse all motifs prior to sampling, then create levenshtein distance graphs to do connectivity analysis

# TODO possible create sequence logos based comparing expanded seuqence against all sequences, and compute empirical p-value



# ------------------------------------------------------------
# Command-line interface
# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Build motif-level graph from TCR sequence network")
    parser.add_argument("--graph", required=True, help="Input sequence-level graph (.gpickle)")
    parser.add_argument("--clusters", required=True, help="CSV file with sequence → cluster mapping")
    parser.add_argument("--motifs", required=True, help="CSV file with motif information")
    parser.add_argument("--out", required=True, help="Output motif-level graph (.gexf)")

    args = parser.parse_args()

    build_motif_graph(
        graph_file=args.graph,
        cluster_file=args.clusters,
        motif_file=args.motifs,
        out_file=args.out
    )


if __name__ == "__main__":
    main()
