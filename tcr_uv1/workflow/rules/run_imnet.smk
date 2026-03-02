import os
import pandas as pd


# ---------------------------------------------------------
# Params
# ---------------------------------------------------------

# Define expanded clones
LFC_CUT = 2
PERSIST = 0.5

# Sampling clone sizes
N_SAMPLES     = 100
SAMPLING_SIZE = 5000

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

# Base Path
BASE = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1"

# Compute Clones
#COUNTS_VST = f"{BASE}/results/deseq_all/counts_vst.csv"
#META = f"{BASE}/workflow/utils/sample_metadata.csv"
#LONE_STATS_ALL_FILE = f"{BASE}/clone_expansion/all_patients_clone_stats.csv"
#CLONE_STATS_EXP_FILE = f"{BASE}/clone_expansion/expanded_clone_stats.csv"


# ClusTCR Input
COUNTS_VST   = f"{BASE}/results/deseq_all/counts_vst.csv"
META         = f"{BASE}/workflow/utils/sample_metadata.csv"
CLUSTERS_CSV = f"{BASE}/results/TCR_clusters.csv"
MOTIFS_CSV   = f"{BASE}/results/TCR_motifs.csv"

# Merged clonal expansion stats
EXPANDED_CLONES_ALL = f"{BASE}/merged_expansion/clone_expansion/expanded_merged_clone_stats.csv"
RANDOM_CLONES_ALL   = f"{BASE}/merged_expansion/clone_expansion/random_merged_clone_stats.csv"

# Merged uncollapsed TCRseqs
SAMPLED_EXPANDED_UNCOLLAPSED = f"{BASE}/merged_expansion/uncollapsed_tcrs/expanded_merged_clone_stats_uncollapsed.csv"
SAMPLED_RANDOM_UNCOLLAPSED   = f"{BASE}/merged_expansion/uncollapsed_tcrs/random_merged_clone_stats_uncollapsed.csv"

# Merged sampled TCRseqs
SAMPLED_EXPANDED = os.path.join(BASE, "merged_expansion/sampled_expanded_tcrs/sampled_{n}_expanded_clone_stats_uncollapsed.csv")
SAMPLED_RANDOM   = os.path.join(BASE, "merged_expansion/sampled_random_tcrs/sampled_{n}_random_clone_stats_uncollapsed.csv")

# Merged sampled TCRseq graphs
SAMPLED_EXPANDED_GRAPH = os.path.join(BASE, "merged_expansion/sampled_expanded_graphs/sampled_{n}_expanded_graph.gpickle")
SAMPLED_RANDOM_GRAPH   = os.path.join(BASE, "merged_expansion/sampled_random_graphs/sampled_{n}_random_graph.gpickle")


# Per patient clonal expansion stats


# SCRIPT PATHS
COMPUTE_EXPANSION = f"{BASE}/workflow/scripts/tcr_longitudinal_expansion.py"
UNCOLLAPSE_SCRIPT = f"{BASE}/workflow/scripts/uncollapse_motifs_to_seqs.py"
IMNET_SCRIPT      = f"{BASE}/workflow/scripts/run_imnet.py"


# ImNet Paths
IMNET_ENV_PATH = f"{BASE}/imnet_env"
SPARK_PATH     = f"{BASE}/opt/spark"

# ---------------------------------------------------------
# Rules
# ---------------------------------------------------------
rule all:
    input:
        EXPANDED_CLONES_ALL,
        RANDOM_CLONES_ALL,
        SAMPLED_EXPANDED_UNCOLLAPSED,
        SAMPLED_RANDOM_UNCOLLAPSED,
        expand(SAMPLED_EXPANDED, n=list(range(N_SAMPLES))),
        expand(SAMPLED_RANDOM, n=list(range(N_SAMPLES))),
        expand(SAMPLED_EXPANDED_GRAPH, n=list(range(N_SAMPLES))),
        expand(SAMPLED_RANDOM_GRAPH, n=list(range(N_SAMPLES)))


rule compute_clone_expansion:
    input:
        vst=COUNTS_VST,
        meta=META

    output:
        expanded=EXPANDED_CLONES_ALL,
        random=RANDOM_CLONES_ALL
    
    params:
        python_script=COMPUTE_EXPANSION,
        lfc=LFC_CUT,
        per=PERSIST
    
    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate
        
        python {params.python_script} \
            --vst {input.vst} \
            --meta {input.meta} \
            --out-exp-all {output.expanded} \
            --out-ran-all {output.random} \
            --expansion-threshold {params.lfc} \
            --presence-ratio {params.per}
        """


rule uncollapse_clustcr_motifs_to_tcr_seqs:
    input:
        expanded=EXPANDED_CLONES_ALL,
        random=RANDOM_CLONES_ALL,
        motifs=MOTIFS_CSV,
        clusters=CLUSTERS_CSV

    output:
        expanded=SAMPLED_EXPANDED_UNCOLLAPSED,
        random=SAMPLED_RANDOM_UNCOLLAPSED

    params:
        script=UNCOLLAPSE_SCRIPT
    
    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        python {params.script} \
            --input-clones {input.expanded} \
            --motifs-csv {input.motifs} \
            --clusters-csv {input.clusters} \
            --out-file {output.expanded}
        
        python {params.script} \
            --input-clones {input.random} \
            --motifs-csv {input.motifs} \
            --clusters-csv {input.clusters} \
            --out-file {output.random}
        """

rule sample_tcrs:
    input:
        expanded=SAMPLED_EXPANDED_UNCOLLAPSED,
        random=SAMPLED_RANDOM_UNCOLLAPSED

    output:
        expanded=SAMPLED_EXPANDED,
        random=SAMPLED_RANDOM

    params:
        n=SAMPLING_SIZE
    
    shell:
        r"""
        # expanded clones: Keep header, shuffle body
        (head -n 1 {input.expanded}; tail -n +2 {input.expanded} | shuf -n {params.n}) > {output.expanded}

        # All clones
        (head -n 1 {input.random}; tail -n +2 {input.random} | shuf -n {params.n}) > {output.random}
        """


rule run_imnet_sampled:
    input:
        expanded=SAMPLED_EXPANDED,
        random=SAMPLED_RANDOM

    output:
        expanded=SAMPLED_EXPANDED_GRAPH,
        random=SAMPLED_RANDOM_GRAPH

    params:
        python_script=IMNET_SCRIPT
    
    threads: 4

    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        {SPARK_PATH}/bin/spark-submit --master local[{threads}] {params.python_script} {input.expanded} {output.expanded}
        {SPARK_PATH}/bin/spark-submit --master local[{threads}] {params.python_script} {input.random} {output.random}

        """

