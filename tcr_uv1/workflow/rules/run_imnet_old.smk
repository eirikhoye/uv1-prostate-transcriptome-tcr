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
SAMPLING_SIZE = 1000

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

# Base Path
BASE = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1"

# Compute Clones
COUNTS_VST            = f"{BASE}/results/deseq_all/counts_vst.csv"
META                  = f"{BASE}/workflow/utils/sample_metadata.csv"
CLONE_STATS_ALL_FILE  = f"{BASE}/clone_expansion/all_patients_clone_stats.csv"
CLONE_STATS_EXP_FILE  = f"{BASE}/clone_expansion/expanded_clone_stats.csv"
COMPUTE_EXPANSION_SCRIPT = f"{BASE}/workflow/scripts/tcr_longitudinal_expansion.py"

# Patient IDS:
meta_df = pd.read_csv(META)
PATIENT_ID = meta_df.patient_id.unique()
CLONE_STATS_PATIENT_EXP_FILES = os.path.join(BASE,"clone_expansion/per_patient/{patient_id}_expanded_clone_stats.csv")

#print(expand(CLONE_STATS_PATIENT_EXP_FILES, patient_id=PATIENT_ID))
# Sampled
SAMPLED_EXPANDED = os.path.join(BASE, "sampled_expanded/sampled_{n}_expanded_clone_stats.csv")
SAMPLED_RANDOM   = os.path.join(BASE, "sampled_random/sampled_{n}_random_clone_stats.csv")

# Sampled uncollapsed
SAMPLED_EXPANDED_UNCOLLAPSED = os.path.join(BASE, "sampled_expanded_uncollapsed/sampled_{n}_expanded_clone_stats_uncollapsed.csv")
SAMPLED_RANDOM_UNCOLLAPSED   = os.path.join(BASE, "sampled_random_uncollapsed/sampled_{n}_random_clone_stats_uncollapsed.csv")

# Sampled GRAPHS
SAMPLED_EXP_GRAPH = os.path.join(BASE, "sampled_exp_graph/sampled_{n}_expanded_graph.gpickle")
SAMPLED_RAN_GRAPH = os.path.join(BASE, "sampled_ran_graph/sampled_{n}_random_graph.gpickle")

# Sampled MOTIF GRAPHS
OUTPUT_MOTIF_GRAPH_EXPANDED = os.path.join(BASE, "sampled_exp_graph_motif/sampled_{n}_expanded_motif.gpickle")
OUTPUT_MOTIF_GRAPH_RANDOM   = os.path.join(BASE, "sampled_ran_graph_motif/sampled_{n}_random_motif.gpickle")

# ImNet Paths
INPUT_TCRs            = f"{BASE}/tcrseqs_for_imnet/expanded_clones_raw_seqs.txt"
IMNET_GRAPH           = f"{BASE}/tcrseqs_graphs/expanded_clones.gpickle"
SPARK_PATH            = f"{BASE}/opt/spark"
IMNET_ENV_PATH        = f"{BASE}/imnet_env"
IMNET_SCRIPT          = f"{BASE}/workflow/scripts/run_imnet.py"
CLUSTERS_CSV          = f"{BASE}/results/TCR_clusters.csv"
MOTIFS_CSV            = f"{BASE}/results/TCR_motifs.csv"
OUTPUT_MOTIF_GRAPH    = f"{BASE}/tcrseqs_graphs/motif_network.gexf"
BUILD_MOTIF_SCRIPT    = f"{BASE}/workflow/scripts/aggregate_to_motif_graph.py"
UNCOLLAPSE_SCRIPT     = f"{BASE}/workflow/scripts/uncollapse_motifs_to_seqs.py"


# ---------------------------------------------------------
# Rules
# ---------------------------------------------------------

rule all:
    input:
        IMNET_GRAPH,
        OUTPUT_MOTIF_GRAPH,
        CLONE_STATS_ALL_FILE,
        CLONE_STATS_EXP_FILE,
        expand(SAMPLED_EXPANDED, n=list(range(N_SAMPLES))),
        expand(SAMPLED_RANDOM, n=list(range(N_SAMPLES))),
        expand(SAMPLED_EXPANDED_UNCOLLAPSED, n=list(range(N_SAMPLES))),
        expand(SAMPLED_RANDOM_UNCOLLAPSED, n=list(range(N_SAMPLES))),
        expand(SAMPLED_EXP_GRAPH, n=list(range(N_SAMPLES))),
        expand(SAMPLED_RAN_GRAPH, n=list(range(N_SAMPLES))),
        expand(OUTPUT_MOTIF_GRAPH_EXPANDED, n=list(range(N_SAMPLES))),
        expand(OUTPUT_MOTIF_GRAPH_RANDOM, n=list(range(N_SAMPLES))),
        expand(CLONE_STATS_PATIENT_EXP_FILES, patient_id=PATIENT_ID)


rule compute_clone_expansion:
    input:
        vst=COUNTS_VST,
        meta=META
    output:
        exp_path=CLONE_STATS_EXP_FILE,
        ran_path=CLONE_STATS_ALL_FILE,


    params:
        python_script=COMPUTE_EXPANSION_SCRIPT,
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
            --out-exp-all {output.exp_path} \
            --out-ran-all {output.ran_path} \
            --expansion-threshold {params.lfc} \
            --presence-ratio {params.per} \

        """


rule compute_clone_expansion_per_patient:
    input:
        vst=COUNTS_VST,
        meta=META
    output:
        exp_path = CLONE_STATS_PATIENT_EXP_FILES,
        ran_path = CLONE_STATS_PATIENT_RAN_FILES
    
    params:
        python_script=COMPUTE_EXPANSION_SCRIPT,
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
            --out-exp-pat {output.exp_path} \
            --out-ran-pat {output.ran_path} \
            --expansion-threshold {params.lfc} \
            --presence-ratio {params.per}      

        """


rule sample_tcrs:
    input:
        stats_all=CLONE_STATS_ALL_FILE,
        stats_exp=CLONE_STATS_EXP_FILE

    output:
        exp=SAMPLED_EXPANDED,
        ran=SAMPLED_RANDOM

    params:
        n=SAMPLING_SIZE

    shell:
        r"""
        # expanded clones: Keep header, shuffle body
        (head -n 1 {input.stats_exp}; tail -n +2 {input.stats_exp} | shuf -n {params.n}) > {output.exp}

        # All clones
        (head -n 1 {input.stats_all}; tail -n +2 {input.stats_all} | shuf -n {params.n}) > {output.ran}
        """


rule uncollapse_clones:
    input:
        exp=SAMPLED_EXPANDED,
        ran=SAMPLED_RANDOM,
        motifs=MOTIFS_CSV,
        clusters=CLUSTERS_CSV
    output:
        exp=SAMPLED_EXPANDED_UNCOLLAPSED,
        ran=SAMPLED_RANDOM_UNCOLLAPSED
    params:
        script=UNCOLLAPSE_SCRIPT
    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        python {params.script} \
            --input-clones {input.exp} \
            --motifs-csv {input.motifs} \
            --clusters-csv {input.clusters} \
            --out-file {output.exp}

        python {params.script} \
            --input-clones {input.ran} \
            --motifs-csv {input.motifs} \
            --clusters-csv {input.clusters} \
            --out-file {output.ran}
        """


rule run_imnet_sampled:
    input:
        exp=SAMPLED_EXPANDED_UNCOLLAPSED,
        ran=SAMPLED_RANDOM_UNCOLLAPSED

    output:
        exp=SAMPLED_EXP_GRAPH,
        ran=SAMPLED_RAN_GRAPH

    params:
        python_script=IMNET_SCRIPT,

    threads: 4

    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        {SPARK_PATH}/bin/spark-submit --master local[{threads}] {params.python_script} {input.exp} {output.exp}
        {SPARK_PATH}/bin/spark-submit --master local[{threads}] {params.python_script} {input.ran} {output.ran}
        
        """


rule aggregate_to_motif_graph_sampled:
    input:
        exp = SAMPLED_EXP_GRAPH,
        ran = SAMPLED_RAN_GRAPH,
        clusters = CLUSTERS_CSV,
        motifs = MOTIFS_CSV
    output:
        exp = OUTPUT_MOTIF_GRAPH_EXPANDED,
        ran = OUTPUT_MOTIF_GRAPH_RANDOM
    params:
        python_script = BUILD_MOTIF_SCRIPT
    threads: 1

    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        python {params.python_script} \
            --graph {input.exp} \
            --clusters {input.clusters} \
            --motifs {input.motifs} \
            --out {output.exp}

        python {params.python_script} \
            --graph {input.ran} \
            --clusters {input.clusters} \
            --motifs {input.motifs} \
            --out {output.ran}
        """


rule run_imnet:
    input:
        tcrs=INPUT_TCRs

    output:
        graph=IMNET_GRAPH

    params:
        python_script=IMNET_SCRIPT,

    threads: 64

    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        {SPARK_PATH}/bin/spark-submit --master local[{threads}] {params.python_script} {input.tcrs} {output.graph}
        
        """


rule aggregate_to_motif_graph:
    input:
        graph = IMNET_GRAPH,
        clusters = CLUSTERS_CSV,
        motifs = MOTIFS_CSV
    output:
        motif_graph = OUTPUT_MOTIF_GRAPH
    params:
        python_script = BUILD_MOTIF_SCRIPT
    threads: 1
    shell:
        """
        unset PYTHONPATH
        unset PYTHONHOME

        source {IMNET_ENV_PATH}/bin/activate

        python {params.python_script} \
            --graph {input.graph} \
            --clusters {input.clusters} \
            --motifs {input.motifs} \
            --out {output.motif_graph}
        """