# Add this rule to the base dir snakefile, and merge the rule all with the rest of them so there is just one.


rule all:
    input:
        expand("results/{patient}_clone_stats.csv", patient=["P1","P2","P3"]),
        "results/all_patients_clone_stats.csv"

rule find_expanded_TCR_motifs:
    input:
        vst="/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_all/counts_vst.csv",
        meta="/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata.csv"
    output:
        directory("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/")
    script:
        "scripts/find_expanded_motifs.py --vst {input.vst} --meta {input.meta}" # Need to define the file outputs here, and add a params with the path to output dir
