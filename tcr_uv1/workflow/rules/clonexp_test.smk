# edit the input names



rule analyze_pbmc_expansion:
    input:
        motif_file="data/motif_counts.csv",
        metadata_file="data/sample_metadata.csv"
    output:
        summary="results/patient_summary.csv",
        expanded_clonotypes=expand("results/patient_{patient}_expanded_clonotypes.csv", patient=[822, 823, 824])  # adjust patient IDs or generate dynamically
    params:
        output_prefix="results/patient"
    conda:
        "envs/pbmc_expansion.yaml"  # optional: specify environment
    shell:
        """
        python analyze_pbmc_expansion.py \
            --motif_file {input.motif_file} \
            --metadata_file {input.metadata_file} \
            --output_prefix {params.output_prefix} \
            --patients 822 823 824
        """
