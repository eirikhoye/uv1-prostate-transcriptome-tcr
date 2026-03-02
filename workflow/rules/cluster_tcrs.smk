rule cluster_tcrs:
    input: MERGED_COUNTS
    output: 
        clusters = CLUSTERS,
        motifs = MOTIFS
    conda: CONDA_ENV
    params:
        aa_colname = CDR3_COLUMN,
        python_script = CLUSTCR_SCRIPT_PATH
    threads: 64
    
    shell:
        """
        python {params.python_script} --input {input} --clusters {output.clusters} --motifs {output.motifs} --threads {threads} --aa_colname {params.aa_colname}
        """

rule merge_motif_counts:
    input:
        counts=MERGED_COUNTS,
        annotation=CLUSTERS,
        motifs=MOTIFS
    output:
        MOTIF_COUNTS
    params:
        python_script = "workflow/scripts/merge_motif_counts.py"
    shell:
        """
        python {params.python_script} \
            -c {input.counts} \
            -a {input.annotation} \
            -m {input.motifs} \
            -o {output}
        """

