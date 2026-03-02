PATH_MIXCR = config['path_mixcr']

rule mixcr_align:
    input:
        fastq=FASTQ_FILES,
        #fasta = FASTA_FILES
    
    output:
        vdjca=VDJCA_FILES,
        report=REPORT_FILES
    
    params:
        species="hs", # human
        preserve_reads="true"
    
    threads: 16
    shell:
        """
        mixcr align -f -s {params.species} -O saveOriginalReads={params.preserve_reads} {input.fastq} {output.vdjca} --report {output.report}
        """

rule mixcr_assemble:
    input:
        vdjca=VDJCA_FILES
    
    output:
        clna=ASSEMBLE_FILES,
        report=ASSEMBLE_REPORT

    threads: 16
    
    shell:
        """
        mixcr assemble -f --write-alignments {input.vdjca} {output.clna} --report {output.report}
        """

rule export_clones:
    input:
        clna=ASSEMBLE_FILES
    
    output:
        clones=CLONE_FILES,
    
    threads: 16

    shell:
        """
        mixcr exportClones -f {input.clna} {output.clones}
        """

rule export_alignments:
    input:
        clna=ASSEMBLE_FILES
    
    output:
        alignment=ALIGNMENT_FILES
    
    threads: 16

    shell:
        """
        mixcr exportAlignments -f -cloneIdWithMappingType -cloneId -readIds -descrsR1 {input.clna} {output.alignment}
        """
rule spatial_counts:
    input:
        alignment=ALIGNMENT_FILES,
        clones=CLONE_FILES
    
    output:
        csv=SPATIAL_FILES
    
    shell:
        "python workflow/scripts/mixcr_spatial_counts.py -a {input.alignment} -c {input.clones} -o {output.csv}"

# TODO make this output different files per sample
rule plot_spatial:
    input:
        csv=SPATIAL_FILES
    output:
        tb=TB_CELL_PLOT,
    conda:
        WORKFLOW_ENV
    shell:
        "python workflow/scripts/plot_spatial_counts.py -i {input.csv} -o {output.tb}"
