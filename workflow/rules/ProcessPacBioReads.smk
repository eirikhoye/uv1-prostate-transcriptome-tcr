PYTHON_SCRIPT = "workflow/scripts/process_spatial_reads.py "

rule process_pacbio_reads:
    input:
        bamfiles = BAM_FILES

    output:
        fastqfiles = FASTQ_FILES

    params:
        script = PYTHON_SCRIPT

    conda:
        WORKFLOW_ENV

    shell:
        """
        python {params.script} -i {input.bamfiles} -o {output.fastqfiles}
        """

