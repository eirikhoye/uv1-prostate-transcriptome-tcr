import pandas as pd


# Output files for selected modules
# -----------------------------------------------------
def get_final_output(ENABLED_MODULES):
    outputs = []

    if "ProcessPacBioReads" in ENABLED_MODULES:
        outputs += expand(BAM_FILES, sample=SAMPLES),
        outputs += expand(FASTQ_FILES, sample=SAMPLES)
    
    if "mixcr" in ENABLED_MODULES:
        outputs += expand(VDJCA_FILES, sample=SAMPLES)
        outputs += expand(REPORT_FILES, sample=SAMPLES)
        outputs += expand(ASSEMBLE_FILES, sample=SAMPLES)
        outputs += expand(ASSEMBLE_REPORT, sample=SAMPLES)
        outputs += expand(CLONE_FILES, sample=SAMPLES)
        outputs += expand(ALIGNMENT_FILES, sample=SAMPLES)
        outputs += expand(SPATIAL_FILES, sample=SAMPLES)
        outputs += expand(TB_CELL_PLOT, sample=SAMPLES)

    return outputs