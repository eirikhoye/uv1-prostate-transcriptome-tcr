import os
import csv
import pandas


# all
COUNTS       = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts.csv"
UV1_CLINICAL = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata.csv"
CONDA_ENV    = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/envs/deseq2.yaml"

# run_DESeq2
DESEQ_SCRIPT_PATH_ALL = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2.R"
RESULTS_ALL           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_all/motifs_deseq_res.csv"
HEATMAP_ALL           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_all/motifs_heatmap.png"
PCA_ALL               = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_all/motifs_pca.png"
VST_ALL               = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_all/counts_vst.csv"

# run_DESeq2_bpmc_tissue
DESEQ_SCRIPT_PATH_PBMC_TISSUE = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_pbmc_tissue.R"
RESULTS_PBMC_TISSUE           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_pbmc_tissue/motifs_deseq_res.csv"
HEATMAP_PBMC_TISSUE           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_pbmc_tissue/motifs_heatmap.png"
PCA_PBMC_TISSUE               = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_pbmc_tissue/motifs_pca.png"

# merge_motifs_per_patient
MERGED_COUNTS_PER_PATIENT        = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts_per_patient.csv"
UV1_CLINICAL_SUBSET              = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_subset.csv"
DESEQ_SCRIPT_PATH_MERGED_PATIENT = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_merged_patient.R"
RESULTS_MERGED_PATIENT           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_merged_patient/motifs_deseq_res.csv"
HEATMAP_MERGED_PATIENT           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_merged_patient/motifs_heatmap.png"
PCA_MERGED_PATIENT               = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_merged_patient/motifs_pca.png"

# select_baseline_sample
BASELINE_COUNTS_PER_PATIENT        = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/baseline_counts_per_patient.csv"
UV1_CLINICAL_BASELINE              = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_baseline.csv"
DESEQ_SCRIPT_PATH_BASELINE_PATIENT = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_baseline_patient.R"
RESULTS_BASELINE_PATIENT           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_baseline_patient/motifs_deseq_res.csv"
HEATMAP_BASELINE_PATIENT           = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_baseline_patient/motifs_heatmap.png"
PCA_BASELINE_PATIENT               = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_baseline_patient/motifs_pca.png"

# run_DESeq2_biopsy_prior
COUNTS_BIOPSIES_PRIOR       = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts_biopsies_prior.csv"
UV1_CLINICAL_BIOPSIES_PRIOR = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_biopsies_prior.csv"
DESEQ_SCRIPT_PATH_BIOPSIES  = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_biopsy_write_counts_vst.R"
RESULTS_BIOPSIES_PRIOR      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior/motifs_deseq_res.csv"
HEATMAP_BIOPSIES_PRIOR      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior/motifs_heatmap.png"
PCA_BIOPSIES_PRIOR          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior/motifs_pca.png"
DDS_BIOPSIES_PRIOR          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior/counts_dds_biopsy.csv"
VST_BIOPSIES_PRIOR          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior/counts_vst_biopsy.csv"

# run_DESeq2_biopsy_prior_with_TEX
COUNTS_BIOPSIES_PRIOR_TEX       = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts_biopsies_prior_TEX.csv"
UV1_CLINICAL_BIOPSIES_PRIOR_TEX = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_biopsies_prior_TEX.csv"
TEX_SIGS                        = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/metadata/signature_scores/Activity_scores_250825.xlsx"
TEX_SAMPLE                      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/metadata/signature_scores/Samplesheet.xlsx"
DESEQ_SCRIPT_PATH_BIOPSIES_TEX  = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_biopsy_write_counts_vst_add_TEX.R"
RESULTS_BIOPSIES_PRIOR_TEX      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/motifs_deseq_res.csv"
HEATMAP_BIOPSIES_PRIOR_TEX      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/motifs_heatmap.png"
HEATMAP_BIOPSIES_PRIOR_POSTER_TEX = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/motifs_heatmap_poster.png"
HEATMAP_BIOPSIES_PRIOR_MANUSCRIPT_TEX = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/motifs_heatmap_manuscript.pdf"
PCA_BIOPSIES_PRIOR_TEX          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/motifs_pca.png"
DDS_BIOPSIES_PRIOR_TEX          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/counts_dds_biopsy.csv"
VST_BIOPSIES_PRIOR_TEX          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_prior_TEX/counts_vst_biopsy.csv"


# run_DESeq2_biopsy_post
COUNTS_BIOPSIES_POST       = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts_biopsies_post.csv"
UV1_CLINICAL_BIOPSIES_POST = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_biopsies_post.csv"
DESEQ_SCRIPT_PATH_BIOPSIES = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_biopsy_write_counts_vst.R"
RESULTS_BIOPSIES_POST      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_post/motifs_deseq_res.csv"
HEATMAP_BIOPSIES_POST      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_post/motifs_heatmap.png"
PCA_BIOPSIES_POST          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_post/motifs_pca.png"
DDS_BIOPSIES_POST          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_post/counts_dds_biopsy.csv"
VST_BIOPSIES_POST          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_post/counts_vst_biopsy.csv"

# run_DESeq2_early_pre_post_biopsy
COUNTS_BIOPSIES_EARLY       = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts_biopsies_early.csv"
UV1_CLINICAL_BIOPSIES_EARLY = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata_biopsies_early.csv"
DESEQ_SCRIPT_PATH_EARLY     = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/scripts/run_deseq2_early_pre_vs_post.R"
RESULTS_BIOPSIES_EARLY      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_early/motifs_deseq_res.csv"
HEATMAP_BIOPSIES_EARLY      = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_early/motifs_heatmap.png"
PCA_BIOPSIES_EARLY          = "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/deseq_biopsy_early/motifs_pca.png"


rule all:
    input:
        RESULTS_ALL,
        HEATMAP_ALL,
        PCA_ALL,
        VST_ALL,
        RESULTS_PBMC_TISSUE,
        HEATMAP_PBMC_TISSUE,
        PCA_PBMC_TISSUE,
        MERGED_COUNTS_PER_PATIENT,
        UV1_CLINICAL_SUBSET,
        RESULTS_MERGED_PATIENT,
        HEATMAP_MERGED_PATIENT,
        PCA_MERGED_PATIENT,
        BASELINE_COUNTS_PER_PATIENT,
        UV1_CLINICAL_BASELINE,
        RESULTS_BASELINE_PATIENT,
        HEATMAP_BASELINE_PATIENT,
        PCA_BASELINE_PATIENT,
        COUNTS_BIOPSIES_PRIOR,
        UV1_CLINICAL_BIOPSIES_PRIOR,
        RESULTS_BIOPSIES_PRIOR,
        HEATMAP_BIOPSIES_PRIOR,
        PCA_BIOPSIES_PRIOR,
        COUNTS_BIOPSIES_POST,
        UV1_CLINICAL_BIOPSIES_POST,
        RESULTS_BIOPSIES_POST,
        HEATMAP_BIOPSIES_POST,
        PCA_BIOPSIES_POST,
        COUNTS_BIOPSIES_EARLY,
        DDS_BIOPSIES_PRIOR,
        VST_BIOPSIES_PRIOR,
        DDS_BIOPSIES_POST,
        VST_BIOPSIES_POST,
        RESULTS_BIOPSIES_PRIOR_TEX,
        HEATMAP_BIOPSIES_PRIOR_TEX,
        PCA_BIOPSIES_PRIOR_TEX,
        DDS_BIOPSIES_PRIOR_TEX,
        VST_BIOPSIES_PRIOR_TEX




rule run_DESeq2:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    output: 
        results = RESULTS_ALL,
        heatmap = HEATMAP_ALL,
        pca = PCA_ALL,
        vst = VST_ALL
    params: 
        script_path = DESEQ_SCRIPT_PATH_ALL
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule run_DESeq2_bpmc_tissue:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    output: 
        results = RESULTS_PBMC_TISSUE,
        heatmap = HEATMAP_PBMC_TISSUE,
        pca = PCA_PBMC_TISSUE
    params: 
        script_path = DESEQ_SCRIPT_PATH_PBMC_TISSUE
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule merge_motifs_per_patient:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    
    output:
        counts_sub = MERGED_COUNTS_PER_PATIENT,
        clinical_sub = UV1_CLINICAL_SUBSET
    
    run:
        import pandas as pd
        import re

        df = pd.read_csv(input.counts)
        clinical = pd.read_csv(input.clinical)

        # Remove 'Biopsy' and 'TCC' columns, map remaining columns to patient IDs using regex
        data_cols = [col for col in df.columns if 'Biopsy' not in col and 'TCC' not in col and col not in ['motif', 'sharing_level']]
        patient_map = {col: re.search(r'patient_([A-Za-z0-9]+)', col).group(1) for col in data_cols if 'patient_' in col}
        
        # Create new DataFrame with just the 'motif' column
        df_patients = df[['motif']].copy()

        # Sum values per patient
        for patient in set(patient_map.values()):
            cols_for_patient = [col for col, pat in patient_map.items() if pat == patient]
            df_patients[patient] = df[cols_for_patient].sum(axis=1)
    
        # Subset clinical data to only those patients
        clinical_subset = clinical[clinical['sample_type'].isin(['PBMC', 'Bulk'])][[
        'patient_id', 'PFS', 'PFS_time', 'OS', 'OS_time', 
        'cohort', 'Response'
        ]].drop_duplicates().sort_values(by='patient_id').reset_index(drop=True)

        # Make the patient id same as counts df
        clinical_subset['patient_id'] = 'Pat' + clinical_subset['patient_id'].astype(str)

        # Write subset counts and clinical data to file
        df_patients.to_csv(output.counts_sub, index=False)
        clinical_subset.to_csv(output.clinical_sub, index=False)


rule run_DESeq2_merged_patients:
    input:
        counts = MERGED_COUNTS_PER_PATIENT,
        clinical = UV1_CLINICAL_SUBSET
    output: 
        results = RESULTS_MERGED_PATIENT,
        heatmap = HEATMAP_MERGED_PATIENT,
        pca = PCA_MERGED_PATIENT
    params: 
        script_path = DESEQ_SCRIPT_PATH_MERGED_PATIENT
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule select_baseline_sample:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    
    output:
        counts_sub = BASELINE_COUNTS_PER_PATIENT,
        clinical_sub = UV1_CLINICAL_BASELINE
    
    run:
        import pandas as pd

        df = pd.read_csv(input.counts)
        clinical = pd.read_csv(input.clinical)

        non_sample_cols = ['motif', 'sharing_level']

        def extract_day(col):
            try:
                return int(col.split('_')[0])
            except ValueError:
                return None

        baseline_cols = [col for col in df.columns
                        if col not in non_sample_cols
                        and extract_day(col) is not None
                        and extract_day(col) <= 0
                        ]

        columns_to_keep = ['motif'] + baseline_cols
        
        df_baseline = df[columns_to_keep]
        clinical_baseline = clinical[clinical['sample'].isin(baseline_cols)].reset_index(drop=True)

        # Write subset counts and clinical data to file
        df_baseline.to_csv(output.counts_sub, index=False)
        clinical_baseline.to_csv(output.clinical_sub, index=False)


rule run_DESeq2_baseline_patients:
    input:
        counts = BASELINE_COUNTS_PER_PATIENT,
        clinical = UV1_CLINICAL_BASELINE
    output: 
        results = RESULTS_BASELINE_PATIENT,
        heatmap = HEATMAP_BASELINE_PATIENT,
        pca = PCA_BASELINE_PATIENT
    params: 
        script_path = DESEQ_SCRIPT_PATH_BASELINE_PATIENT
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule select_biopsies_prior:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    
    output:
        counts_sub = COUNTS_BIOPSIES_PRIOR,
        clinical_sub = UV1_CLINICAL_BIOPSIES_PRIOR
    
    run:
        import pandas as pd
        import re

        # Load data        
        counts_df = pd.read_csv(input.counts)
        clinical_df = pd.read_csv(input.clinical)

        # Define columns to exclude from biopsy selection
        biopsy_cols = [col for col in counts_df.columns if re.match(r"^-\d+_Biopsy", col)]
        selected_cols = ['motif'] + biopsy_cols

        counts_biopsy = counts_df[selected_cols]
        clinical_biopsy = clinical_df[clinical_df['sample'].isin(biopsy_cols)].reset_index(drop=True)

        # Write output files
        counts_biopsy.to_csv(output.counts_sub, index=False)
        clinical_biopsy.to_csv(output.clinical_sub, index=False)


rule run_DESeq2_biopsy_prior:
    input:
        counts   = COUNTS_BIOPSIES_PRIOR,
        clinical = UV1_CLINICAL_BIOPSIES_PRIOR

    output:
        results     = RESULTS_BIOPSIES_PRIOR,
        heatmap     = HEATMAP_BIOPSIES_PRIOR,
        pca         = PCA_BIOPSIES_PRIOR,
        norm_counts = DDS_BIOPSIES_PRIOR,
        vst_counts  = VST_BIOPSIES_PRIOR

    params:
        script_path = DESEQ_SCRIPT_PATH_BIOPSIES
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule select_biopsies_prior_tex:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL,
        signature = TEX_SIGS,
        sig_sample = TEX_SAMPLE
    output:
        counts_sub = COUNTS_BIOPSIES_PRIOR_TEX,
        clinical_sub = UV1_CLINICAL_BIOPSIES_PRIOR_TEX
    conda: "/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/envs/make_tex_file.yaml"
    shell:
        "python workflow/scripts/select_biopsies_prior_tex.py "
        "--counts {input.counts} "
        "--clinical {input.clinical} "
        "--signature {input.signature} "
        "--sig_sample {input.sig_sample} "
        "--counts_out {output.counts_sub} "
        "--clinical_out {output.clinical_sub}"


rule run_DESeq2_biopsy_prior_tex:
    input:
        counts   = COUNTS_BIOPSIES_PRIOR_TEX,
        clinical = UV1_CLINICAL_BIOPSIES_PRIOR_TEX

    output:
        results            = RESULTS_BIOPSIES_PRIOR_TEX,
        heatmap            = HEATMAP_BIOPSIES_PRIOR_TEX,
        heatmap_poster     = HEATMAP_BIOPSIES_PRIOR_POSTER_TEX,
        heatmap_manuscript = HEATMAP_BIOPSIES_PRIOR_MANUSCRIPT_TEX,
        pca                = PCA_BIOPSIES_PRIOR_TEX,
        norm_counts        = DDS_BIOPSIES_PRIOR_TEX,
        vst_counts         = VST_BIOPSIES_PRIOR_TEX

    params:
        script_path = DESEQ_SCRIPT_PATH_BIOPSIES_TEX
    conda: CONDA_ENV
    script:
        "{params.script_path}"


rule select_biopsies_post:
    input:
        counts = COUNTS,
        clinical = UV1_CLINICAL
    
    output:
        counts_sub = COUNTS_BIOPSIES_POST,
        clinical_sub = UV1_CLINICAL_BIOPSIES_POST
    
    run:
        import pandas as pd
        import re

        # Load data        
        counts_df = pd.read_csv(input.counts)
        clinical_df = pd.read_csv(input.clinical)

        # Define columns to exclude from biopsy selection
        nonneg_biopsy_cols = [col for col in counts_df.columns if re.match(r"^\d+_Biopsy", col, re.IGNORECASE)]
        selected_cols = ['motif'] + nonneg_biopsy_cols

        counts_biopsy = counts_df[selected_cols]
        clinical_biopsy = clinical_df[clinical_df['sample'].isin(nonneg_biopsy_cols)].reset_index(drop=True)

        # Write output files
        counts_biopsy.to_csv(output.counts_sub, index=False)
        clinical_biopsy.to_csv(output.clinical_sub, index=False)


rule run_DESeq2_biopsy_post:
    input:
        counts = COUNTS_BIOPSIES_POST,
        clinical = UV1_CLINICAL_BIOPSIES_POST

    output:
        results     = RESULTS_BIOPSIES_POST,
        heatmap     = HEATMAP_BIOPSIES_POST,
        pca         = PCA_BIOPSIES_POST,
        norm_counts = DDS_BIOPSIES_POST,
        vst_counts  = VST_BIOPSIES_POST

    params:
        script_path = DESEQ_SCRIPT_PATH_BIOPSIES
    conda: CONDA_ENV
    script:
        "{params.script_path}"






# TODO make it filter for biopsies, and only early responders
#rule select_biopsies_responders:
#    input:
#        counts = COUNTS,
#        clinical = UV1_CLINICAL
#    output:
#        counts_sub = COUNTS_BIOPSIES_EARLY,
#        clinical_sub = UV1_CLINICAL_BIOPSIES_EARLY
#    run:
#        import pandas as pd
#        import re

        # Load data
#        counts_df = pd.read_csv(input.counts)
#        clinical_df = pd.read_csv(input.clinical)

        # Define columns to exclude from biopsy selection
#        nonneg_biopsy_cols = [col for col in counts_df.columns if re.match(r"^\d+_Biopsy", col, re.IGNORECASE)]
#        selected_cols = ['motif'] + nonneg_biopsy_cols

#        counts_biopsy = counts_df[selected_cols]
#        clinical_biopsy = clinical_df[clinical_df['sample'].isin(nonneg_biopsy_cols)].reset_index(drop=True)

        # Write output files
#        counts_biopsy.to_csv(output.counts_sub, index=False)
#        clinical_biopsy.to_csv(output.clinical_sub, index=False)

# TODO Need to modify the linked script so it runs the correct diffexp

#rule run_DESeq2_biopsy_responders:
#    input:
#        counts = COUNTS_BIOPSIES_EARLY,
#        clinical = UV1_CLINICAL_BIOPSIES_EARLY
#
#    output:
#        results = RESULTS_BIOPSIES_EARLY,
#        heatmap = HEATMAP_BIOPSIES_EARLY,
#        pca = PCA_BIOPSIES_EARLY
#
#    params:
#        script_path = DESEQ_SCRIPT_PATH_EARLY
#    conda: CONDA_ENV
#    script:
#        "{params.script_path}"

