rule split_motif_counts:
    input:
        MOTIF_COUNTS
    output:
        SPLIT_MOTIFS
    params:
        python_script = "workflow/scripts/split_motif_counts.py",
        outdir = SPLIT_MOTIFS_DIR
    shell:
        """
        python {params.python_script} -i {input} -o {params.outdir}
        """


rule run_hill_profiles:
    input:
        SPLIT_MOTIFS
    output:
        MOTIF_HILL_PROFILES
    params:
        r_script = "workflow/scripts/hill_diversity_single.R"
    conda:
        ALAKAZAM_ENV
    shell:
        "Rscript {params.r_script} {input} {output}"
    
rule merge_hill_profiles:
    input:
        expand(MOTIF_HILL_PROFILES, sample=SAMPLE_NAMES)
    output:
        MERGED_MOTIF_HILL
    run:
        import pandas as pd

        dfs = []
        for f in input:
            df = pd.read_csv(f, sep="\t")
            # Optionally add a column to track the sample
            sample_name = os.path.basename(f).replace(".csv", "")
            df["sample"] = sample_name
            dfs.append(df)

        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv(output[0], sep="\t", index=False)
        print(f"Merged {len(dfs)} files into {output[0]}")