#!/usr/bin/env python

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------------------------

def check_patient_samples(patient_adata):
    """Ensure exactly one baseline and ≥1 follow-up sample."""
    baseline_mask = patient_adata.obs["day_after_treatment"] == 0
    n_baseline = baseline_mask.sum()
    n_followup = (~baseline_mask).sum()

    if n_baseline != 1:
        raise ValueError(
            f"Patient {patient_adata.obs['patient_id'].iloc[0]} "
            f"has {n_baseline} baseline samples — expected exactly 1."
        )
    if n_followup < 1:
        raise ValueError(
            f"Patient {patient_adata.obs['patient_id'].iloc[0]} has 0 follow-up samples."
        )
    return baseline_mask


def compute_clone_stats(patient_adata):
    """Compute delta, rank, expanded, and presence_score for a patient."""
    baseline_mask = check_patient_samples(patient_adata)
    baseline_vals = patient_adata.X[baseline_mask].A[0] if hasattr(patient_adata.X, "A") else patient_adata.X[baseline_mask][0]

    # first follow-up
    followups = patient_adata.obs[~baseline_mask]
    first_day = followups["day_after_treatment"].min()
    followup_vals = patient_adata.X[patient_adata.obs["day_after_treatment"] == first_day]
    followup_vals = followup_vals.A[0] if hasattr(followup_vals, "A") else followup_vals[0]

    delta = followup_vals - baseline_vals

    var_df = patient_adata.var.copy()
    var_df["delta_VST"] = delta
    var_df["rank_delta"] = var_df["delta_VST"].rank(ascending=False)
    var_df["expanded"] = np.where(var_df["delta_VST"] >= 2, "yes", "no")

    # presence_score = number of follow-ups where clone is > 0
    followup_vals_all = patient_adata.X[~baseline_mask]
    followup_vals_all = followup_vals_all.A if hasattr(followup_vals_all, "A") else followup_vals_all
    var_df["presence_score"] = (followup_vals_all > 0).sum(axis=0)

    return var_df


def plot_patient(patient_adata, var_df, outpath):
    """Plot clone trajectories with color-coded expansion."""
    df = pd.DataFrame(
        patient_adata.X.A if hasattr(patient_adata.X, "A") else patient_adata.X,
        columns=patient_adata.var_names,
        index=patient_adata.obs_names
    )
    df["day_after_treatment"] = patient_adata.obs["day_after_treatment"].values

    long = df.melt(id_vars="day_after_treatment", var_name="clone", value_name="VST")\
             .merge(var_df[["expanded", "presence_score"]], left_on="clone", right_index=True, how="left")

    long["color_value"] = np.where(long["expanded"] == "yes", long["presence_score"], -1)
    max_score = var_df["presence_score"].max()
    palette = [(0.7, 0.7, 0.7)] + sns.color_palette("Oranges", n_colors=max_score + 1)
    long["color_index"] = long["color_value"] + 1

    plt.figure(figsize=(13, 6))
    sns.lineplot(
        data=long,
        x="day_after_treatment",
        y="VST",
        hue="color_index",
        units="clone",
        estimator=None,
        palette=palette,
        alpha=0.7,
        lw=1,
        legend=False
    )
    pid = patient_adata.obs["patient_id"].iloc[0]
    plt.title(f"Clone VST trajectories — Patient {pid}")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


# ---------------------------------------------------------------------
# MAIN FUNCTION
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Identify expanded PBMC TCR motifs")
    parser.add_argument("--vst", required=True, help="Path to VST CSV file")
    parser.add_argument("--meta", required=True, help="Path to metadata CSV file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    OUTDIR = Path(args.outdir)
    OUTDIR.mkdir(exist_ok=True, parents=True)

    vst_df = pd.read_csv(args.vst, index_col=0)
    meta = pd.read_csv(args.meta)
    meta["Response"] = meta["Response"].replace({"Late": "Late/No", "No": "Late/No"})
    meta = meta[meta["sample"].isin(vst_df.columns)].set_index("sample").loc[vst_df.columns]

    adata = anndata.AnnData(
        X=vst_df.T,
        obs=meta,
        var=pd.DataFrame(index=vst_df.index)
    )

    pbmc = adata[adata.obs["sample_type"] == "PBMC"]

    all_results = []
    for patient_id in sorted(pbmc.obs["patient_id"].unique()):
        patient_adata = pbmc[pbmc.obs["patient_id"] == patient_id]
        try:
            var_df = compute_clone_stats(patient_adata)
        except ValueError as e:
            print(f"[ERROR] {e}")
            continue

        out_csv = OUTDIR / f"patient_{patient_id}_clone_stats.csv"
        out_fig = OUTDIR / f"patient_{patient_id}_trajectories.png"

        var_df.to_csv(out_csv, index=True)
        plot_patient(patient_adata, var_df, out_fig)

        var_df["patient_id"] = patient_id
        all_results.append(var_df)

    combined = pd.concat(all_results, axis=0)
    combined.to_csv(OUTDIR / "all_patients_clone_stats.csv")
    print("Finished.")


if __name__ == "__main__":
    main()
