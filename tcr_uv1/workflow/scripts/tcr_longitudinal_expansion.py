#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata
import os

# ---------------------------------------------------------
# CHECKS
# ---------------------------------------------------------

def check_patient_samples(patient_adata):
    """Ensure exactly one baseline and at least one follow-up sample."""
    baseline_mask = patient_adata.obs["day_after_treatment"] == 0
    n_baseline = baseline_mask.sum()
    n_followup = (~baseline_mask).sum()

    if n_baseline != 1:
        raise ValueError(
            f"Patient {patient_adata.obs['patient_id'].iloc[0]} "
            f"has {n_baseline} baseline samples (expected 1)."
        )
    if n_followup < 1:
        raise ValueError(
            f"Patient {patient_adata.obs['patient_id'].iloc[0]} has no follow-up samples."
        )

    return baseline_mask


# ---------------------------------------------------------
# CORE COMPUTATION
# ---------------------------------------------------------

def compute_clone_stats(patient_adata, expansion_threshold):
    """
    Compute clonal expansion metrics for a single patient.
    expansion_threshold = minimum delta_VST to call a clone "expanded".
    """
    baseline_mask = check_patient_samples(patient_adata)

    baseline_vals = (
        patient_adata.X[baseline_mask].A[0]
        if hasattr(patient_adata.X, "A")
        else patient_adata.X[baseline_mask][0]
    )

    followups = patient_adata.obs[~baseline_mask]
    first_day = followups["day_after_treatment"].min()

    followup_vals = (
        patient_adata.X[patient_adata.obs["day_after_treatment"] == first_day].A[0]
        if hasattr(patient_adata.X, "A")
        else patient_adata.X[patient_adata.obs["day_after_treatment"] == first_day][0]
    )

    delta = followup_vals - baseline_vals

    var_df = patient_adata.var.copy()
    var_df["delta_VST"] = delta
    var_df["rank_delta"] = var_df["delta_VST"].rank(ascending=False)
    var_df["expanded"] = np.where(
        var_df["delta_VST"] >= expansion_threshold, "yes", "no"
    )

    # Presence score across all follow-ups
    followup_vals_all = (
        patient_adata.X[~baseline_mask].A
        if hasattr(patient_adata.X, "A")
        else patient_adata.X[~baseline_mask]
    )

    var_df["presence_score"] = (followup_vals_all > 0).sum(axis=0)

    # Compute scalar max
    max_presence = var_df["presence_score"].max()
    var_df["max_presence"] = max_presence

    # Use scalar max to compute proportion
    if max_presence > 0:
        var_df["proportion_presence"] = var_df["presence_score"] / max_presence
    else:
        var_df["proportion_presence"] = 0
    
    return var_df


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute clonal expansion statistics across all patients."
    )

    parser.add_argument("--vst", required=True, help="Path to VST-transformed counts CSV.")
    parser.add_argument("--meta", required=True, help="Sample metadata CSV.")
    parser.add_argument("--out-exp-all", help="Output CSV file for all expanded patients.")
    parser.add_argument("--out-ran-all", help="Output CSV file for all random patients.")
    parser.add_argument("--out-exp-pat", help="Output CSV file for expanded single patient")
    parser.add_argument("--out-ran-pat", help="Output CSV file for random single patient")
    parser.add_argument("--expansion-threshold", type=float, default=2.0,
                        help="Delta_VST threshold for calling clones expanded.")
    parser.add_argument("--presence-ratio", type=float, default=0.5,
                        help="Fraction of follow up PBMC where cloen is detected")
    parser.add_argument("--sample-type", default="PBMC", help="Sample type to subset (default PBMC).")

    args = parser.parse_args()
    outpath = Path(args.out_ran_all)
    out_exp = Path(args.out_exp_all)

    # -------------------------------
    # Load data
    # -------------------------------
    vst_df = pd.read_csv(args.vst, index_col=0)
    meta = pd.read_csv(args.meta)

    meta["Response"] = meta["Response"].replace({"Late": "Late/No", "No": "Late/No"})
    meta = meta[meta["sample"].isin(vst_df.columns)].set_index("sample").loc[vst_df.columns]

    adata = anndata.AnnData(
        X=vst_df.T,
        obs=meta,
        var=pd.DataFrame(index=vst_df.index)
    )

    adata = adata[adata.obs["sample_type"] == args.sample_type]

    all_results = []


    # -------------------------------
    # Per-patient processing
    # -------------------------------

    if args.out_exp_all is not None:
        for patient_id in sorted(adata.obs["patient_id"].unique()):
            patient_adata = adata[adata.obs["patient_id"] == patient_id]

            try:
                var_df = compute_clone_stats(patient_adata, args.expansion_threshold)
            except ValueError as e:
                print(f"[WARNING] {e}")
                continue

            var_df["patient_id"] = patient_id
            all_results.append(var_df)


        # -------------------------------
        # Combine all and write single file
        # -------------------------------
        if len(all_results) > 0:
            combined = pd.concat(all_results, axis=0)
            combined.reset_index().rename(columns={"index": "amino_acid"}).to_csv(outpath, index=False)

            exp_clones = combined[
                (combined['proportion_presence'] > args.presence_ratio) &
                (combined['expanded'] == 'yes')
            ]
            exp_clones.reset_index().rename(columns={"index": "amino_acid"}).to_csv(out_exp, index=False)


        print(f"Finished. Output written to:\n{outpath}\n{out_exp}")
    
    # -------------------------------
    # Individual patient processing
    # -------------------------------

#    if args.out_exp_pat is not None:
#        path = args.out_exp_pat
#        patient_id = os.path.basename(path).split("_", 1)[0]
#        patient_adata = adata[adata.obs["patient_id"] == patient_id]

#        try
        






if __name__ == "__main__":
    main()
