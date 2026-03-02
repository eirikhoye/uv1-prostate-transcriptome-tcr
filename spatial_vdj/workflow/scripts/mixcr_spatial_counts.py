#!/usr/bin/env python3
import pandas as pd
import re
import argparse

def parse_coords(descr):
    """Extract x,y coordinates from descrsR1 string."""
    match = re.search(r"x_coor=([\d\.]+)_y_coor=([\d\.]+)", descr)
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None

def classify_receptor(vhits):
    """Classify receptor type based on V gene."""
    if pd.isna(vhits):
        return "Other"
    genes = [g.split("*")[0] for g in vhits.split(",")]
    if any(g.startswith(("IGH", "IGK", "IGL")) for g in genes):
        return "B"
    if any(g.startswith("TRA") for g in genes):
        return "TRA"
    if any(g.startswith("TRB") for g in genes):
        return "TRB"
    if any(g.startswith("TRG") for g in genes):
        return "TRG"
    if any(g.startswith("TRD") for g in genes):
        return "TRD"
    return "Other"

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate MiXCR alignment + clone files into spatial counts."
    )
    parser.add_argument(
        "-a", "--alignment", required=True, help="MiXCR alignment file (TSV)"
    )
    parser.add_argument(
        "-c", "--clones", required=True, help="MiXCR clones file (TSV)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output CSV file with spatial counts"
    )
    args = parser.parse_args()

    # --- Load data ---
    align_df = pd.read_csv(args.alignment, sep="\t", comment="#")
    clones_df = pd.read_csv(args.clones, sep="\t", comment="#")

    # --- Extract coordinates ---
    align_df[["x", "y"]] = align_df["descrsR1"].apply(
        lambda d: pd.Series(parse_coords(str(d)))
    )

    # Keep only valid coords and valid clones
    align_df = align_df[
        (align_df["x"].notna()) & (align_df["y"].notna()) & (align_df["cloneId"] >= 0)
    ]

    # --- Assign receptor type ---
    clones_df["Receptor"] = clones_df["allVHitsWithScore"].apply(classify_receptor)

    # --- Join alignment with clones ---
    merged = align_df.merge(clones_df[["cloneId", "Receptor"]], on="cloneId", how="left")

    # --- Aggregate counts ---
    counts = (
        merged.groupby(["x", "y", "Receptor"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    # Ensure consistent columns
    for col in ["B", "TRA", "TRB", "TRG", "TRD", "Other"]:
        if col not in counts.columns:
            counts[col] = 0

    # --- Save ---
    counts.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
