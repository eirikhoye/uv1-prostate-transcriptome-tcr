#!/usr/bin/env python3
import pandas as pd
import os
import re
import argparse

def main(input_file, outdir):
    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Load CSV
    df = pd.read_csv(input_file)

    # Drop 'sharing_level' if present
    if "sharing_level" in df.columns:
        df = df.drop(columns=["sharing_level"])

    # Iterate over all sample columns (everything except 'motif')
    for col in df.columns:
        if col == "motif":
            continue

        # Extract patient ID from column name
        match = re.search(r"(Pat\d+)", col)
        patient_id = match.group(1) if match else "Unknown"

        # Extract motif + this sample
        sub = df[["motif", col]].copy()
        
        # Filter out zero counts
        sub = sub[sub[col] != 0]

        # Add patient_id column
        sub["patient_id"] = patient_id

        # Rename count column to "counts"
        sub = sub.rename(columns={col: "counts"})

        # Reorder columns: motif, patient_id, counts
        sub = sub[["motif", "patient_id", "counts"]]

        # Output path
        outpath = os.path.join(outdir, f"{col}.csv")

        # Save
        sub.to_csv(outpath, index=False)
        print(f"Saved {outpath} with {len(sub)} rows")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split TCR counts file into per-sample CSVs with patient ID.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for per-sample CSVs")
    args = parser.parse_args()
    
    main(args.input, args.outdir)
