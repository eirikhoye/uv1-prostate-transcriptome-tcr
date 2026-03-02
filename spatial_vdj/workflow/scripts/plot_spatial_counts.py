#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def plot_spatial(input_csv, output_file):
    # --- Load counts ---
    df = pd.read_csv(input_csv)

    # Add a total T-cell count column
    df["T"] = df[["TRA", "TRB", "TRG", "TRD"]].sum(axis=1)

    # --- Create side-by-side plots ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)

    # Left: B cells
    sc1 = axes[0].scatter(df["x"], df["y"], c=df["B"], cmap="Blues", s=100, marker="s")
    axes[0].set_title("B cells")
    axes[0].set_xlabel("X coordinate")
    axes[0].set_ylabel("Y coordinate")
    axes[0].invert_yaxis()
    cbar1 = fig.colorbar(sc1, ax=axes[0])
    cbar1.set_label("B cell count")

    # Right: T cells
    sc2 = axes[1].scatter(df["x"], df["y"], c=df["T"], cmap="Reds", s=100, marker="s")
    axes[1].set_title("T cells")
    axes[1].set_xlabel("X coordinate")
    axes[1].invert_yaxis()
    cbar2 = fig.colorbar(sc2, ax=axes[1])
    cbar2.set_label("T cell count")

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Plot spatial distribution of B and T cells side by side."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input CSV file with spatial counts"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output image file (e.g. PNG, PDF)"
    )
    args = parser.parse_args()

    plot_spatial(args.input, args.output)

if __name__ == "__main__":
    main()
