#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import sys

def parse_mixcr_report(report_path):
    """
    Parse a MiXCR alignment report and return a dict with relevant counts.
    """
    report_path = Path(report_path)
    if not report_path.is_file():
        raise FileNotFoundError(f"Report file not found: {report_path}")

    data = {"sample": report_path.stem}
    chains = ["TRA", "TRB", "TRG", "IGH", "IGK", "IGL"]
    
    with open(report_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Total sequencing reads:"):
                data["total_reads"] = int(line.split(":")[1].split()[0].replace(",", ""))
            elif line.startswith("Successfully aligned reads:"):
                data["aligned_reads"] = int(line.split(":")[1].split("(")[0].replace(",", ""))
            elif any(line.startswith(chain + " chains:") for chain in chains):
                chain_name = line.split()[0]
                count = int(line.split()[2].replace(",", ""))
                data[chain_name] = count

    # Compute unaligned reads
    data["unaligned_reads"] = data["total_reads"] - data.get("aligned_reads", 0)
    # Compute aligned but not assigned to known chains
    assigned_sum = sum(data.get(c, 0) for c in chains)
    data["aligned_other"] = data.get("aligned_reads", 0) - assigned_sum
    return data

def plot_stacked_bar(df, output_file):
    """
    Create and save a stacked barplot of MiXCR read counts.
    """
    df = df.set_index("sample").fillna(0)
    chains = ["TRA", "TRB", "TRG", "IGH", "IGK", "IGL", "aligned_other", "unaligned_reads"]
    colors = {
        "TRA": "#1f77b4",
        "TRB": "#ff7f0e",
        "TRG": "#2ca02c",
        "IGH": "#d62728",
        "IGK": "#9467bd",
        "IGL": "#8c564b",
        "aligned_other": "#e377c2",
        "unaligned_reads": "#7f7f7f"
    }

    ax = df[chains].plot(kind="bar", stacked=True, color=[colors[c] for c in chains], figsize=(8,12))
    ax.set_ylabel("Read counts")
    ax.set_xlabel("Sample")
    ax.set_xticklabels(df.index, rotation=90, ha="right")
    ax.legend(title="Chain / Read category", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate stacked barplot from MiXCR alignment reports."
    )
    parser.add_argument(
        "reports", nargs="+",
        help="Paths to MiXCR _IGH-align.report files."
    )
    parser.add_argument(
        "-o", "--output", default="mixcr_reads_barplot.png",
        help="Output filename for the plot (default: mixcr_reads_barplot.png)."
    )
    args = parser.parse_args()

    # Parse all reports
    data_list = []
    for report in args.reports:
        try:
            data_list.append(parse_mixcr_report(report))
        except FileNotFoundError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    df = pd.DataFrame(data_list)
    plot_stacked_bar(df, args.output)

if __name__ == "__main__":
    main()
