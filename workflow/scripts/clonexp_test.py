import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# ------------------------------------ Functions ------------------------------------

def analyze_expansion_per_patient(merged_df, patient_id, min_reads=3, pseudo=1):
    """
    Identify significantly expanded clonotypes for a given patient using one-sided Fisher's test.
    """

    df = merged_df.query("patient_id == @patient_id and sample_type == 'PBMC'").copy()
    if df.empty:
        raise ValueError(f"No PBMC samples found for patient {patient_id}.")

    # Determine baseline and first post-treatment day
    days = np.sort(df['day_after_treatment'].unique())
    baseline_day = days[0]
    post_day = days[days > baseline_day][0] if len(days) > 1 else None
    if post_day is None:
        raise ValueError(f"Patient {patient_id} has no post-treatment samples.")

    # Run Fisher’s test per motif
    results = []
    for motif, g in df.groupby('motif'):
        counts = g.groupby('day_after_treatment')['count'].sum()
        totals = g.groupby('day_after_treatment')['total_per_sample'].first()

        a = counts.get(baseline_day, 0)
        b = counts.get(post_day, 0)

        if (a < min_reads) and (b < min_reads):
            continue

        a += pseudo
        b += pseudo
        c = totals.get(baseline_day, 0) - a + pseudo
        d = totals.get(post_day, 0) - b + pseudo

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='less')

        log2_fc = np.log2((b / totals.get(post_day, 1) + 1e-9) /
                          (a / totals.get(baseline_day, 1) + 1e-9))

        results.append({
            'patient_id': patient_id,
            'motif': motif,
            'baseline_count': a,
            'post_count': b,
            'p_value': p_value,
            'log2_fc': log2_fc
        })

    res_df = pd.DataFrame(results)
    if res_df.empty:
        return pd.DataFrame(), 0

    res_df['p_adj'] = multipletests(res_df['p_value'], method='fdr_by')[1]
    res_df['significant'] = res_df['p_adj'] < 0.05
    n_expanded = res_df['significant'].sum()

    return res_df, n_expanded


def load_and_prepare_data(motif_file, metadata_file):
    """
    Load motif counts and sample metadata, merge into long-format DataFrame.
    """

    motif_df = pd.read_csv(motif_file, sep=",")
    samples_metadata = pd.read_csv(metadata_file)
    samples_metadata['Response'] = samples_metadata['Response'].replace({'Late': 'Late/No', 'No': 'Late/No'})
    samples_metadata = samples_metadata[samples_metadata['sample'].isin(motif_df.columns)]

    long_df = motif_df.melt(
        id_vars=['motif', 'sharing_level'],
        var_name='sample',
        value_name='count'
    )

    merged_df = long_df.merge(samples_metadata, on='sample', how='left')
    merged_df = merged_df.dropna(subset=['patient_id'])
    merged_df['count'] = merged_df['count'].astype(float)

    sample_totals = merged_df.groupby('sample')['count'].sum()
    merged_df['total_per_sample'] = merged_df['sample'].map(sample_totals)
    merged_df['count_cpm'] = (merged_df['count'] / merged_df['total_per_sample']) * 1_000_000

    return merged_df


# ------------------------------------ Main ------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Analyze PBMC clonotype expansion per patient.")
    parser.add_argument('--motif_file', required=True, help='CSV file with motif counts')
    parser.add_argument('--metadata_file', required=True, help='CSV file with sample metadata')
    parser.add_argument('--output_prefix', default='patient', help='Prefix for output CSV files')
    parser.add_argument('--patients', nargs='+', type=int, default=None, help='Patient IDs to analyze (default: all)')
    args = parser.parse_args()

    merged_df = load_and_prepare_data(args.motif_file, args.metadata_file)

    if args.patients is None:
        patient_ids = merged_df['patient_id'].unique()
    else:
        patient_ids = args.patients

    summary_list = []
    for pid in patient_ids:
        res_df, n_expanded = analyze_expansion_per_patient(merged_df, pid)
        summary_list.append({'patient_id': pid, 'n_significant_expanded': n_expanded})
        res_df.to_csv(f"{args.output_prefix}_{pid}_expanded_clonotypes.csv", index=False)

    summary_df = pd.DataFrame(summary_list)
    summary_df.to_csv(f"{args.output_prefix}_summary.csv", index=False)


if __name__ == "__main__":
    main()
