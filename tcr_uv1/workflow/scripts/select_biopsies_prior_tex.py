#!/usr/bin/env python

import pandas as pd
import numpy as np
import re
import argparse

def main(counts_path, clinical_path, signature_path, sig_sample_path, counts_out, clinical_out):
    # Load data
    counts_df = pd.read_csv(counts_path)
    clinical_df = pd.read_csv(clinical_path)
    tex_df = pd.read_excel(signature_path, engine='openpyxl')
    tex_sample_df = pd.read_excel(sig_sample_path, engine='openpyxl')

    # Process TEX
    tex_values = tex_df.set_index('Unnamed: 0').loc['TEX']
    tex_transposed = tex_values.to_frame(name='TEX').reset_index()
    tex_transposed.rename(columns={'index': 'Original_Column'}, inplace=True)
    tex_transposed[['Prefix', 'Aribtrary tube number']] = tex_transposed['Original_Column'].str.split('-', n=1, expand=True)

    tex_sample_df['Original_Column'] = tex_sample_df['Prefix'].astype(str) + tex_sample_df['Aribtrary tube number'].astype(str)

    merged_df = tex_transposed.merge(
        tex_sample_df[['Prefix', 'Aribtrary tube number', 'Correspond to patient ID number', 'Original_Column']],
        on='Original_Column',
        how='left'
    )

    # Merge TEX into clinical
    merged_for_clinical = pd.merge(
        clinical_df, 
        merged_df[['TEX', 'Prefix_x', 'Correspond to patient ID number']].rename(columns={'Prefix':'Prefix_x'}), 
        left_on='patient_id', 
        right_on='Correspond to patient ID number', 
        how='left'
    )

    # Mask conditions
    mask = (
        (merged_for_clinical['sample_type'] == 'Biopsy') &
        (
            ((merged_for_clinical['Prefix_x'] == '1stB') & (merged_for_clinical['day_after_treatment'] < 0)) |
            ((merged_for_clinical['Prefix_x'] == '2stB') & (merged_for_clinical['day_after_treatment'] > 0))
        )
    )
    merged_for_clinical['TEX_assigned'] = np.where(mask, merged_for_clinical['TEX'], np.nan)

    # Select biopsy columns
    biopsy_cols = [col for col in counts_df.columns if re.match(r"^-\d+_Biopsy", col)]
    counts_biopsy = counts_df[['motif'] + biopsy_cols]

    result_df = merged_for_clinical[merged_for_clinical['sample'].isin(biopsy_cols)].reset_index(drop=True)
    result_df_sorted = result_df.sort_values(by=['patient_id', 'TEX_assigned'], ascending=[True, False])
    result_df_unique = result_df_sorted.groupby('patient_id', as_index=False).first().reset_index(drop=True)

    # Write output
    counts_biopsy.to_csv(counts_out, index=False)
    result_df_unique.to_csv(clinical_out, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Select biopsy samples and merge TEX scores")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--sig_sample", required=True)
    parser.add_argument("--counts_out", required=True)
    parser.add_argument("--clinical_out", required=True)
    args = parser.parse_args()

    main(args.counts, args.clinical, args.signature, args.sig_sample, args.counts_out, args.clinical_out)
