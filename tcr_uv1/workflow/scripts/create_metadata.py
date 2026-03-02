import pandas as pd
import re

# Load or define uv1_clinical as a DataFrame
# Import and process files
surv_uv1 = pd.read_excel("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/data/Survival_table_UV1_RNA_patients_UPDATED.xlsx", sheet_name=0)
ir_uv1 = pd.read_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/metadata/patient_data.csv")
ir_uv1 = ir_uv1[ir_uv1['UV1'] == 'x']
ir_uv1 = ir_uv1[['Correspond to patient ID number', 'Response']].drop_duplicates().rename(columns={"Correspond to patient ID number":"medid"})
ir_uv1['medid'] = ir_uv1['medid'].replace('811-UV1', '811')
surv_uv1['medid'] = surv_uv1['medid'].astype(str)
ir_uv1['medid'] = ir_uv1['medid'].astype(str)

# Merge files
uv1_clinical = pd.merge(surv_uv1, ir_uv1, on='medid', how='left')



# Sample header list (replace this with your actual header extraction)
with open("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/results/motif_counts.csv") as f:
    header = f.readline().strip().split(',')

# Drop the first column name ('motif' or similar) and the last if it's 'sharing_level'
sample_columns = header[1:]
if sample_columns[-1].lower() == "sharing_level":
    sample_columns = sample_columns[:-1]

# Parse metadata from column names
sample_meta = []
pattern = re.compile(r"(?P<day>-?\d+)_?(?P<type>PBMC|TCC|Biopsy|Bulk)_patient_Pat(?P<patient>\d+)")

for sample in sample_columns:
    match = pattern.match(sample)
    if match:
        meta = match.groupdict()
        meta["sample"] = sample
        meta["day_after_treatment"] = int(meta.pop("day"))
        meta["sample_type"] = meta.pop("type")
        meta["patient_id"] = int(meta.pop("patient"))
        sample_meta.append(meta)
    else:
        print(f"Warning: Could not parse sample name: {sample}")

# Create metadata DataFrame
meta_df = pd.DataFrame(sample_meta)

# Ensure both are string type for merge
meta_df["patient_id"] = meta_df["patient_id"].astype(str)
uv1_clinical["medid"] = uv1_clinical["medid"].astype(str)


# Merge with clinical metadata
meta_df = meta_df.merge(uv1_clinical, how="left", left_on="patient_id", right_on="medid")

# Optional: reorder columns
columns_order = [
    "sample", "patient_id", "day_after_treatment", "sample_type",
    "PFS", "PFS_time", "OS", "OS_time", "cohort", "Response"
]
meta_df = meta_df[columns_order]

# Save output
meta_df.to_csv("/storage/mathelierarea/processed/eirikhoy/tcr_uv1/workflow/utils/sample_metadata.csv", index=False)
print("✅ Metadata written to sample_metadata.csv")