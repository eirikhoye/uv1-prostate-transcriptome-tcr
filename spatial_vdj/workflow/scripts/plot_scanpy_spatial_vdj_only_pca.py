import scanpy as sc
import pandas as pd
from scipy import io
import squidpy as sq
import gzip
import json
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# Define your samples with file paths
samples = [
    {
        "name": "PCa1",
        "matrix": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_matrix.mtx.gz",
        "barcodes": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_barcodes.tsv.gz",
        "features": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_features.tsv.gz",
        "positions": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_tissue_positions_list.csv.gz",
        "scalefactors": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_scalefactors_json.json.gz",
        "hires": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557997_NEADT_1_tissue_hires_image.png.gz",
        "vdj": "/storage/mathelierarea/processed/eirikhoy/vdj_spatial/mixcr_spatial_counts/m84212_250410_123503_s3.hifi_reads.bcAd1031T_spatial_conts.txt",
        "vmax_B": 0.8,
        "vmax_TRB": 0.95
    },
    {
        "name": "PCa2",
        "matrix": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_matrix.mtx.gz",
        "barcodes": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_barcodes.tsv.gz",
        "features": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_features.tsv.gz",
        "positions": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_tissue_positions_list.csv.gz",
        "scalefactors": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_scalefactors_json.json.gz",
        "hires": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557983_TRNA_4_tissue_hires_image.png.gz",
        "vdj": "/storage/mathelierarea/processed/eirikhoy/vdj_spatial/mixcr_spatial_counts/m84212_250410_123503_s3.hifi_reads.bcAd1032T_spatial_conts.txt",
        "vmax_B": 0.8,
        "vmax_TRB": 0.95
    },
#    {
#        "name": "CRPC",
#        "matrix": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_matrix.mtx.gz",
#        "barcodes": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_barcodes.tsv.gz",
#        "features": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_features.tsv.gz",
#        "positions": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_tissue_positions_list.csv.gz",
#        "scalefactors": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_scalefactors_json.json.gz",
#        "hires": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8558023_CRPC_5_tissue_hires_image.png.gz",
#        "vdj": "/storage/mathelierarea/processed/eirikhoy/vdj_spatial/mixcr_spatial_counts/m84212_250410_123503_s3.hifi_reads.bcAd1033T_spatial_conts.txt",
#        "vmax_B": 0.98,
#        "vmax_TRB": 10
#    },
#    {
#        "name": "BPH",
#        "matrix": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_matrix.mtx.gz",
#        "barcodes": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_barcodes.tsv.gz",
#        "features": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_features.tsv.gz",
#        "positions": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_tissue_positions_list.csv.gz",
#        "scalefactors": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_scalefactors_json.json.gz",
#        "hires": "/storage/mathelierarea/raw/urbanucci_vdj_spatial/data/GSE278936_RAW/GSM8557976_BPH_1_tissue_hires_image.png.gz",
#        "vdj": "/storage/mathelierarea/processed/eirikhoy/vdj_spatial/mixcr_spatial_counts/m84212_250410_123503_s3.hifi_reads.bcAd1034T_spatial_conts.txt",
#        "vmax_B": 0.8,
#        "vmax_TRB": 0.95
#    },
]

def load_adata(sample):
    # Load count matrix
    X = io.mmread(sample["matrix"]).T.tocsr()
    obs = pd.read_csv(sample["barcodes"], header=None, names=["barcode"])
    var = pd.read_csv(sample["features"], header=None, sep="\t")
    var.columns = ["gene_id","gene_name","feature_type"]
    adata = sc.AnnData(X=X, obs=obs, var=var)

    # Load spatial positions
    pos = pd.read_csv(sample["positions"], header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    adata.obs = adata.obs.set_index("barcode").join(pos.set_index("barcode"))
    adata.obsm["spatial"] = adata.obs[["pxl_col","pxl_row"]].to_numpy()

    # Load images and scalefactors
    with gzip.open(sample["scalefactors"], "rt") as f:
        scalefactors = json.load(f)
    with gzip.open(sample["hires"], "rb") as f:
        hires = np.array(Image.open(f))
    adata.uns["spatial"] = {"sample": {"scalefactors": scalefactors, "images": {"hires": hires}}}

    # Load VDJ counts
    vdj = pd.read_csv(sample["vdj"])
    vdj = vdj.rename(columns={"x":"array_col","y":"array_row"})
    adata.obs = adata.obs.merge(vdj, on=["array_row","array_col"], how="left")
    for col in ["B","TRA","TRB","TRD","TRG","Other"]:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].fillna(0)
    return adata

# Prepare grid: 4 samples x 3 columns (H&E, B, TRB)
fig, axes = plt.subplots(len(samples), 3, figsize=(15, 10))
plt.subplots_adjust(wspace=0.05, hspace=0.05)

for i, sample in enumerate(samples):
    adata = load_adata(sample)
    
    # Column 1: H&E
    sq.pl.spatial_scatter(
        adata,
        img="hires",
        size=1.2,
        legend_loc=None,
        ax=axes[i,0],
        #show=False
    )
    axes[i,0].set_title(f"{sample['name']} H&E", fontsize=12)
    axes[i,0].axis("off")
    
    # Column 2: B cells
    vmax_B = adata.obs["B"].quantile(sample["vmax_B"])
    sq.pl.spatial_scatter(
        adata,
        color=["B"],
        img="hires",
        img_alpha=0.25,
        size=1.2,
        cmap="viridis",
        vmin=0,
        vmax=vmax_B,
        legend_loc=None,
        ax=axes[i,1],
        colorbar=False
    )
    axes[i,1].set_title("BCRs", fontsize=12)
    axes[i,1].axis("off")
    
    # Column 3: TRB T cells
    vmax_TRB = sample["vmax_TRB"]
    sq.pl.spatial_scatter(
        adata,
        color=["TRB"],
        img="hires",
        img_alpha=0.25,
        size=1.2,
        cmap="Oranges",
        vmin=0,
        vmax=vmax_TRB,
        legend_loc=None,
        ax=axes[i,2],
        colorbar=False
    )
    axes[i,2].set_title("TCRs", fontsize=12)
    axes[i,2].axis("off")

plt.tight_layout()
plt.savefig("spatial_grid_only_pca.png", dpi=1200)
plt.savefig("spatial_grid_only_pca.pdf", dpi=1200)
