import scanpy as sc
import pandas as pd
from scipy import io
import squidpy as sq
import gzip, json
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# --- Utility: Load hires and force white background ---
def load_hires_white(path_gz):
    with gzip.open(path_gz, "rb") as f:
        img = np.array(Image.open(f).convert("RGBA"))
    # Define mask of tissue (non-bright pixels)
    tissue_mask = img[..., :3].sum(axis=2) < 740  # tweak threshold if needed
    img_bg_white = np.ones_like(img) * 255
    img_bg_white[tissue_mask] = img[tissue_mask]
    return img_bg_white

# --- Define samples (your 4 dicts go here) ---
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
    }
]

def load_adata(sample):
    X = io.mmread(sample["matrix"]).T.tocsr()
    obs = pd.read_csv(sample["barcodes"], header=None, names=["barcode"])
    var = pd.read_csv(sample["features"], header=None, sep="\t")
    var.columns = ["gene_id","gene_name","feature_type"]
    adata = sc.AnnData(X=X, obs=obs, var=var)

    pos = pd.read_csv(sample["positions"], header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    adata.obs = adata.obs.set_index("barcode").join(pos.set_index("barcode"))
    adata.obsm["spatial"] = adata.obs[["pxl_col","pxl_row"]].to_numpy()

    hires_white = load_hires_white(sample["hires"])
    with gzip.open(sample["scalefactors"], "rt") as f:
        scalefactors = json.load(f)
    adata.uns["spatial"] = {"sample": {"scalefactors": scalefactors, "images": {"hires": hires_white}}}

    vdj = pd.read_csv(sample["vdj"])
    vdj = vdj.rename(columns={"x":"array_col","y":"array_row"})
    adata.obs = adata.obs.merge(vdj, on=["array_row","array_col"], how="left")
    for col in ["B","TRA","TRB","TRD","TRG","Other"]:
        if col in adata.obs:
            adata.obs[col] = adata.obs[col].fillna(0)
    return adata

# --- Collect global vmax values for consistency ---
all_B = []
all_TRB = []
for sample in samples:
    adata = load_adata(sample)
    all_B.extend(adata.obs["B"])
    all_TRB.extend(adata.obs["TRB"])
global_vmax_B = np.percentile(all_B, 99)
global_vmax_TRB = np.percentile(all_TRB, 99)

# --- Prepare grid ---
fig, axes = plt.subplots(len(samples), 3, figsize=(15, 20))
plt.subplots_adjust(wspace=0.05, hspace=0.05)

for i, sample in enumerate(samples):
    adata = load_adata(sample)

    # Column 1: H&E (lightened background only)
    sq.pl.spatial_scatter(
        adata,
        img="hires",
        size=1.0,
        legend_loc=None,
        ax=axes[i,0],
        img_alpha=1.0
    )
    axes[i,0].set_title(f"{sample['name']} H&E", fontsize=12)
    axes[i,0].axis("off")

    # Column 2: B cells (blue overlay)
    sq.pl.spatial_scatter(
        adata,
        color="B",
        img="hires",
        img_alpha=0.2,
        size=1.0,
        cmap="Blues",
        vmin=0,
        vmax=global_vmax_B,
        legend_loc=None,
        ax=axes[i,1],
    )
    if i == 0: axes[i,1].set_title("B cells", fontsize=12)
    axes[i,1].axis("off")

    # Column 3: TRB (orange overlay)
    sq.pl.spatial_scatter(
        adata,
        color="TRB",
        img="hires",
        img_alpha=0.2,
        size=1.0,
        cmap="Oranges",
        vmin=0,
        vmax=global_vmax_TRB,
        legend_loc=None,
        ax=axes[i,2],
    )
    if i == 0: axes[i,2].set_title("TRB T cells", fontsize=12)
    axes[i,2].axis("off")

plt.savefig("spatial_grid_clean_nocbar.png", dpi=600, bbox_inches="tight")
plt.show()
