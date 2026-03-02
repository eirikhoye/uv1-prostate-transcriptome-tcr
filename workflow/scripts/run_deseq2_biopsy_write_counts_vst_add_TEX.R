suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library('viridis'))
suppressPackageStartupMessages(library('RColorBrewer'))

counts <- read.csv(snakemake@input[["counts"]], row.names=1)
colnames(counts) <- gsub("\\.", "-", gsub("^X", "", colnames(counts)))
counts$sharing_level <- NULL

coldata <- read.csv(snakemake@input[["clinical"]], row.names=2)
coldata$Response <- ifelse(coldata$Response %in% c("Late", "No"), "Late_No", "Early")
coldata$Response <- factor(coldata$Response, levels = c("Late_No", "Early"))

stopifnot(all(colnames(counts) %in% rownames(coldata)))
coldata <- coldata[colnames(counts), ]

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=coldata, 
                              design=~Response)

dds <- DESeq(dds, sfType = 'poscounts')
res <- results(dds)

write.csv(as.data.frame(res), snakemake@output[["results"]])

# Filter out motif sequences where there are less than 3 samples with 1 count
keep <- rowSums(counts(dds) >= 1) >= 3
dds_filtered <- dds[keep, ]

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds_filtered, blind=FALSE)

# --- Save counts to CSV for Python ---
# Normalized counts
norm_counts <- counts(dds_filtered, normalized=TRUE)
write.csv(as.data.frame(norm_counts), snakemake@output[["norm_counts"]])

# VST counts
vst_counts <- assay(vsd)
write.csv(as.data.frame(vst_counts), snakemake@output[["vst_counts"]])
# -------------------------------------

# Select top variable genes for heatmap
select <- order(apply(norm_counts, 1, var), decreasing=TRUE)[1:100]

#df <- as.data.frame(colData(dds_filtered)[,c("day_after_treatment","sample_type",
#                                    "PFS", "PFS_time",
#                                    "OS", "OS_time", "Response", "TEX_assigned")])

df <- as.data.frame(colData(dds_filtered)[,c("day_after_treatment", "OS_time", "Response", "TEX_assigned")])

colnames(df) <- c("Days prior to vaccination", "OS", "Response", "TEX")

annotation_colors <- list(
    TEX = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    Response = c(
        "Late_No" = "#B0B0B0",   # gray for non-responder
        "Early"   = "#E41A1C"    # red for responder
    ),
    OS = colorRampPalette(c("white", "seagreen3", "darkgreen"))(100),
    "Days prior to vaccination" = colorRampPalette(c("navy", "white", "goldenrod1"))(100)        
)


# --- Standard heatmap for analysis ---
png(filename = snakemake@output[["heatmap"]], width=800, height=800)
pheatmap(
    assay(vsd)[select,], 
    cluster_rows=TRUE, 
    show_rownames=FALSE,
    cluster_cols=TRUE, 
    annotation_col=df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
    )

dev.off()

# --- High-resolution poster heatmap ---
png(filename = snakemake@output[["heatmap_poster"]], width=10000, height=12000, res=600)
pheatmap(
    assay(vsd)[select,],
    cluster_rows=TRUE,
    show_rownames=FALSE,
    show_colnames=FALSE,
    cluster_cols=TRUE,
    annotation_col=df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(200),
    fontsize = 28,              # main text
    fontsize_col = 32,          # column labels
    fontsize_row = 22,          # row labels (if shown)
    fontsize_number = 20,       # for numbers if added
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 80,
    treeheight_col = 80,
    border_color = NA
)
dev.off()

# --- High-resolution poster heatmap ---
pdf(
    file = snakemake@output[["heatmap_manuscript"]],
    width = 16.7,   # roughly 10000 px / 600 dpi
    height = 14,    # roughly 12000 px / 600 dpi
    useDingbats = FALSE
)
pheatmap(
    assay(vsd)[select,],
    cluster_rows=TRUE,
    show_rownames=FALSE,
    show_colnames=FALSE,
    cluster_cols=TRUE,
    annotation_col=df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(200),
    fontsize = 28,              # main text
    fontsize_col = 32,          # column labels
    fontsize_row = 22,          # row labels (if shown)
    fontsize_number = 20,       # for numbers if added
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 80,
    treeheight_col = 80,
    border_color = NA
)
dev.off()

png(filename = snakemake@output[["pca"]], width = 800, height=800)
plotPCA(vsd, intgroup=c("Response"))
dev.off()
