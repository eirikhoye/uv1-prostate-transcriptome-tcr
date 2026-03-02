suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("pheatmap"))

counts <- read.csv(snakemake@input[["counts"]], row.names=1)
colnames(counts) <- gsub("\\.", "-", gsub("^X", "", colnames(counts)))
counts$sharing_level <- NULL

coldata <- read.csv(snakemake@input[["clinical"]], row.names=1)
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
print(dds_filtered)

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

df <- as.data.frame(colData(dds_filtered)[,c("day_after_treatment","sample_type",
                                    "PFS", "PFS_time",
                                    "OS", "OS_time", "Response")])

png(filename = snakemake@output[["heatmap"]], width=800, height=800)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

png(filename = snakemake@output[["pca"]], width = 800, height=800)
plotPCA(vsd, intgroup=c("Response"))
dev.off()
