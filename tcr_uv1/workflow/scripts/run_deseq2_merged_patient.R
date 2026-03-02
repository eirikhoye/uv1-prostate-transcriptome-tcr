suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("pheatmap"))


counts <- read.csv(snakemake@input[["counts"]], row.names=1)
colnames(counts) <- gsub("\\.", "-", gsub("^X", "", colnames(counts)))
counts$sharing_level <- NULL
coldata <- read.csv(snakemake@input[["clinical"]], row.names=1)

coldata$Response <- ifelse(coldata$Response %in% c("Late", "No"), "Late_No", "Early")
coldata$Response <- factor(coldata$Response, levels = c("Early", "Late_No"))

stopifnot(all(colnames(counts) %in% rownames(coldata)))
coldata <- coldata[colnames(counts), ]

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=coldata, 
                              design=~Response)
print(dds)
dds <- DESeq(dds, sfType = 'poscounts')
print(dds)
res <- results(dds)
print(res)
write.csv(as.data.frame(res), snakemake@output[["results"]])

vsd <- vst(dds, blind=FALSE)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("PFS_time", "OS_time", "Response")])
png(filename = snakemake@output[["heatmap"]], width=800, height=800)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

png(filename = snakemake@output[["pca"]], width = 800, height=800)
plotPCA(vsd, intgroup=c("Response"))
dev.off()
