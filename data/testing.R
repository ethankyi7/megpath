library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)

poolOneData <- ReadMtx(
    mtx = "genePool1/GSE213902_MQpool1_matrix.mtx.gz",
    cells = "genePool1/GSE213902_MQpool1_barcodes.tsv.gz",
    features = "genePool1/GSE213902_MQpool1_features.tsv.gz",
)

poolTwoData <- ReadMtx(
    mtx = "genePool2/GSE213902_MQpool2_matrix.mtx.gz",
    cells = "genePool2/GSE213902_MQpool2_barcodes.tsv.gz",
    features = "genePool2/GSE213902_MQpool2_features.tsv.gz",
)

poolOneHTO <- ReadMtx(
    mtx =  "sra/poolOne/umi_count/matrix.mtx.gz",
    cells = "sra/poolOne/umi_count/barcodes.tsv.gz",
    features = "sra/poolOne/umi_count/features.tsv.gz",
    feature.column = 1
)

poolTwoHTO <- ReadMtx(
    mtx =  "sra/poolTwo/umi_count/matrix.mtx.gz",
    cells = "sra/poolTwo/umi_count/barcodes.tsv.gz",
    features = "sra/poolTwo/umi_count/features.tsv.gz",
    feature.column = 1
)

poolOneHTO <- poolOneHTO[rownames(poolOneHTO) != "unmapped", ]
poolTwoHTO <- poolTwoHTO[rownames(poolTwoHTO) != "unmapped", ]

colnames(poolOneData) <- gsub("-1$", "", colnames(poolOneData))
colnames(poolTwoData) <- gsub("-1$", "", colnames(poolTwoData))

#find overlapping barcodes in hto and genome expression (common cells)
poolOneBarcodes <- intersect(colnames(poolOneData), colnames(poolOneHTO))
poolTwoBarcodes <- intersect(colnames(poolTwoData), colnames(poolTwoHTO))
#check for how many were removed
removed_one <- length(colnames(poolOneData)) - length(poolOneBarcodes)
removed_two <- length(colnames(poolTwoData)) - length(poolTwoBarcodes)
print(paste("removed from pool one: ", removed_one))
print(paste("removed from pool two: ", removed_two))

poolOneData <- poolOneData[, poolOneBarcodes]
poolOneHTO <- as.matrix(poolOneHTO[, poolOneBarcodes])
poolTwoData <- poolTwoData[, poolTwoBarcodes]
poolTwoHTO <- as.matrix(poolTwoHTO[, poolTwoBarcodes])

poolOne <- CreateSeuratObject(counts = poolOneData)
poolOne[["percent.mt"]] <- PercentageFeatureSet(poolOne, pattern = "^MT-")
poolOne[["HTO"]] <- CreateAssayObject(counts = poolOneHTO)
poolTwo <- CreateSeuratObject(counts = poolTwoData)
poolTwo[["percent.mt"]] <- PercentageFeatureSet(poolTwo, pattern = "^MT-")
poolTwo[['HTO']] <- CreateAssayObject(counts = poolTwoHTO)

poolOne <- NormalizeData(poolOne, assay = "HTO", normalization.method = "CLR")
poolTwo <- NormalizeData(poolTwo, assay = "HTO", normalization.method = "CLR")

poolOne <- HTODemux(
    poolOne,
    assay = "HTO",
    positive.quantile = 0.995,
    kfunc = 'clara'
)

poolTwo <- HTODemux(
    poolTwo,
    assay = "HTO",
    positive.quantile = 0.995,
    kfunc = 'clara'
)

pooled <- merge(x = poolOne, y = poolTwo, add.cell.ids = c("P1", "P2"))
print(pooled)
# #FeatureScatter(pooled, feature2 = "percent.mt", feature1 = "nFeature_RNA")

pooled <- JoinLayers(pooled, assay = "RNA", layer = "counts")
counts <- GetAssayData(pooled, assay = "RNA", layer = "counts")

# writeMM(counts, "counts.mtx")
# write.table(rownames(counts), "genes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
# write.table(colnames(counts), "barcodes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)

scrub <- read.csv("scrublet_results.csv")
#print(scrub[match(colnames(pooled), scrub$barcode), ])
#print(scrub$doublet_score)
pooled[["doublet_score"]] <- scrub$doublet_score
#summary(pooled$doublet_score)

# df <- pooled@meta.data
# df <- mutate (
#         df, 
#         log10_features = log10(nFeature_RNA),
#         log10_umi = log10(nCount_RNA)
#     )

# ggplot(df, aes(x = log10_umi, y = log10_features)) +
#     geom_bin2d(bins = 150) +
#     scale_fill_gradient(low = "lightblue", high = "darkblue") +
#     theme_classic() +
#     labs(
#         x = "log10(# of counts)",
#         y = "log10(# of features)"
#     ) +
#     geom_hline(yintercept = log10(300))

# ggplot(df, aes(x = log10_features, y = percent.mt)) +
#     geom_bin2d(bins = 150) +
#     scale_fill_gradient(low = "lightblue", high = "darkblue") +
#     theme_classic() +
#     labs(
#         x = "log10(# of Features)",
#         y = "percent.mt"
#     ) + 
#     geom_vline(xintercept = log10(300)) +
#     geom_hline(yintercept = 7.5)

# ggplot(df, aes(x = log10_features, y = doublet_score)) +
#     geom_bin2d(bins = 150) +
#     scale_fill_gradient(low = "lightblue", high = "darkblue") +
#     theme_classic() +
#     labs(
#         x = "log10(# of features)",
#         y = "doublet score"
#     ) +
#     geom_hline(yintercept = 0.2) +
#     geom_vline(xintercept = log10(300))

# pooled <- subset(pooled, subset = nFeature_RNA >= 300 & nFeature_RNA < 7500 & doublet_score < 0.2 & percent.mt <= 7.5)
# print(pooled)
# table(pooled$HTO_classification.global)
# table(pooled$hash.ID)
# table(pooled$HTO_maxID)

# pooled <- NormalizeData(pooled)
# pooled <- FindVariableFeatures(pooled)
# pooled <- ScaleData(pooled)
# pooled <- RunPCA(pooled)
# pooled <- RunUMAP(pooled, dims = 1:30)
# print(Reductions(pooled))

# FeaturePlot(
#   pooled,
#   features = "nCount_RNA",
#   cols = c("lightgrey", "red"),
#   pt.size = 0.4,
#   ncol = 1
# )

# FeaturePlot(
#   pooled,
#   features = "percent.mt",
#   cols = c("lightgrey", "red"),
#   pt.size = 0.4,
#   ncol = 1
# )

# FeaturePlot(
#   pooled,
#   features = "doublet_score",
#   cols = c("lightgrey", "red"),
#   pt.size = 0.4,
#   ncol = 1
# )

# DimPlot(
#   pooled,
#   group.by = "HTO_classification.global",
#   pt.size = 0.4
# )

# dims_to_use <- 30
# pooled <- FindNeighbors(pooled, dims = 1:dims_to_use)
# pooled <- FindClusters(pooled, resolution = 0.8)
# pooled <- RunUMAP(pooled, dims = 1:dims_to_use)

# DimPlot(
#     pooled, 
#     reduction = "umap", 
#     label = TRUE
# )



#utilization of post published data, passed qc filtering
published_data = readRDS("Ito_MQ_scRNAseq_sct.rds")
print(published_data)
print(table(published_data$HTO_maxID))

#might need to perfrom additional gene filtering for genes expressed in at >= 3 cells to support relevance
#dimensionality reduction to get a resulting 19k x 4 matrix
#either reduce by average of each gene quanity in cell per tiepoint or reduce to percentage increase

barcodes_pooled <- Cells(pooled)
barcodes_published <- Cells(published_data)

sub1 <- gsub("^P[0-9]_", "", barcodes_pooled)
sub2 <- gsub("^.*_([A-Z]+)-1$", "\\1", barcodes_published)
print(paste("overlapping barcodes from manual analysis vs. paper filter:", length(intersect(sub1, sub2))))


annotations_p1 <- read.csv("GSE213902_MQpool1_filtered_contig_annotations.csv") %>%
    select(barcode, cdr3)
annotations_p2 <- read.csv("GSE213902_MQpool2_filtered_contig_annotations.csv") %>%
    select(barcode, cdr3)

annotations_p1$barcode <- paste0("MQpool1_RS-03742653_", annotations_p1$barcode)
annotations_p2$barcode <- paste0("MQpool2_RS-03742653_", annotations_p2$barcode)

merged_annotations <- bind_rows(annotations_p1, annotations_p2) %>%
    distinct(barcode, .keep_all = TRUE)
common_cells <- intersect(Cells(published_data), merged_annotations$barcode)
print(paste("Match found for", length(common_cells), "cells"))

published_data$clonotype <- merged_annotations$cdr3[match(Cells(published_data), merged_annotations$barcode)]

clones <- sort(table(published_data$clonotype), decreasing = TRUE)
print(head(clones, 10))

total_tcr_barcodes <- nrow(merged_annotations)
total_tcr_barcodes

dim_reduction <- AverageExpression(published_data, group.by = "HTO_maxID", assays = "RNA", slot = "data", return.seurat = TRUE)

#print(head(dim_reduction[["RNA"]]$counts))

dim_reduction <- NormalizeData(dim_reduction, normalization.method = "LogNormalize", scale.factor = 10000)
final <- LayerData(dim_reduction, assay = "RNA", layer = "data")
final <- as.matrix(final)

#write.table(final, file = "dense_dataLT.csv", sep = ",", row.names = FALSE, col.names = FALSE)