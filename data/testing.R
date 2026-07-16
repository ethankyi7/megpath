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

# # writeMM(counts, "counts.mtx")
# # write.table(rownames(counts), "genes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)
# # write.table(colnames(counts), "barcodes.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)

# scrub <- read.csv("scrublet_results.csv")
# #print(scrub[match(colnames(pooled), scrub$barcode), ])
# #print(scrub$doublet_score)
# pooled[["doublet_score"]] <- scrub$doublet_score
# #summary(pooled$doublet_score)

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

clonotype_map <- read.csv("clonotype_map.csv")

annotations_p1$barcode <- paste0("MQpool1_RS-03742653_", annotations_p1$barcode)
annotations_p2$barcode <- paste0("MQpool2_RS-03742653_", annotations_p2$barcode)

merged_annotations <- bind_rows(annotations_p1, annotations_p2) %>%
    distinct(barcode, .keep_all = TRUE)
common_cells <- intersect(Cells(published_data), merged_annotations$barcode)
print(paste("Match found for", length(common_cells), "cells"))

published_data$clonotype <- merged_annotations$cdr3[match(Cells(published_data), merged_annotations$barcode)]

published_data$cluster <- clonotype_map$cluster[match(published_data$clonotype, clonotype_map$clonotype)]
gene_list <- rownames(published_data[["RNA"]])
gene_list
write.table(gene_list, file = "gene_list.csv", sep = ",", row.names = FALSE, col.names = FALSE)

clones <- sort(table(published_data$clonotype), decreasing = TRUE)
print(table(published_data$cluster, published_data$HTO_maxID))


cell_data <- data.frame(
  Cell_Barcode = rownames(published_data@meta.data),
  Cluster = published_data$cluster # or my_seurat_object$cluster depending on your exact column name
)

cell_data

#filtering out housekeeping genes
hk_url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HOUNKPE_HOUSEKEEPING_GENES&fileType=grp"
hounkpe <- read.table(hk_url, skip = 2, header = FALSE)$V1
ribo_cyto <- grep("^RPL|^RPS", rownames(published_data), value = TRUE)
mito      <- grep("^MT-", rownames(published_data), value = TRUE)

all_hk <- union(hounkpe, union(ribo_cyto, mito))

cat("Hounkpe genes in data:          ", length(intersect(rownames(published_data), hounkpe)), "\n")
cat("Cytosolic ribo genes added:     ", length(intersect(rownames(published_data), ribo_cyto)), "\n")
cat("MT genes added:                 ", length(intersect(rownames(published_data), mito)), "\n")
cat("Total genes removed:            ", length(intersect(rownames(published_data), all_hk)), "\n")

published_data <- published_data[!rownames(published_data) %in% all_hk, ]
cat("Genes remaining for NMF:        ", nrow(published_data), "\n")
# write.table(row.names(published_data), "gene_list_filtered.csv", row.names=FALSE, col.names=FALSE)

#<-------- CLUSTER STUFF ---------->
cluster1 <- subset(published_data, subset = cluster == 1)
cluster2 <- subset(published_data, subset = cluster == 2)
cluster3 <- subset(published_data, subset = cluster == 3)
cluster4 <- subset(published_data, subset = cluster == 4)


#<-------- CLUSTER ONE ------------>
avg_exp1 <- AverageExpression(cluster1, group.by = "HTO_maxID", assays = "RNA", slot = "data")$RNA
print("range", paste(range(avg_exp1)))
lfc_cluster1 <- (avg_exp1[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp1[, "HashTag1"]) / log(2)

lfc_cluster1[lfc_cluster1 > 5] <- 5
lfc_cluster1[lfc_cluster1 < -5] <- -5
lfc_cluster1[is.na(lfc_cluster1) | is.infinite(lfc_cluster1)] <- 0

g_min <- min(lfc_cluster1)
g_max <- max(lfc_cluster1)
lfc_scaled <- (lfc_cluster1 - g_min) / (g_max - g_min)
#write.table(lfc_cluster1, "cluster_oneLFC.csv", sep=",", row.names=TRUE, col.names=FALSE)
# write.table(lfc_scaled, "cluster_oneMM.csv", sep=",", row.names=FALSE, col.names=FALSE)

raw_tcr <- c(1.45, 12.0, 5.21, 6.51)
lfc_tcr <- log2(raw_tcr / raw_tcr[1]) # [0, 3.05, 1.85, 2.17]
print(paste(lfc_tcr))
scaled_tcr <- (lfc_tcr - g_min) / (g_max - g_min)
scaled_tcr

#<-------- CLUSTER TWO ------------>

avg_exp2 <- AverageExpression(cluster2, group.by = "HTO_maxID", assays = "RNA", slot = "data")$RNA
lfc_cluster2 <- (avg_exp2[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp2[, "HashTag1"]) / log(2)

lfc_cluster2[lfc_cluster2 > 5] <- 5
lfc_cluster2[lfc_cluster2 < -5] <- -5
lfc_cluster2[is.na(lfc_cluster2) | is.infinite(lfc_cluster2)] <- 0

g_min <- min(lfc_cluster2)
g_max <- max(lfc_cluster2)
lfc_scaled <- (lfc_cluster2 - g_min) / (g_max - g_min)
#write.table(lfc_cluster2, "cluster_twoLFC.csv", sep=",", row.names=TRUE, col.names=FALSE)
# write.table(lfc_scaled, "cluster_twoMM.csv", sep=",", row.names=FALSE, col.names=FALSE)

raw_tcr <- c(0.371,0.803,1.22,1.78)
lfc_tcr <- log2(raw_tcr / raw_tcr[1])
print(paste(lfc_tcr))
scaled_tcr <- (lfc_tcr - g_min) / (g_max - g_min)
scaled_tcr

#<-------- CLUSTER THREE ---------->

avg_exp3 <- AverageExpression(cluster3, group.by = "HTO_maxID", assays = "RNA", slot = "data")$RNA
lfc_cluster3 <- (avg_exp3[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp3[, "HashTag1"]) / log(2)

lfc_cluster3[lfc_cluster3 > 5] <- 5
lfc_cluster3[lfc_cluster3 < -5] <- -5
lfc_cluster3[is.na(lfc_cluster3) | is.infinite(lfc_cluster3)] <- 0

g_min <- min(lfc_cluster3)
g_max <- max(lfc_cluster3)
lfc_scaled <- (lfc_cluster3 - g_min) / (g_max - g_min)
#write.table(lfc_cluster3, "cluster_threeLFC.csv", sep=",", row.names=TRUE, col.names=FALSE)
# write.table(lfc_scaled, "cluster_threeMM.csv", sep=",", row.names=FALSE, col.names=FALSE)

raw_tcr <- c(0.0405,0,0.0369,0.0353)
lfc_tcr <- log2(raw_tcr / raw_tcr[1])
print(paste(lfc_tcr))
scaled_tcr <- (lfc_tcr - g_min) / (g_max - g_min)
scaled_tcr

#<-------- CLUSTER FOUR ----------->

avg_exp4 <- AverageExpression(cluster4, group.by = "HTO_maxID", assays = "RNA", slot = "data")$RNA
lfc_cluster4 <- (avg_exp4[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp4[, "HashTag1"]) / log(2)

lfc_cluster4[lfc_cluster4 > 5] <- 5
lfc_cluster4[lfc_cluster4 < -5] <- -5
lfc_cluster4[is.na(lfc_cluster4) | is.infinite(lfc_cluster4)] <- 0

g_min <- min(lfc_cluster4)
g_max <- max(lfc_cluster4)
lfc_scaled <- (lfc_cluster4 - g_min) / (g_max - g_min)
#write.table(lfc_cluster4, "cluster_fourLFC.csv", sep=",", row.names=TRUE, col.names=FALSE)
# write.table(lfc_scaled, "cluster_fourMM.csv", sep=",", row.names=FALSE, col.names=FALSE)

raw_tcr <- c(0.00913,0.638,0.0178,0.0243)
lfc_tcr <- log2(raw_tcr / raw_tcr[1])
print(paste(lfc_tcr))
scaled_tcr <- (lfc_tcr - g_min) / (g_max - g_min)
scaled_tcr

print("<---QUANTILES--->")
quantile(lfc_cluster1, probs = c(0.01, 0.025, 0.975, 0.99))
quantile(lfc_cluster2, probs = c(0.01, 0.025, 0.975, 0.99))
quantile(lfc_cluster3, probs = c(0.01, 0.025, 0.975, 0.99))
quantile(lfc_cluster4, probs = c(0.01, 0.025, 0.975, 0.99))
print("<--------------->")

# CDF Pipeline
lfc1 <- (avg_exp1[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp1[, "HashTag1"]) / log(2)
lfc2 <- (avg_exp2[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp1[, "HashTag1"]) / log(2)
lfc3 <- (avg_exp3[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp1[, "HashTag1"]) / log(2)
lfc4 <- (avg_exp1[, c("HashTag1", "HashTag2", "HashTag3", "HashTag4")] - avg_exp1[, "HashTag1"]) / log(2)

normalize_cluster <- function(lfc_cluster) {
    lfc_cluster <- as.matrix(lfc_cluster)
    avg <- mean(lfc_cluster[is.finite(lfc_cluster)])
    sdev <- sd(lfc_cluster[is.finite(lfc_cluster)])

    out <- matrix(
        pnorm(lfc_cluster, mean=avg, sd=sdev),
        nrow = nrow(lfc_cluster),
        dimnames = dimnames(lfc_cluster)
    )
    return(out)
}

normalize_tcr <- function(raw_tcr, lfc_cluster) {
  avg <- mean(lfc_cluster[is.finite(lfc_cluster)])
  sdev <- sd(lfc_cluster[is.finite(lfc_cluster)])
  
  lfc_tcr <- log2(raw_tcr / raw_tcr[1])
  lfc_tcr[is.infinite(lfc_tcr)] <- min(lfc_cluster[is.finite(lfc_cluster)])
  pnorm(lfc_tcr, mean = avg, sd = sdev)
}

scaled_tcr1 <- normalize_tcr(c(1.45, 12.0, 5.21, 6.51),   lfc1)
scaled_tcr2 <- normalize_tcr(c(0.371, 0.803, 1.22, 1.78),  lfc2)
scaled_tcr3 <- normalize_tcr(c(0.0405, 0, 0.0369, 0.0353), lfc3)
scaled_tcr4 <- normalize_tcr(c(0.00913, 0.638, 0.0178, 0.0243), lfc4)

lfc_cdf1 <- normalize_cluster(lfc1)
lfc_cdf2 <- normalize_cluster(lfc2)
lfc_cdf3 <- normalize_cluster(lfc3)
lfc_cdf4 <- normalize_cluster(lfc4)

print("<-------->")
scaled_tcr1
scaled_tcr2
scaled_tcr3
scaled_tcr4

# write.table(lfc_cdf1, "clusterOneCDF_annotated.csv", sep=',', row.names=TRUE, col.names=FALSE)
# write.table(lfc_cdf2, "clusterTwoCDF_annotated.csv", sep=',', row.names=TRUE, col.names=FALSE)
# write.table(lfc_cdf3, "clusterThreeCDF_annotated.csv", sep=',', row.names=TRUE, col.names=FALSE)
# write.table(lfc_cdf4, "clusterFourCDF_annotated.csv", sep=',', row.names=TRUE, col.names=FALSE)

#inverse logit pipeline


dim_reduction <- AverageExpression(published_data, group.by = "HTO_maxID", assays = "RNA", slot = "counts", return.seurat = TRUE)

#print(head(dim_reduction[["RNA"]]$counts))

dim_reduction <- NormalizeData(dim_reduction, normalization.method = "LogNormalize", scale.factor = 10000)
final <- LayerData(dim_reduction, assay = "RNA", layer = "counts")
final <- as.matrix(final)

#<------ raw averaged counts for each cluster ----->
counts1 <- AverageExpression(cluster1, group.by = "HTO_maxID", assays = "RNA", layer = "counts", return.seurat=TRUE)
counts_matrix1 <- LayerData(counts1, assay = "RNA", layer = "data")

counts2 <- AverageExpression(cluster2, group.by = "HTO_maxID", assays = "RNA", layer = "counts", return.seurat=TRUE)
counts_matrix2 <- LayerData(counts2, assay = "RNA", layer = "data")

counts3 <- AverageExpression(cluster3, group.by = "HTO_maxID", assays = "RNA", layer = "counts", return.seurat=TRUE)
counts_matrix3 <- LayerData(counts3, assay = "RNA", layer = "data")

counts4 <- AverageExpression(cluster4, group.by = "HTO_maxID", assays = "RNA", layer = "counts", return.seurat=TRUE)
counts_matrix4 <- LayerData(counts4, assay = "RNA", layer = "data")

#<--- EXPRESSIONS --->
print(table(cluster1$HTO_maxID))
gene_sum <- sum(GetAssayData(cluster1, layer = "counts")["FAM41C", ])
gene_sum
raw_avg1 <- AverageExpression(cluster1, group.by = "HTO_maxID", assays = "RNA", layer = "counts")$RNA
raw_avg2 <- AverageExpression(cluster2, group.by = "HTO_maxID", assays = "RNA", layer = "counts")$RNA
raw_avg3 <- AverageExpression(cluster3, group.by = "HTO_maxID", assays = "RNA", layer = "counts")$RNA
raw_avg4 <- AverageExpression(cluster4, group.by = "HTO_maxID", assays = "RNA", layer = "counts")$RNA


#write.table(final, file = "dense_dataLT.csv", sep = ",", row.names = FALSE, col.names = FALSE)
# write.table(as.matrix(counts_matrix1), file = "clusterOne_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(counts_matrix2), file = "clusterTwo_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(counts_matrix3), file = "clusterThree_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(counts_matrix4), file = "clusterFour_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)

# write.table(as.matrix(raw_avg1), file = "clusterOne_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(raw_avg2), file = "clusterTwo_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(raw_avg3), file = "clusterThree_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)
# write.table(as.matrix(raw_avg4), file = "clusterFour_annotated.csv", sep = ',', row.names = TRUE, col.names = FALSE)