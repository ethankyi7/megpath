library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)


#utilization of post published data, passed qc filtering
published_data = readRDS("Ito_MQ_scRNAseq_sct.rds")
print(published_data)
print(table(published_data$HTO_maxID))

#might need to perfrom additional gene filtering for genes expressed in at >= 3 cells to support relevance
#dimensionality reduction to get a resulting 19k x 4 matrix
#either reduce by average of each gene quanity in cell per tiepoint or reduce to percentage increase

# barcodes_pooled <- Cells(pooled)
# barcodes_published <- Cells(published_data)

# sub1 <- gsub("^P[0-9]_", "", barcodes_pooled)
# sub2 <- gsub("^.*_([A-Z]+)-1$", "\\1", barcodes_published)
# print(paste("overlapping barcodes from manual analysis vs. paper filter:", length(intersect(sub1, sub2))))


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


tcr_list <- list(
  "1" = c(1.45, 12.0, 5.21, 6.51),
  "2" = c(0.0371, 0.803, 1.22, 1.78),
  "3" = c(0.0405, 0, 0.0369, 0.0353),
  "4" = c(0.00913, 0.638, 0.0178, 0.0243)
)

# Loop through each cluster and process
for (cl in names(tcr_list)) {
  cat("\n--- Processing Cluster", cl, "---\n")
  
  cluster_sub <- subset(published_data, subset = cluster == cl)
  
  avg_exp <- AverageExpression(cluster_sub, group.by = "HTO_maxID", assays = "RNA", slot = "data")$RNA

  avg_exp_log2 <- log2(avg_exp + 1)
  
  hashtags <- c("HashTag1", "HashTag2", "HashTag3", "HashTag4")
  lfc_cluster <- avg_exp_log2[, hashtags] - avg_exp_log2[, "HashTag1"]
  
  lfc_cluster[is.na(lfc_cluster) | is.infinite(lfc_cluster)] <- 0
  
  g_min <- min(lfc_cluster)
  g_max <- max(lfc_cluster)
  lfc_scaled <- (lfc_cluster - g_min) / (g_max - g_min)
  
  raw_tcr <- tcr_list[[cl]]
  scaled_tcr <- (raw_tcr - min(raw_tcr)) / (max(raw_tcr) - min(raw_tcr))

  cat("Scaled TCR values: ", paste(round(scaled_tcr, 5), collapse=", "), "\n")
  
  #write.table(lfc_scaled, sep=",", paste0("cluster_", cl, "_MM.csv"), row.names = FALSE, col.names = FALSE)
  #write.table(lfc_scaled, sep=",", paste0("cluster_", cl, "_MM_annotated.csv"), row.names = TRUE, col.names = FALSE)
}

# print("<---QUANTILES--->")
# quantile(lfc_cluster1, probs = c(0.01, 0.025, 0.975, 0.99))
# quantile(lfc_cluster2, probs = c(0.01, 0.025, 0.975, 0.99))
# quantile(lfc_cluster3, probs = c(0.01, 0.025, 0.975, 0.99))
# quantile(lfc_cluster4, probs = c(0.01, 0.025, 0.975, 0.99))
# print("<--------------->")

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