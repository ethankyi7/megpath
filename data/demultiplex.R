library(Seurat)
library(reticulate)

dataPath <- "genePool1/"
poolOnePBMC.umis <- ReadMtx(
    mtx = "genePool1/GSE213902_MQpool1_matrix.mtx.gz",
    cells = "genePool1/GSE213902_MQpool1_barcodes.tsv.gz",
    features = "genePool1/GSE213902_MQpool1_features.tsv.gz"
)

poolOnePBMC.htos <- ReadMtx(
    mtx =  "sra/poolOne/umi_count/matrix.mtx.gz",
    cells = "sra/poolOne/umi_count/barcodes.tsv.gz",
    features = "sra/poolOne/umi_count/features.tsv.gz",
    feature.column = 1
)

poolTwoPBMC.umis <- ReadMtx(
    mtx = "genePool2/GSE213902_MQpool2_matrix.mtx.gz",
    cells = "genePool2/GSE213902_MQpool2_barcodes.tsv.gz",
    features = "genePool2/GSE213902_MQpool2_features.tsv.gz",
)

poolTwoPBMC.htos <- ReadMtx(
    mtx =  "sra/poolTwo/umi_count/matrix.mtx.gz",
    cells = "sra/poolTwo/umi_count/barcodes.tsv.gz",
    features = "sra/poolTwo/umi_count/features.tsv.gz",
    feature.column = 1
)

colnames(poolOnePBMC.umis) <- gsub("-1$", "", colnames(poolOnePBMC.umis))
colnames(poolTwoPBMC.umis) <- gsub("-1$", "", colnames(poolTwoPBMC.umis))

joint.bcs <- intersect(colnames(poolOnePBMC.umis), colnames(poolOnePBMC.htos))
joint2.bcs <- intersect(colnames(poolTwoPBMC.umis), colnames(poolTwoPBMC.htos))

poolOnePBMC.umis <- poolOnePBMC.umis[, joint.bcs]
poolOnePBMC.htos <- as.matrix(poolOnePBMC.htos[, joint.bcs])

poolTwoPBMC.umis <- poolTwoPBMC.umis[, joint2.bcs]
poolTwoPBMC.htos <- as.matrix(poolTwoPBMC.htos[, joint2.bcs])

poolOnePBMC.hashtags <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(poolOnePBMC.umis), sparse = T))
poolTwoPBMC.hashtags <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(poolTwoPBMC.umis), sparse = T))

# poolOnePBMC.umis[["percent.mt"]] <- PercentageFeatureSet(merged_pbmc, pattern = "^MT-")
# poolOnePBMC_filterd <- subset(
#     poolOnePBMC.umis,
#     subset = nFeature_RNA >= 300 & nFeature_RNA < 7500, min_cells >= 3, percent.mt > 7.5
# )


poolOnePBMC.hashtags[["HTO"]] <- CreateAssayObject(counts = poolOnePBMC.htos)
poolTwoPBMC.hashtags[["HTO"]] <- CreateAssayObject(counts = poolTwoPBMC.htos)

poolOnePBMC.hashtags <- RenameCells(poolOnePBMC.hashtags, add.cell.id = c("P1"))
poolTwoPBMC.hashtags <- RenameCells(poolTwoPBMC.hashtags, add.cell.id = c("P2"))

merged_pbmc <- merge(x = poolOnePBMC.hashtags, y = poolTwoPBMC.hashtags)

merged_pbmc <- NormalizeData(merged_pbmc, assay = "HTO", normalization.method = "CLR")

merged_pbmc <- HTODemux(
    merged_pbmc, 
    assay = "HTO", 
    positive.quantile = 0.99,
    kfunc = 'clara'   
)

table(merged_pbmc$HTO_classification.global)
table(merged_pbmc$hash.ID)
#RidgePlot(merged_pbmc, assay = "HTO", features = rownames(merged_pbmc[["HTO"]])[1:2], ncol = 2)

doublet_cells <- subset(merged_pbmc, subset = HTO_classification.global == "Doublet")
max_id <- doublet_cells$HTO_maxID
second_id <- doublet_cells$HTO_secondID

doublet_combinations <- sapply(1:length(max_id), function(i) {
    tags <- sort(c(max_id[i], second_id[i]))
    paste(tags, collapse = " + ")
})

doublet_cells$doublet_origin <- doublet_combinations
table(doublet_cells$doublet_origin)

merged_pbmc[["percent.mt"]] <- PercentageFeatureSet(merged_pbmc, pattern = "^MT-")
VlnPlot(merged_pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

merged_pbmc_filtered <- subset(
    merged_pbmc, 
    subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 7.5
)

table(merged_pbmc_filtered$HTO_classification.global)
#scrublet time
use_virtualenv(Sys.getenv("ENV_PATH"), required = TRUE)
scrublet <- import("scrublet")

merged_pbmc_filtered <- JoinLayers(merged_pbmc_filtered, assay = "RNA")
counts_matrix <- LayerData(merged_pbmc_filtered, assay = "RNA", layer = "counts")
counts_matrix_t <- t(as.matrix(counts_matrix))

python_counts_matrix <- r_to_py(counts_matrix_t)

scrub <- scrublet$Scrublet(
    python_counts_matrix,
    expected_doublet_rate = 0.06
)

results <- scrub$scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30L)
doublet_scores <- py_to_r(results[[1]]) 
predicted_doublets <- py_to_r(results[[2]])

merged_pbmc_filtered <- AddMetaData(
    merged_pbmc_filtered, 
    metadata = doublet_scores, 
    col.name = "scrublet_score"
)

final_cells_filtered <- subset(merged_pbmc_filtered, subset = scrublet_score < 0.2)

cells_removed_scrublet <- ncol(merged_pbmc_filtered) - ncol(final_cells_filtered)

print(paste("Cells before Scrublet filter:", ncol(merged_pbmc_filtered)))
print(paste("Cells removed by Scrublet (score >= 0.2):", cells_removed_scrublet))
print(paste("Final high-quality cells remaining:", ncol(final_cells_filtered)))

# final_cells_filtered <- NormalizeData(final_cells_filtered, assay = "HTO", normalization.method = "CLR")
# final_cells_filtered <- HTODemux(
#     final_cells_filtered, 
#     assay = "HTO", 
#     positive.quantile = 0.99,
#     kfunc = 'clara'
# )

# table(final_cells_filtered$HTO_classification.global)
# table(final_cells_filtered$hash.ID)

final_cells_filtered$timepoint <- final_cells_filtered$HTO_maxID

tag_to_sample_map <- c(
    "HTO1-GTCAACTCTTTAGCG" = "W0",
    "HTO2-TGATGGCCTATTGGG" = "W3",
    "HTO3-TTCCGCCTCTCTTTG" = "W6",
    "HTO4-AGTAAGTTCAGCGTA" = "W9"
)

final_cells_filtered$timepoint <- factor(
    final_cells_filtered$timepoint, 
    levels = names(tag_to_sample_map), 
    labels = unname(tag_to_sample_map)
)

print(table(final_cells_filtered$timepoint))