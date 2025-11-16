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

poolOnePBMC <- CreateSeuratObject(counts = poolOnePBMC.umis)
poolTwoPBMC <- CreateSeuratObject(counts = poolTwoPBMC.umis)

poolOnePBMC[["HTO"]] <- CreateAssayObject(counts = poolOnePBMC.htos)
poolOnePBMC[["percent.mt"]] <- PercentageFeatureSet(poolOnePBMC, pattern = "^MT-")
poolOnePBMC_filtered <- subset(
    poolOnePBMC,
    subset = nFeature_RNA >= 300 & nFeature_RNA < 7500 & percent.mt < 7.5
)

poolTwoPBMC[["HTO"]] <- CreateAssayObject(counts = poolTwoPBMC.htos)
poolTwoPBMC[["percent.mt"]] <- PercentageFeatureSet(poolTwoPBMC, pattern = "^MT-")
poolTwoPBMC_filtered <- subset(
    poolTwoPBMC,
    subset = nFeature_RNA >= 300 & nFeature_RNA < 7500 & percent.mt < 7.5
)

poolOnePBMC_filtered <- NormalizeData(poolOnePBMC_filtered, assay = "HTO", normalization.method = "CLR")
poolTwoPBMC_filtered <- NormalizeData(poolTwoPBMC_filtered, assay = "HTO", normalization.method = "CLR")

poolOnePBMC_filtered <- HTODemux(
    poolOnePBMC_filtered,
    assay = "HTO",
    positive.quantile = 0.99,
    kfunc = 'clara'
)

poolTwoPBMC_filtered <- HTODemux(
    poolTwoPBMC_filtered,
    assay = "HTO",
    positive.quantile = 0.99,
    kfunc = 'clara'
)

print("Pool One")
table(poolOnePBMC_filtered$HTO_classification.global)
table(poolOnePBMC_filtered$hash.ID)

print("Pool Two")
table(poolTwoPBMC_filtered$HTO_classification.global)
table(poolTwoPBMC_filtered$hash.ID)

use_virtualenv(Sys.getenv("ENV_PATH"), required = TRUE)
scr <- import("scrublet")

doublet_cells <- subset(poolOnePBMC_filtered, subset = HTO_classification.global == "Doublet")
doublet_cells1 <- subset(poolTwoPBMC_filtered, subset = HTO_classification.global == "Doublet")

poolOnePBMC_filtered <- JoinLayers(poolOnePBMC_filtered, assay = "RNA")
counts_matrix <- LayerData(poolOnePBMC_filtered, assay = "RNA", layer = "counts")
counts_matrix_t <- t(as.matrix(counts_matrix))

python_counts_matrix <- r_to_py(counts_matrix_t)

scrub <- scr$Scrublet(
    python_counts_matrix,
    expected_doublet_rate = 0.06
)

results <- scrub$scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30L)
doublet_scores <- py_to_r(results[[1]]) 
predicted_doublets <- py_to_r(results[[2]])

poolOnePBMC_filtered <- AddMetaData(
    poolOnePBMC_filtered, 
    metadata = doublet_scores, 
    col.name = "scrublet_score"
)

final_poolOne <- subset(poolOnePBMC_filtered, subset = scrublet_score < 0.2)
cells_removed_scrublet <- ncol(poolOnePBMC_filtered) - ncol(final_poolOne)

print("Pool One")
print(paste("Cells before Scrublet filter:", ncol(poolOnePBMC_filtered)))
print(paste("Cells removed by Scrublet (score >= 0.2):", cells_removed_scrublet))
print(paste("Final high-quality cells remaining:", ncol(final_poolOne)))


poolTwoPBMC_filtered <- JoinLayers(poolTwoPBMC_filtered, assay = "RNA")
count_matrix1 <- LayerData(poolTwoPBMC_filtered, assay = "RNA", layer = "counts")
count_matrix1_t <- t(as.matrix(count_matrix1))

python_counts_matrix1 <- r_to_py(count_matrix1_t)
scrub1 <- scr$Scrublet(
    python_counts_matrix1,
    expected_doublet_rate = 0.06
)

results1 <- scrub1$scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 85, n_prin_comps = 30L)
doublet_scores <- py_to_r(results1[[1]]) 
predicted_doublets <- py_to_r(results1[[2]])

poolTwoPBMC_filtered <- AddMetaData(
    poolTwoPBMC_filtered, 
    metadata = doublet_scores, 
    col.name = "scrublet_score"
)

final_poolTwo <- subset(poolTwoPBMC_filtered, subset = scrublet_score < 0.2)
cells_removed_scrublet <- ncol(poolTwoPBMC_filtered) - ncol(final_poolTwo)

print("Pool Two")
print(paste("Cells before Scrublet filter:", ncol(poolTwoPBMC_filtered)))
print(paste("Cells removed by Scrublet (score >= 0.2):", cells_removed_scrublet))
print(paste("Final high-quality cells remaining:", ncol(final_poolTwo)))


final_poolOne$timepoint <- final_poolOne$HTO_maxID
final_poolTwo$timepoint <- final_poolTwo$HTO_maxID

tag_to_sample_map <- c(
    "HTO1-GTCAACTCTTTAGCG" = "W0",
    "HTO2-TGATGGCCTATTGGG" = "W3",
    "HTO3-TTCCGCCTCTCTTTG" = "W6",
    "HTO4-AGTAAGTTCAGCGTA" = "W9"
)

final_poolOne$timepoint <- factor(
    final_poolOne$timepoint, 
    levels = names(tag_to_sample_map), 
    labels = unname(tag_to_sample_map)
)

print("pool One")
print(table(final_poolOne$timepoint))


final_poolTwo$timepoint <- factor(
    final_poolTwo$timepoint,
    levels = names(tag_to_sample_map),
    labels = unname(tag_to_sample_map)
)
print("pool Two")
print(table(final_poolTwo$timepoint))

print(nrow(LayerData(final_poolOne, assay = "RNA", layer = "counts")))
print(nrow(LayerData(final_poolTwo, assay = "RNA", layer = "counts")))

final_poolOne <- RenameCells(final_poolOne, add.cell.id = c("P1"))
final_poolTwo <- RenameCells(final_poolTwo, add.cell.id = c("P2"))

merged <- merge(x = final_poolOne, y = final_poolTwo)
merged <- JoinLayers(merged, assay = "RNA")

genes.use <- rownames(merged)[rowSums(GetAssayData(merged, assay = "RNA", layer = "counts") > 0) >= 3]
merged <- subset(merged, features = genes.use)
print(table(merged$timepoint))
print(nrow(LayerData(merged, assay = "RNA", layer = "counts")))
print(table(merged$HTO_classification.global))

# merged <- NormalizeData(merged)
# merged <- ScaleData(merged)
# merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 10000)
# merged <- RunPCA(merged, features = VariableFeatures(object = merged))
# merged <- RunUMAP(merged, dims = 1:10)
# DimPlot(merged, reduction = "umap", group.by = "HTO_classification.global")