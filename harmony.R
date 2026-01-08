
#Redo the human integration with harmony without SCT

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(presto)

setwd("/r_workspace/cristal_data/tabula_sapiens_singlecell/harmony_integration")



# -----------------------------
# 1) Load your human RDS objects
# -----------------------------
base_dir <- "/projects/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples"

ts14_ant_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts14_ant_qc_sct30k.rds")
ts14_post_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts14_post_qc_sct30k.rds")
ts25_ant_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts25_ant_qc_sct30k.rds")
ts25_post_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts25_post_qc_sct30k.rds")
ts27_ant_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts27_ant_qc_sct30k.rds")
ts27_post_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts27_post_qc_sct30k.rds")
ts7_ant_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts7_ant_qc_sct30k.rds")
ts7_post_qc_sct30k <- readRDS("/r_workspace/cristal_data/tabula_sapiens_singlecell/rds_samples/ts7_post_qc_sct30k.rds")


objs <- list(
  ts7_ant  = ts7_ant_qc_sct30k,
  ts7_post = ts7_post_qc_sct30k,
  ts14_ant  = ts14_ant_qc_sct30k,
  ts14_post = ts14_post_qc_sct30k,
  ts25_ant  = ts25_ant_qc_sct30k,
  ts25_post = ts25_post_qc_sct30k,
  ts27_ant  = ts27_ant_qc_sct30k,
  ts27_post = ts27_post_qc_sct30k
)

objs <- lapply(names(objs), function(nm) {
  obj <- objs[[nm]]
  DefaultAssay(obj) <- "RNA"
  obj$sample_id <- nm
  obj
})
names(objs) <- names(objs)
# -----------------------------
# 2) Merge
# -----------------------------
combined <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = names(objs),
  project = "human_ts_harmony_RNA"
)

DefaultAssay(combined) <- "RNA"

# -----------------------------
# 3) Standard RNA workflow (no SCT)
# -----------------------------
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 8000, verbose = FALSE)

# If you have huge cell counts, you can restrict scaling to HVGs:
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)

combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 50, verbose = FALSE)

# -----------------------------
# 4) Harmony (batch = sample_id)
# -----------------------------
combined <- IntegrateLayers(
  object = combined,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by = "sample_id"   # <- batch variable
)

# -----------------------------
# 5) UMAP/Neighbors/Clusters on Harmony embedding
# -----------------------------
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:30
)

#Run clustering at multiple resolutions
resolutions <- seq(0.2, 0.6, by = 0.1)

for (r in resolutions) {
  combined <- FindClusters(
    combined,
    resolution = r,
    algorithm = 1,   # Louvain (default, stable)
    verbose = FALSE
  )
  
  # Rename the active cluster column so it’s preserved
  colnames(combined@meta.data)[
    ncol(combined@meta.data)
  ] <- paste0("harmony_res_", r)
}

library(patchwork)

plots <- lapply(resolutions, function(r) {
  DimPlot(
    combined,
    reduction = "umap",
    group.by = paste0("harmony_res_", r),
    label = TRUE,
    repel = TRUE,
    pt.size = 0.05
  ) + ggtitle(paste("Resolution", r))
})

wrap_plots(plots)


resolutions <- seq(0.2, 0.6, by = 0.1)
cluster_cols <- paste0("harmony_res_", resolutions)

# sanity check
cluster_cols %in% colnames(combined@meta.data)
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}
library(dplyr)

jaccard_results <- list()

for (i in 1:(length(cluster_cols) - 1)) {
  r1 <- cluster_cols[i]
  r2 <- cluster_cols[i + 1]
  
  clust1 <- combined@meta.data[[r1]]
  clust2 <- combined@meta.data[[r2]]
  
  for (c1 in unique(clust1)) {
    cells1 <- WhichCells(combined, expression = get(r1) == c1)
    
    best_j <- 0
    best_c2 <- NA
    
    for (c2 in unique(clust2)) {
      cells2 <- WhichCells(combined, expression = get(r2) == c2)
      j <- jaccard(cells1, cells2)
      
      if (j > best_j) {
        best_j <- j
        best_c2 <- c2
      }
    }
    
    jaccard_results[[length(jaccard_results) + 1]] <-
      data.frame(
        res_from = r1,
        res_to   = r2,
        cluster_from = c1,
        best_match_cluster = best_c2,
        jaccard = best_j
      )
  }
}

jaccard_df <- bind_rows(jaccard_results)
jaccard_summary <- jaccard_df %>%
  group_by(res_from) %>%
  summarise(
    mean_jaccard = mean(jaccard),
    median_jaccard = median(jaccard),
    min_jaccard = min(jaccard)
  )

jaccard_summary

#####################################
#########Interpretation#############
#####################################
# Mean/median ≥0.6 → stable resolution
# Sharp drop as resolution increases → over-splitting
#########################################################





#Optional: Set one resolution as the default
Idents(combined) <- combined$harmony_res_0.2


# -----------------------------
# 6) Quick plots
# -----------------------------
p1 <- DimPlot(combined, reduction = "umap", group.by = "sample_id", pt.size = 0.05) + ggtitle("Harmony (RNA) by sample_id")
p2 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.05) + ggtitle("Harmony (RNA) clusters")
p3 <- DimPlot(combined, reduction = "umap", group.by = "region", pt.size = 0.05) + ggtitle("Harmony (RNA) by region")
p4 <- DimPlot(combined, reduction = "umap", group.by = "donor", pt.size = 0.05) + ggtitle("Harmony (RNA) by donor")

print(p1); print(p2); print(p3); print(p4)


# -----------------------------
# 7) Find markers
# -----------------------------

DefaultAssay(combined) <- "RNA"
Idents(combined) <- combined$harmony_res_0.2
combined <- JoinLayers(combined, assay = "RNA")
markers_all <- FindAllMarkers(
  object = combined,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

write.csv(markers_all, "markers_FindAllMarkers_res0.2.csv", row.names = FALSE)


DefaultAssay(combined) <- "RNA"
Idents(combined) <- combined$harmony_res_0.5
combined <- JoinLayers(combined, assay = "RNA")
markers_all <- FindAllMarkers(
  object = combined,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)

write.csv(markers_all, "markers_FindAllMarkers_res0.5.csv", row.names = FALSE)


#Markers res=0.5

library(dplyr)
top10 <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top10, "markers_top10_res0.2.csv", row.names = FALSE)

markers_all_clean <- markers_all %>%
  filter(!grepl("^MT-", gene)) %>%
  filter(!grepl("^RPL|^RPS", gene))

top10_clean <- markers_all_clean %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(markers_all_clean, "markers_FindAllMarkers_res0.2_noMTRibo.csv", row.names = FALSE)
write.csv(top10_clean, "markers_top10_res0.5_noMTRibo.csv", row.names = FALSE)
DotPlot(combined, features = unique(top10_clean$gene)) + RotatedAxis()
DimPlot(combined, label = TRUE)



# Make sure identities are seurat clusters
Idents(combined) <- "RNA_snn_res.0.5"

# Define mapping (example — adjust names as needed)
new.cluster.ids <- c(
  "0"  = "Epithelial",
  "1"  = "Epithelial",
  "2"  = "Epithelial",
  "3"  = "Immune",
  "4"  = "Epithelial",
  "5"  = "Epithelial",
  "6"  = "Endothelial",
  "7"  = "Immune",
  "8"  = "Fibroblast",
  "9"  = "Immune",
  "10"  = "Doublet",
  "11"  = "Epithelial",
  "12"  = "Epithelial",
  "13"  = "Epithelial",
  "14"  = "Epithelial",
  "15"  = "Doublet",
  "16"  = "Pericyte",
  "17"  = "Epithelial",
  "18"  = "Schwann cells",
  "19"  = "Immune",
  "20"  = "Immune",
  "21"  = "Immune",
  "22"  = "Endothelial",
  "23"  = "Doublet",
  "24"  = "Doublet",
  "25"  = "Doublet"

)

# Rename
combined <- RenameIdents(combined, new.cluster.ids)
# replace 'celltype_clean' with your column name
combined <- subset(combined, subset = cell_type_harmonized != "Doublet")

# (Optional but recommended) store as metadata
combined$cell_type_harmonized <- Idents(combined)
DimPlot(combined, group.by = "cell_type_harmonized")
DimPlot(combined, group.by = "cell_type_harmonized", split.by = "region")
DimPlot(combined, group.by = "sample_id")

# highlight cluster 12
DimPlot(combined, cells.highlight = WhichCells(combined, idents = "12"),
        sizes.highlight = 0.2) + NoLegend()

DimPlot(combined, cells.highlight = WhichCells(combined, idents = "13"),
        sizes.highlight = 0.2) + NoLegend()


DimPlot(combined, cells.highlight = WhichCells(combined, idents = "14"),
        sizes.highlight = 0.2) + NoLegend()

DimPlot(combined, cells.highlight = WhichCells(combined, idents = "19"),
        sizes.highlight = 0.2) + NoLegend()

DimPlot(combined, cells.highlight = WhichCells(combined, idents = "25"),
        sizes.highlight = 0.2) + NoLegend()

# markers to check
epi <- c("KRT14","KRT5","TP63","KRT15","COL17A1","EPCAM")
neu <- c("TUBB3","SYT1","SNAP25","RBFOX3","PRPH")
immune <- c("PTPRC","LST1","TYROBP", "MRC1", "CD3D")
endo <- c("PECAM1","VWF","KDR", "NOTCH3", "MYH11", "ACTA2")
fibro <- c("COL1A1","DCN","LUM", "MGP", "LAMA2")
schwann <- c("SOX10","S100B","MPZ", "PLP1", "MBP")

FeaturePlot(combined, features = c(epi, neu, immune), ncol = 4, pt.size = 0.05)
FeaturePlot(combined, features = c(endo, fibro, schwann), ncol = 4, pt.size = 0.05)

VlnPlot(combined, features = c("nCount_RNA","nFeature_RNA"),  pt.size = 0)


# co-expression check (very informative)
FeatureScatter(combined, feature1 = "KRT14", feature2 = "TUBB3")
FeatureScatter(combined, feature1 = "TP63", feature2 = "SYT1")




# Rename
combined <- RenameIdents(combined, new.cluster.ids)

# (Optional but recommended) store as metadata
combined$human_major_type <- Idents(combined)
DimPlot(combined, group.by = "human_major_type")




# Table: clusters x sample
cluster_sample_counts <- table(
  Cluster = Idents(combined),
  Sample  = combined$sample_id
)

# Convert to data.frame
cluster_sample_counts_df <- as.data.frame(cluster_sample_counts)

head(cluster_sample_counts_df)

library(pheatmap)

# Convert to matrix
prop_mat <- xtabs(Freq ~ Cluster + Sample, data = cluster_sample_props_df)

pheatmap(
  prop_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = colorRampPalette(c("white","orange","red"))(100),
  main = "Sample contribution per cluster"
)

library(ggplot2)

ggplot(cluster_sample_counts_df,
       aes(x = Cluster, y = Freq, fill = Sample)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(
    title = "Cell counts per cluster by sample",
    y = "Number of cells"
  )

ggplot(cluster_sample_props_df,
       aes(x = Cluster, y = Freq, fill = Sample)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  labs(
    title = "Sample contribution per cluster",
    y = "Proportion of cells"
  )


##clusters where >80% of cells come from a single sample
library(dplyr)

sample_specific_clusters <- cluster_sample_props_df %>%
  group_by(Cluster) %>%
  filter(Freq > 0.8) %>%
  arrange(desc(Freq))

sample_specific_clusters

DimPlot(
  combined,
  reduction = "umap",
  split.by = "sample_id",
  group.by = "RNA_snn_res.0.5",
  pt.size = 0.2
)

##Perform deg analysis to subtype
Idents(combined) <- "cell_type_harmonized"


# -----------------------------
# 8) Save
# -----------------------------
saveRDS(combined,  "human_integrated_majortypes.rds")
