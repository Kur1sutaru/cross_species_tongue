library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)

species_vec <- c("Human", "Marmoset", "Rat", "Mouse")

imm_list <- Splitimmune_humanect(immune, split.by = "species")
imm_list <- imm_list[names(imm_list) %in% species_vec]
process_species_immune_human <- function(immune_human,
                                nfeatures = 3000,
                                npcs = 30,
                                resolution = 0.5,
                                dims = 1:30) {
  
  DefaultAssay(immune_human) <- "RNA"
  
  immune_human <- NormalizeData(immune_human, verbose = FALSE)
  immune_human <- FindVariableFeatures(immune_human, nfeatures = nfeatures, verbose = FALSE)
  immune_human <- ScaleData(immune_human, verbose = FALSE)
  immune_human <- RunPCA(immune_human, npcs = npcs, verbose = FALSE)
  immune_human <- FindNeighbors(immune_human, dims = dims, verbose = FALSE)
  immune_human <- FindClusters(immune_human, resolution = resolution, verbose = FALSE)
  immune_human <- RunUMAP(immune_human, dims = dims, verbose = FALSE)
  
  immune_human
}



imm_proc_list <- map(
  imm_list,
  ~ process_species_immune_human(.x, nfeatures = 3000, npcs = 30, resolution = 0.5)
)

all_markers<- FindAllMarkers(immune_human,
                   only.pos = TRUE,
                   min.pct = 0.25,
                   logfc.threshold = 0.25,
                   test.use = "wilcox")


write.csv(all_markers, "all_markershuman.csv")

feature_candidates <- list(
  "Naive/Central Memory T cells" = c("IL7R","LEF1","TCF7","CCR7","CD3D","BCL11B"),
  "γδ T cells"                   = c("IL23R","TRDC","ID2","PTPRCAP","LTB","DOCK3"),
  "NK/Cytotoxic lymphocytes"     = c("NKG7","TBX21","CTSW","GZMB","GNLY","CXCR6"),
  "B cells"                      = c("MS4A1","CD79B","PAX5","BANK1","CD79A","EBF1"),
  "Plasma cells"                 = c("JCHAIN","XBP1","DERL3","SSR4","TXNDC5","MZB1"),
  "Proliferating immune cells"   = c("MKI67","TOP2A","CENPF","RRM2","STMN1","HMGB2"),
  "Tissue-resident macrophages"  = c("MRC1","MSR1","GPNMB","ACP5","APOC1","CD163"),
  "Monocyte-derived macrophages" = c("CCR2","ADGRE1","PLD4","NRROS","SLAMF9","PID1"),
  "Dendritic cells (AXL+)"       = c("AXL","FCER1A","ZBTB46","CLEC5A","CD1C","IRF4"),
  "Neutrophils"                  = c("CSF3R","CXCR2","FPR2","MME","S100A8","S100A9"),
  "Activated neutrophils"        = c("NOS2","SLPI","CAMP","LRG1","IFIT3","PGLYRP1"),
  "Mast cells (progenitor-like)" = c("KIT","GATA2","HPGDS","MS4A2","SLC18A2","CLU"),
  "Mast cells (mature MCTC)"     = c("CPA3","CMA1","HDC","IL1RL1","TPSAB1","TPSB2")
)



pdf_markers <- c(
  "IL7R","LEF1","IL23R","NKG7","TBX21","CD79B","MS4A1","JCHAIN","XBP1",
  "MKI67","TOP2A","MRC1","MSR1","ADGRE1","CCR2","AXL","FCER1A",
  "CSF3R","CXCR2","NOS2","SLPI","GATA2","KIT","CMA1","CPA3"
)

feature_candidates[["PDF_markers"]] <- pdf_markers

immune_markers <- lapply(feature_candidates, unique)


annotate_clusters_from_markers <- function(markers_df,
                                           marker_signatures,
                                           pval_thresh = 0.05,
                                           logfc_thresh = 0.25) {
  
  sig <- markers_df %>%
    filter(p_val_adj < pval_thresh,
           avg_log2FC > logfc_thresh)
  
  cluster_ids <- sort(unique(sig$cluster))
  
  map_dfr(cluster_ids, function(cl) {
    genes_cl <- sig %>% filter(cluster == cl) %>% pull(gene)
    
    scores <- map_dbl(marker_signatures, ~ length(intersect(genes_cl, .x)))
    
    tibble(
      cluster = cl,
      best_label = names(scores)[which.max(scores)],
      best_score = max(scores)
    )
  })
}


cluster_annotations_list <- map(markers_list, ~ annotate_clusters_from_markers(immune_human, immune_markers))

imm_proc_list <- map2(
  immune_human,
  cluster_annotations_list,
  function(immune_human, ann) {
    
    ann$cluster <- as.character(ann$cluster)
    immune_human$seurat_clusters <- as.character(immune_human$seurat_clusters)
    
    meta <- immune_human@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      left_join(ann, by = c("seurat_clusters" = "cluster")) %>%
      tibble::column_to_rownames("cell")
    
    immune_human@meta.data <- meta
    immune_human$immune_subtype_annot <- immune_human$best_label
    immune_human
  }
)


DefaultAssay(imm_proc_list[["Rat"]]) <- "RNA"

all_genes <- unique(unlist(immune_markers))

p <- DotPlot(
  imm_proc_list[["Human"]],
  features = all_genes,
  group.by = "immune_subtype_annot"
) +
  scale_color_gradient(
    limits = c(0, 2),
    breaks = c(0, 1, 2),
    low = "lightgrey",
    high = "royalblue"
  ) +
  RotatedAxis() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

p


##Human dotplot
DefaultAssay(imm_proc_list[["Human"]]) <- "RNA"

p <- DotPlot(
  imm_proc_list[["Human"]],
  features = gene_df$gene,
  group.by = "seurat_clusters"     # ← HERE is the change
) +
  scale_color_gradient(
    limits = c(0, 2),
    breaks = c(0, 1, 2),
    low = "lightgrey",
    high = "royalblue"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # ← horizontal labels
    panel.grid = element_blank()
  )

group_boundaries <- cumsum(lengths(pdf_groups)) + 0.5

p <- p +
  geom_vline(
    xintercept = group_boundaries,
    color = "grey60",
    linewidth = 0.4
  )

p <- p +
  annotate(
    "text",
    x = tapply(seq_along(gene_order), gene_group, mean),
    y = length(unique(imm_proc_list[["Human"]]$seurat_clusters)) + 0.7,
    label = names(pdf_groups),
    angle = 0,      # ← no rotation
    hjust = 0.5,
    size = 3.2
  )

p


cluster_to_subtype <- c(
  "0"  = "Naive/CM T cells",
  "1"  = "B cells",
  "2"  = "Masts Cells - mature",
  "3"  = "Tissue-resident macrophages",
  "4"  = "AXL+ dendritic cells",
  "5"  = "γδ T cells",
  "6"  = "Activated neutrophils",
  "7"  = "Mast cell progenitors",
  "8"  = "Proliferating immune cells",
  "9"  = "NK / Cytotoxic lymphocytes",
  "10" = "Monocyte‑derived macrophages",
  "11" = "Naive B cells",
  "12" = "Plasma cells",
  "13" = "B cells (mature)",
  "14" = "Neutrophils"
)

immune_human$immune_subtype <- cluster_to_subtype[immune_human$seurat_clusters]

# Colors
immune_colors <- setNames(scales::hue_pal()(length(levels(immune$immune_subtype))), 
                          levels(immune$immune_subtype))

# ============================================================================
# CELL PROPORTION BARPLOT
# ============================================================================

prop_df <- immune@meta.data %>%
  group_by(species, immune_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

p_prop <- ggplot(prop_df, aes(x = species, y = proportion, fill = immune_subtype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = ifelse(n > 5, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = immune_colors) +
  labs(x = NULL, y = "Proportion (%)", fill = "Cell type",
       title = "Immune cell type proportions by species") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        legend.text = element_text(size = 8))

print(p_prop)

# ============================================================================
# DOTPLOT
# ============================================================================

DefaultAssay(immune) <- "RNA"

p_dotplot <- DotPlot(immune, features = unique(gene_groups$features.plot),
                     group.by = "immune_subtype",
                     dot.scale = 6) +
  scale_color_gradientn(colors = c("lightgrey", "#4393C3", "#2166AC", "#053061"),
                        limits = c(0, 2), oob = scales::squish,
                        name = "Avg\nExpression")

p_dotplot$data <- left_join(p_dotplot$data, gene_groups, by = "features.plot")

p_dotplot <- p_dotplot +
  facet_grid(~ gene_group, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size = 7, face = "bold"),
        strip.background = element_rect(fill = "white", color = "grey70"),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "right") +
  ggtitle("Canonical markers per immune subtype") +
  ylab("Immune Subtype") + xlab(NULL)

print(p_dotplot)