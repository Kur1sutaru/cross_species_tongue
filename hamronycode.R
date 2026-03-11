
# ══════════════════════════════════════════════════════════════════
# PREPARE — convert mouse genes to uppercase
# ══════════════════════════════════════════════════════════════════

mouse_upper <- mouse
mouse_upper@assays$RNA <- RenameFeatures(
  mouse_upper@assays$RNA,
  newnames = toupper(rownames(mouse_upper@assays$RNA))
)

# If RenameFeatures doesn't work, use this instead:
counts_mouse <- GetAssayData(mouse, assay = "RNA", layer = "counts")
rownames(counts_mouse) <- toupper(rownames(counts_mouse))
mouse_upper <- CreateSeuratObject(counts  = counts_mouse,
                                  meta.data = mouse@meta.data)

# ══════════════════════════════════════════════════════════════════
# ADD SPECIES LABEL + SUBSET TO COMMON GENES
# ══════════════════════════════════════════════════════════════════

human$species  <- "Human"
marmoset$species    <- "Marmoset"
mouse_upper$species <- "Mouse"
rat$species         <- "Rat"

# Common genes (already computed = 12,761)
common_genes <- Reduce(intersect, list(
  rownames(human),
  rownames(marmoset),
  rownames(mouse_upper),
  rownames(rat)
))
cat("Common genes:", length(common_genes), "\n")

obj_human    <- human[common_genes, ]
obj_marmoset <- marmoset[common_genes, ]
obj_mouse    <- mouse_upper[common_genes, ]
obj_rat      <- rat[common_genes, ]

# ══════════════════════════════════════════════════════════════════
# MERGE INTO ONE SEURAT OBJECT
# SeuratWrappers LIGER requires a single merged object
# with dataset split by a metadata variable
# ══════════════════════════════════════════════════════════════════

tongue_merged <- merge(
  obj_human,
  y     = list(obj_marmoset, obj_mouse, obj_rat),
  add.cell.ids = c("Human", "Marmoset", "Mouse", "Rat"),
  merge.data   = TRUE
)

# Check
dim(tongue_merged)
table(tongue_merged$species)

# ══════════════════════════════════════════════════════════════════
# NORMALIZE + VARIABLE FEATURES per dataset
# LIGER requires this done per dataset, not globally
# ══════════════════════════════════════════════════════════════════

# Split by species for per-dataset processing
tongue_list <- SplitObject(tongue_merged, split.by = "species")

tongue_list <- lapply(tongue_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  x
})

# Merge back
tongue_merged <- merge(
  tongue_list[["Human"]],
  y            = list(tongue_list[["Marmoset"]],
                      tongue_list[["Mouse"]],
                      tongue_list[["Rat"]]),
  add.cell.ids = NULL,
  merge.data   = TRUE
)

# ══════════════════════════════════════════════════════════════════
# HARMONY INTEGRATION — works natively with Seurat v5
# ══════════════════════════════════════════════════════════════════

install.packages("harmony")
library(harmony)

# Join layers first for PCA
tongue_merged <- JoinLayers(tongue_merged)

# Standard preprocessing
tongue_merged <- NormalizeData(tongue_merged)
tongue_merged <- FindVariableFeatures(tongue_merged, nfeatures = 3000)
tongue_merged <- ScaleData(tongue_merged, do.center = TRUE)
tongue_merged <- RunPCA(tongue_merged, npcs = 50)

ElbowPlot(tongue_merged, ndims = 50)

# Run Harmony — integrates by species
tongue_merged <- RunHarmony(tongue_merged,
                            group.by.vars  = "species",
                            reduction      = "pca",
                            reduction.save = "harmony",
                            dims.use       = 1:30)

# Cluster + embeddings using harmony reduction
tongue_merged <- FindNeighbors(tongue_merged,
                               reduction = "harmony", dims = 1:30)
tongue_merged <- FindClusters(tongue_merged, resolution = 0.4)
tongue_merged <- RunUMAP(tongue_merged,
                         reduction = "harmony", dims = 1:30, seed.use = 42)
tongue_merged <- RunTSNE(tongue_merged,
                         reduction  = "harmony", dims = 1:30,
                         seed.use   = 42, perplexity = 50)

cat("Harmony integration done!\n")
names(tongue_merged@reductions)



# ══════════════════════════════════════════════════════════════════
# VISUALIZE
# ══════════════════════════════════════════════════════════════════

species_colors <- c(
  "Human"    = "#C0392B",
  "Marmoset" = "#2980B9",
  "Mouse"    = "#27AE60",
  "Rat"      = "#E67E22"
)

major_colors <- c(
  "Epithelial Cells"     = "#F08080",
  "Immune Cells"         = "#808000",
  "Endothelial Cells"    = "#228B22",
  "Fibroblasts"          = "#00CED1",
  "Muscle Cells"         = "#6495ED",
  "Neural/Schwann Cells" = "#FF69B4"
)

p1 <- DimPlot(tongue_merged,
              reduction = "umap",
              group.by  = "species",
              cols      = species_colors,
              pt.size   = 0.1,
              raster    = FALSE) +
  ggtitle("Cross-Species Atlas — by Species") +
  theme_classic() +
  theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank())

p2 <- DimPlot(tongue_merged,
              reduction = "umap",
              group.by  = "major_celltype",
              cols      = major_colors,
              pt.size   = 0.1,
              raster    = FALSE) +
  ggtitle("Cross-Species Atlas — by Cell Type") +
  theme_classic() +
  theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_blank())

p1 + p2
ggsave("tongue_atlas_harmony_umap.pdf",
       width = 20, height = 8, dpi = 300)

# Split by species
DimPlot(tongue_merged,
        reduction = "umap",
        group.by  = "major_celltype",
        cols      = major_colors,
        split.by  = "species",
        pt.size   = 0.1,
        raster    = FALSE,
        ncol      = 2) + theme_classic()
ggsave("tongue_atlas_harmony_umap_split.pdf",
       width = 18, height = 14, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# IF MOUSE GENES ARE UPPERCASE (after integration)
# use same gene list for all species
# ══════════════════════════════════════════════════════════════════

genes_all <- genes_hmr  # same for all 4 species

# ── Helper function
make_species_dotplot <- function(obj, genes, species_name,
                                 show_y = TRUE, show_legend = FALSE) {
  
  genes_avail <- genes[genes %in% rownames(obj)]
  cat(species_name, "- genes available:", length(genes_avail), "\n")
  
  Idents(obj) <- obj$major_celltype
  
  dot_data <- DotPlot(obj,
                      features = genes_avail,
                      group.by = "major_celltype")$data
  
  dot_data$id <- factor(dot_data$id, levels = rev(celltype_order))
  
  p <- ggplot(dot_data, aes(x = features.plot, y = id)) +
    geom_point(aes(size  = pct.exp,
                   fill  = id,
                   alpha = avg.exp),
               shape  = 21,
               color  = "grey30",
               stroke = 0.4) +
    scale_fill_manual(values = major_colors, name = "Cell Type") +
    scale_alpha_continuous(range  = c(0.1, 1),
                           limits = c(0.5, 2),
                           oob    = scales::squish,
                           name   = "Avg\nExpression") +
    scale_size(range = c(0.5, 7), name = "% Expressed") +
    scale_x_discrete(limits = genes_avail) +
    scale_y_discrete(limits = rev(celltype_order)) +
    labs(title = species_name, x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1,
                                  size  = 8,  face = "italic"),
      axis.text.y  = if (show_y) element_text(size = 9, face = "bold")
      else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = if (show_legend) "right" else "none",
      plot.margin  = margin(5, 5, 5, if (show_y) 5 else 2)
    )
  
  return(p)
}

# ══════════════════════════════════════════════════════════════════
# BUILD PLOTS
# ══════════════════════════════════════════════════════════════════

celltype_order <- c("Epithelial Cells", "Immune Cells", "Endothelial Cells",
                    "Fibroblasts", "Muscle Cells", "Neural/Schwann Cells")

major_colors <- c(
  "Epithelial Cells"     = "#F08080",
  "Immune Cells"         = "#808000",
  "Endothelial Cells"    = "#228B22",
  "Fibroblasts"          = "#00CED1",
  "Muscle Cells"         = "#6495ED",
  "Neural/Schwann Cells" = "#FF69B4"
)

p_human    <- make_species_dotplot(tongue_list[["Human"]],
                                   genes_all, "Human",
                                   show_y = TRUE,  show_legend = FALSE)

p_marmoset <- make_species_dotplot(tongue_list[["Marmoset"]],
                                   genes_all, "Marmoset",
                                   show_y = FALSE, show_legend = FALSE)

p_mouse    <- make_species_dotplot(tongue_list[["Mouse"]],
                                   genes_all, "Mouse",
                                   show_y = FALSE, show_legend = FALSE)

p_rat      <- make_species_dotplot(tongue_list[["Rat"]],
                                   genes_all, "Rat",
                                   show_y = FALSE, show_legend = TRUE)

# ══════════════════════════════════════════════════════════════════
# COMBINE
# ══════════════════════════════════════════════════════════════════

p_combined <- p_human + p_marmoset + p_mouse + p_rat +
  plot_layout(nrow   = 1,
              widths = c(1.4, 1, 1, 1.3),
              guides = "collect") &
  theme(legend.position = "right")

p_combined <- p_combined +
  plot_annotation(
    title = "Canonical Markers Across Species — Major Cell Types",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face  = "bold",
                                            size  = 13))
  )

print(p_combined)
ggsave("dotplot_all_species_major_celltypes.pdf",
       p_combined, width = 28, height = 7, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 1. CELL PROPORTIONS ACROSS SPECIES
# ══════════════════════════════════════════════════════════════════

df_counts <- as.data.frame(table(
  species        = tongue_merged$species,
  major_celltype = tongue_merged$major_celltype
))
colnames(df_counts)[3] <- "n"
df_counts <- df_counts[df_counts$n > 0, ]
df_counts$major_celltype <- factor(df_counts$major_celltype,
                                   levels = names(major_colors))
df_counts$species <- factor(df_counts$species,
                            levels = c("Human", "Marmoset", "Mouse", "Rat"))

# ── Proportional
p_prop <- ggplot(df_counts,
                 aes(x = species, y = n, fill = major_celltype)) +
  geom_bar(stat = "identity", width = 0.7, position = "fill") +
  scale_fill_manual(values = major_colors, name = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cell Type Proportions Across Species",
       x = NULL, y = "Proportion (%)") +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 11, face = "bold.italic"),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 11, face = "bold"),
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.text  = element_text(size  = 10),
        legend.title = element_text(size  = 11, face = "bold"))

# ── Absolute counts
p_abs <- ggplot(df_counts,
                aes(x = species, y = n, fill = major_celltype)) +
  geom_bar(stat = "identity", width = 0.7, position = "stack") +
  scale_fill_manual(values = major_colors, name = "Cell Type") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cell Type Counts Across Species",
       x = NULL, y = "Number of Cells") +
  theme_classic() +
  theme(axis.text.x  = element_text(size = 11, face = "bold.italic"),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_text(size = 11, face = "bold"),
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.text  = element_text(size  = 10),
        legend.title = element_text(size  = 11, face = "bold"))

p_prop + p_abs + plot_layout(guides = "collect")
ggsave("tongue_atlas_proportions_by_species.pdf",
       width = 14, height = 6, dpi = 300)

# ══════════════════════════════════════════════════════════════════
# 2. FEATURE PLOTS — one canonical marker per major cell type
# common to all 4 species (all uppercase in integrated obj)
# ══════════════════════════════════════════════════════════════════

# One representative marker per class — well conserved across species
feature_genes <- c(
  "KRT5",   # Epithelial
  "PTPRC",  # Immune
  "PECAM1", # Endothelial
  "COL1A1", # Fibroblasts
  "ACTA2",  # Muscle
  "S100B"   # Neural/Schwann
)

# Labels for titles
feature_labels <- c(
  "KRT5"   = "Epithelial Cells",
  "PTPRC"  = "Immune Cells",
  "PECAM1" = "Endothelial Cells",
  "COL1A1" = "Fibroblasts",
  "ACTA2"  = "Muscle Cells",
  "S100B"  = "Neural/Schwann Cells"
)

# Check all genes present
cat("Genes available in integrated object:\n")
print(feature_genes[feature_genes %in% rownames(tongue_merged)])
cat("Missing:\n")
print(feature_genes[!feature_genes %in% rownames(tongue_merged)])

# ── Feature plots — all species together
fp_list <- lapply(feature_genes, function(gene) {
  FeaturePlot(tongue_merged,
              features   = gene,
              reduction  = "umap",
              pt.size    = 0.05,
              raster     = FALSE,
              order      = TRUE) +
    scale_color_gradientn(
      colors = c("lightgrey", "#FEE08B", "#FC8D59", "#D53E4F"),
      name   = "Expression"
    ) +
    labs(title    = feature_labels[gene],
         subtitle = paste0("(", gene, ")")) +
    theme_classic() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold",   size = 10),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 9,
                                   color = "#555555"),
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 8),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank()
    )
})

# Combine 6 feature plots
p_features <- wrap_plots(fp_list, ncol = 6) +
  plot_annotation(
    title = "Canonical Marker Expression — Cross-Species Atlas",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face  = "bold",
                                            size  = 13))
  )

print(p_features)
ggsave("tongue_atlas_featureplots.pdf",
       p_features, width = 24, height = 5, dpi = 300)

# ── Feature plots split by species (4 rows x 6 genes = 24 panels)
fp_split_list <- lapply(feature_genes, function(gene) {
  FeaturePlot(tongue_merged,
              features   = gene,
              reduction  = "umap",
              split.by   = "species",
              pt.size    = 0.05,
              raster     = FALSE,
              order      = TRUE) &
    scale_color_gradientn(
      colors = c("lightgrey", "#FEE08B", "#FC8D59", "#D53E4F"),
      name   = "Expression"
    ) &
    theme_classic() &
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
})

# Save each gene split plot separately
for (i in seq_along(feature_genes)) {
  gene <- feature_genes[i]
  ggsave(paste0("tongue_atlas_featureplot_", gene, "_split.pdf"),
         fp_split_list[[i]],
         width = 20, height = 5, dpi = 300)
  
  
  
}

# ══════════════════════════════════════════════════════════════════
# FEATURE PLOTS — each gene split by species (4 panels per gene)
# ══════════════════════════════════════════════════════════════════

feature_genes <- c(
  "KRT5",   # Epithelial
  "PTPRC",  # Immune
  "PECAM1", # Endothelial
  "COL1A1", # Fibroblasts
  "ACTA2",  # Muscle
  "S100B"   # Neural/Schwann
)

feature_labels <- c(
  "KRT5"   = "Epithelial Cells",
  "PTPRC"  = "Immune Cells",
  "PECAM1" = "Endothelial Cells",
  "COL1A1" = "Fibroblasts",
  "ACTA2"  = "Muscle Cells",
  "S100B"  = "Neural/Schwann Cells"
)

# Make sure species order is consistent
tongue_merged$species <- factor(tongue_merged$species,
                                levels = c("Human", "Marmoset",
                                           "Mouse", "Rat"))

# ══════════════════════════════════════════════════════════════════
# LOOP — one PDF per gene, 4 panels (one per species)
# ══════════════════════════════════════════════════════════════════

for (gene in feature_genes) {
  
  cat("Plotting", gene, "...\n")
  
  # Split by species manually for full control
  p_list <- lapply(c("Human", "Marmoset", "Mouse", "Rat"), function(sp) {
    
    # Subset to species
    obj_sp <- subset(tongue_merged, subset = species == sp)
    
    FeaturePlot(obj_sp,
                features  = gene,
                reduction = "umap",
                pt.size   = 0.1,
                raster    = FALSE,
                order     = TRUE) +
      scale_color_gradientn(
        colors = c("lightgrey", "#FEE08B", "#FC8D59", "#D53E4F"),
        limits = c(0, NA),
        name   = "Expression"
      ) +
      labs(title = sp) +
      theme_classic() +
      theme(
        plot.title   = element_text(hjust = 0.5, face = "bold.italic",
                                    size  = 11),
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        axis.title   = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9)
      )
  })
  
  # Combine 4 panels
  p_combined <- wrap_plots(p_list, nrow = 1) +
    plot_annotation(
      title    = feature_labels[gene],
      subtitle = paste0("Gene: ", gene),
      theme    = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, face = "italic",
                                     size  = 10, color = "#555555")
      )
    )
  
  print(p_combined)
  
  ggsave(paste0("featureplot_", gene, "_by_species.pdf"),
         p_combined,
         width  = 20,
         height = 6,
         dpi    = 300)
}

cat("All feature plots saved!\n")

# ══════════════════════════════════════════════════════════════════
# OPTIONAL — all 6 genes in one figure (6 rows x 4 species = 24 panels)
# ══════════════════════════════════════════════════════════════════

all_panels <- lapply(feature_genes, function(gene) {
  
  lapply(c("Human", "Marmoset", "Mouse", "Rat"), function(sp) {
    
    obj_sp <- subset(tongue_merged, subset = species == sp)
    
    FeaturePlot(obj_sp,
                features  = gene,
                reduction = "umap",
                pt.size   = 0.05,
                raster    = FALSE,
                order     = TRUE) +
      scale_color_gradientn(
        colors = c("lightgrey", "#FEE08B", "#FC8D59", "#D53E4F"),
        limits = c(0, NA),
        name   = "Expr"
      ) +
      labs(title = if (gene == feature_genes[1]) sp else NULL,
           subtitle = if (sp == "Human") paste0(feature_labels[gene],
                                                "\n(", gene, ")")
           else NULL) +
      theme_classic() +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold.italic",
                                     size  = 10),
        plot.subtitle = element_text(hjust = 0,   face = "bold",
                                     size  = 8),
        axis.text     = element_blank(),
        axis.ticks    = element_blank(),
        axis.title    = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text   = element_text(size = 7),
        legend.title  = element_text(size = 7),
        plot.margin   = margin(2, 2, 2, 2)
      )
  })
}) %>% unlist(recursive = FALSE)

p_all <- wrap_plots(all_panels, nrow = 6, ncol = 4) +
  plot_annotation(
    title = "Canonical Markers — Cross-Species Expression",
    theme = theme(plot.title = element_text(hjust = 0.5,
                                            face  = "bold",
                                            size  = 14))
  )

ggsave("featureplots_all_genes_by_species.pdf",
       p_all,
       width  = 20,
       height = 26,
       dpi    = 300)

cat("Combined figure saved!\n")


#### edits on figures
# =============================================================================
# Cross-Species Tongue Atlas — Violin & Feature Plots
# Matching scale of canonical marker dotplot (scaled avg expression 0–2)
# Species: Human, Marmoset, Mouse, Rat
# Markers: 24 canonical markers grouped by cell type
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

# =============================================================================
# 0. USER SETTINGS — adjust these paths/column names to your object
# =============================================================================

RDS_PATH       <- "your_integrated_seurat.rds"   # path to your Seurat object
CELL_TYPE_COL  <- "cell_type"                     # metadata column: Epi, Immune, Endo, etc.
SPECIES_COL    <- "species"                       # metadata column: Human, Marmoset, Mouse, Rat
OUTPUT_DIR     <- "."                             # where to save PDFs

# Species levels and colors (match dotplot row label colors)
SPECIES_LEVELS  <- c("Human", "Marmoset", "Mouse", "Rat")
SPECIES_COLORS  <- c(Human   = "#E41A1C",   # red
                     Marmoset = "#FF7F00",  # orange
                     Mouse   = "#4DAF4A",   # green
                     Rat     = "#984EA3")   # purple

# Cell type levels and colors (match dotplot row label colors)
CELLTYPE_LEVELS <- c("Epi", "Immune", "Endo", "Fibro", "Muscle", "Schwann")
CELLTYPE_COLORS <- c(Epi     = "#E41A1C",
                     Immune  = "#FF7F00",
                     Endo    = "#4DAF4A",
                     Fibro   = "#00BFC4",
                     Muscle  = "#9966FF",
                     Schwann = "#FF69B4")

# Markers grouped by cell type (same order as dotplot columns)
MARKER_GROUPS <- list(
  Epi     = c("KRT5", "TP63", "KRT14", "EPCAM"),
  Immune  = c("PTPRC", "CD68", "AIF1", "CD3D"),
  Endo    = c("PECAM1", "CDH5", "VWF", "FLT1"),
  Fibro   = c("COL1A1", "DCN", "LUM", "PDGFRA"),
  Muscle  = c("ACTA2", "MYH11", "CNN1", "DES"),
  Schwann = c("S100B", "SOX10", "MPZ", "MBP")
)

ALL_MARKERS <- unlist(MARKER_GROUPS, use.names = FALSE)

# Color scale matching dotplot (white → dark red, 0–2 scaled expression)
EXPR_PALETTE <- c("grey90", "#FCBBA1", "#FC9272", "#FB6A4A",
                  "#EF3B2C", "#CB181D", "#99000D")
EXPR_LIMITS  <- c(0, 2)

# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================

cat("Loading Seurat object...\n")
seu <- readRDS(RDS_PATH)

# Ensure factors with correct levels
seu@meta.data[[SPECIES_COL]]   <- factor(seu@meta.data[[SPECIES_COL]],
                                         levels = SPECIES_LEVELS)
seu@meta.data[[CELL_TYPE_COL]] <- factor(seu@meta.data[[CELL_TYPE_COL]],
                                         levels = CELLTYPE_LEVELS)

# Filter to markers present in the object
available_markers <- ALL_MARKERS[ALL_MARKERS %in% rownames(seu)]
missing <- setdiff(ALL_MARKERS, available_markers)
if (length(missing) > 0) {
  cat("⚠ Markers not found in object (skipped):", paste(missing, collapse = ", "), "\n")
}

# Use scaled data if available, otherwise scale on the fly
if ("scale.data" %in% names(seu@assays$RNA) &&
    length(rownames(seu@assays$RNA@scale.data)) > 0) {
  cat("Using existing scale.data slot.\n")
  slot_use <- "scale.data"
} else {
  cat("Scaling data for marker genes...\n")
  seu <- ScaleData(seu, features = available_markers, verbose = FALSE)
  slot_use <- "scale.data"
}

# Helper: clamp scaled values to 0–2 (matching dotplot display range)
clamp <- function(x, lo = 0, hi = 2) pmax(lo, pmin(hi, x))

# =============================================================================
# 2. SHARED THEME
# =============================================================================

base_theme <- theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 7, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y      = element_text(size = 7),
    axis.title       = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40")
  )

# =============================================================================
# 3. VIOLIN PLOTS — per cell type, split by species, y-axis fixed 0–2
# =============================================================================

cat("Generating violin plots...\n")

make_violin_panel <- function(celltype, markers, seu, species_col, celltype_col) {
  
  markers <- markers[markers %in% rownames(seu)]
  if (length(markers) == 0) return(NULL)
  
  cells_use <- colnames(seu)[seu@meta.data[[celltype_col]] == celltype]
  if (length(cells_use) < 10) return(NULL)
  
  seu_sub <- subset(seu, cells = cells_use)
  
  # Extract scaled expression, clamp to 0–2
  expr_mat <- GetAssayData(seu_sub, slot = slot_use)[markers, , drop = FALSE]
  expr_mat  <- clamp(expr_mat)
  
  df <- as.data.frame(t(as.matrix(expr_mat)))
  df$species <- seu_sub@meta.data[[species_col]]
  df <- tidyr::pivot_longer(df, cols = all_of(markers),
                            names_to = "gene", values_to = "expression")
  df$gene    <- factor(df$gene, levels = markers)
  df$species <- factor(df$species, levels = SPECIES_LEVELS)
  
  ggplot(df, aes(x = species, y = expression, fill = species)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, linewidth = 0.3,
                 fill = "white", alpha = 0.7) +
    facet_wrap(~gene, nrow = 1, scales = "fixed") +
    scale_fill_manual(values = SPECIES_COLORS, name = "Species") +
    scale_y_continuous(limits = EXPR_LIMITS, breaks = c(0, 0.5, 1, 1.5, 2)) +
    labs(title = celltype,
         x = NULL,
         y = "Scaled Expression (0–2)") +
    base_theme +
    theme(legend.position = "none",
          plot.title = element_text(color = CELLTYPE_COLORS[[celltype]]))
}

# Build one page per cell type
violin_plots <- lapply(names(MARKER_GROUPS), function(ct) {
  make_violin_panel(ct, MARKER_GROUPS[[ct]], seu, SPECIES_COL, CELL_TYPE_COL)
})
names(violin_plots) <- names(MARKER_GROUPS)
violin_plots <- Filter(Negate(is.null), violin_plots)

# Legend panel
legend_df <- data.frame(species = factor(SPECIES_LEVELS, levels = SPECIES_LEVELS), x = 1, y = 1)
legend_plot <- ggplot(legend_df, aes(x, y, fill = species)) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = SPECIES_COLORS, name = "Species") +
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"))

# Combine all cell types into one figure
violin_combined <- wrap_plots(violin_plots, ncol = 1) +
  plot_annotation(
    title    = "Canonical Marker Expression — Cross-Species Tongue Atlas",
    subtitle = "Violin plots | Scaled expression (0–2) | Split by species",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"))
  )

# Save violin plots
violin_path <- file.path(OUTPUT_DIR, "violin_crossspecies_canonical_markers.pdf")
ggsave(violin_path,
       plot   = violin_combined,
       width  = 18,
       height = length(violin_plots) * 3.5,
       units  = "in",
       device = "pdf")
cat("✓ Violin plots saved:", violin_path, "\n")

# =============================================================================
# 4. FEATURE PLOTS — Layout A: 4 panels per gene (one per species)
# =============================================================================

cat("Generating feature plots (4-panel per gene)...\n")

# Check UMAP exists
if (!"umap" %in% names(seu@reductions) && !"UMAP" %in% names(seu@reductions)) {
  stop("No UMAP reduction found. Run RunUMAP() first or check reduction name.")
}
umap_key <- if ("umap" %in% names(seu@reductions)) "umap" else "UMAP"

make_feature_4panel <- function(gene, seu, species_col, umap_key) {
  
  if (!gene %in% rownames(seu)) return(NULL)
  
  umap_coords <- as.data.frame(seu@reductions[[umap_key]]@cell.embeddings[, 1:2])
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords$species <- seu@meta.data[[species_col]]
  
  expr_vals <- clamp(GetAssayData(seu, slot = slot_use)[gene, ])
  umap_coords$expression <- expr_vals
  
  panels <- lapply(SPECIES_LEVELS, function(sp) {
    df_sp <- umap_coords[umap_coords$species == sp, ]
    # plot non-expressing cells first (grey), then expressing on top
    df_sp <- df_sp[order(df_sp$expression), ]
    
    ggplot(df_sp, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
      geom_point(size = 0.15, alpha = 0.7, stroke = 0) +
      scale_color_gradientn(
        colors  = EXPR_PALETTE,
        limits  = EXPR_LIMITS,
        oob     = squish,
        name    = "Scaled\nExpr (0–2)"
      ) +
      labs(title = sp, x = NULL, y = NULL) +
      base_theme +
      theme(
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = SPECIES_COLORS[[sp]], linewidth = 1),
        plot.title = element_text(color = SPECIES_COLORS[[sp]], size = 9),
        legend.position = "none"
      )
  })
  
  # Shared color legend
  legend_g <- ggplot(data.frame(x = 1, y = seq(0, 2, 0.1)),
                     aes(x, y, color = y)) +
    geom_point() +
    scale_color_gradientn(colors = EXPR_PALETTE, limits = EXPR_LIMITS,
                          name = "Scaled\nExpr") +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text  = element_text(size = 6),
          legend.key.height = unit(0.8, "cm"))
  legend_grob <- cowplot::get_legend(legend_g)
  
  row_plot <- wrap_plots(panels, nrow = 1) +
    plot_annotation(title = gene,
                    theme = theme(plot.title = element_text(
                      size = 11, face = "bold.italic", hjust = 0.5))) |
    legend_grob
  
  row_plot
}

# Save one PDF with all genes, 4-panel layout
# Group by cell type, one page per cell type
feat4_path <- file.path(OUTPUT_DIR, "featureplot_4panel_crossspecies_canonical_markers.pdf")
pdf(feat4_path, width = 16, height = 3.5)

for (ct in names(MARKER_GROUPS)) {
  markers <- MARKER_GROUPS[[ct]][MARKER_GROUPS[[ct]] %in% rownames(seu)]
  for (gene in markers) {
    p <- make_feature_4panel(gene, seu, SPECIES_COL, umap_key)
    if (!is.null(p)) print(p)
  }
}
dev.off()
cat("✓ Feature plots (4-panel) saved:", feat4_path, "\n")

# =============================================================================
# 5. FEATURE PLOTS — Layout B: merged UMAP (all species), gene-faceted grid
# =============================================================================

cat("Generating feature plots (merged UMAP, faceted by cell type group)...\n")

umap_coords_all <- as.data.frame(seu@reductions[[umap_key]]@cell.embeddings[, 1:2])
colnames(umap_coords_all) <- c("UMAP_1", "UMAP_2")
umap_coords_all$species   <- seu@meta.data[[SPECIES_COL]]
umap_coords_all$cell_type <- seu@meta.data[[CELL_TYPE_COL]]
umap_coords_all$cell       <- colnames(seu)

feat_merged_path <- file.path(OUTPUT_DIR, "featureplot_merged_crossspecies_canonical_markers.pdf")
pdf(feat_merged_path, width = 18, height = ceiling(length(available_markers) / 6) * 3 + 1)

for (ct in names(MARKER_GROUPS)) {
  markers <- MARKER_GROUPS[[ct]][MARKER_GROUPS[[ct]] %in% rownames(seu)]
  if (length(markers) == 0) next
  
  expr_mat <- clamp(
    as.matrix(GetAssayData(seu, slot = slot_use)[markers, , drop = FALSE])
  )
  
  df_long <- as.data.frame(t(expr_mat))
  df_long$cell <- rownames(df_long)
  df_long <- merge(df_long, umap_coords_all, by = "cell")
  df_long <- tidyr::pivot_longer(df_long, cols = all_of(markers),
                                 names_to = "gene", values_to = "expression")
  df_long$gene <- factor(df_long$gene, levels = markers)
  df_long <- df_long[order(df_long$expression), ]  # plot high expr on top
  
  p <- ggplot(df_long, aes(x = UMAP_1, y = UMAP_2, color = expression)) +
    geom_point(size = 0.1, alpha = 0.6, stroke = 0) +
    facet_wrap(~gene, nrow = 1) +
    scale_color_gradientn(
      colors  = EXPR_PALETTE,
      limits  = EXPR_LIMITS,
      oob     = squish,
      name    = "Scaled\nExpr (0–2)"
    ) +
    labs(title    = paste0(ct, " Markers — Merged Cross-Species UMAP"),
         subtitle = "All 4 species combined | Scaled expression (0–2)",
         x = "UMAP 1", y = "UMAP 2") +
    base_theme +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
  
  print(p)
}
dev.off()
cat("✓ Feature plots (merged) saved:", feat_merged_path, "\n")

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n=== All outputs saved ===\n")
cat("1.", violin_path, "\n")
cat("2.", feat4_path, "\n")
cat("3.", feat_merged_path, "\n")
cat("\nRemember to set RDS_PATH, CELL_TYPE_COL, and SPECIES_COL at the top of the script.\n")





cat("All plots saved!\n")

saveRDS(tongue_merged, "tongue_atlas_harmony_seurat.rds")
cat("Saved!\n")
