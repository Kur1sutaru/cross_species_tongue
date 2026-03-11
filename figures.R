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

RDS_PATH       <- "tongue_merged.rds"             # path to your Seurat object
CELL_TYPE_COL  <- "major_celltype"                # metadata column: Epi, Immune, Endo, etc.
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