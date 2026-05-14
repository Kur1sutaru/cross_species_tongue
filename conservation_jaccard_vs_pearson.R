# ============================================================
# Conservation Analysis: Jaccard Index vs Pearson Correlation
# Step 3 replacement/comparison
# Input: conserved_markers output from FindConservedMarkers
#        already filtered (1,342 genes, all 4 species)
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(corrplot)

# ============================================================
# SECTION 1 — Compute mean expression per gene per species
#             (one value per cell type x gene x species)
# ============================================================

species_list  <- c("Human", "Mouse", "Marmoset", "Rat")
cell_types    <- c("Epithelial Cells", "Immune Cells", "Endothelial Cells",
                   "Fibroblasts", "Muscle Cells", "Neural/Schwann Cells")

CT_COL <- "major_celltype"   # adjust if needed
SP_COL <- "species"

# Build a long data frame: gene x species x cell_type -> mean log-norm expression
mean_expr_list <- lapply(species_list, function(sp) {
  
  cells <- colnames(atlas)[atlas[[SP_COL, drop = TRUE]] == sp]
  obj   <- subset(atlas, cells = cells)
  DefaultAssay(obj) <- "RNA"
  obj   <- NormalizeData(obj, verbose = FALSE)
  
  # Average expression per cell type
  Idents(obj) <- CT_COL
  avg <- AverageExpression(obj, assays = "RNA", return.seurat = FALSE,
                           verbose = FALSE)$RNA
  
  # Keep only conserved genes present in this species
  genes_keep <- rownames(avg)[rownames(avg) %in% conserved_genes]
  avg_df <- as.data.frame(avg[genes_keep, , drop = FALSE])
  avg_df$gene    <- rownames(avg_df)
  avg_df$species <- sp
  
  pivot_longer(avg_df, -c(gene, species),
               names_to = "cell_type", values_to = "mean_expr")
})

mean_expr <- bind_rows(mean_expr_list)

# ============================================================
# SECTION 2 — Pearson Correlation (Human as reference)
# ============================================================

pearson_results <- lapply(cell_types, function(ct) {
  
  df_ct <- mean_expr %>%
    filter(cell_type == ct) %>%
    pivot_wider(names_from = species, values_from = mean_expr) %>%
    drop_na()
  
  if (nrow(df_ct) < 10) return(NULL)
  
  comparisons <- list(
    c("Human", "Mouse"),
    c("Human", "Marmoset"),
    c("Human", "Rat")
  )
  
  lapply(comparisons, function(pair) {
    x <- df_ct[[pair[1]]]
    y <- df_ct[[pair[2]]]
    r <- cor(x, y, method = "pearson", use = "complete.obs")
    data.frame(
      cell_type  = ct,
      comparison = paste(pair[1], "vs", pair[2]),
      metric     = "Pearson r",
      value      = r,
      n_genes    = nrow(df_ct)
    )
  }) %>% bind_rows()
  
}) %>% bind_rows()

# ============================================================
# SECTION 3 — Jaccard Index
# Uses binary presence/absence: gene expressed if mean > threshold
# ============================================================

EXPR_THRESHOLD <- 0.1   # log-norm mean expression cutoff for "expressed"

jaccard_index <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

jaccard_results <- lapply(cell_types, function(ct) {
  
  df_ct <- mean_expr %>%
    filter(cell_type == ct) %>%
    pivot_wider(names_from = species, values_from = mean_expr) %>%
    drop_na()
  
  if (nrow(df_ct) < 10) return(NULL)
  
  # Binary gene sets per species
  gene_sets <- lapply(species_list, function(sp) {
    df_ct$gene[df_ct[[sp]] > EXPR_THRESHOLD]
  })
  names(gene_sets) <- species_list
  
  comparisons <- list(
    c("Human", "Mouse"),
    c("Human", "Marmoset"),
    c("Human", "Rat")
  )
  
  lapply(comparisons, function(pair) {
    j <- jaccard_index(gene_sets[[pair[1]]], gene_sets[[pair[2]]])
    data.frame(
      cell_type  = ct,
      comparison = paste(pair[1], "vs", pair[2]),
      metric     = "Jaccard Index",
      value      = j,
      n_genes    = nrow(df_ct)
    )
  }) %>% bind_rows()
  
}) %>% bind_rows()

# ============================================================
# SECTION 4 — Combine and summarize
# ============================================================

combined <- bind_rows(pearson_results, jaccard_results)

# Print summary table
summary_table <- combined %>%
  group_by(metric, comparison) %>%
  summarise(
    mean_value = round(mean(value, na.rm = TRUE), 3),
    sd_value   = round(sd(value,   na.rm = TRUE), 3),
    .groups = "drop"
  )

print(summary_table)

# Per-cell-type table
per_ct_table <- combined %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(cell_type, comparison)

print(per_ct_table)

# ============================================================
# SECTION 5 — Plots
# ============================================================

cell_type_colors <- c(
  "Epithelial Cells"    = "#B5A44E",
  "Immune Cells"        = "#00BFC4",
  "Endothelial Cells"   = "#FA7F6F",
  "Fibroblasts"         = "#44B54A",
  "Muscle Cells"        = "#619CFF",
  "Neural/Schwann Cells"= "#F564E3"
)

# ---- Plot 1: Grouped bar — Pearson vs Jaccard per cell type ----
p1 <- combined %>%
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors))) %>%
  ggplot(aes(x = comparison, y = value, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = cell_type_colors, name = "Cell Type") +
  labs(
    title = "Conservation across species: Pearson r vs Jaccard Index",
    x     = NULL,
    y     = "Value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.text      = element_text(face = "bold")
  )

# ---- Plot 2: Scatter — Pearson vs Jaccard per observation ----
scatter_df <- combined %>%
  pivot_wider(names_from = metric, values_from = value)

p2 <- scatter_df %>%
  mutate(cell_type = factor(cell_type, levels = names(cell_type_colors))) %>%
  ggplot(aes(x = `Jaccard Index`, y = `Pearson r`,
             color = cell_type, shape = comparison)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = cell_type_colors, name = "Cell Type") +
  scale_shape_manual(values = c(16, 17, 15), name = "Comparison") +
  labs(
    title = "Pearson r vs Jaccard Index per cell type and comparison",
    x     = "Jaccard Index",
    y     = "Pearson r"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# Pearson between the two metrics across all observations
r_metrics <- cor(scatter_df$`Pearson r`, scatter_df$`Jaccard Index`,
                 use = "complete.obs", method = "pearson")
message(sprintf("Correlation between Pearson r and Jaccard Index: r = %.3f", r_metrics))

# ---- Plot 3: Heatmap-style tile — value per cell type x comparison ----
make_tile <- function(df, metric_name, fill_low, fill_high) {
  df %>%
    filter(metric == metric_name) %>%
    mutate(cell_type = factor(cell_type, levels = rev(names(cell_type_colors)))) %>%
    ggplot(aes(x = comparison, y = cell_type, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(value, 2)), size = 3.2) +
    scale_fill_gradient(low = fill_low, high = fill_high,
                        name = metric_name, limits = c(0, 1)) +
    labs(title = metric_name, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

p3a <- make_tile(combined, "Pearson r",    "#DEEBF7", "#08519C")
p3b <- make_tile(combined, "Jaccard Index","#E5F5E0", "#238B45")

p3 <- p3a | p3b

# ============================================================
# SECTION 6 — Save
# ============================================================

ggsave("conservation_barplot_pearson_vs_jaccard.pdf",
       p1, width = 10, height = 8, device = cairo_pdf)
ggsave("conservation_barplot_pearson_vs_jaccard.png",
       p1, width = 10, height = 8, dpi = 300)

ggsave("conservation_scatter_pearson_vs_jaccard.pdf",
       p2, width = 8, height = 6, device = cairo_pdf)
ggsave("conservation_scatter_pearson_vs_jaccard.png",
       p2, width = 8, height = 6, dpi = 300)

ggsave("conservation_heatmap_pearson_vs_jaccard.pdf",
       p3, width = 12, height = 5, device = cairo_pdf)
ggsave("conservation_heatmap_pearson_vs_jaccard.png",
       p3, width = 12, height = 5, dpi = 300)

message("All done! 3 figure sets saved.")

# ============================================================
# SECTION 7 — Which metric is better? Decision helper
# ============================================================
# 
# Pearson r   → measures whether expression LEVELS are correlated
#               across species (quantitative). Best for asking:
#               "Do the same genes show high expression in both species?"
#
# Jaccard     → measures binary OVERLAP of expressed gene sets.
#               Best for asking: "Do the same genes turn on/off?"
#               Sensitive to the threshold you choose (EXPR_THRESHOLD).
#
# Recommendation:
#   - Use PEARSON as the primary metric (as in your flowchart step 3),
#     because it captures graded expression differences and is threshold-free.
#   - Use JACCARD as a complementary validation: if Pearson is high but
#     Jaccard is low, it means a few highly-expressed genes drive the
#     correlation while many others differ in presence/absence.
#   - If Pearson ~ Jaccard correlation (r_metrics above) is > 0.8,
#     both metrics agree and Pearson is sufficient.
# ============================================================
