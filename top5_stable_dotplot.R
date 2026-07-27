################################################################################
# FIGURE 3 (new panel) — Top 5 Highly-Stable Genes per Major Cell Type
# Same visual style as reference Panel F: dot size = % expressing,
# dot color = avg expression, dashed blocks per cell type, blank outside
# each gene's own cell type
################################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cairoDevice)   # for cairo_pdf, per your convention

# ---- USER SETTINGS ----
CELLTYPE_COL <- "major_celltype"
SPECIES_COL  <- "species"

CELLTYPE_ORDER <- c("Epithelial Cells", "Immune Cells", "Endothelial Cells",
                     "Fibroblasts", "Muscle Cells", "Neural/Schwann Cells")
CT_SHORT <- c("Epithelial Cells" = "Epithelial", "Immune Cells" = "Immune",
              "Endothelial Cells" = "Endothelial", "Fibroblasts" = "Fibroblasts",
              "Muscle Cells" = "Muscle", "Neural/Schwann Cells" = "Schwann")

# Block order top-to-bottom in the plot (reverse of CELLTYPE_ORDER, matching
# the reference image where Schwann is the top block and Epithelial the bottom)
BLOCK_ORDER_TOP_TO_BOTTOM <- rev(CELLTYPE_ORDER)
SPECIES_ROW_ORDER <- c("Rat", "Marmoset", "Mouse", "Human")  # top-to-bottom within each block

N_TOP_GENES <- 5

# ---- STEP 1: pull top 5 highly-stable genes per cell type from Step4 output ----
# rank_df comes from run_step4() in your conservation pipeline (step4_rank_cv.csv),
# or pass the in-memory rank_df object directly if running right after the pipeline.
rank_df <- read.csv("step4_rank_cv.csv")

top5_stable <- rank_df %>%
  filter(sd_stability == "Highly stable") %>%
  group_by(celltype) %>%
  arrange(norm_sd) %>%
  slice_head(n = N_TOP_GENES) %>%
  ungroup()

for (ct in CELLTYPE_ORDER) {
  n <- sum(top5_stable$celltype == ct)
  if (n < N_TOP_GENES) {
    message("WARNING: ", ct, " has only ", n, " 'Highly stable' genes (<", N_TOP_GENES,
            "). Consider allowing 'Moderately stable' genes to fill the remaining slots.")
  }
}

gene_celltype_map <- setNames(top5_stable$celltype, top5_stable$gene)

# Gene order along x-axis: grouped by cell type (matching BLOCK order), ranked by norm_sd within block
gene_order <- top5_stable %>%
  mutate(celltype = factor(celltype, levels = CELLTYPE_ORDER)) %>%
  arrange(celltype, norm_sd) %>%
  pull(gene)

# ---- STEP 2: compute per-gene x per-celltype x per-species pct.exp / avg.exp ----
DefaultAssay(atlas) <- "RNA"
atlas <- JoinLayers(atlas)

genes_use <- unique(top5_stable$gene)
expr_mat  <- FetchData(atlas, vars = genes_use, layer = "data")
expr_mat$celltype <- atlas@meta.data[[CELLTYPE_COL]]
expr_mat$species  <- atlas@meta.data[[SPECIES_COL]]

dot_df <- expr_mat %>%
  pivot_longer(cols = all_of(genes_use), names_to = "gene", values_to = "expr") %>%
  # keep only rows where this gene belongs to THIS cell type (blank elsewhere, like the reference)
  filter(celltype == gene_celltype_map[gene]) %>%
  group_by(celltype, species, gene) %>%
  summarise(
    pct.exp = mean(expr > 0) * 100,
    avg.exp = mean(expm1(expr)),   # linear-scale mean, matches Seurat DotPlot convention
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(species %in% SPECIES_ROW_ORDER)

# ---- STEP 3: build combined row factor (celltype x species) in display order ----
row_levels <- c()
for (ct in BLOCK_ORDER_TOP_TO_BOTTOM) {
  row_levels <- c(row_levels, paste(ct, SPECIES_ROW_ORDER))
}
# reverse for ggplot (first factor level = bottom of plot)
row_levels <- rev(row_levels)

dot_df <- dot_df %>%
  mutate(
    row_id = paste(celltype, species),
    row_id = factor(row_id, levels = row_levels),
    gene   = factor(gene, levels = gene_order)
  )

# ---- STEP 4: block label positions + dashed separator positions ----
gene_positions <- tibble(gene = gene_order, x = seq_along(gene_order)) %>%
  left_join(top5_stable %>% select(gene, celltype), by = "gene") %>%
  mutate(celltype = factor(celltype, levels = CELLTYPE_ORDER))

block_labels <- gene_positions %>%
  group_by(celltype) %>%
  summarise(x_mid = mean(x), .groups = "drop")

# y position for block label = middle row of that celltype's 4-species block
row_df <- tibble(row_id = row_levels, y = seq_along(row_levels))
block_y <- tibble(celltype = BLOCK_ORDER_TOP_TO_BOTTOM) %>%
  mutate(rows = purrr::map(celltype, ~ paste(.x, SPECIES_ROW_ORDER))) %>%
  tidyr::unnest(rows) %>%
  left_join(row_df, by = c("rows" = "row_id")) %>%
  group_by(celltype) %>%
  summarise(y_mid = max(y) + 1.5, .groups = "drop")  # place label just above its block

# dashed vertical lines between cell type gene-blocks
vsep <- gene_positions %>%
  group_by(celltype) %>%
  summarise(x_max = max(x), .groups = "drop") %>%
  filter(celltype != tail(CELLTYPE_ORDER, 1)) %>%
  pull(x_max)
vsep <- vsep + 0.5

# ---- STEP 5: plot ----
p <- ggplot(dot_df, aes(x = gene, y = row_id)) +
  geom_vline(xintercept = vsep, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, color = avg.exp)) +
  scale_size(range = c(0.5, 6), breaks = c(25, 50, 75), name = "% Expressing",
             limits = c(0, 100)) +
  scale_color_gradient(low = "white", high = "#B2182B", name = "Avg. Expression") +
  geom_text(data = block_labels,
            aes(x = x_mid, y = Inf, label = CT_SHORT[as.character(celltype)]),
            vjust = -0.5, inherit.aes = FALSE, fontface = "plain", size = 4) +
  labs(title = "Suggested Conserved Markers per Cell Type Across Species",
       subtitle = paste0("Top ", N_TOP_GENES, " highly-stable genes per major cell type"),
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
  ) +
  coord_cartesian(clip = "off")

ggsave("panelF_top5_stable_dotplot.pdf", p, width = 14, height = 8, dpi = 300, device = cairo_pdf)
message("Saved: panelF_top5_stable_dotplot.pdf")
p
