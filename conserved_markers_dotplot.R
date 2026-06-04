# ============================================================
# Conserved markers dotplot across species
# Reproduces figure G style from the tongue atlas manuscript
# Input: step4_rank_cv.csv
# ============================================================

library(tidyverse)
library(ggplot2)

# ── 1. Load data ─────────────────────────────────────────────
df <- read.csv("step4_rank_cv.csv", stringsAsFactors = FALSE)

# ── 2. Select top 5 genes per cell type (lowest sd_rank = most stable) ──
top5 <- df %>%
  group_by(celltype) %>%
  slice_min(order_by = sd_rank, n = 5) %>%
  ungroup()

# ── 3. Pivot to long format (one row per gene × species) ──────
expr_long <- top5 %>%
  select(gene, celltype, Human_expr, Mouse_expr, Marmoset_expr, Rat_expr) %>%
  pivot_longer(
    cols = c(Human_expr, Mouse_expr, Marmoset_expr, Rat_expr),
    names_to  = "species",
    values_to = "avg_expr"
  ) %>%
  mutate(species = str_remove(species, "_expr"))

pct_long <- top5 %>%
  select(gene, celltype, Human_pct, Mouse_pct, Marmoset_pct, Rat_pct) %>%
  pivot_longer(
    cols = c(Human_pct, Mouse_pct, Marmoset_pct, Rat_pct),
    names_to  = "species",
    values_to = "pct_expr"
  ) %>%
  mutate(species = str_remove(species, "_pct"))

plot_df <- left_join(expr_long, pct_long, by = c("gene", "celltype", "species"))

# ── 4. Factor ordering ────────────────────────────────────────
# Cell types: match original figure order (bottom to top in ggplot y-axis)
celltype_order <- c(
  "Epithelial Cells",
  "Immune Cells",
  "Endothelial Cells",
  "Fibroblasts",
  "Muscle Cells",
  "Neural/Schwann Cells"
)

# Species: Rat top → Human bottom (as in original)
species_order <- c("Rat", "Marmoset", "Mouse", "Human")

# Gene order: per cell type, sorted by sd_rank ascending; then grouped by cell type
gene_order_df <- top5 %>%
  arrange(
    factor(celltype, levels = celltype_order),
    sd_rank
  )
gene_order <- unique(gene_order_df$gene)

plot_df <- plot_df %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    species  = factor(species,  levels = species_order),
    gene     = factor(gene,     levels = gene_order)
  )

# ── 5. Build the plot ─────────────────────────────────────────
p <- ggplot(plot_df, aes(x = gene, y = species)) +

  # Dots: size = % expressing, fill = avg expression (log2FC)
  geom_point(aes(size = pct_expr, fill = avg_expr),
             shape = 21, color = "white", stroke = 0.2) +

  # Color scale: white → deep red (matching original)
  # Use oob = scales::squish so values above max clamp to darkest red, not grey
  scale_fill_gradient(
    low      = "white",
    high     = "#8B1A1A",
    name     = "Avg. Expression\n(log\u2082FC)",
    limits   = c(0, 5),
    breaks   = c(0, 1, 2, 3, 4, 5),
    oob      = scales::squish,
    guide    = guide_colorbar(
      barwidth  = 0.7,
      barheight = 4,
      ticks     = FALSE,
      title.position = "top"
    )
  ) +

  # Size scale: % expressing
  scale_size_continuous(
    name   = "% Expressing",
    range  = c(0.5, 6),
    breaks = c(25, 50, 75, 100),
    labels = c("25%", "50%", "75%", "100%"),
    guide  = guide_legend(
      override.aes = list(fill = "grey50", color = "white"),
      title.position = "top"
    )
  ) +

  # Facet by cell type (columns), free x so each panel shows only its genes
  facet_grid(
    . ~ celltype,
    scales = "free_x",
    space  = "free_x"
  ) +

  # Labels
  labs(
    title = "Suggested conserved markers per cell type across species",
    x     = NULL,
    y     = NULL
  ) +

  # Theme
  theme_bw(base_size = 10) +
  theme(
    # Title
    plot.title         = element_text(hjust = 0.5, size = 11, face = "plain",
                                      margin = margin(b = 8)),

    # Facet strips
    strip.text         = element_text(size = 9, face = "plain"),
    strip.background   = element_rect(fill = "white", color = NA),

    # Axes
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y        = element_text(size = 8),
    axis.ticks         = element_line(linewidth = 0.3),

    # Panel
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.15, "cm"),

    # Legend
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.4, "cm"),
    legend.box         = "vertical",
    legend.position    = "right",
    legend.box.spacing = unit(0.1, "cm"),

    # Margins
    plot.margin = margin(6, 6, 6, 6)
  )

# ── 6. Save ───────────────────────────────────────────────────
ggsave(
  "conserved_markers_dotplot.pdf",
  plot   = p,
  width  = 12,
  height = 4,
  device = cairo_pdf
)

ggsave(
  "conserved_markers_dotplot.png",
  plot   = p,
  width  = 12,
  height = 4,
  dpi    = 300
)

message("Done! Files saved: conserved_markers_dotplot.pdf / .png")
