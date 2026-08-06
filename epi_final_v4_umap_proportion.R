#############################################
# Epithelial Figure — Panel A (UMAP) + Panel B (Species Proportion)
# Matches Fig 4 style: labeled UMAP + species-split UMAP + stacked
# proportion bar with cell counts
#############################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

dir.create("DEG_comparisons", showWarnings = FALSE)

## ---- Parameters ----
INPUT_RDS   <- "epi_clean_v4_final.rds"
CLUSTER_COL <- "epi_subtype_v4"
SPECIES_COL <- "species"
SPECIES_ORDER <- c("Human", "Marmoset", "Mouse", "Rat")   # edit if names differ

## ---- Load ----
epi_clean_final <- readRDS(INPUT_RDS)
DefaultAssay(epi_clean_final) <- "RNA"
epi_clean_final <- JoinLayers(epi_clean_final)

## ---- Final confirmed identity order (11 categories) ----
identity_order <- c(
  "Basal Keratinocytes",
  "Basal-to-Suprabasal Transitional",
  "Granular-layer Keratinocytes",
  "Keratinized/Barrier Keratinocytes",
  "Late Differentiating",
  "Mature Non-keratinized",
  "Minor Salivary Gland Cells",
  "Proliferating Keratinocytes",
  "Sensory-associated epithelial cells",
  "Stress Responsive",
  "Suprabasal"
)

raw_labels <- as.character(epi_clean_final@meta.data[[CLUSTER_COL]])
mismatched <- setdiff(unique(raw_labels), identity_order)
if (length(mismatched) > 0) {
  cat("ERROR - labels present in data but NOT in identity_order:\n"); print(mismatched)
  stop("Fix label mismatch before proceeding.")
}
epi_clean_final@meta.data[[CLUSTER_COL]] <- factor(raw_labels, levels = identity_order)
Idents(epi_clean_final) <- CLUSTER_COL

identity_colors <- setNames(scales::hue_pal()(length(identity_order)), identity_order)

#############################################
## PANEL A — Integrated UMAP, labeled directly + legend
#############################################

p_umap <- DimPlot(
  epi_clean_final, reduction = "umap", group.by = CLUSTER_COL,
  cols = identity_colors, label = TRUE, label.size = 3, repel = TRUE
) +
  ggtitle("Epithelial Atlas — Final Annotation") +
  theme(legend.position = "right")

cairo_pdf("epi_fig_panelA_umap.pdf", width = 11, height = 7)
print(p_umap)
dev.off()

## ---- Species-split UMAP ----
p_umap_species <- DimPlot(
  epi_clean_final, reduction = "umap", group.by = CLUSTER_COL,
  cols = identity_colors, split.by = SPECIES_COL, label = FALSE
) +
  ggtitle("Epithelial Atlas — Final Annotation, split by species") +
  theme(legend.position = "bottom")

cairo_pdf("epi_fig_panelA_umap_by_species.pdf", width = 20, height = 6)
print(p_umap_species)
dev.off()

#############################################
## PANEL B — Species proportion stacked bar (with cell counts)
#############################################

prop_df <- epi_clean_final@meta.data %>%
  count(!!sym(SPECIES_COL), !!sym(CLUSTER_COL)) %>%
  rename(species = !!sym(SPECIES_COL), subtype = !!sym(CLUSTER_COL)) %>%
  mutate(species = factor(species, levels = SPECIES_ORDER),
         subtype = factor(subtype, levels = rev(identity_order))) %>%   # reversed so top-of-stack = first in legend
  group_by(species) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

write.csv(prop_df, "DEG_comparisons/epi_fig_panelB_proportion_table.csv", row.names = FALSE)

p_prop <- ggplot(prop_df, aes(x = species, y = pct, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.5, color = "black") +
  scale_fill_manual(values = identity_colors, breaks = identity_order, name = NULL) +
  labs(title = "Epithelial Subtype Proportion by Species", x = NULL, y = "Cell proportion (%)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 9)
  )

cairo_pdf("epi_fig_panelB_proportion.pdf", width = 9, height = 8)
print(p_prop)
dev.off()

cat("\nDone. Outputs:\n")
cat(" - epi_fig_panelA_umap.pdf\n")
cat(" - epi_fig_panelA_umap_by_species.pdf\n")
cat(" - epi_fig_panelB_proportion.pdf\n")
cat(" - DEG_comparisons/epi_fig_panelB_proportion_table.csv\n")
