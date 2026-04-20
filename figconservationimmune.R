library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

DefaultAssay(immune) <- "RNA"

# ============================================================================
# 1. Read conserved gene counts from your files
# ============================================================================

conserved_all <- immune_all_conserved_markers

conserved_counts <- conserved_all %>%
  group_by(subtype) %>%
  summarise(n_conserved = n(), .groups = "drop")

# Add subtypes with 0 conserved genes
all_subtypes <- levels(immune$immune_subtype)
missing <- setdiff(all_subtypes, conserved_counts$subtype)
if (length(missing) > 0) {
  conserved_counts <- bind_rows(
    conserved_counts,
    data.frame(subtype = missing, n_conserved = 0)
  )
}

cat("=== Conserved genes per subtype ===\n")
print(conserved_counts)

# ============================================================================
# 2. Run FindAllMarkers per species to get total DEGs per species per subtype
# ============================================================================

cat("\n=== Finding DEGs per species ===\n")

# Map immune_celltype to immune_subtype
celltype_to_subtype <- setNames(
  levels(immune$immune_subtype),
  levels(immune$immune_celltype)
)

species_degs <- list()

for (sp in c("Human", "Mouse", "Marmoset", "Rat")) {
  cat("\nProcessing:", sp, "\n")
  
  sp_obj <- subset(immune, species == sp)
  Idents(sp_obj) <- "immune_subtype"
  
  # Drop subtypes with too few cells
  keep_ct <- names(table(sp_obj$immune_subtype))[table(sp_obj$immune_subtype) >= 5]
  
  if (length(keep_ct) < 2) {
    cat("  Too few subtypes, skipping\n")
    next
  }
  
  sp_obj <- subset(sp_obj, immune_subtype %in% keep_ct)
  sp_obj$immune_subtype <- droplevels(sp_obj$immune_subtype)
  
  tryCatch({
    sp_markers <- FindAllMarkers(sp_obj, only.pos = TRUE, min.pct = 0.1,
                                 logfc.threshold = 0.25, verbose = FALSE)
    
    sp_summary <- sp_markers %>%
      filter(p_val_adj < 0.05) %>%
      group_by(cluster) %>%
      summarise(n_total = n(), .groups = "drop") %>%
      mutate(species = sp) %>%
      rename(subtype = cluster)
    
    species_degs[[sp]] <- sp_summary
    cat("  Found DEGs for", nrow(sp_summary), "subtypes\n")
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

species_deg_df <- bind_rows(species_degs)

cat("\n=== DEGs per species per subtype ===\n")
print(species_deg_df)

# ============================================================================
# 3. Build plot data: conserved (same for all species) + species-specific
# ============================================================================

plot_df <- species_deg_df %>%
  left_join(conserved_counts, by = "subtype") %>%
  mutate(
    n_conserved = ifelse(is.na(n_conserved), 0, n_conserved),
    n_specific = pmax(n_total - n_conserved, 0),
    total = n_total
  )

plot_df$subtype <- factor(plot_df$subtype, levels = levels(immune$immune_subtype))
plot_df$species <- factor(plot_df$species, levels = c("Human", "Mouse", "Marmoset", "Rat"))

# Pivot for stacking
plot_long <- plot_df %>%
  dplyr::select(subtype, species, n_conserved, n_specific) %>%
  pivot_longer(cols = c(n_conserved, n_specific),
               names_to = "category", values_to = "count") %>%
  mutate(category = recode(category,
                           "n_conserved" = "Conserved (all species)",
                           "n_specific"  = "Species-specific"
  ),
  category = factor(category, levels = c("Species-specific", "Conserved (all species)")))

# ============================================================================
# 4. Plot — faceted by subtype (matching your example figure)
# ============================================================================

p <- ggplot(plot_long, aes(x = species, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  # Conserved count label inside blue bar
  geom_text(data = plot_df, aes(x = species, y = n_conserved / 2, 
                                label = n_conserved, fill = NULL),
            size = 2.5, color = "white", fontface = "bold") +
  # Total label on top
  geom_text(data = plot_df, aes(x = species, y = total,
                                label = total, fill = NULL),
            vjust = -0.5, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Conserved (all species)" = "#6BAED6",
    "Species-specific"        = "#FC8D62"
  )) +
  facet_wrap(~ subtype, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "No. of Genes", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.5, "lines"))

pdf("Immune_Conserved_vs_Specific_PerSpecies.pdf", width = 16, height = 14)
print(p)
dev.off()

print(p)

cat("\nDone!\n")