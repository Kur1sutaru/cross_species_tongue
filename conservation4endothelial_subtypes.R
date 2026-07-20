# =============================================================================
# Full Conservation Analysis — Endothelial (EC) Subtypes
# Adapted from conservation4immunesubtypes_june3rd.R
#
# KEY DIFFERENCE FROM IMMUNE PIPELINE: Rat is fully excluded (not just
# filtered post-hoc) — species_list has 3 entries, "Conserved" now means
# significant in ALL 3 species (Human, Marmoset, Mouse), and every
# Rat-specific column/comparison/palette entry has been removed.
# =============================================================================

# =============================================================================
# CONFIG — check/edit these four lines before running
# =============================================================================
setwd("/r_workspace/cristal_data/cross_species/subtypes_cross_species/endothelial/endo_reintegration")
source("utils.R")   # loads theme_pub, outdir, palettes (shared across compartments)

OBJ_NAME    <- "endo_clean"     # your EC_CCA_annotated_v4_no_rat.rds object, loaded into this name
SUBTYPE_COL <- "cca_subtypes"   # metadata column holding the 11 EC subtype labels — confirm this
                                 # matches what you used for EC (e.g. colnames(endo_clean@meta.data))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(tibble)
})

LOG2FC      <- 0.6
PVAL_CUTOFF <- 0.05
MIN_PCT     <- 0.20

# ── Rat excluded here — 3-species pipeline throughout ───────────────────────
species_list <- c("Human", "Marmoset", "Mouse")
fc_cols      <- paste0(species_list, "_avg_log2FC")
pv_cols      <- paste0(species_list, "_p_val_adj")

endo_clean <- get(OBJ_NAME)

# Defensive: drop any stray Rat cells even though the object is already
# supposed to be Rat-free (EC_CCA_annotated_v4_no_rat.rds)
if ("Rat" %in% unique(endo_clean$species)) {
  cat("NOTE: Rat cells found in object despite _no_rat naming — removing them.\n")
  endo_clean <- subset(endo_clean, subset = species != "Rat")
  endo_clean$species <- droplevels(factor(endo_clean$species))
}

Idents(endo_clean) <- SUBTYPE_COL
subtypes <- levels(factor(endo_clean[[SUBTYPE_COL]][[1]]))
subtypes <- subtypes[!is.na(subtypes)]
cat(sprintf("EC subtypes found (%d): %s\n", length(subtypes), paste(subtypes, collapse = ", ")))

# EC-specific display order for figures. Edit this vector to match your
# established 11-subtype naming/order (EC_Arterial, EC_Arterial_transitional,
# EC_Activated_Capillary, etc.) — falls back to data-derived order if you
# haven't customized it, but a fixed order keeps figures reproducible.
if (!exists("subtype_order_ec")) {
  subtype_order_ec <- subtypes
  cat("NOTE: subtype_order_ec not defined in utils.R — using data-derived order.\n")
}

# EC subtype color palette — edit hex codes / names to match your 11 subtypes
subtype_colors_ec <- setNames(
  scales::hue_pal()(length(subtype_order_ec)),
  subtype_order_ec
)

# =============================================================================
# STEP 1 — FindConservedMarkers per subtype (Human/Marmoset/Mouse only)
# =============================================================================
cat(sprintf("\n── Step 1: FindConservedMarkers (min.pct = %.2f, no Rat) ──\n", MIN_PCT))

all_conserved <- list()

for (ct in subtypes) {
  cat(sprintf("  Processing: %s\n", ct))

  cell_counts <- table(
    endo_clean$species[endo_clean[[SUBTYPE_COL]][[1]] == ct]
  )
  cell_counts <- cell_counts[species_list]                 # enforce order, drop any non-target species
  cell_counts[is.na(cell_counts)] <- 0
  valid_species <- names(cell_counts)[cell_counts >= 3]

  if (length(valid_species) < 2) {
    cat("    Skipping — fewer than 2 species with >=3 cells\n")
    next
  }

  tryCatch({
    markers <- FindConservedMarkers(
      endo_clean,
      ident.1         = ct,
      grouping.var    = "species",
      only.pos        = TRUE,
      min.pct         = MIN_PCT,
      logfc.threshold = LOG2FC,
      verbose         = FALSE
    )

    if (nrow(markers) > 0) {
      markers$gene    <- rownames(markers)
      markers$subtype <- ct
      all_conserved[[ct]] <- markers
      cat(sprintf("    Found %d candidate markers\n", nrow(markers)))
    }
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
  })
}

conserved_raw <- bind_rows(all_conserved)

write.csv(conserved_raw,
          file.path(outdir, "ec_step1_pct20_conserved_markers.csv"),
          row.names = FALSE)
cat(sprintf("\nStep 1: %d total candidate markers\n", nrow(conserved_raw)))

# =============================================================================
# STEP 2 — Conservation filter (3-species logic)
# =============================================================================
cat("\n── Step 2: Conservation filter ──\n")

fc_cols_present <- intersect(fc_cols, colnames(conserved_raw))
pv_cols_present <- intersect(pv_cols, colnames(conserved_raw))
sp_names        <- gsub("_avg_log2FC", "", fc_cols_present)

conserved_raw$n_species_pass <- rowSums(
  sapply(seq_along(fc_cols_present), function(i) {
    fc <- conserved_raw[[fc_cols_present[i]]]
    pv <- conserved_raw[[pv_cols_present[i]]]
    (!is.na(fc) & !is.na(pv) & fc > LOG2FC & pv < PVAL_CUTOFF)
  })
)

n_sp <- length(sp_names)   # 3 for EC (vs. 4 for immune)

conserved_raw$conservation <- dplyr::case_when(
  conserved_raw$n_species_pass == n_sp ~ "Conserved",
  conserved_raw$n_species_pass == 1    ~ "Species-specific",
  conserved_raw$n_species_pass >= 2    ~ "Partially conserved",
  TRUE                                  ~ "Not significant"
)

for (sp in sp_names) {
  conserved_raw[[paste0(sp, "_pass")]] <- (
    !is.na(conserved_raw[[paste0(sp, "_avg_log2FC")]]) &
      !is.na(conserved_raw[[paste0(sp, "_p_val_adj")]]) &
      conserved_raw[[paste0(sp, "_avg_log2FC")]] > LOG2FC &
      conserved_raw[[paste0(sp, "_p_val_adj")]]  < PVAL_CUTOFF
  )
}

conserved_raw$specific_species <- apply(conserved_raw, 1, function(row) {
  if (row[["conservation"]] != "Species-specific") return(NA_character_)
  pass_cols <- paste0(sp_names, "_pass")
  sp_pass   <- sp_names[unlist(row[pass_cols]) == "TRUE"]
  if (length(sp_pass) == 1) return(sp_pass) else return(NA_character_)
})

conserved_strict <- conserved_raw %>% filter(conservation == "Conserved")

write.csv(conserved_raw,
          file.path(outdir, "ec_step2_pct20_classified.csv"),
          row.names = FALSE)

n_conserved <- nrow(conserved_strict)
n_step1     <- nrow(conserved_raw)

cat(sprintf("Conserved (all %d species):  %d\n", n_sp, n_conserved))
cat(sprintf("Partially conserved:         %d\n",
            sum(conserved_raw$conservation == "Partially conserved")))
cat(sprintf("Species-specific:            %d\n",
            sum(conserved_raw$conservation == "Species-specific")))

# =============================================================================
# STEP 3 — Pearson r + Jaccard index (Human vs Marmoset, Human vs Mouse only)
# =============================================================================
cat("\n── Step 3: Pearson r + Jaccard index ──\n")

genes_conserved <- intersect(conserved_strict$gene, rownames(endo_clean))
cat("Conserved genes in object:", length(genes_conserved), "\n")

avg_list <- lapply(species_list, function(sp) {
  cells_sp <- colnames(endo_clean)[endo_clean$species == sp]
  sub_obj  <- subset(endo_clean, cells = cells_sp)
  Idents(sub_obj) <- SUBTYPE_COL
  AverageExpression(sub_obj, features = genes_conserved,
                    assay = "RNA", layer = "data",
                    verbose = FALSE)$RNA
})
names(avg_list) <- species_list

# Rat dropped from comparisons entirely
comparisons <- list(
  "Human vs Marmoset" = c("Human", "Marmoset"),
  "Human vs Mouse"    = c("Human", "Mouse")
)

# Pearson r
pearson_results <- lapply(subtypes, function(ct) {
  res <- sapply(names(comparisons), function(comp_name) {
    sps <- comparisons[[comp_name]]
    sp1 <- sps[1]; sp2 <- sps[2]
    if (!ct %in% colnames(avg_list[[sp1]]) ||
        !ct %in% colnames(avg_list[[sp2]])) return(NA)
    cor(avg_list[[sp1]][, ct], avg_list[[sp2]][, ct],
        method = "pearson", use = "complete.obs")
  })
  data.frame(subtype = ct, comparison = names(comparisons),
             pearson_r = res, stringsAsFactors = FALSE)
}) %>% bind_rows()

# Jaccard index
jaccard_results <- lapply(subtypes, function(ct) {
  markers_per_sp <- lapply(sp_names, function(sp) {
    conserved_raw %>%
      filter(subtype == ct,
             .data[[paste0(sp, "_pass")]] == TRUE) %>%
      pull(gene)
  })
  names(markers_per_sp) <- sp_names

  res <- sapply(names(comparisons), function(comp_name) {
    sps <- comparisons[[comp_name]]
    g1  <- markers_per_sp[[sps[1]]]
    g2  <- markers_per_sp[[sps[2]]]
    if (length(g1) == 0 || length(g2) == 0) return(NA)
    length(intersect(g1, g2)) / length(union(g1, g2))
  })
  data.frame(subtype = ct, comparison = names(comparisons),
             jaccard = res, stringsAsFactors = FALSE)
}) %>% bind_rows()

write.csv(pearson_results,
          file.path(outdir, "ec_step3_pct20_pearson_r.csv"),
          row.names = FALSE)
write.csv(jaccard_results,
          file.path(outdir, "ec_step3_pct20_jaccard.csv"),
          row.names = FALSE)

cat("Pearson r:\n")
print(pearson_results %>%
        pivot_wider(names_from = comparison, values_from = pearson_r) %>%
        mutate(across(where(is.numeric), \(x) round(x, 2))))

cat("\nJaccard:\n")
print(jaccard_results %>%
        pivot_wider(names_from = comparison, values_from = jaccard) %>%
        mutate(across(where(is.numeric), \(x) round(x, 2))))

# =============================================================================
# STEP 4 — Stability classification (quartile-based, on 3 species' log2FC SD)
# =============================================================================
cat("\n── Step 4: Stability classification ──\n")

conserved_strict$sd_log2fc <- apply(
  conserved_strict[, fc_cols_present, drop = FALSE], 1,
  function(x) sd(as.numeric(x), na.rm = TRUE)
)

T1 <- quantile(conserved_strict$sd_log2fc, 0.25, na.rm = TRUE)
T2 <- quantile(conserved_strict$sd_log2fc, 0.50, na.rm = TRUE)
T3 <- quantile(conserved_strict$sd_log2fc, 0.75, na.rm = TRUE)

cat(sprintf("Quartile thresholds: Q1=%.3f Q2=%.3f Q3=%.3f\n", T1, T2, T3))

conserved_strict$stability <- dplyr::case_when(
  conserved_strict$sd_log2fc <= T1 ~ "Highly stable",
  conserved_strict$sd_log2fc <= T2 ~ "Moderately stable",
  conserved_strict$sd_log2fc <= T3 ~ "Variable",
  TRUE                              ~ "Highly variable"
)

stability_summary <- conserved_strict %>%
  count(subtype, stability) %>%
  mutate(
    subtype   = factor(subtype,   levels = subtype_order_ec),
    stability = factor(stability, levels = c("Highly stable",
                                             "Moderately stable",
                                             "Variable",
                                             "Highly variable"))
  ) %>%
  filter(!is.na(subtype))

write.csv(conserved_strict,
          file.path(outdir, "ec_step4_pct20_stability.csv"),
          row.names = FALSE)

cat("Stability distribution:\n")
print(table(conserved_strict$subtype, conserved_strict$stability))

# =============================================================================
# STEP 5 — Top conserved markers per subtype (for dotplot panel)
# =============================================================================
cat("\n── Step 5: Top conserved markers ──\n")

top_markers <- conserved_strict %>%
  mutate(mean_fc = rowMeans(
    select(., all_of(fc_cols_present)), na.rm = TRUE)
  ) %>%
  group_by(subtype) %>%
  slice_max(order_by = mean_fc, n = 3) %>%
  ungroup() %>%
  arrange(factor(subtype, levels = subtype_order_ec), desc(mean_fc))

cat("Top 3 conserved markers per EC subtype:\n")
print(top_markers %>% select(subtype, gene, mean_fc, all_of(fc_cols_present)))

write.csv(top_markers,
          file.path(outdir, "ec_step5_pct20_top_markers.csv"),
          row.names = FALSE)

# =============================================================================
# PANELS
# =============================================================================
cat("\n── Building panels ──\n")

# ── Panel B: conserved gene count per subtype ────────────────────────────────
conserved_counts <- conserved_strict %>%
  count(subtype, name = "n") %>%
  mutate(subtype = factor(subtype, levels = subtype_order_ec)) %>%
  filter(!is.na(subtype))

p_B <- ggplot(conserved_counts,
              aes(x = subtype, y = n, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = subtype_colors_ec) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Number of Genes",
       title = "Conserved genes per endothelial subtype") +
  theme_pub(base_size = 9) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1,
                                          size = 7, face = "bold"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25))

# ── Panel C: Pearson r heatmap (2 comparison columns, no Rat) ───────────────
pearson_wide <- pearson_results %>%
  mutate(subtype = factor(subtype, levels = rev(subtype_order_ec))) %>%
  filter(!is.na(pearson_r))

p_C <- ggplot(pearson_wide,
              aes(x = comparison, y = subtype, fill = pearson_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(pearson_r, 2)), size = 2.8, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#9ECAE1", high = "#08519C",
                       midpoint = 0.5, limits = c(0, 1), name = "Pearson r") +
  scale_x_discrete(labels = c("Human vs\nMarmoset", "Human vs\nMouse")) +
  labs(x = NULL, y = NULL, title = "Pearson r") +
  theme_pub(base_size = 9) +
  theme(axis.text.x     = element_text(size = 8),
        axis.text.y     = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"),
        panel.border    = element_rect(color = "grey80",
                                       fill = NA, linewidth = 0.3))

# ── Panel D: Jaccard heatmap (2 comparison columns, no Rat) ─────────────────
jaccard_wide <- jaccard_results %>%
  mutate(subtype = factor(subtype, levels = rev(subtype_order_ec))) %>%
  filter(!is.na(jaccard))

p_D <- ggplot(jaccard_wide,
              aes(x = comparison, y = subtype, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(jaccard, 2)), size = 2.8, color = "black") +
  scale_fill_gradient2(low = "white", mid = "#A1D99B", high = "#006D2C",
                       midpoint = 0.5, limits = c(0, 1), name = "Jaccard\nIndex") +
  scale_x_discrete(labels = c("Human vs\nMarmoset", "Human vs\nMouse")) +
  labs(x = NULL, y = NULL, title = "Jaccard Index") +
  theme_pub(base_size = 9) +
  theme(axis.text.x     = element_text(size = 8),
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        legend.key.size = unit(0.35, "cm"),
        panel.border    = element_rect(color = "grey80",
                                       fill = NA, linewidth = 0.3))

# ── Panel E: Species-specific vs conserved per subtype per species ──────────
df_conserved_sp <- conserved_raw %>%
  filter(conservation == "Conserved") %>%
  count(subtype, name = "conserved")

pass_per_sp <- lapply(sp_names, function(sp) {
  conserved_raw %>%
    filter(.data[[paste0(sp, "_pass")]] == TRUE) %>%
    count(subtype, name = "n") %>%
    mutate(species = sp)
}) %>% bind_rows()

pass_per_sp <- pass_per_sp %>%
  left_join(df_conserved_sp, by = "subtype") %>%
  mutate(
    species_specific = n - conserved,
    species_specific = pmax(species_specific, 0),
    subtype = factor(subtype, levels = subtype_order_ec),
    species = factor(species, levels = species_list)   # Human, Marmoset, Mouse — no Rat
  ) %>%
  filter(!is.na(subtype))

bar_palette_ec <- c(
  "Conserved (all 3 species)" = "#2B83BA",
  "Species-specific"          = "#F46D43"
)

n_panels_e <- min(6, length(levels(pass_per_sp$subtype)))
p_E_list <- lapply(levels(pass_per_sp$subtype)[1:n_panels_e],
                   function(st) {
                     df_st <- pass_per_sp %>% filter(subtype == st)
                     ggplot(df_st) +
                       geom_bar(aes(x = species, y = n, fill = "Species-specific"),
                                stat = "identity", width = 0.7) +
                       geom_bar(aes(x = species, y = conserved, fill = "Conserved (all 3 species)"),
                                stat = "identity", width = 0.7) +
                       geom_text(aes(x = species, y = n, label = n),
                                 vjust = -0.4, size = 2.2, fontface = "bold") +
                       geom_text(aes(x = species, y = conserved / 2, label = conserved),
                                 size = 2.2, color = "white", fontface = "bold") +
                       scale_fill_manual(values = bar_palette_ec, name = NULL) +
                       scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
                       labs(x = NULL, y = "No. of Genes", title = st) +
                       theme_pub(base_size = 8) +
                       theme(legend.position  = "none",
                             plot.title       = element_text(size = 7.5, face = "bold"),
                             axis.text.x      = element_text(size = 7, angle = 30, hjust = 1),
                             panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2))
                   })

p_E <- wrap_plots(p_E_list, ncol = 3) +
  plot_annotation(title = "Species-specific vs conserved genes per EC subtype")

# ── Panel F: Stability stacked barplot ───────────────────────────────────────
stability_palette <- c(
  "Highly stable"     = "#2166AC",
  "Moderately stable" = "#92C5DE",
  "Variable"          = "#FDDBC7",
  "Highly variable"   = "#B2182B"
)

p_F <- ggplot(stability_summary,
              aes(x = subtype, y = n, fill = stability)) +
  geom_bar(stat = "identity", width = 0.72,
           color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(n >= 2, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, color = "black", fontface = "bold") +
  scale_fill_manual(values = stability_palette, name = "Stability") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Number of Genes",
       title    = "Stability of expression classification",
       subtitle = paste0("Q1=", round(T1, 2),
                         " | Q2=", round(T2, 2),
                         " | Q3=", round(T3, 2))) +
  theme_pub(base_size = 9) +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1,
                                          size = 7, face = "bold"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
        legend.text        = element_text(size = 8),
        legend.key.size    = unit(0.4, "cm"),
        plot.subtitle      = element_text(size = 7.5, color = "grey40")) +
  guides(fill = guide_legend(ncol = 1))

# ── Panel G: Top conserved markers dotplot (3 species, no Rat row) ──────────
top_genes    <- unique(top_markers$gene)
top_genes_ok <- intersect(top_genes, rownames(endo_clean))

Idents(endo_clean) <- SUBTYPE_COL

dot_data <- lapply(species_list, function(sp) {
  cells_sp <- colnames(endo_clean)[endo_clean$species == sp]
  sub_obj  <- subset(endo_clean, cells = cells_sp)
  Idents(sub_obj) <- SUBTYPE_COL
  suppressWarnings(
    DotPlot(sub_obj, features = top_genes_ok, scale = FALSE)$data
  ) %>%
    rename(gene = features.plot, subtype_y = id) %>%
    mutate(species = sp)
}) %>% bind_rows() %>%
  mutate(
    gene      = factor(gene, levels = top_genes_ok),
    subtype_y = factor(subtype_y, levels = subtype_order_ec),
    species   = factor(species, levels = rev(species_list))
  )

p_G <- ggplot(dot_data,
              aes(x = gene, y = interaction(species, subtype_y),
                  size = pct.exp, color = avg.exp)) +
  geom_point() +
  scale_color_gradient(low = "#FFE5D9", high = "#B22222",
                       name = "Avg. Expression\n(log\u2082 FC)") +
  scale_size_continuous(name = "% Expressing",
                        range = c(0.5, 5),
                        breaks = c(25, 50, 75)) +
  scale_y_discrete(labels = function(x) {
    sp  <- gsub("\\..*$", "", x)
    sub <- gsub("^[^.]*\\.", "", x)
    ifelse(sp == "Human", paste0(sub, " — ", sp), sp)
  }) +
  labs(x = NULL, y = NULL,
       title = "Suggested conserved markers per endothelial subtype across species") +
  theme_pub(base_size = 8) +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1,
                                        size = 7, face = "italic"),
        axis.text.y      = element_text(size = 6.5),
        legend.key.size  = unit(0.35, "cm"),
        legend.text      = element_text(size = 7),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.2))

# =============================================================================
# COMBINE
# =============================================================================
cat("\n── Combining panels ──\n")

row1 <- p_B | p_C | p_D +
  plot_layout(widths = c(1.2, 1, 0.8))

row2 <- p_E | p_F +
  plot_layout(widths = c(1.5, 1))

row3 <- p_G

fig_combined <- row1 / row2 / row3 +
  plot_layout(heights = c(0.8, 1, 1.2)) +
  plot_annotation(
    title    = "Conservation analysis — endothelial subtypes (Rat excluded)",
    subtitle = paste0(
      "Step 1: ", n_step1, " candidate genes  |  ",
      "Step 2: ", n_conserved, " conserved genes across ", n_sp, " species (Human, Marmoset, Mouse)  |  ",
      "log\u2082FC > ", LOG2FC, " | Bonferroni p < ", PVAL_CUTOFF
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
    )
  )

# =============================================================================
# SAVE
# =============================================================================
cairo_pdf(file.path(outdir, "ec_fig4_endothelial_pct20.pdf"),
          width = 22, height = 28)
print(fig_combined)
dev.off()

cairo_pdf(file.path(outdir, "ecB_conserved_counts.pdf"),    width = 10, height = 5); print(p_B); dev.off()
cairo_pdf(file.path(outdir, "ecC_pearson_heatmap.pdf"),     width = 6,  height = 7); print(p_C); dev.off()
cairo_pdf(file.path(outdir, "ecD_jaccard_heatmap.pdf"),     width = 6,  height = 7); print(p_D); dev.off()
cairo_pdf(file.path(outdir, "ecE_species_specific.pdf"),    width = 16, height = 10); print(p_E); dev.off()
cairo_pdf(file.path(outdir, "ecF_stability.pdf"),           width = 12, height = 6); print(p_F); dev.off()
cairo_pdf(file.path(outdir, "ecG_top_markers_dotplot.pdf"), width = 20, height = 16); print(p_G); dev.off()

cat("\n── All outputs saved ────────────────────────────────────────────────\n")
cat("  ec_step1_pct20_conserved_markers.csv\n")
cat("  ec_step2_pct20_classified.csv\n")
cat("  ec_step3_pct20_pearson_r.csv\n")
cat("  ec_step3_pct20_jaccard.csv\n")
cat("  ec_step4_pct20_stability.csv\n")
cat("  ec_step5_pct20_top_markers.csv\n")
cat("  ec_fig4_endothelial_pct20.pdf  (combined)\n")
cat("  ecB–ecG individual panel PDFs\n")
