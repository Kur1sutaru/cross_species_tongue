## ============================================================
## Integration Quality Assessment — Cross-Species Epithelial Atlas
## Cristal Villalba Silva | Ruparel Lab | UT Health San Antonio
## ============================================================
## Evaluates 7 integration strategies from IntegrationMethods_epithelial.pptx
## Metrics: iLISI, cLISI, ASW (batch + bio), per-cluster species coverage,
##           silhouette scores, and a ranked summary table.
## ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cluster)    # silhouette()
library(lisi)       # compute_lisi() — install: devtools::install_github("immunogenomics/lisi")

## ─────────────────────────────────────────────────────────────
## 0. SETUP — Load your integration objects
## ─────────────────────────────────────────────────────────────
## Each object should already be:
##   • Scaled + PCA + UMAP computed
##   • Metadata columns: $species, $orig.ident, $cell_type
##   • Clusters stored in $seurat_clusters (or re-cluster below)
##
## Replace these paths / object names with your actual ones.

# OPTION A: Load saved RDS files
# obj_list <- list(
#   "1_Original"        = readRDS("path/to/obj1.rds"),
#   "2_SpeciesSplit_noNorm" = readRDS("path/to/obj2.rds"),
#   "3_13ind_1step_CCA" = readRDS("path/to/obj3.rds"),
#   "4_Hierarchical_CCA"= readRDS("path/to/obj4.rds"),
#   "5_Hier_CCA_code"   = readRDS("path/to/obj5.rds"),  # same as 4 — skip if identical
#   "6_Hier_CCA_259mk"  = readRDS("path/to/obj6.rds"),
#   "7_Harmony"         = readRDS("path/to/obj7.rds")
# )

# OPTION B: Objects already in your environment
obj_list <- list(
  "1_Original"           = obj_original,
  "2_SpeciesSplit_noNorm"= obj_no_norm,
  "3_13ind_1step"        = combined_1step,
  "4_Hierarchical_CCA"   = combined_hier,
  "6_Hier_CCA_259mk"     = combined_hier_markers,
  "7_Harmony"            = obj_harmony
)

## Which embedding to use for distance/LISI calculations
REDUCTION  <- "umap"   # or "pca", "harmony"
N_DIMS     <- 2        # number of dims from that reduction to use
BATCH_COL  <- "orig.ident"   # individual-level batch
SPECIES_COL<- "species"      # species column
CELLTYPE_COL <- "cell_type"  # annotated cell type column

## ─────────────────────────────────────────────────────────────
## 1. HELPER FUNCTIONS
## ─────────────────────────────────────────────────────────────

# Extract embedding matrix from a Seurat object
get_embedding <- function(obj, reduction = REDUCTION, n_dims = N_DIMS) {
  emb <- Embeddings(obj, reduction = reduction)
  emb[, 1:min(n_dims, ncol(emb))]
}

# Compute LISI scores (iLISI for batch, cLISI for cell type)
compute_lisi_scores <- function(obj, reduction = REDUCTION) {
  emb  <- get_embedding(obj, reduction, n_dims = 2)
  meta <- obj@meta.data[, c(BATCH_COL, SPECIES_COL, CELLTYPE_COL), drop = FALSE]
  lisi_res <- compute_lisi(emb, meta, c(BATCH_COL, SPECIES_COL, CELLTYPE_COL))
  list(
    iLISI_individual = mean(lisi_res[[BATCH_COL]], na.rm = TRUE),
    iLISI_species    = mean(lisi_res[[SPECIES_COL]], na.rm = TRUE),
    cLISI_celltype   = mean(lisi_res[[CELLTYPE_COL]], na.rm = TRUE),
    lisi_df          = lisi_res
  )
}

# Silhouette score by a given grouping (cell type or species)
compute_silhouette <- function(obj, group_col, reduction = REDUCTION, n_dims = 30) {
  emb    <- get_embedding(obj, reduction, n_dims = n_dims)
  labels <- as.numeric(as.factor(obj@meta.data[[group_col]]))
  if (length(unique(labels)) < 2) return(NA)
  d   <- dist(emb)
  sil <- silhouette(labels, d)
  mean(sil[, 3], na.rm = TRUE)
}

# Per-cluster species coverage: fraction of clusters with ≥ 3 species
cluster_species_coverage <- function(obj) {
  tbl <- table(obj$seurat_clusters, obj@meta.data[[SPECIES_COL]])
  n_species <- apply(tbl, 1, function(x) sum(x > 0))
  list(
    mean_species_per_cluster = mean(n_species),
    pct_clusters_3plus       = mean(n_species >= 3) * 100,
    cluster_table            = tbl
  )
}

## ─────────────────────────────────────────────────────────────
## 2. RUN ALL METRICS ACROSS INTEGRATION METHODS
## ─────────────────────────────────────────────────────────────

results <- list()

for (nm in names(obj_list)) {
  cat("\n─── Evaluating:", nm, "───\n")
  obj <- obj_list[[nm]]

  # Cluster if not done yet
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    obj <- FindNeighbors(obj, reduction = REDUCTION, dims = 1:N_DIMS)
    obj <- FindClusters(obj, resolution = 0.5)
    obj_list[[nm]] <- obj
  }

  # LISI
  cat("  → Computing LISI...\n")
  lisi_out <- tryCatch(compute_lisi_scores(obj), error = function(e) {
    cat("  ⚠ LISI failed:", conditionMessage(e), "\n")
    list(iLISI_individual = NA, iLISI_species = NA, cLISI_celltype = NA, lisi_df = NULL)
  })

  # Silhouette (use PCA dims for distance — more informative than 2D UMAP)
  cat("  → Computing silhouette...\n")
  n_pca <- min(30, ncol(Embeddings(obj, "pca")))
  asw_bio   <- tryCatch(compute_silhouette(obj, CELLTYPE_COL, "pca", n_pca), error = function(e) NA)
  asw_batch <- tryCatch(compute_silhouette(obj, SPECIES_COL,  "pca", n_pca), error = function(e) NA)

  # Cluster species coverage
  cat("  → Computing cluster species coverage...\n")
  cov_out <- tryCatch(cluster_species_coverage(obj), error = function(e) {
    list(mean_species_per_cluster = NA, pct_clusters_3plus = NA, cluster_table = NULL)
  })

  results[[nm]] <- data.frame(
    Method                   = nm,
    iLISI_individual         = round(lisi_out$iLISI_individual, 3),
    iLISI_species            = round(lisi_out$iLISI_species, 3),
    cLISI_celltype           = round(lisi_out$cLISI_celltype, 3),
    ASW_bio_celltype         = round(asw_bio, 3),
    ASW_batch_species        = round(asw_batch, 3),
    mean_species_per_cluster = round(cov_out$mean_species_per_cluster, 2),
    pct_clusters_3plus_spp   = round(cov_out$pct_clusters_3plus, 1),
    lisi_df_ref              = nm   # pointer for downstream plots
  )
}

# Combine into summary table
summary_df <- bind_rows(results)
cat("\n\n══════════════════════════════════════════\n")
cat("INTEGRATION QUALITY SUMMARY\n")
cat("══════════════════════════════════════════\n")
print(summary_df, row.names = FALSE)

## ─────────────────────────────────────────────────────────────
## 3. COMPOSITE SCORE (rank-based, higher = better)
## ─────────────────────────────────────────────────────────────
## Goal: high iLISI (species mixed) + high ASW_bio (cell types preserved)
##       + high cLISI (cell types well-separated) + high % clusters with ≥3 spp
##       + LOW ASW_batch (species NOT artificially separated)

score_df <- summary_df %>%
  mutate(
    rank_iLISI_sp    = rank(iLISI_species,           ties.method = "min"),
    rank_cLISI_ct    = rank(cLISI_celltype,           ties.method = "min"),
    rank_ASW_bio     = rank(ASW_bio_celltype,          ties.method = "min"),
    rank_ASW_batch_inv = rank(-ASW_batch_species,      ties.method = "min"),  # lower batch = better
    rank_coverage    = rank(pct_clusters_3plus_spp,   ties.method = "min"),
    composite_score  = rank_iLISI_sp + rank_cLISI_ct + rank_ASW_bio +
                       rank_ASW_batch_inv + rank_coverage
  ) %>%
  arrange(desc(composite_score))

cat("\n\n── RANKED BY COMPOSITE SCORE ──\n")
print(score_df[, c("Method", "composite_score",
                   "iLISI_species", "cLISI_celltype",
                   "ASW_bio_celltype", "ASW_batch_species",
                   "pct_clusters_3plus_spp")], row.names = FALSE)

cat("\n★ RECOMMENDED METHOD:", score_df$Method[1], "\n")

## ─────────────────────────────────────────────────────────────
## 4. PLOTS
## ─────────────────────────────────────────────────────────────

out_dir <- "integration_quality_plots"
dir.create(out_dir, showWarnings = FALSE)

## 4a. Summary heatmap-style dot plot of all metrics
long_df <- summary_df %>%
  select(Method, iLISI_species, cLISI_celltype, ASW_bio_celltype,
         ASW_batch_species, pct_clusters_3plus_spp) %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Value")

p_dot <- ggplot(long_df, aes(x = Metric, y = Method, size = Value, color = Value)) +
  geom_point() +
  scale_color_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                        midpoint = median(long_df$Value, na.rm = TRUE)) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal(base_size = 12) +
  labs(title = "Integration Quality Metrics — Cross-Species Epithelial Atlas",
       subtitle = "Larger/greener = better for iLISI, cLISI, ASW_bio, % coverage\n(ASW_batch lower = better batch removal)",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "01_metrics_dotplot.pdf"), p_dot,
       width = 9, height = 5, useDingbats = FALSE)

## 4b. Composite score bar chart
p_bar <- ggplot(score_df, aes(x = reorder(Method, composite_score),
                               y = composite_score, fill = composite_score)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = composite_score), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac") +
  theme_minimal(base_size = 12) +
  labs(title = "Composite Integration Score (rank-based)",
       subtitle = "Higher = better overall integration quality",
       x = NULL, y = "Composite Score") +
  theme(legend.position = "none") +
  expand_limits(y = max(score_df$composite_score) * 1.15)

ggsave(file.path(out_dir, "02_composite_score.pdf"), p_bar,
       width = 8, height = 4, useDingbats = FALSE)

## 4c. LISI violin plots per method (if lisi_df objects available)
lisi_plot_list <- list()
for (nm in names(obj_list)) {
  lisi_out <- tryCatch(compute_lisi_scores(obj_list[[nm]]), error = function(e) NULL)
  if (!is.null(lisi_out$lisi_df)) {
    df_tmp <- lisi_out$lisi_df %>%
      mutate(Method = nm) %>%
      pivot_longer(cols = c(all_of(SPECIES_COL), all_of(CELLTYPE_COL)),
                   names_to = "Metric", values_to = "LISI")
    lisi_plot_list[[nm]] <- df_tmp
  }
}

if (length(lisi_plot_list) > 0) {
  lisi_all <- bind_rows(lisi_plot_list)
  p_lisi <- ggplot(lisi_all, aes(x = Method, y = LISI, fill = Metric)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.3, position = position_dodge(0.9)) +
    facet_wrap(~Metric, scales = "free_y") +
    scale_fill_manual(values = c("#4393c3", "#d6604d")) +
    theme_minimal(base_size = 11) +
    labs(title = "LISI Score Distributions by Integration Method",
         subtitle = "iLISI (species): higher = better mixing | cLISI (cell type): higher = better separation",
         x = NULL, y = "LISI Score") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          legend.position = "none")
  ggsave(file.path(out_dir, "03_lisi_violins.pdf"), p_lisi,
         width = 12, height = 5, useDingbats = FALSE)
}

## 4d. UMAP grid colored by species — one panel per method
umap_plots <- lapply(names(obj_list), function(nm) {
  obj <- obj_list[[nm]]
  df  <- as.data.frame(Embeddings(obj, REDUCTION))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$species <- obj@meta.data[[SPECIES_COL]]
  ggplot(df, aes(UMAP_1, UMAP_2, color = species)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_manual(values = c(human   = "#E41A1C",
                                  marmoset= "#377EB8",
                                  mouse   = "#4DAF4A",
                                  rat     = "#FF7F00")) +
    theme_void(base_size = 9) +
    labs(title = nm) +
    guides(color = guide_legend(override.aes = list(size = 2)))
})

p_umap_grid <- wrap_plots(umap_plots, ncol = 3) +
  plot_annotation(title = "UMAP by Species — All Integration Methods",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(out_dir, "04_umap_species_grid.pdf"), p_umap_grid,
       width = 15, height = 10, useDingbats = FALSE)

## 4e. Per-cluster species composition bar plot — best method
best_obj <- obj_list[[score_df$Method[1]]]
cov      <- cluster_species_coverage(best_obj)

prop_df <- as.data.frame(prop.table(cov$cluster_table, margin = 1)) %>%
  rename(Cluster = Var1, Species = Var2, Proportion = Freq)

p_stack <- ggplot(prop_df, aes(x = Cluster, y = Proportion, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = c(human   = "#E41A1C",
                               marmoset= "#377EB8",
                               mouse   = "#4DAF4A",
                               rat     = "#FF7F00")) +
  theme_minimal(base_size = 11) +
  labs(title = paste("Per-Cluster Species Composition —", score_df$Method[1]),
       subtitle = "Ideal: each cluster has cells from all 4 species",
       x = "Cluster", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "05_cluster_species_composition.pdf"), p_stack,
       width = 10, height = 5, useDingbats = FALSE)

## ─────────────────────────────────────────────────────────────
## 5. SAVE SUMMARY TABLE AS CSV
## ─────────────────────────────────────────────────────────────
write.csv(score_df[, -ncol(score_df)],   # drop internal rank cols if desired
          file.path(out_dir, "integration_quality_summary.csv"),
          row.names = FALSE)

cat("\n✅ Done. Outputs saved to:", out_dir, "\n")
cat("   01_metrics_dotplot.pdf\n")
cat("   02_composite_score.pdf\n")
cat("   03_lisi_violins.pdf\n")
cat("   04_umap_species_grid.pdf\n")
cat("   05_cluster_species_composition.pdf\n")
cat("   integration_quality_summary.csv\n\n")

## ─────────────────────────────────────────────────────────────
## 6. QUICK DECISION GUIDE (printed to console)
## ─────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════\n")
cat("HOW TO INTERPRET YOUR RESULTS\n")
cat("══════════════════════════════════════════\n")
cat("
Metric              Want    Meaning
─────────────────── ─────── ─────────────────────────────────────────────────
iLISI (species)     HIGH    Cells from all species well-mixed in neighborhood
cLISI (cell type)   HIGH    Cell type labels well-separated (bio preserved)
ASW_bio (cell type) HIGH    Cell types form tight, distinct clusters
ASW_batch (species) LOW     Species NOT clustering together (batch removed)
% clusters ≥3 spp   HIGH    Most clusters contain cells from ≥3 species

Composite score = sum of ranks across all metrics (higher = better overall).

Typical red flags:
  • iLISI_species < 1.5  → species still strongly separated
  • ASW_batch_species > 0.3 → over-correction, species blobs remain
  • cLISI_celltype < 1.2  → cell types collapsed/blended
  • % clusters ≥3 spp < 50% → many species-specific clusters remain
")
