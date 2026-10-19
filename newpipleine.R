################################################################################
##  Cross-species SIMILARITY / DIFFERENCES pipeline  --  immune subtypes
##  (implements "For Similarity/Differences" document, panels A / B / C)
##
##  Per cell type / subtype, using 1:1 orthologs shared by Human, Marmoset, Mouse:
##
##   1) GENE REPERTOIRE SIMILARITY  (Panel A)
##      pseudobulk per biological sample -> CPM -> gene "expressed" in a species
##      if CPM > 1 in >= 50 % of that species' samples -> Jaccard index for
##      Human-Mouse, Human-Marmoset, Marmoset-Mouse (+ bar of n genes, 3 Venns)
##
##   2) FUNCTIONAL / TRANSCRIPTIONAL PROGRAM SIMILARITY  (Panel B)
##      GSVA (or ssGSEA) of the 50 MSigDB Hallmark pathways on every pseudobulk
##      sample -> mean activity per species -> pathway-wise z-scores ->
##      species effect (one-way ANOVA) + pairwise tests (H-M, H-Mar, Mar-M) ->
##      rank most CONSERVED / most DIVERGENT -> top-10 heatmaps
##
##   3) SIMILARITY OF PATHWAY PROGRAMS TO HUMAN  (Panel C)
##      Pearson r of the per-species mean pathway scores (all pathways):
##      Human vs Marmoset, Human vs Mouse (bootstrap 95 % CI) + scatter plots
##
##  Outputs: <OUT_DIR>/tables/*.csv  and  <OUT_DIR>/figures/*.pdf|png
##           one composite A|B|C figure per subtype + cross-subtype summaries
##
##  Standalone script -- all parameters inlined below, nothing sourced.
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(GSVA)
  library(msigdbr)
})

## =============================================================================
## 0. PARAMETERS  -------------------------------------------------------------
## =============================================================================

immune_PATH      <- "immune_subtypesannotmarch30th.rds"   # <- immune_clean immuneect (Rat already excluded, or filtered below)
SUBTYPE_COL   <- "immune_subtype_v2"                    # final immune subtype labels
SPECIES_COL   <- "species"                              # Human / Marmoset / Mouse
SAMPLE_COL    <- "samples"                              # biological sample (created below from orig.ident)
SPECIES_LEVELS <- c("Human", "Marmoset", "Mouse")       # order used in all figures

## T-cell subtypes are analysed separately (no merge).  Set TRUE only if a
## pooled "T cells" run is wanted in addition.
MERGE_T_CELLS <- FALSE
T_CELL_LABELS <- c("Naive/CM T cells", "Th17-like T cells", "CD8+ Cytotoxic T cells")

## Subtypes to run (NULL = every subtype present after merging).
## A subtype is skipped only when a species has < MIN_CELLS_PER_SPECIES cells
## (e.g. Plasma cells = Human-only, Scavenger Mph = no Marmoset).
SUBTYPES_TO_RUN <- NULL

## Ortholog universe.  The integrated immuneect was built on 1:1 orthologs shared by
## Human, Marmoset and Mouse, so its gene space IS the ortholog set: leave NULL.
## (Optional: a CSV with columns human_symbol, marmoset_symbol, mouse_symbol,
##  one2one can further restrict it -- see build_ortholog_table_biomart() at the end.)
ORTHOLOG_TABLE <- NULL

## Pseudobulk design (species level is the PRIMARY unit, per PI request):
##  * SPECIES pseudobulk  = all cells of the species for this subtype, summed. This is
##    what Panel A (expressed genes), Panel B (pathway scores / z) and Panel C (r) use.
##    A subtype is skipped only if a species has < MIN_CELLS_PER_SPECIES cells.
##  * SAMPLE pseudobulks  = per biological sample (>= MIN_CELLS_PER_PSEUDOBULK cells).
##    Used ONLY to add replicate-level statistics (ANOVA / t-tests / bootstrap CI)
##    when every species has >= MIN_SAMPLES_FOR_STATS of them; otherwise those are NA.
MIN_CELLS_PER_SPECIES    <- 10   # species total cells needed to include the subtype
MIN_CELLS_PER_PSEUDOBULK <- 5    # per-sample pseudobulk (statistics only)
MIN_SAMPLES_FOR_STATS    <- 2
CPM_THRESHOLD            <- 15   # PI choice (Oct 19 2026); sensitivity curve in SUMMARY_cpm_cutoff_sensitivity
FRAC_SAMPLES_EXPRESSED   <- 0.50 # only used when REPERTOIRE_MODE == "per_sample"
## Sensitivity of Panel A to the CPM cutoff: n expressed genes + Jaccard per pair are
## recomputed at each of these cutoffs -> SUMMARY_cpm_cutoff_sensitivity.csv / .pdf.
## The CPM distribution per species (quantiles, max, top genes) is also written per
## subtype so the cutoff can be judged against the real dynamic range.
CPM_SENSITIVITY          <- c(1, 2, 5, 10, 15, 20, 50, 100)

## Pathway scoring.  "ssgsea" (default) gives ABSOLUTE per-sample activity, which is
## required for Panel C: GSVA scores are centred across samples per pathway, so with
## GSVA the species means are anti-correlated by construction (r < 0 is an artefact).
GSVA_METHOD <- "ssgsea"          # "ssgsea" or "gsva"
HALLMARK_SPECIES <- "Homo sapiens"   # gene symbols in the immuneect are human-style
GSVA_MIN_SET <- 5
GSVA_MAX_SET <- 500

## z-score mode for the heatmaps
##   "within_species": z of the species' mean scores ACROSS the 50 pathways (per species
##                     column). A conserved pathway then shows the same colour in all
##                     three columns, a divergent one does not -- the mock-up's look.
##   "species_means" : z across the 3 species means (literal reading of the doc). NOTE:
##                     z of 3 values is always ~+/-1 whatever their spread, so conserved
##                     rows look as contrasted as divergent ones.
##   "samples"       : z per pathway across all pseudobulk samples, then mean per species
Z_MODE  <- "within_species"
N_TOP   <- 10                    # top-N conserved / divergent pathways shown
DIVERGENT_PADJ <- 0.05           # "statistically supported differences"

## How conserved / divergent are RANKED (must match what the heatmap shows):
##   "z_range" (default): effect size on the displayed scale = max(z) - min(z) across
##             the 3 species. Divergent = largest range; conserved = smallest range.
##             ANOVA BH p is kept as a secondary "statistically supported" flag.
##   "F"     : ANOVA F statistic (old behaviour). NOTE: F rewards tight replicates,
##             so pathways with tiny but reproducible differences rank as divergent
##             even when their z-scores look identical.
RANK_BY <- "z_range"
## Conserved pathways are drawn from pathways ACTIVE in all three species (z > 0 in
## each) so the list shows shared programs (IFN, TNF...) rather than flat, low ones.
CONSERVED_REQUIRE_ACTIVE <- TRUE
## Divergent pathways are drawn from pathways active (z > 0) in AT LEAST ONE species,
## so the list is not dominated by programs that are low everywhere (bile acid,
## pancreas beta cells, angiogenesis...) whose rank differences are mostly noise.
DIVERGENT_REQUIRE_ACTIVE_ANY <- TRUE

## How "expressed" is defined for Panel A:
##   "species_pooled" (default): all cells of the species summed into one pseudobulk,
##             expressed = CPM > CPM_THRESHOLD in that profile (works for every subtype
##             where each species has cells; the per-sample fraction is still reported)
##   "per_sample": CPM > CPM_THRESHOLD in >= FRAC_SAMPLES_EXPRESSED of the species' samples
REPERTOIRE_MODE <- "species_pooled"

## Genes used for the pathway scores (steps 2 and 3) = the Panel A expressed genes:
##   "union"  (default): genes expressed (CPM > CPM_THRESHOLD) in AT LEAST ONE species.
##             One common list for all three species -- required for ssGSEA ranks to be
##             comparable; a gene expressed only in human still ranks at the bottom in
##             mouse and lowers that pathway's mouse score (that IS the species difference).
##   "shared":  genes expressed in ALL three species (the Venn intersection). Stricter;
##             species differences then come only from expression level, not presence.
##   "per_species": EACH species is scored on ITS OWN expressed-gene list (e.g. 10,852
##             human / 9,670 marmoset / 8,616 mouse genes), one ssGSEA run per species.
##             Scores are comparable through the within-species z-scores (Z_MODE =
##             "within_species") and Pearson r, both of which are scale-free; raw scores
##             of different species should not be compared directly in this mode.
PATHWAY_GENE_SET <- "per_species"

## Expected Hallmark programs per subtype (PI list). For each subtype the pipeline
## reports where these sit: activity rank within each species (1 = highest of 50),
## z-scores, active in all species, conserved / divergent rank -> Panel D + CSV.
## Names must match the subtype labels in SUBTYPE_COL (after T-cell merge, if any).
EXPECTED_PATHWAYS <- list(
  "Naive/CM T cells" = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM",
                         "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING"),
  "Th17-like T cells" = c("HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                          "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS", "HALLMARK_IL2_STAT5_SIGNALING"),
  "CD8+ Cytotoxic T cells" = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL2_STAT5_SIGNALING",
                               "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS", "HALLMARK_ALLOGRAFT_REJECTION",
                               "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  "Proliferating T cells" = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
                              "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_DNA_REPAIR", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_GLYCOLYSIS"),
  "Tissue-resident Macrophages" = c("HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                                    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                                    "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  "B cells" = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING"),
  "Mature Mast Cells (MCTC)" = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                                 "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                 "HALLMARK_HYPOXIA"),
  "Classical Neutrophils" = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_COMPLEMENT",
                              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
                              "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_APOPTOSIS")
)

## Correlation bootstrap (Panel C error bars)
BOOT_N  <- 1000
SEED    <- 1234

## Species colours -- manuscript palette (ShinyCell / Fig 2 convention).
## To reproduce the mock-up colours instead use:
##   c(Human = "#4C8DF5", Marmoset = "#F28B82", Mouse = "#57BB6A")
SPECIES_COLS <- c(Human = "#1D9E75", Marmoset = "#377EB8", Mouse = "#E41A1C")

OUT_DIR <- "similarity_differences_immune"

## =============================================================================
## 1. SET-UP  -----------------------------------------------------------------
## =============================================================================

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)
## path helpers -- (re)create the folders on every call, so the script can be run
## section by section or after the output folder was deleted
TAB <- function(...) { d <- file.path(OUT_DIR, "tables");  dir.create(d, showWarnings = FALSE, recursive = TRUE); file.path(d, ...) }
FIG <- function(...) { d <- file.path(OUT_DIR, "figures"); dir.create(d, showWarnings = FALSE, recursive = TRUE); file.path(d, ...) }

save_fig <- function(p, name, w, h) {
  cairo_pdf(FIG(paste0(name, ".pdf")), width = w, height = h); print(p); dev.off()
  ggsave(FIG(paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x); x <- tolower(gsub("_", " ", x))
  x <- paste0(toupper(substring(x, 1, 1)), substring(x, 2))
  ## keep common acronyms readable
  for (a in c("Il2", "Il6", "Stat3", "Stat5", "Tnfa", "Nfkb", "Kras", "Mtorc1", "Pi3k",
              "Akt", "Mtor", "Tgf", "Wnt", "Uv", "Dna", "G2m", "E2f", "Myc", "Jak"))
    x <- gsub(paste0("\\b", a, "\\b"), toupper(a), x)
  x <- gsub("Tnfa signaling via NFKB", "TNF signaling via NFKB", x)
  x
}
theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(axis.text = element_text(colour = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = base + 1),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", size = base - 1),
          legend.key.size = unit(0.35, "cm"))
}

message(">>> Loading immuneect: ", immune_PATH)
immune <- readRDS(immune_PATH)
DefaultAssay(immune) <- "RNA"
immune <- JoinLayers(immune)

## biological-sample column: human libraries tsXX_ant / tsXX_post collapse to the
## donor (ts14, ts25, ts27, ts7); every other orig.ident is kept as is
if (!SAMPLE_COL %in% colnames(immune@meta.data)) {
  immune[[SAMPLE_COL]] <- unname(sub("^(ts\\d+)_(ant|post)$", "\\1", as.character(immune$orig.ident)))
}
print(table(immune[[SAMPLE_COL, drop = TRUE]], immune[[SPECIES_COL, drop = TRUE]]))

## keep the 3 species only (drops Rat if still present)
immune <- subset(immune, cells = colnames(immune)[immune[[SPECIES_COL, drop = TRUE]] %in% SPECIES_LEVELS])

## T-cell merge (Oct 15 decision)
subtype_vec <- as.character(immune[[SUBTYPE_COL, drop = TRUE]])
if (MERGE_T_CELLS) subtype_vec[subtype_vec %in% T_CELL_LABELS] <- "T cells"
immune$subtype_run <- unname(subtype_vec)

md <- immune@meta.data %>%
  transmute(cell = rownames(.),
            species = factor(.data[[SPECIES_COL]], levels = SPECIES_LEVELS),
            sample  = as.character(.data[[SAMPLE_COL]]),
            subtype = subtype_run) %>%
  filter(!is.na(subtype), !is.na(sample), !is.na(species))

## quick sanity print: samples per species per subtype
samp_tab <- md %>% distinct(subtype, species, sample) %>% count(subtype, species, name = "n_samples") %>%
  pivot_wider(names_from = species, values_from = n_samples, values_fill = 0)
print(samp_tab)
write.csv(samp_tab, TAB("samples_per_species_per_subtype.csv"), row.names = FALSE)

## cells per sample per subtype (BEFORE any filter) -- explains why a subtype is skipped
cells_tab <- md %>% count(subtype, species, sample, name = "n_cells") %>%
  arrange(subtype, species, sample)
write.csv(cells_tab, TAB("cells_per_sample_per_subtype.csv"), row.names = FALSE)
## how many samples per species would survive MIN_CELLS_PER_PSEUDOBULK
surv_tab <- cells_tab %>% filter(n_cells >= MIN_CELLS_PER_PSEUDOBULK) %>%
  count(subtype, species, name = "n_samples_ok") %>%
  pivot_wider(names_from = species, values_from = n_samples_ok, values_fill = 0)
message(">>> samples with >= ", MIN_CELLS_PER_PSEUDOBULK, " cells, per species:"); print(surv_tab)
write.csv(surv_tab, TAB("samples_passing_min_cells_per_subtype.csv"), row.names = FALSE)

## ---- 1:1 orthologs -----------------------------------------------------------
all_genes <- rownames(immune)
if (!is.null(ORTHOLOG_TABLE)) {
  ortho <- read.csv(ORTHOLOG_TABLE, stringsAsFactors = FALSE)
  ortho <- ortho %>%
    filter(one2one %in% c(TRUE, "TRUE", "true", 1),
           !is.na(human_symbol), human_symbol != "",
           !is.na(marmoset_symbol), marmoset_symbol != "",
           !is.na(mouse_symbol), mouse_symbol != "")
  ortho_genes <- intersect(unique(ortho$human_symbol), all_genes)
  message(">>> 1:1 orthologs in immuneect: ", length(ortho_genes), " / ", nrow(ortho), " in table")
} else {
  ortho_genes <- all_genes
  message(">>> Ortholog universe = all ", length(ortho_genes),
          " genes of the immuneect (immuneect integrated on 1:1 orthologs)")
}
writeLines(ortho_genes, TAB("ortholog_gene_universe.txt"))

## ---- Hallmark gene sets --------------------------------------------------------
hm <- tryCatch(msigdbr(species = HALLMARK_SPECIES, collection = "H"),   # msigdbr >= 10
               error = function(e) msigdbr(species = HALLMARK_SPECIES, category = "H"))
hallmark_sets <- split(hm$gene_symbol, hm$gs_name)
hallmark_sets <- lapply(hallmark_sets, function(g) intersect(unique(g), ortho_genes))
message(">>> Hallmark sets: ", length(hallmark_sets), "; median genes/set after ortholog restriction: ",
        median(lengths(hallmark_sets)))

## ---- raw counts ------------------------------------------------------------------
counts_all <- GetAssayData(immune, assay = "RNA", layer = "counts")

## =============================================================================
## 2. HELPER FUNCTIONS  -------------------------------------------------------
## =============================================================================

## Pseudobulks (sum of raw counts) within a subtype:
##   species level  : one column per species, ALL its cells           (primary unit)
##   sample level   : one column per sample with >= MIN_CELLS_PER_PSEUDOBULK cells (stats only)
## CPM is computed on ALL genes (library size of the full pseudobulk), then the
## matrices are restricted to the ortholog universe.
sum_cells <- function(cell_ids, group) {              # sparse aggregation helper
  grp <- factor(group, levels = unique(group))
  design <- sparse.model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  as.matrix(counts_all[, cell_ids, drop = FALSE] %*% design)
}
to_cpm <- function(m) sweep(m, 2, colSums(m), "/") * 1e6

make_pseudobulk <- function(cells_md) {
  ## -- species level
  n_sp <- table(factor(cells_md$species, levels = SPECIES_LEVELS))
  if (any(n_sp < MIN_CELLS_PER_SPECIES)) {
    return(list(ok = FALSE, n_species_cells = n_sp,
                reason = paste0("< ", MIN_CELLS_PER_SPECIES, " cells in ",
                                paste(names(n_sp)[n_sp < MIN_CELLS_PER_SPECIES], collapse = ", "))))
  }
  sp_pb  <- sum_cells(cells_md$cell, as.character(cells_md$species))[, SPECIES_LEVELS, drop = FALSE]
  sp_cpm <- to_cpm(sp_pb)
  ## -- sample level (statistics only)
  cells_md <- cells_md %>% mutate(pb_id = paste(species, sample, sep = "||"))
  keep_ids <- cells_md %>% count(pb_id) %>% filter(n >= MIN_CELLS_PER_PSEUDOBULK) %>% pull(pb_id)
  smp <- cells_md %>% filter(pb_id %in% keep_ids)
  if (nrow(smp) > 0) {
    s_pb  <- sum_cells(smp$cell, smp$pb_id)
    s_cpm <- to_cpm(s_pb)
    info <- tibble(pb_id = colnames(s_pb)) %>%
      separate(pb_id, into = c("species", "sample"), sep = "\\|\\|", remove = FALSE) %>%
      mutate(species = factor(species, levels = SPECIES_LEVELS),
             n_cells = as.integer(table(smp$pb_id)[pb_id]),
             lib_size = colSums(s_pb)[pb_id])
  } else {
    s_pb <- s_cpm <- matrix(numeric(0), nrow = nrow(sp_pb), ncol = 0, dimnames = list(rownames(sp_pb), NULL))
    info <- tibble(pb_id = character(), species = factor(character(), levels = SPECIES_LEVELS),
                   sample = character(), n_cells = integer(), lib_size = numeric())
  }
  list(ok = TRUE,
       n_species_cells = n_sp,
       species_counts = sp_pb[ortho_genes, , drop = FALSE],
       species_cpm    = sp_cpm[ortho_genes, , drop = FALSE],
       species_lib    = colSums(sp_pb),
       counts = s_pb[ortho_genes, , drop = FALSE],
       cpm    = s_cpm[ortho_genes, , drop = FALSE],
       info   = info)
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## GSVA wrapper compatible with old (<1.50) and new (>=1.50) GSVA APIs
run_gsva <- function(expr, gsets, method = GSVA_METHOD) {
  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    prm <- if (method == "gsva") {
      GSVA::gsvaParam(expr, gsets, kcdf = "Gaussian", minSize = GSVA_MIN_SET, maxSize = GSVA_MAX_SET)
    } else {
      GSVA::ssgseaParam(expr, gsets, minSize = GSVA_MIN_SET, maxSize = GSVA_MAX_SET, normalize = TRUE)
    }
    GSVA::gsva(prm, verbose = FALSE)
  } else {
    GSVA::gsva(expr, gsets, method = method, kcdf = "Gaussian",
               min.sz = GSVA_MIN_SET, max.sz = GSVA_MAX_SET, verbose = FALSE)
  }
}

## two-circle Venn drawn with plain ggplot (no extra dependency)
circle_df <- function(cx, cy, r, n = 200) {
  t <- seq(0, 2 * pi, length.out = n); data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
}
plot_venn2 <- function(set1, set2, name1, name2, col1, col2, jac) {
  only1 <- length(setdiff(set1, set2)); both <- length(intersect(set1, set2)); only2 <- length(setdiff(set2, set1))
  c1 <- circle_df(-0.55, 0, 1); c2 <- circle_df(0.55, 0, 1)
  ggplot() +
    geom_polygon(data = c1, aes(x, y), fill = col1, alpha = 0.55, colour = "grey30", linewidth = 0.3) +
    geom_polygon(data = c2, aes(x, y), fill = col2, alpha = 0.55, colour = "grey30", linewidth = 0.3) +
    annotate("text", x = -1.05, y = 0, label = format(only1, big.mark = ","), size = 2.8) +
    annotate("text", x = 0,     y = 0, label = format(both,  big.mark = ","), size = 2.8, fontface = "bold") +
    annotate("text", x = 1.05,  y = 0, label = format(only2, big.mark = ","), size = 2.8) +
    annotate("text", x = 0, y = -1.35, label = sprintf("Jaccard = %.2f", jac), size = 2.9) +
    coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.6, 1.3)) +
    labs(title = paste(name1, "vs", name2)) +
    theme_void(base_size = 9) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8.5))
}

## pathway labels on scatter plots: ggrepel if installed, otherwise plain geom_text
ggrepel_or_text <- function(df) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = df, aes(label = label), size = 2.2, max.overlaps = 12,
                             segment.size = 0.2, colour = "grey25")
  } else {
    geom_text(data = df, aes(label = label), size = 2.0, vjust = -0.7, colour = "grey25", check_overlap = TRUE)
  }
}

## heatmap of species z-scores for a set of pathways (rows in the given order)
plot_zheat <- function(zdf, title, subtitle, zlim) {
  ggplot(zdf, aes(x = species, y = pathway_label, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
                         limits = c(-zlim, zlim), oob = scales::squish,
                         name = "Pathway activity\n(z-score)") +
    scale_y_discrete(limits = rev) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 0),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5),
          plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 8),
          legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"))
}

## =============================================================================
## 3. MAIN LOOP OVER SUBTYPES  ------------------------------------------------
## =============================================================================

subtypes <- if (is.null(SUBTYPES_TO_RUN)) sort(unique(md$subtype)) else SUBTYPES_TO_RUN

`%||%` <- function(a, b) if (is.null(a)) b else a
jaccard_all   <- list()   # per subtype: Jaccard per pair
expected_all  <- list()   # per subtype: expected-pathway check
cpm_dist_all  <- list()   # per subtype: CPM dynamic range per species
sens_all      <- list()   # per subtype: Jaccard vs CPM cutoff
ngenes_all    <- list()   # per subtype: n expressed genes per species
pathway_all   <- list()   # per subtype: full pathway table
corr_all      <- list()   # per subtype: r to human
skipped       <- character()

for (st in subtypes) {
  message("\n==================== ", st, " ====================")
  st_md <- md %>% filter(subtype == st)
  pb <- make_pseudobulk(st_md)
  message("  cells per species: ", paste(names(pb$n_species_cells), pb$n_species_cells, sep = "=", collapse = ", "))
  if (!pb$ok) { skipped <- c(skipped, st); message("  skipped: ", pb$reason); next }
  
  n_per_sp <- table(factor(pb$info$species, levels = SPECIES_LEVELS))   # per-sample pseudobulks (stats)
  message("  per-sample pseudobulks (>= ", MIN_CELLS_PER_PSEUDOBULK, " cells): ",
          paste(names(n_per_sp), n_per_sp, sep = "=", collapse = ", "))
  tag <- safe_name(st)
  write.csv(bind_rows(tibble(pb_id = paste(SPECIES_LEVELS, "ALL", sep = "||"), species = SPECIES_LEVELS, sample = "ALL cells",
                             n_cells = as.integer(pb$n_species_cells[SPECIES_LEVELS]), lib_size = unname(pb$species_lib[SPECIES_LEVELS])),
                      pb$info %>% mutate(species = as.character(species))),
            TAB(paste0(tag, "_pseudobulk_samples.csv")), row.names = FALSE)
  stats_ok <- all(n_per_sp >= MIN_SAMPLES_FOR_STATS)
  if (!stats_ok) message("  NOTE: < ", MIN_SAMPLES_FOR_STATS, " per-sample pseudobulks in ",
                         paste(names(n_per_sp)[n_per_sp < MIN_SAMPLES_FOR_STATS], collapse = ", "),
                         " -> ANOVA / t-tests / bootstrap CI reported as NA")
  
  ## ---------------------------------------------------------------------------
  ## 3.1  GENE REPERTOIRE  (Panel A)  -- species-level pseudobulk
  ## ---------------------------------------------------------------------------
  frac_expr <- sapply(SPECIES_LEVELS, function(sp) {       # fraction of samples with CPM > thr (NA if none)
    ids <- pb$info$pb_id[pb$info$species == sp]
    if (!length(ids)) return(rep(NA_real_, length(ortho_genes)))
    rowMeans(pb$cpm[, ids, drop = FALSE] > CPM_THRESHOLD)
  })
  expressed <- lapply(SPECIES_LEVELS, function(sp) {
    if (REPERTOIRE_MODE == "species_pooled") ortho_genes[pb$species_cpm[, sp] > CPM_THRESHOLD]
    else ortho_genes[!is.na(frac_expr[, sp]) & frac_expr[, sp] >= FRAC_SAMPLES_EXPRESSED]
  }); names(expressed) <- SPECIES_LEVELS
  
  ## per-gene presence table (useful for species-specific gene lists later)
  pres <- tibble(gene = ortho_genes,
                 Human = ortho_genes %in% expressed$Human,
                 Marmoset = ortho_genes %in% expressed$Marmoset,
                 Mouse = ortho_genes %in% expressed$Mouse,
                 speciesCPM_Human = round(pb$species_cpm[, "Human"], 2),
                 speciesCPM_Marmoset = round(pb$species_cpm[, "Marmoset"], 2),
                 speciesCPM_Mouse = round(pb$species_cpm[, "Mouse"], 2),
                 fracSamples_Human = round(frac_expr[, "Human"], 2),
                 fracSamples_Marmoset = round(frac_expr[, "Marmoset"], 2),
                 fracSamples_Mouse = round(frac_expr[, "Mouse"], 2)) %>%
    mutate(pattern = case_when(Human & Marmoset & Mouse ~ "All three",
                               Human & Mouse ~ "Human+Mouse", Human & Marmoset ~ "Human+Marmoset",
                               Marmoset & Mouse ~ "Marmoset+Mouse",
                               Human ~ "Human only", Marmoset ~ "Marmoset only", Mouse ~ "Mouse only",
                               TRUE ~ "None"))
  write.csv(pres, TAB(paste0(tag, "_expressed_gene_presence.csv")), row.names = FALSE)
  
  pairs <- list(c("Human", "Mouse"), c("Human", "Marmoset"), c("Marmoset", "Mouse"))
  
  ## -- CPM distribution per species (dynamic range) + cutoff sensitivity ----------
  cpm_dist <- bind_rows(lapply(SPECIES_LEVELS, function(sp) {
    v <- pb$species_cpm[, sp]; q <- quantile(v[v > 0], c(.5, .9, .99))
    tibble(subtype = st, species = sp, n_cells = as.integer(pb$n_species_cells[sp]),
           library_counts = unname(pb$species_lib[sp]),
           cpm_of_one_count = 1e6 / unname(pb$species_lib[sp]),
           genes_detected = sum(v > 0), median_cpm_detected = q[1], p90_cpm = q[2], p99_cpm = q[3],
           max_cpm = max(v), top5_genes = paste(names(sort(v, decreasing = TRUE))[1:5], collapse = ";"),
           !!!setNames(lapply(CPM_SENSITIVITY, function(ct) sum(v > ct)), paste0("n_genes_cpm_gt_", CPM_SENSITIVITY)))
  }))
  write.csv(cpm_dist, TAB(paste0(tag, "_cpm_distribution.csv")), row.names = FALSE)
  cpm_dist_all[[st]] <- cpm_dist
  sens <- bind_rows(lapply(CPM_SENSITIVITY, function(ct) {
    sets <- lapply(SPECIES_LEVELS, function(sp) ortho_genes[pb$species_cpm[, sp] > ct]); names(sets) <- SPECIES_LEVELS
    bind_rows(lapply(pairs, function(p) tibble(
      subtype = st, cpm_cutoff = ct, pair = paste(p[1], "vs", p[2]),
      n1 = length(sets[[p[1]]]), n2 = length(sets[[p[2]]]),
      shared = length(intersect(sets[[p[1]]], sets[[p[2]]])), jaccard = jaccard(sets[[p[1]]], sets[[p[2]]]))))
  }))
  sens_all[[st]] <- sens
  
  jac_df <- bind_rows(lapply(pairs, function(p) tibble(
    subtype = st, species1 = p[1], species2 = p[2],
    n1 = length(expressed[[p[1]]]), n2 = length(expressed[[p[2]]]),
    shared = length(intersect(expressed[[p[1]]], expressed[[p[2]]])),
    union  = length(union(expressed[[p[1]]], expressed[[p[2]]])),
    jaccard = jaccard(expressed[[p[1]]], expressed[[p[2]]]))))
  jac_df$pair <- paste(jac_df$species1, "vs", jac_df$species2)
  jaccard_all[[st]] <- jac_df
  ng <- tibble(subtype = st, species = factor(SPECIES_LEVELS, levels = SPECIES_LEVELS),
               n_expressed = lengths(expressed)[SPECIES_LEVELS])
  ngenes_all[[st]] <- ng
  
  pA_bar <- ggplot(ng, aes(species, n_expressed, fill = species)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = format(n_expressed, big.mark = ",")), vjust = -0.4, size = 2.9) +
    scale_fill_manual(values = SPECIES_COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
    labs(title = "Number of expressed\northologous genes", x = NULL, y = "Number of genes") +
    theme_pub()
  venns <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    plot_venn2(expressed[[p[1]]], expressed[[p[2]]], p[1], p[2],
               SPECIES_COLS[p[1]], SPECIES_COLS[p[2]], jac_df$jaccard[i])
  })
  pA_venn <- wrap_plots(venns, nrow = 1) +
    plot_annotation(title = "Overlap of expressed genes", subtitle = "(Venn diagrams)",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5),
                                  plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 8)))
  pA <- (pA_bar / wrap_elements(full = pA_venn)) + plot_layout(heights = c(1.1, 1))
  
  ## ---------------------------------------------------------------------------
  ## 3.2  PATHWAY ACTIVITY  (Panel B)
  ##   Species-level pseudobulks ("<species>||ALL") are scored TOGETHER with the
  ##   per-sample pseudobulks in one ssGSEA call so all scores share one scale.
  ##   Species score (means / z / r)  = score of the species pseudobulk.
  ##   Per-sample scores              = used for ANOVA, t-tests and bootstrap CI only.
  ## ---------------------------------------------------------------------------
  sp_ids <- paste(SPECIES_LEVELS, "ALL", sep = "||")
  cpm_all <- cbind(`colnames<-`(pb$species_cpm, sp_ids), pb$cpm)
  ## gene list for scoring = Panel A expressed genes (CPM > CPM_THRESHOLD, species level)
  if (PATHWAY_GENE_SET == "per_species") {
    ## one ssGSEA run per species on its own expressed genes (species profile + its samples)
    per_sp <- lapply(SPECIES_LEVELS, function(sp) {
      cols <- c(paste(sp, "ALL", sep = "||"), pb$info$pb_id[pb$info$species == sp])
      lc <- log2(cpm_all[expressed[[sp]], cols, drop = FALSE] + 1)
      if (ncol(lc) > 1) lc <- lc[apply(lc, 1, sd) > 0, , drop = FALSE]
      writeLines(rownames(lc), TAB(paste0(tag, "_pathway_gene_list_", sp, ".txt")))
      message("  ", toupper(GSVA_METHOD), " ", sp, ": ", nrow(lc), " expressed genes (CPM > ", CPM_THRESHOLD, ") x ", ncol(lc), " profiles")
      as.matrix(run_gsva(lc, hallmark_sets))
    })
    common_pw <- Reduce(intersect, lapply(per_sp, rownames))
    dropped <- setdiff(unique(unlist(lapply(per_sp, rownames))), common_pw)
    if (length(dropped)) message("  pathways not scorable in every species (dropped): ", paste(dropped, collapse = ", "))
    gs_all <- do.call(cbind, lapply(per_sp, function(m) m[common_pw, , drop = FALSE]))
    logcpm <- cpm_all[Reduce(union, expressed), , drop = FALSE]      # only used for n_genes_scored below
  } else {
    pathway_genes <- if (PATHWAY_GENE_SET == "shared") Reduce(intersect, expressed) else Reduce(union, expressed)
    logcpm <- log2(cpm_all[pathway_genes, , drop = FALSE] + 1)
    logcpm <- logcpm[apply(logcpm, 1, sd) > 0, , drop = FALSE]
    writeLines(rownames(logcpm), TAB(paste0(tag, "_pathway_gene_list.txt")))
    message("  ", toupper(GSVA_METHOD), " input: ", nrow(logcpm), " genes (", PATHWAY_GENE_SET,
            " of species expressed sets, CPM > ", CPM_THRESHOLD, "; same list for all species) x ", ncol(logcpm), " profiles")
    gs_all <- as.matrix(run_gsva(logcpm, hallmark_sets))   # pathways x (3 species + samples)
  }
  write.csv(as.data.frame(gs_all) %>% rownames_to_column("pathway"),
            TAB(paste0(tag, "_", GSVA_METHOD, "_scores_per_pseudobulk.csv")), row.names = FALSE)
  gs_species <- gs_all[, sp_ids, drop = FALSE]; colnames(gs_species) <- SPECIES_LEVELS
  gs_mat     <- gs_all[, setdiff(colnames(gs_all), sp_ids), drop = FALSE]   # per-sample
  
  long <- if (ncol(gs_mat) > 0) {
    as.data.frame(gs_mat) %>% rownames_to_column("pathway") %>%
      pivot_longer(-pathway, names_to = "pb_id", values_to = "score")
  } else tibble(pathway = character(), pb_id = character(), score = numeric())
  long <- long %>% left_join(pb$info %>% mutate(pb_id = as.character(pb_id)), by = "pb_id")
  
  ## species scores (primary) + per-sample summaries (for reference)
  means <- as.data.frame(gs_species) %>% rownames_to_column("pathway") %>%
    pivot_longer(-pathway, names_to = "species", values_to = "mean_score") %>%
    mutate(species = factor(species, levels = SPECIES_LEVELS)) %>%
    left_join(long %>% group_by(pathway, species) %>%
                summarise(sample_mean = mean(score), sample_sd = sd(score), n_samples = n(), .groups = "drop"),
              by = c("pathway", "species")) %>%
    mutate(n_samples = ifelse(is.na(n_samples), 0L, n_samples))
  mean_wide <- means %>% select(pathway, species, mean_score) %>%
    pivot_wider(names_from = species, values_from = mean_score)
  
  ## z-scores
  if (Z_MODE == "within_species") {
    z_long <- means %>% group_by(species) %>%
      mutate(z = as.numeric(scale(mean_score))) %>% ungroup() %>% select(pathway, species, z)
  } else if (Z_MODE == "species_means") {
    z_long <- means %>% group_by(pathway) %>%
      mutate(z = as.numeric(scale(mean_score))) %>% ungroup() %>% select(pathway, species, z)
  } else {
    z_long <- as.data.frame(gs_all) %>% rownames_to_column("pathway") %>%
      pivot_longer(-pathway, names_to = "pb_id", values_to = "score") %>%
      group_by(pathway) %>% mutate(zs = as.numeric(scale(score))) %>% ungroup() %>%
      filter(pb_id %in% sp_ids) %>%
      transmute(pathway, species = factor(sub("\\|\\|ALL$", "", pb_id), levels = SPECIES_LEVELS), z = zs)
  }
  z_wide <- z_long %>% pivot_wider(names_from = species, values_from = z, names_prefix = "z_")
  
  ## species effect: one-way ANOVA on per-sample scores (+ Kruskal as robustness)
  if (stats_ok) {
    species_effect <- long %>% group_by(pathway) %>%
      group_modify(~{
        a <- summary(aov(score ~ species, data = .x))[[1]]
        k <- tryCatch(kruskal.test(score ~ species, data = .x)$p.value, error = function(e) NA_real_)
        ss <- a[["Sum Sq"]]
        tibble(F_species = a[["F value"]][1], p_anova = a[["Pr(>F)"]][1],
               eta_sq = ss[1] / sum(ss), p_kruskal = k)
      }) %>% ungroup() %>%
      mutate(padj_anova = p.adjust(p_anova, "BH"))
  } else {
    species_effect <- tibble(pathway = mean_wide$pathway, F_species = NA_real_, p_anova = NA_real_,
                             eta_sq = NA_real_, p_kruskal = NA_real_, padj_anova = NA_real_)
  }
  
  ## pairwise comparisons: difference of SPECIES scores (always), Welch t-test on the
  ## per-sample scores when both sides have >= 2 samples (else NA), BH within each pair
  pw <- bind_rows(lapply(pairs, function(pr) {
    pair_name <- paste(pr[1], "vs", pr[2])      # computed OUTSIDE mutate (the p column would shadow it)
    x_ids <- pb$info$pb_id[pb$info$species == pr[1]]
    y_ids <- pb$info$pb_id[pb$info$species == pr[2]]
    pvals <- if (length(x_ids) >= 2 && length(y_ids) >= 2) {
      sapply(rownames(gs_species), function(pw_) tryCatch(
        t.test(gs_mat[pw_, x_ids], gs_mat[pw_, y_ids])$p.value, error = function(e) NA_real_))
    } else rep(NA_real_, nrow(gs_species))
    tibble(pathway = rownames(gs_species),
           diff = unname(gs_species[, pr[1]] - gs_species[, pr[2]]),
           p = unname(pvals), pair = pair_name) %>%
      mutate(padj = p.adjust(p, "BH"))
  }))
  pw_wide <- pw %>% mutate(pair = gsub(" ", "", pair)) %>%
    pivot_wider(id_cols = pathway, names_from = pair, values_from = c(diff, p, padj))
  
  path_tab <- mean_wide %>%
    left_join(z_wide, by = "pathway") %>%
    left_join(species_effect, by = "pathway") %>%
    left_join(pw_wide, by = "pathway") %>%
    mutate(subtype = st,
           pathway_label = clean_pathway(pathway),
           range_means = pmax(Human, Marmoset, Mouse) - pmin(Human, Marmoset, Mouse),
           grand_mean  = (Human + Marmoset + Mouse) / 3,
           z_range     = pmax(z_Human, z_Marmoset, z_Mouse) - pmin(z_Human, z_Marmoset, z_Mouse),
           z_mean      = (z_Human + z_Marmoset + z_Mouse) / 3,
           active_all  = z_Human > 0 & z_Marmoset > 0 & z_Mouse > 0,
           divergent_supported = !is.na(padj_anova) & padj_anova < DIVERGENT_PADJ,
           ## divergence key: effect size on the displayed scale (z range), or F
           active_any  = pmax(z_Human, z_Marmoset, z_Mouse) > 0,
           div_key = if (RANK_BY == "F" && stats_ok) F_species else z_range,
           rank_divergent = rank(-(div_key - if (DIVERGENT_REQUIRE_ACTIVE_ANY) ifelse(active_any, 0, 1e6) else 0),
                                 ties.method = "first"),
           ## conserved key: smallest divergence; active-in-all-species pathways first
           cons_key = div_key + if (CONSERVED_REQUIRE_ACTIVE) ifelse(active_all, 0, 1e6) else 0,
           rank_conserved = rank(cons_key, ties.method = "first")) %>%
    select(-div_key, -cons_key) %>%
    arrange(rank_conserved) %>%
    relocate(subtype, pathway, pathway_label)
  ## activity rank within each species (1 = most active of the 50 pathways)
  path_tab <- path_tab %>%
    mutate(activity_rank_Human    = rank(-z_Human,    ties.method = "first"),
           activity_rank_Marmoset = rank(-z_Marmoset, ties.method = "first"),
           activity_rank_Mouse    = rank(-z_Mouse,    ties.method = "first"),
           expected = pathway %in% (EXPECTED_PATHWAYS[[st]] %||% character()))
  if (PATHWAY_GENE_SET == "per_species") {
    for (sp in SPECIES_LEVELS) path_tab[[paste0("n_genes_scored_", sp)]] <- length(expressed[[sp]])
  } else path_tab$n_genes_scored <- nrow(logcpm)
  write.csv(path_tab, TAB(paste0(tag, "_hallmark_species_comparison.csv")), row.names = FALSE)
  pathway_all[[st]] <- path_tab
  
  ## ---------------------------------------------------------------------------
  ## 3.2b  EXPECTED PROGRAMS CHECK  (Panel D)  -- where the PI's expected pathways sit
  ## ---------------------------------------------------------------------------
  exp_set <- EXPECTED_PATHWAYS[[st]]
  if (!is.null(exp_set)) {
    missing_exp <- setdiff(exp_set, path_tab$pathway)
    if (length(missing_exp)) message("  expected pathways not in Hallmark/results: ", paste(missing_exp, collapse = ", "))
    chk <- path_tab %>% filter(expected) %>%
      arrange(desc(z_mean)) %>%
      transmute(subtype, pathway, pathway_label,
                z_Human, z_Marmoset, z_Mouse, z_mean, z_range,
                activity_rank_Human, activity_rank_Marmoset, activity_rank_Mouse,
                active_all, rank_conserved, rank_divergent, padj_anova,
                in_top_conserved = rank_conserved <= N_TOP, in_top_divergent = rank_divergent <= N_TOP,
                verdict = case_when(active_all & z_range <= 0.5 ~ "high & conserved",
                                    active_all ~ "high in all, species differ in degree",
                                    z_mean > 0 ~ "high in some species only",
                                    TRUE ~ "not among the active programs"))
    write.csv(chk, TAB(paste0(tag, "_expected_pathways_check.csv")), row.names = FALSE)
    expected_all[[st]] <- chk
    
    chk_long <- chk %>%
      select(pathway_label, z_Human, z_Marmoset, z_Mouse,
             activity_rank_Human, activity_rank_Marmoset, activity_rank_Mouse) %>%
      pivot_longer(-pathway_label, names_to = c(".value", "species"), names_pattern = "(z|activity_rank)_(.*)") %>%
      mutate(species = factor(species, levels = SPECIES_LEVELS),
             pathway_label = factor(pathway_label, levels = rev(chk$pathway_label)))
    pD <- ggplot(chk_long, aes(species, pathway_label, fill = z)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = paste0("#", activity_rank)), size = 2.7) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
                           limits = c(-2, 2), oob = scales::squish, name = "Pathway activity\n(z-score)") +
      labs(title = paste0(st, ": expected programs"),
           subtitle = "cell label = activity rank within species (1 = most active of 50 Hallmark pathways)",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
            plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 7.5))
    save_fig(pD, paste0(tag, "_D_expected_pathways"), 4.6, 0.32 * nrow(chk) + 1.8)
  }
  
  top_cons <- path_tab %>% arrange(rank_conserved) %>% slice_head(n = N_TOP)
  top_div  <- path_tab %>% arrange(rank_divergent) %>% slice_head(n = N_TOP)
  n_sig    <- sum(top_div$divergent_supported, na.rm = TRUE)
  
  to_zlong <- function(tab) tab %>% select(pathway, pathway_label, z_Human, z_Marmoset, z_Mouse) %>%
    pivot_longer(starts_with("z_"), names_to = "species", values_to = "z", names_prefix = "z_") %>%
    mutate(species = factor(species, levels = SPECIES_LEVELS),
           pathway_label = factor(pathway_label, levels = tab$pathway_label))
  zlim <- if (Z_MODE == "species_means") 1.2 else 2
  z_sub <- c(within_species = "(z across pathways, per species)",
             species_means  = "(z across species, per pathway)",
             samples        = "(z across pseudobulk samples, per pathway)")[[Z_MODE]]
  pB_cons <- plot_zheat(to_zlong(top_cons), sprintf("Top %d conserved pathways", N_TOP), "(similar across species)", zlim)
  pB_div  <- plot_zheat(to_zlong(top_div),  sprintf("Top %d divergent pathways", N_TOP),
                        if (stats_ok) sprintf("(species-enriched; %d/%d BH p < %.2f)", n_sig, N_TOP, DIVERGENT_PADJ)
                        else "(species-enriched; ranked by range, no test possible)", zlim)
  pB <- (pB_cons + theme(legend.position = "none")) + pB_div
  
  ## supplementary: all pathways ordered conserved -> divergent, with significance
  full_z <- to_zlong(path_tab %>% arrange(rank_conserved)) %>%
    left_join(path_tab %>% transmute(pathway, sig = case_when(padj_anova < 0.001 ~ "***", padj_anova < 0.01 ~ "**",
                                                              padj_anova < 0.05 ~ "*", TRUE ~ "")), by = "pathway")
  p_full <- plot_zheat(full_z, paste0(st, ": Hallmark pathway activity"),
                       "rows ordered most conserved (top) -> most divergent (bottom); * BH-adjusted ANOVA p", zlim) +
    geom_text(data = full_z %>% filter(species == "Mouse"), aes(label = sig), nudge_x = 0.75, size = 2.5) +
    coord_cartesian(clip = "off") + theme(plot.margin = margin(5, 20, 5, 5))
  save_fig(p_full, paste0(tag, "_B_supp_all_hallmark_zscores"), 5.5, 10)
  
  ## supplementary: mean-score heatmap (not z) for the same ordering
  p_mean <- path_tab %>% arrange(rank_conserved) %>%
    select(pathway_label, Human, Marmoset, Mouse) %>%
    pivot_longer(-pathway_label, names_to = "species", values_to = "mean_score") %>%
    mutate(species = factor(species, levels = SPECIES_LEVELS),
           pathway_label = factor(pathway_label, levels = unique(pathway_label))) %>%
    ggplot(aes(species, pathway_label, fill = mean_score)) + geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
                         name = paste0("Mean ", toupper(GSVA_METHOD), "\nscore")) +
    scale_y_discrete(limits = rev) + labs(x = NULL, y = NULL, title = paste0(st, ": mean pathway scores")) +
    theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
                                         plot.title = element_text(face = "bold", hjust = 0.5))
  save_fig(p_mean, paste0(tag, "_B_supp_all_hallmark_meanscores"), 5, 10)
  
  ## ---------------------------------------------------------------------------
  ## 3.3  CORRELATION OF PATHWAY PROGRAMS WITH HUMAN  (Panel C)
  ## ---------------------------------------------------------------------------
  mm <- as.matrix(mean_wide[, SPECIES_LEVELS]); rownames(mm) <- mean_wide$pathway
  r_obs <- c("Marmoset" = cor(mm[, "Human"], mm[, "Marmoset"]),
             "Mouse"    = cor(mm[, "Human"], mm[, "Mouse"]),
             "Marmoset_vs_Mouse" = cor(mm[, "Marmoset"], mm[, "Mouse"]))
  r_p <- c(cor.test(mm[, "Human"], mm[, "Marmoset"])$p.value,
           cor.test(mm[, "Human"], mm[, "Mouse"])$p.value,
           cor.test(mm[, "Marmoset"], mm[, "Mouse"])$p.value)
  
  ## bootstrap CI (replicate level): resample per-sample pseudobulks within each species,
  ## recompute species means and r. Only possible when every species has >= 2 samples.
  if (stats_ok) {
    samp_ids <- split(as.character(pb$info$pb_id), pb$info$species)
    boot <- replicate(BOOT_N, {
      bm <- sapply(SPECIES_LEVELS, function(sp) {
        ids <- sample(samp_ids[[sp]], length(samp_ids[[sp]]), replace = TRUE)
        rowMeans(gs_mat[, ids, drop = FALSE])
      })
      c(cor(bm[, "Human"], bm[, "Marmoset"]), cor(bm[, "Human"], bm[, "Mouse"]), cor(bm[, "Marmoset"], bm[, "Mouse"]))
    })
    ci <- apply(boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  } else {
    ci <- matrix(NA_real_, nrow = 2, ncol = 3)
  }
  corr_df <- tibble(subtype = st,
                    comparison = c("Human vs Marmoset", "Human vs Mouse", "Marmoset vs Mouse"),
                    other = c("Marmoset", "Mouse", "Marmoset vs Mouse"),
                    r = unname(r_obs), p = r_p, ci_low = ci[1, ], ci_high = ci[2, ],
                    n_pathways = nrow(mm))
  corr_all[[st]] <- corr_df
  
  pC_df <- corr_df %>% filter(other %in% c("Marmoset", "Mouse")) %>%
    mutate(other = factor(other, levels = c("Marmoset", "Mouse")))
  pC <- ggplot(pC_df, aes(other, r, fill = other)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.4, na.rm = TRUE) +
    geom_text(aes(y = pmax(ci_high, r, na.rm = TRUE), label = sprintf("%.2f", r)), vjust = -0.6, size = 3) +
    scale_fill_manual(values = SPECIES_COLS, guide = "none") +
    scale_y_continuous(limits = c(min(0, min(c(pC_df$ci_low, pC_df$r), na.rm = TRUE) - 0.05), 1.08),
                       breaks = seq(-1, 1, 0.2), expand = c(0, 0)) +
    labs(title = "Pathway profile similarity",
         subtitle = sprintf("(correlation with human, %d pathways)", nrow(mm)),
         x = NULL, y = "Correlation (r)") +
    theme_pub()
  
  ## scatter plots of mean scores (Human vs Marmoset, Human vs Mouse)
  scat <- bind_rows(
    tibble(pathway = rownames(mm), x = mm[, "Human"], y = mm[, "Marmoset"], other = "Marmoset"),
    tibble(pathway = rownames(mm), x = mm[, "Human"], y = mm[, "Mouse"],    other = "Mouse")) %>%
    mutate(other = factor(other, levels = c("Marmoset", "Mouse")), label = clean_pathway(pathway))
  lab_df <- corr_df %>% filter(other %in% c("Marmoset", "Mouse")) %>%
    mutate(other = factor(other, levels = c("Marmoset", "Mouse")),
           txt = sprintf("r = %.2f\np = %.2g", r, p))
  p_scat <- ggplot(scat, aes(x, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5) +
    geom_point(aes(colour = other), size = 1.6, alpha = 0.85) +
    ggrepel_or_text(scat) +
    geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = txt), hjust = -0.15, vjust = 1.3, size = 3) +
    facet_wrap(~ other, nrow = 1, labeller = labeller(other = function(x) paste("Human vs", x))) +
    scale_colour_manual(values = SPECIES_COLS, guide = "none") +
    labs(x = paste("Human mean", toupper(GSVA_METHOD), "score"), y = paste("Mean", toupper(GSVA_METHOD), "score"),
         title = paste0(st, ": pathway program similarity to human")) +
    theme_pub() + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_fig(p_scat, paste0(tag, "_C_supp_scatter_mean_scores"), 7.5, 3.8)
  
  ## ---------------------------------------------------------------------------
  ## 3.4  COMPOSITE FIGURE  (A | B | C)
  ## ---------------------------------------------------------------------------
  pA_t <- wrap_elements(full = pA + plot_annotation(
    title = "A.  Gene repertoire identity",
    subtitle = if (REPERTOIRE_MODE == "species_pooled")
      sprintf("(orthologous genes with CPM > %g in the species-pooled profile)", CPM_THRESHOLD)
    else sprintf("(orthologous genes with CPM > %g in ≥%d%% samples)", CPM_THRESHOLD, round(FRAC_SAMPLES_EXPRESSED * 100)),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(face = "italic", size = 8.5))))
  pB_t <- wrap_elements(full = pB + plot_annotation(
    title = "B.  Pathway activity comparison",
    subtitle = sprintf("(%s using MSigDB Hallmark pathways)", toupper(GSVA_METHOD)),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(face = "italic", size = 8.5))))
  pC_t <- wrap_elements(full = pC + plot_annotation(
    title = "C.  Similarity of pathway programs to human",
    theme = theme(plot.title = element_text(face = "bold", size = 11))))
  ncell_txt <- paste0("cells: ", paste(SPECIES_LEVELS, pb$n_species_cells[SPECIES_LEVELS], sep = " ", collapse = ", "))
  st_title <- if (!stats_ok) {
    paste0(st, "   [", ncell_txt, "; < ", MIN_SAMPLES_FOR_STATS, " replicate samples in ",
           paste(names(n_per_sp)[n_per_sp < MIN_SAMPLES_FOR_STATS], collapse = ", "), " -> no test / CI]")
  } else paste0(st, "   [", ncell_txt, "]")
  composite <- (pA_t | pB_t | pC_t) + plot_layout(widths = c(1.15, 1.9, 0.8)) +
    plot_annotation(title = st_title, theme = theme(plot.title = element_text(face = "bold", size = 13)))
  save_fig(composite, paste0("FIG_", tag, "_similarity_differences"), 16, 6)
  message("  done: ", tag)
}

## =============================================================================
## 4. CROSS-SUBTYPE SUMMARIES  ------------------------------------------------
## =============================================================================

if (length(skipped)) writeLines(skipped, TAB("skipped_subtypes.txt"))

jac_tab  <- bind_rows(jaccard_all);  write.csv(jac_tab,  TAB("SUMMARY_jaccard_by_subtype.csv"), row.names = FALSE)
ng_tab   <- bind_rows(ngenes_all);   write.csv(ng_tab,   TAB("SUMMARY_n_expressed_genes_by_subtype.csv"), row.names = FALSE)
corr_tab <- bind_rows(corr_all);     write.csv(corr_tab, TAB("SUMMARY_pathway_correlation_by_subtype.csv"), row.names = FALSE)
path_tab_all <- bind_rows(pathway_all)
write.csv(path_tab_all, TAB("SUMMARY_hallmark_species_comparison_all_subtypes.csv"), row.names = FALSE)
write.csv(path_tab_all %>% group_by(subtype) %>% slice_min(rank_conserved, n = N_TOP) %>% ungroup(),
          TAB(paste0("SUMMARY_top", N_TOP, "_conserved_pathways_by_subtype.csv")), row.names = FALSE)
write.csv(path_tab_all %>% group_by(subtype) %>% slice_min(rank_divergent, n = N_TOP) %>% ungroup(),
          TAB(paste0("SUMMARY_top", N_TOP, "_divergent_pathways_by_subtype.csv")), row.names = FALSE)

subtype_order <- unique(ng_tab$subtype)

## CPM dynamic range + cutoff sensitivity across subtypes
cpm_dist_tab <- bind_rows(cpm_dist_all); write.csv(cpm_dist_tab, TAB("SUMMARY_cpm_distribution.csv"), row.names = FALSE)
sens_tab <- bind_rows(sens_all);         write.csv(sens_tab, TAB("SUMMARY_cpm_cutoff_sensitivity.csv"), row.names = FALSE)
pair_cols <- c("Human vs Marmoset" = "#377EB8", "Human vs Mouse" = "#E41A1C", "Marmoset vs Mouse" = "#984EA3")
p_sens <- ggplot(sens_tab %>% mutate(subtype = factor(subtype, levels = subtype_order)),
                 aes(cpm_cutoff, jaccard, colour = pair)) +
  geom_line() + geom_point(size = 1.6) +
  geom_vline(xintercept = CPM_THRESHOLD, linetype = 2, colour = "grey50") +
  scale_x_log10(breaks = CPM_SENSITIVITY) + scale_colour_manual(values = pair_cols, name = NULL) +
  facet_wrap(~ subtype, ncol = 4) +
  labs(x = "CPM cutoff for 'expressed' (log scale)", y = "Jaccard index",
       title = "Gene-repertoire similarity vs CPM cutoff", subtitle = "dashed line = cutoff used in this run") +
  theme_pub(9) + theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p_sens, "SUMMARY_cpm_cutoff_sensitivity", 11, 2.6 * ceiling(length(sens_all) / 4) + 1.4)
p_ngs <- ggplot(sens_tab %>% filter(pair == "Human vs Marmoset") %>% select(subtype, cpm_cutoff, Human = n1, Marmoset = n2) %>%
                  left_join(sens_tab %>% filter(pair == "Human vs Mouse") %>% select(subtype, cpm_cutoff, Mouse = n2), by = c("subtype", "cpm_cutoff")) %>%
                  pivot_longer(c(Human, Marmoset, Mouse), names_to = "species", values_to = "n_genes") %>%
                  mutate(species = factor(species, levels = SPECIES_LEVELS), subtype = factor(subtype, levels = subtype_order)),
                aes(cpm_cutoff, n_genes, colour = species)) +
  geom_line() + geom_point(size = 1.6) +
  geom_vline(xintercept = CPM_THRESHOLD, linetype = 2, colour = "grey50") +
  scale_x_log10(breaks = CPM_SENSITIVITY) + scale_y_continuous(labels = scales::comma) +
  scale_colour_manual(values = SPECIES_COLS, name = NULL) + facet_wrap(~ subtype, ncol = 4) +
  labs(x = "CPM cutoff for 'expressed' (log scale)", y = "Expressed orthologous genes",
       title = "Number of expressed genes vs CPM cutoff") +
  theme_pub(9) + theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p_ngs, "SUMMARY_cpm_cutoff_ngenes", 11, 2.6 * ceiling(length(sens_all) / 4) + 1.4)

## expected-programs summary across subtypes (CSV + one faceted heatmap)
if (length(expected_all)) {
  exp_tab <- bind_rows(expected_all)
  write.csv(exp_tab, TAB("SUMMARY_expected_pathways_check.csv"), row.names = FALSE)
  exp_long <- exp_tab %>%
    select(subtype, pathway_label, z_Human, z_Marmoset, z_Mouse,
           activity_rank_Human, activity_rank_Marmoset, activity_rank_Mouse) %>%
    pivot_longer(-c(subtype, pathway_label), names_to = c(".value", "species"), names_pattern = "(z|activity_rank)_(.*)") %>%
    mutate(species = factor(species, levels = SPECIES_LEVELS), subtype = factor(subtype, levels = subtype_order))
  p_exp <- ggplot(exp_long, aes(species, pathway_label, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = activity_rank), size = 2.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
                         limits = c(-2, 2), oob = scales::squish, name = "z-score") +
    scale_y_discrete(limits = rev) +
    facet_wrap(~ subtype, scales = "free_y", ncol = 3) +
    labs(title = "Expected Hallmark programs per subtype",
         subtitle = "colour = within-species z-score of the species pseudobulk; number = activity rank within species (1 = highest of 50)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 8.5) +
    theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 7.5))
  save_fig(p_exp, "SUMMARY_expected_pathways_heatmap", 11, 2.6 * ceiling(length(expected_all) / 3) + 1.5)
}
blue_scale <- scale_fill_gradient(low = "#F7FBFF", high = "#08306B", limits = c(0, 1), name = NULL)

## Jaccard heatmap (subtype x species pair) -- Fig-4-Panel-F style
p_jac <- jac_tab %>% mutate(subtype = factor(subtype, levels = rev(subtype_order)),
                            pair = factor(pair, levels = c("Human vs Marmoset", "Human vs Mouse", "Marmoset vs Mouse"))) %>%
  ggplot(aes(pair, subtype, fill = jaccard)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", jaccard), colour = jaccard > 0.6), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  blue_scale + labs(title = "Gene repertoire similarity (Jaccard)", x = NULL, y = NULL, fill = "Jaccard") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
                                        plot.title = element_text(face = "bold", hjust = 0.5))
save_fig(p_jac, "SUMMARY_jaccard_heatmap_by_subtype", 5, 0.45 * length(subtype_order) + 1.5)

## Pearson r heatmap (subtype x Human-vs-X) -- Fig-4-Panel-F style
p_cor <- corr_tab %>% filter(other %in% c("Marmoset", "Mouse")) %>%
  mutate(subtype = factor(subtype, levels = rev(subtype_order)), other = factor(other, levels = c("Marmoset", "Mouse"))) %>%
  ggplot(aes(other, subtype, fill = r)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r), colour = r > 0.6), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  blue_scale + labs(title = "Pathway-program similarity to human (Pearson r)", x = "Human vs.", y = NULL) +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), axis.text = element_text(colour = "black"),
                                        plot.title = element_text(face = "bold", hjust = 0.5))
save_fig(p_cor, "SUMMARY_pathway_correlation_heatmap_by_subtype", 4.2, 0.45 * length(subtype_order) + 1.5)

## grouped bar: n expressed genes per species per subtype
p_ng <- ggplot(ng_tab %>% mutate(subtype = factor(subtype, levels = subtype_order)),
               aes(subtype, n_expressed, fill = species)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = SPECIES_COLS, name = NULL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Expressed orthologous genes", title = "Expressed gene repertoire per species") +
  theme_pub(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_fig(p_ng, "SUMMARY_n_expressed_genes_by_subtype", 1 + 0.9 * length(subtype_order), 4)

## correlation-with-human bars for all subtypes (Panel C across subtypes)
p_corbar <- corr_tab %>% filter(other %in% c("Marmoset", "Mouse")) %>%
  mutate(subtype = factor(subtype, levels = subtype_order), other = factor(other, levels = c("Marmoset", "Mouse"))) %>%
  ggplot(aes(subtype, r, fill = other)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(0.7), width = 0.2, linewidth = 0.4, na.rm = TRUE) +
  scale_fill_manual(values = SPECIES_COLS, name = "Human vs.") +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2)) + geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Correlation with human (r)", title = "Pathway profile similarity to human (bootstrap 95% CI)") +
  theme_pub(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_fig(p_corbar, "SUMMARY_pathway_correlation_bars_by_subtype", 1 + 0.9 * length(subtype_order), 4)

## =============================================================================
## 5. INDEX FILE  (INDEX.md in OUT_DIR -- regenerated on every run)
## =============================================================================

run_subtypes <- names(jaccard_all)
idx <- c(
  "# Cross-species similarity / differences -- immune subtypes",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), "  |  immuneect: `", immune_PATH, "`  |  subtype column: `", SUBTYPE_COL,
         "`  |  sample column: `", SAMPLE_COL, "`"),
  paste0("Method: SPECIES-level pseudobulk (all cells of the species; subtype included when every species has >= ",
         MIN_CELLS_PER_SPECIES, " cells) -> CPM; expressed = CPM > ", CPM_THRESHOLD,
         if (REPERTOIRE_MODE == "species_pooled") " in the species pseudobulk; " else paste0(" in >= ", round(FRAC_SAMPLES_EXPRESSED * 100), "% of samples; "),
         toupper(GSVA_METHOD),
         " on log2(CPM+1) with ", length(hallmark_sets), " Hallmark sets; z-mode = ", Z_MODE,
         "; conserved/divergent ranked by ", if (RANK_BY == "F") "ANOVA F" else "range of species z-scores (effect size)",
         if (CONSERVED_REQUIRE_ACTIVE) ", conserved drawn from pathways active (z > 0) in all species" else "",
         "; species effect = one-way ANOVA (BH) as supporting flag; Pearson r of species mean scores with ", BOOT_N, "x bootstrap CI.",
         " Expressed-gene mode: ", REPERTOIRE_MODE, "."),
  paste0("Ortholog universe: ", length(ortho_genes), " genes."),
  "",
  paste0("## Subtypes analysed (", length(run_subtypes), ")"), "",
  paste0("- ", run_subtypes),
  "",
  paste0("## Subtypes skipped (", length(skipped), ")"), "",
  if (length(skipped)) paste0("- ", skipped, "  (see tables/samples_passing_min_cells_per_subtype.csv)") else "- none",
  "",
  "## Figures (figures/, each as .pdf + .png)", "",
  "### Per subtype", "",
  "| File | Content |", "|---|---|",
  "| `FIG_<subtype>_similarity_differences` | MAIN composite: A gene repertoire (bar + 3 Venns + Jaccard), B top-10 conserved / divergent Hallmark pathways (z-score heatmaps), C Pearson r with human (bootstrap 95% CI) |",
  "| `<subtype>_B_supp_all_hallmark_zscores` | all 50 Hallmark pathways, z-scores, ordered most conserved -> most divergent, stars = BH-adjusted ANOVA p |",
  "| `<subtype>_B_supp_all_hallmark_meanscores` | same ordering, raw mean scores instead of z |",
  "| `<subtype>_C_supp_scatter_mean_scores` | scatter of mean pathway scores, Human vs Marmoset and Human vs Mouse, with r and p |",
  "| `<subtype>_D_expected_pathways` | PI's expected programs for this subtype: z per species + activity rank within species |",
  "| `SUMMARY_expected_pathways_heatmap` | the same for all subtypes in one faceted heatmap (table: SUMMARY_expected_pathways_check.csv) |",
  "", "### Cross-subtype summaries", "",
  "| File | Content |", "|---|---|",
  "| `SUMMARY_jaccard_heatmap_by_subtype` | Jaccard index, subtype x species pair (blue scale) |",
  "| `SUMMARY_pathway_correlation_heatmap_by_subtype` | Pearson r with human, subtype x species (blue scale) |",
  "| `SUMMARY_pathway_correlation_bars_by_subtype` | same r as grouped bars with bootstrap 95% CI |",
  "| `SUMMARY_n_expressed_genes_by_subtype` | number of expressed orthologous genes per species per subtype |",
  "",
  "## Tables (tables/)", "",
  "| File | Content |", "|---|---|",
  "| `samples_per_species_per_subtype.csv` | samples per species per subtype BEFORE the min-cells filter |",
  "| `cells_per_sample_per_subtype.csv` | cells per sample per subtype (diagnostic for skipped subtypes) |",
  "| `samples_passing_min_cells_per_subtype.csv` | samples per species that pass MIN_CELLS_PER_PSEUDOBULK |",
  "| `skipped_subtypes.txt` | subtypes not analysed |",
  "| `ortholog_gene_universe.txt` | gene universe used |",
  "| `<subtype>_pseudobulk_samples.csv` | pseudobulk samples used (cells, library size) |",
  "| `<subtype>_expressed_gene_presence.csv` | per gene: expressed in Human / Marmoset / Mouse + presence pattern |",
  paste0("| `<subtype>_", GSVA_METHOD, "_scores_per_pseudobulk.csv` | pathway x pseudobulk-sample score matrix |"),
  "| `<subtype>_hallmark_species_comparison.csv` | per pathway: species means, z-scores, ANOVA F / p / BH p, eta2, pairwise diffs and p, ranks |",
  "| `SUMMARY_jaccard_by_subtype.csv` | Jaccard, set sizes, shared / union per pair per subtype |",
  "| `SUMMARY_n_expressed_genes_by_subtype.csv` | expressed gene counts |",
  "| `SUMMARY_pathway_correlation_by_subtype.csv` | r, p, bootstrap CI for Human-Marmoset, Human-Mouse, Marmoset-Mouse |",
  "| `SUMMARY_hallmark_species_comparison_all_subtypes.csv` | all pathway tables concatenated |",
  paste0("| `SUMMARY_top", N_TOP, "_conserved_pathways_by_subtype.csv` / `..._divergent_...` | top-", N_TOP, " lists per subtype |"),
  "",
  "## Key numbers", "",
  "| Subtype | Jaccard H-Mar | Jaccard H-M | Jaccard Mar-M | r Human-Marmoset | r Human-Mouse |", "|---|---|---|---|---|---|"
)
for (st in run_subtypes) {
  j <- jaccard_all[[st]]; cc <- corr_all[[st]]
  idx <- c(idx, sprintf("| %s | %.2f | %.2f | %.2f | %.2f | %.2f |", st,
                        j$jaccard[j$pair == "Human vs Marmoset"], j$jaccard[j$pair == "Human vs Mouse"],
                        j$jaccard[j$pair == "Marmoset vs Mouse"],
                        cc$r[cc$other == "Marmoset"], cc$r[cc$other == "Mouse"]))
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
writeLines(idx, file.path(OUT_DIR, "INDEX.md"))

## ---- index.html (relative links; open it from inside OUT_DIR) -----------------
esc <- function(x) { x <- gsub("&", "&amp;", x); x <- gsub("<", "&lt;", x); gsub(">", "&gt;", x) }
fig_links <- function(base, extra = "") sprintf(
  '<p class="links"><a href="figures/%s.pdf">PDF</a><a href="figures/%s.png">PNG</a>%s</p>', base, base, extra)
html <- c(
  '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>Immune similarity index</title><style>',
  ':root{--ink:#1f2328;--muted:#59636e;--line:#d9dee3;--panel:#f6f8fa}*{box-sizing:border-box}',
  'body{margin:0;background:#fff;color:var(--ink);font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}',
  'main{max-width:1100px;margin:0 auto;padding:28px 16px 60px}h1{font-size:26px;margin:0 0 4px}',
  'h2{font-size:19px;margin:36px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}h3{font-size:16px;margin:22px 0 8px}',
  '.meta,.small{color:var(--muted);font-size:13.5px}code{background:var(--panel);padding:1px 5px;border-radius:4px;font-size:13px}',
  'table{border-collapse:collapse;width:100%;margin:8px 0 14px;font-size:14px}th,td{border:1px solid var(--line);padding:6px 9px;text-align:left;vertical-align:top}',
  'th{background:var(--panel)}td.num{text-align:right;font-variant-numeric:tabular-nums}',
  '.chips span{display:inline-block;background:var(--panel);border:1px solid var(--line);border-radius:999px;padding:2px 10px;margin:3px 4px 3px 0;font-size:13px}',
  '.fig{border:1px solid var(--line);border-radius:8px;padding:12px;margin:12px 0 18px;background:#fff}.fig img{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:4px}',
  '.links a{margin-right:12px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px}.grid .fig{margin:0}.grid img{max-height:260px;immuneect-fit:contain}',
  'a{color:#0b57d0}nav a{margin-right:14px;font-size:14px}</style></head><body><main>',
  '<h1>Cross-species similarity / differences &mdash; immune subtypes</h1>',
  sprintf('<p class="meta">Generated %s &middot; immuneect <code>%s</code> &middot; subtype column <code>%s</code> &middot; sample column <code>%s</code> &middot; ortholog universe %d genes</p>',
          format(Sys.time(), "%Y-%m-%d %H:%M"), esc(immune_PATH), SUBTYPE_COL, SAMPLE_COL, length(ortho_genes)),
  '<nav><a href="#settings">Settings</a><a href="#subtypes">Subtypes</a><a href="#main">Main figures</a><a href="#summary">Cross-subtype</a><a href="#supp">Supplementary</a><a href="#cpm">CPM cutoff</a><a href="#expected">Expected programs</a><a href="#tables">Tables</a><a href="#numbers">Key numbers</a></nav>',
  '<h2 id="settings">Settings</h2><table><tr><th style="width:32%">Step</th><th>Value</th></tr>',
  '<tr><td>Pseudobulk</td><td>SPECIES level: all cells of the species summed (used for expressed genes, pathway scores, z and r); SAMPLE level: per biological sample, used only for ANOVA / t-tests / bootstrap CI</td></tr>',
  sprintf('<tr><td>Inclusion</td><td>subtype analysed when every species has &ge; %d cells; replicate-level tests / CIs when every species has &ge; %d per-sample pseudobulks (&ge; %d cells each)</td></tr>', MIN_CELLS_PER_SPECIES, MIN_SAMPLES_FOR_STATS, MIN_CELLS_PER_PSEUDOBULK),
  if (REPERTOIRE_MODE == "species_pooled") sprintf('<tr><td>Expressed gene</td><td>CPM &gt; %g in the species-pooled pseudobulk</td></tr>', CPM_THRESHOLD)
  else sprintf('<tr><td>Expressed gene</td><td>CPM &gt; %g in &ge; %d%% of the species&#39; samples</td></tr>', CPM_THRESHOLD, round(FRAC_SAMPLES_EXPRESSED * 100)),
  sprintf('<tr><td>Conserved / divergent ranking</td><td>%s%s; ANOVA BH p kept as supporting flag</td></tr>',
          if (RANK_BY == "F") "ANOVA F statistic" else "range of species z-scores (effect size on the displayed scale)",
          if (CONSERVED_REQUIRE_ACTIVE) "; conserved drawn from pathways active (z &gt; 0) in all three species" else ""),
  if (PATHWAY_GENE_SET == "per_species")
    sprintf('<tr><td>Pathway scoring</td><td>%s on log2(CPM+1), %d MSigDB Hallmark sets; each species scored on its OWN expressed-gene list (CPM &gt; %g; tables/&lt;subtype&gt;_pathway_gene_list_&lt;species&gt;.txt); species compared through within-species z-scores and Pearson r</td></tr>',
            toupper(GSVA_METHOD), length(hallmark_sets), CPM_THRESHOLD)
  else
    sprintf('<tr><td>Pathway scoring</td><td>%s on log2(CPM+1), %d MSigDB Hallmark sets; gene list = %s of the species expressed-gene sets (CPM &gt; %g), identical for all species (tables/&lt;subtype&gt;_pathway_gene_list.txt)</td></tr>',
            toupper(GSVA_METHOD), length(hallmark_sets), PATHWAY_GENE_SET, CPM_THRESHOLD),
  sprintf('<tr><td>z-score</td><td>%s</td></tr>', Z_MODE),
  '<tr><td>Species effect</td><td>one-way ANOVA on per-sample scores, BH-adjusted; pairwise Welch t-tests</td></tr>',
  sprintf('<tr><td>Similarity to human</td><td>Pearson r of species mean scores over all pathways; %d&times; bootstrap 95%% CI</td></tr></table>', BOOT_N),
  sprintf('<h2 id="subtypes">Subtypes</h2><h3>Analysed (%d)</h3><p class="chips">%s</p>', length(run_subtypes),
          paste0("<span>", esc(run_subtypes), "</span>", collapse = "")),
  sprintf('<h3>Skipped (%d)</h3><p class="chips">%s</p><p class="small">See <a href="tables/samples_passing_min_cells_per_subtype.csv">samples_passing_min_cells_per_subtype.csv</a> and <a href="tables/cells_per_sample_per_subtype.csv">cells_per_sample_per_subtype.csv</a>.</p>',
          length(skipped), if (length(skipped)) paste0("<span>", esc(skipped), "</span>", collapse = "") else "<span>none</span>"),
  '<h2 id="main">Main composite figures (A | B | C)</h2>',
  '<p class="small">A &mdash; gene repertoire (bar + 3 Venns + Jaccard) &middot; B &mdash; top conserved / divergent Hallmark pathways (z-score heatmaps) &middot; C &mdash; Pearson r with human (bootstrap 95% CI)</p>'
)
for (st in run_subtypes) {
  tag <- safe_name(st); base <- paste0("FIG_", tag, "_similarity_differences")
  html <- c(html, sprintf('<div class="fig"><h3>%s</h3><a href="figures/%s.png"><img src="figures/%s.png" alt="%s"></a>%s</div>',
                          esc(st), base, base, esc(st),
                          fig_links(base, sprintf('<a href="tables/%s_hallmark_species_comparison.csv">pathway table</a><a href="tables/%s_pseudobulk_samples.csv">pseudobulk samples</a><a href="tables/%s_expressed_gene_presence.csv">gene presence</a><a href="tables/%s_%s_scores_per_pseudobulk.csv">%s scores</a>',
                                                  tag, tag, tag, tag, GSVA_METHOD, toupper(GSVA_METHOD)))))
}
summ <- list(c("SUMMARY_jaccard_heatmap_by_subtype", "Jaccard by subtype", "SUMMARY_jaccard_by_subtype.csv"),
             c("SUMMARY_pathway_correlation_heatmap_by_subtype", "Pearson r with human (heatmap)", "SUMMARY_pathway_correlation_by_subtype.csv"),
             c("SUMMARY_pathway_correlation_bars_by_subtype", "Pearson r with human (bars + CI)", ""),
             c("SUMMARY_n_expressed_genes_by_subtype", "Expressed orthologs per species", "SUMMARY_n_expressed_genes_by_subtype.csv"))
html <- c(html, '<h2 id="summary">Cross-subtype summaries</h2><div class="grid">',
          sapply(summ, function(s) sprintf('<div class="fig"><h3>%s</h3><a href="figures/%s.png"><img src="figures/%s.png" alt="%s"></a>%s</div>',
                                           s[2], s[1], s[1], s[2],
                                           fig_links(s[1], if (nzchar(s[3])) sprintf('<a href="tables/%s">CSV</a>', s[3]) else ""))),
          '</div>',
          '<h2 id="supp">Supplementary figures per subtype</h2><table><tr><th>Subtype</th><th>All pathways, z-scores</th><th>All pathways, mean scores</th><th>Scatter of mean scores vs human</th></tr>')
for (st in run_subtypes) {
  tag <- safe_name(st)
  cell <- function(b) sprintf('<a href="figures/%s.pdf">PDF</a> &middot; <a href="figures/%s.png">PNG</a>', b, b)
  html <- c(html, sprintf('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>', esc(st),
                          cell(paste0(tag, "_B_supp_all_hallmark_zscores")), cell(paste0(tag, "_B_supp_all_hallmark_meanscores")),
                          cell(paste0(tag, "_C_supp_scatter_mean_scores"))))
}
tab_files <- list.files(TAB(), pattern = "\\.(csv|txt)$")
html <- c(html, '</table>')
html <- c(html, '<h2 id="cpm">CPM cutoff sensitivity</h2>',
          sprintf('<p class="small">Jaccard and expressed-gene counts recomputed at CPM &gt; %s; dashed line = cutoff used (%g). Dynamic range per species in <a href="tables/SUMMARY_cpm_distribution.csv">SUMMARY_cpm_distribution.csv</a> (median / p90 / p99 / max CPM, CPM of one count, top genes).</p>',
                  paste(CPM_SENSITIVITY, collapse = ", "), CPM_THRESHOLD),
          '<div class="grid">',
          sprintf('<div class="fig"><h3>Jaccard vs cutoff</h3><a href="figures/SUMMARY_cpm_cutoff_sensitivity.png"><img src="figures/SUMMARY_cpm_cutoff_sensitivity.png" alt="Jaccard vs CPM cutoff"></a>%s</div>',
                  fig_links("SUMMARY_cpm_cutoff_sensitivity", '<a href="tables/SUMMARY_cpm_cutoff_sensitivity.csv">CSV</a>')),
          sprintf('<div class="fig"><h3>Expressed genes vs cutoff</h3><a href="figures/SUMMARY_cpm_cutoff_ngenes.png"><img src="figures/SUMMARY_cpm_cutoff_ngenes.png" alt="Genes vs CPM cutoff"></a>%s</div>',
                  fig_links("SUMMARY_cpm_cutoff_ngenes")),
          '</div>')
if (length(expected_all)) {
  html <- c(html, '<h2 id="expected">Expected programs check</h2>',
            '<p class="small">Colour = within-species z-score of the species pseudobulk; number = activity rank within species (1 = most active of 50 Hallmark pathways). Per-subtype versions: <code>&lt;subtype&gt;_D_expected_pathways</code>.</p>',
            '<div class="fig"><a href="figures/SUMMARY_expected_pathways_heatmap.png"><img src="figures/SUMMARY_expected_pathways_heatmap.png" alt="Expected programs"></a>',
            fig_links("SUMMARY_expected_pathways_heatmap", '<a href="tables/SUMMARY_expected_pathways_check.csv">CSV</a>'), '</div>')
}
html <- c(html, sprintf('<h2 id="tables">Tables (%d files)</h2><p class="small">Column descriptions are in INDEX.md.</p><ul>', length(tab_files)),
          sprintf('<li><a href="tables/%s">%s</a></li>', tab_files, tab_files), '</ul>',
          '<h2 id="numbers">Key numbers</h2><table><tr><th>Subtype</th><th>Expressed genes H / Mar / M</th><th class="num">Jaccard H&ndash;Mar</th><th class="num">Jaccard H&ndash;M</th><th class="num">Jaccard Mar&ndash;M</th><th class="num">r Human&ndash;Marmoset</th><th class="num">r Human&ndash;Mouse</th></tr>')
for (st in run_subtypes) {
  j <- jaccard_all[[st]]; cc <- corr_all[[st]]; ng <- ngenes_all[[st]]
  html <- c(html, sprintf('<tr><td>%s</td><td>%s</td><td class="num">%.2f</td><td class="num">%.2f</td><td class="num">%.2f</td><td class="num">%.2f</td><td class="num">%.2f</td></tr>',
                          esc(st), paste(format(ng$n_expressed, big.mark = ","), collapse = " / "),
                          j$jaccard[j$pair == "Human vs Marmoset"], j$jaccard[j$pair == "Human vs Mouse"], j$jaccard[j$pair == "Marmoset vs Mouse"],
                          cc$r[cc$other == "Marmoset"], cc$r[cc$other == "Mouse"]))
}
html <- c(html, '</table></main></body></html>')
writeLines(html, file.path(OUT_DIR, "index.html"))

message("\n>>> Finished. Skipped subtypes: ", if (length(skipped)) paste(skipped, collapse = ", ") else "none")
message(">>> Outputs in: ", normalizePath(OUT_DIR), "  (see INDEX.md)")

################################################################################
## APPENDIX -- build the 1:1 ortholog table with biomaRt (run once; the CSV is
## then reused offline).  Produces the CSV expected by ORTHOLOG_TABLE.
##   source this script's bottom part, then:  build_ortholog_table_biomart()
################################################################################
build_ortholog_table_biomart <- function(out_csv = "orthologs_1to1_human_marmoset_mouse.csv") {
  library(biomaRt)
  hs <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mm <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                             "mmusculus_homolog_associated_gene_name", "mmusculus_homolog_orthology_type"),
              mart = hs)
  cj <- getBM(attributes = c("ensembl_gene_id",
                             "cjacchus_homolog_associated_gene_name", "cjacchus_homolog_orthology_type"),
              mart = hs)
  tab <- mm %>% inner_join(cj, by = "ensembl_gene_id") %>%
    transmute(human_ensembl = ensembl_gene_id, human_symbol = external_gene_name,
              mouse_symbol = mmusculus_homolog_associated_gene_name,
              marmoset_symbol = cjacchus_homolog_associated_gene_name,
              one2one = mmusculus_homolog_orthology_type == "ortholog_one2one" &
                cjacchus_homolog_orthology_type == "ortholog_one2one") %>%
    filter(human_symbol != "", mouse_symbol != "", marmoset_symbol != "") %>%
    distinct()
  ## a human gene must map to exactly one mouse and one marmoset gene (and vice-versa)
  dup_h <- tab %>% count(human_symbol) %>% filter(n > 1) %>% pull(human_symbol)
  dup_m <- tab %>% count(mouse_symbol) %>% filter(n > 1) %>% pull(mouse_symbol)
  dup_c <- tab %>% count(marmoset_symbol) %>% filter(n > 1) %>% pull(marmoset_symbol)
  tab <- tab %>% mutate(one2one = one2one & !(human_symbol %in% dup_h) &
                          !(mouse_symbol %in% dup_m) & !(marmoset_symbol %in% dup_c))
  write.csv(tab, out_csv, row.names = FALSE)
  message("1:1 orthologs: ", sum(tab$one2one), " written to ", out_csv)
  invisible(tab)
}