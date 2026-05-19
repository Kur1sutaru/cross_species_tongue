# =============================================================================
# Epithelial Compartment Re-integration & Subtype Annotation
# Input:    `epithelial` Seurat v5 object (4 species)
# Strategy: Human + Marmoset → Harmony | Mouse + Rat → Harmony → annotate
# Markers:  Tabula Sapiens + Tabula Muris + PanglaoDB (oral/tongue epithelium)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

set.seed(42)

# ─── Settings ─────────────────────────────────────────────────────────────────
setwd("/projects/r_workspace/cristal_data/cross_species/subtypes_cross_species/epithelial")

outdir <- "epithelial_reintegration"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

RESOLUTION_PRIMATE <- 0.5
RESOLUTION_RODENT  <- 0.4
N_PCS              <- 30
N_HARMONY_DIMS     <- 25

species_colors <- c(
  Human    = "#1D9E75",
  Marmoset = "#377EB8",
  Mouse    = "#E41A1C",
  Rat      = "#984EA3"
)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "") +
    theme(
      axis.line         = element_line(linewidth = 0.4),
      axis.ticks        = element_line(linewidth = 0.3),
      axis.text         = element_text(color = "black"),
      legend.background = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold", hjust = 0.5)
    )
}

# =============================================================================
# STEP 0 — Confirm object
# =============================================================================
cat("── Step 0: Confirming epithelial object ──\n")
cat("Total cells:", ncol(epithelial), "\n")
cat("Species breakdown:\n")
print(table(epithelial$species))
cat("Metadata columns:\n")
print(colnames(epithelial@meta.data))

# =============================================================================
# STEP 1 — Split into primate and rodent
# =============================================================================
cat("\n── Step 1: Splitting into primate and rodent groups ──\n")

epi_primate <- subset(epithelial, subset = species %in% c("Human", "Marmoset"))
epi_rodent  <- subset(epithelial, subset = species %in% c("Mouse", "Rat"))

cat("Primate — Human:", sum(epi_primate$species == "Human"),
    "| Marmoset:", sum(epi_primate$species == "Marmoset"),
    "| Total:", ncol(epi_primate), "\n")
cat("Rodent  — Mouse:", sum(epi_rodent$species == "Mouse"),
    "| Rat:", sum(epi_rodent$species == "Rat"),
    "| Total:", ncol(epi_rodent), "\n")

# =============================================================================
# STEP 2 — Pre-processing helper (NormalizeData → PCA)
# =============================================================================
preprocess_seurat <- function(seu, group_name, n_pcs = N_PCS) {
  cat(sprintf("\n  [%s] Normalizing...\n", group_name))
  seu <- NormalizeData(seu, verbose = FALSE)
  
  cat(sprintf("  [%s] Finding variable features (split by species)...\n", group_name))
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])           # join first in case already split
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$species)
  seu <- FindVariableFeatures(seu, selection.method = "vst",
                              nfeatures = 3000, verbose = FALSE)
  
  cat(sprintf("  [%s] Scaling...\n", group_name))
  seu <- ScaleData(seu, verbose = FALSE)
  
  cat(sprintf("  [%s] PCA...\n", group_name))
  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  
  return(seu)
}

# =============================================================================
# STEP 3 — Harmony integration helper
# =============================================================================
run_harmony <- function(seu, group_name, resolution,
                        n_harmony_dims = N_HARMONY_DIMS) {
  cat(sprintf("\n  [%s] Running Harmony...\n", group_name))
  seu <- IntegrateLayers(
    object         = seu,
    method         = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction  = "harmony",
    group.by.vars  = "species",
    verbose        = TRUE
  )
  
  cat(sprintf("  [%s] Joining layers...\n", group_name))
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  
  cat(sprintf("  [%s] FindNeighbors...\n", group_name))
  seu <- FindNeighbors(seu, reduction = "harmony",
                       dims = 1:n_harmony_dims, verbose = FALSE)
  
  cat(sprintf("  [%s] FindClusters (res = %s)...\n", group_name, resolution))
  seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
  
  cat(sprintf("  [%s] UMAP...\n", group_name))
  seu <- RunUMAP(seu, reduction = "harmony",
                 dims = 1:n_harmony_dims, verbose = FALSE)
  
  cat(sprintf("  [%s] → %d clusters\n",
              group_name, length(levels(seu$seurat_clusters))))
  print(table(seu$seurat_clusters))
  return(seu)
}

# ── Run pipeline ──────────────────────────────────────────────────────────────
epi_primate <- preprocess_seurat(epi_primate, "Primate")
epi_primate <- run_harmony(epi_primate, "Primate", RESOLUTION_PRIMATE)
saveRDS(epi_primate, file.path(outdir, "epi_primate.rds"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
epi_rodent  <- preprocess_seurat(epi_rodent,  "Rodent")
epi_rodent  <- run_harmony(epi_rodent,  "Rodent",  RESOLUTION_RODENT)
saveRDS(epi_rodent, file.path(outdir, "epi_rodent.rds"))

cat("\n  Both objects integrated and saved.\n")

# =============================================================================
# STEP 4 — Canonical epithelial marker genes
#           Sources: Tabula Sapiens (2022 Science), Tabula Muris (2018 Nature),
#                    PanglaoDB (Franzén et al. 2019), oral/tongue literature
# =============================================================================

# ── Subtype-specific markers ──────────────────────────────────────────────────
subtype_markers <- list(
  
  # ── Basal / stem-like ─────────────────────────────────────────────────────
  # Tabula Sapiens: basal cells (tongue/oral); PanglaoDB: basal epithelial
  "Basal"                   = c("KRT5", "KRT14", "KRT15", "TP63",
                                "ITGA6", "COL17A1", "DST", "LAMA3"),
  
  # ── Suprabasal / spinous ──────────────────────────────────────────────────
  # Tabula Sapiens: suprabasal keratinocytes; oral squamous epithelium
  "Suprabasal/Spinous"      = c("KRT1", "KRT10", "KRT4", "KRT13",
                                "IVL", "SPRR1A", "SPRR1B", "DSG1"),
  
  # ── Superficial / cornified ───────────────────────────────────────────────
  # Tabula Muris: differentiated epithelium; cornification program
  "Superficial/Cornified"   = c("KRT2", "FLG", "LOR", "SPRR3",
                                "CDSN", "CALML5", "KRTDAP", "SBSN"),
  
  # ── Simple / columnar epithelium ──────────────────────────────────────────
  # Tabula Sapiens: simple epithelial (ductal/glandular); Tabula Muris: gland
  "Simple Epithelial"       = c("KRT8", "KRT18", "KRT19", "EPCAM",
                                "CDH1", "CLDN3", "CLDN4", "MUC1"),
  
  # ── Taste bud cells ───────────────────────────────────────────────────────
  # Tongue-specific; Tabula Muris tongue dataset
  "Taste receptor (Type II)" = c("TRPM5", "PLCB2", "GNAT3", "TAS1R3",
                                 "ITPR3", "SNAP25", "NCAM1"),
  
  "Taste receptor (Type III)" = c("PKDL2", "SNAP25", "SYP",
                                  "CACNA1A", "NCAM1", "ADCY3"),
  
  "Taste bud (Type I/Glial)" = c("NTPDase2", "GLAST", "KCNK5",
                                 "ENTPD2", "SLC1A3"),
  
  # ── Filiform / non-taste papillae ─────────────────────────────────────────
  # Tongue dorsal surface specific
  "Filiform papillae"       = c("KRT4", "KRT13", "KRTDAP", "CDSN",
                                "RHCG", "AQP3"),
  
  # ── Proliferating / transit-amplifying ───────────────────────────────────
  # Tabula Sapiens + Tabula Muris: cycling epithelial cells
  "Proliferating"           = c("MKI67", "TOP2A", "PCNA", "STMN1",
                                "CENPF", "UBE2C", "BIRC5", "CDC20"),
  
  # ── Progenitor / stem ─────────────────────────────────────────────────────
  # Tabula Sapiens: epithelial progenitors; tongue stem niche
  "Progenitor/Stem"         = c("SOX2", "TP63", "KRT15", "ITGA6",
                                "LGR5", "LRIG1", "ALDH1A1", "CDH1"),
  
  # ── Suprabasal differentiating ────────────────────────────────────────────
  # Intermediate differentiation state; Tabula Muris oral epithelium
  "Differentiating"         = c("KRT6A", "KRT6B", "KRT16", "KRT17",
                                "S100A2", "S100A8", "S100A9", "SPRR2A"),
  
  # ── Secretory / mucous ────────────────────────────────────────────────────
  # Tabula Sapiens: secretory cells; salivary gland-associated
  "Secretory/Mucous"        = c("MUC5B", "MUC5AC", "MUC4", "BPIFB1",
                                "BPIFB2", "LTF", "PRB1", "TFF1"),
  
  # ── Ionocyte / rare secretory ─────────────────────────────────────────────
  # Tabula Sapiens: rare airway/oral ionocytes
  "Ionocyte"                = c("FOXI1", "CFTR", "ATP6V1B1",
                                "ATP6V0D2", "ASCL3"),
  
  # ── Pan-epithelial ────────────────────────────────────────────────────────
  # Tabula Sapiens + Tabula Muris: broad epithelial markers
  "Pan-epithelial"          = c("EPCAM", "CDH1", "KRT8", "KRT18",
                                "CLDN7", "OCLN")
)

# ── Ordered gene list for DotPlot (lineage order) ─────────────────────────────
genes_ordered <- unique(unlist(subtype_markers))

# ── Subtype order for Y axis ──────────────────────────────────────────────────
subtype_order <- c(
  "Pan-epithelial",
  "Basal",
  "Progenitor/Stem",
  "Proliferating",
  "Suprabasal/Spinous",
  "Differentiating",
  "Superficial/Cornified",
  "Simple Epithelial",
  "Secretory/Mucous",
  "Ionocyte",
  "Filiform papillae",
  "Taste receptor (Type II)",
  "Taste receptor (Type III)",
  "Taste bud (Type I/Glial)"
)

# ── Color palette ─────────────────────────────────────────────────────────────
subtype_palette <- c(
  "Pan-epithelial"            = "#AAAAAA",
  "Basal"                     = "#E41A1C",
  "Progenitor/Stem"           = "#FF7F00",
  "Proliferating"             = "#FFED6F",
  "Suprabasal/Spinous"        = "#4DAF4A",
  "Differentiating"           = "#80B1D3",
  "Superficial/Cornified"     = "#377EB8",
  "Simple Epithelial"         = "#984EA3",
  "Secretory/Mucous"          = "#A65628",
  "Ionocyte"                  = "#F781BF",
  "Filiform papillae"         = "#B2DF8A",
  "Taste receptor (Type II)"  = "#1F78B4",
  "Taste receptor (Type III)" = "#33A02C",
  "Taste bud (Type I/Glial)"  = "#FB9A99"
)

# =============================================================================
# STEP 5 — Module scores per subtype
# =============================================================================
cat("\n── Step 5: Module scores ──\n")

add_module_scores <- function(seu, markers, obj_name) {
  for (nm in names(markers)) {
    genes_present <- intersect(markers[[nm]], rownames(seu))
    if (length(genes_present) == 0) {
      cat(sprintf("  [%s] Skipping '%s' — no genes found\n", obj_name, nm))
      next
    }
    safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
    seu <- AddModuleScore(seu,
                          features = list(genes_present),
                          name     = paste0("score_", safe_nm),
                          ctrl     = 50, seed = 42)
    colnames(seu@meta.data)[ncol(seu@meta.data)] <- paste0("score_", safe_nm)
  }
  return(seu)
}

epi_primate <- add_module_scores(epi_primate, subtype_markers, "Primate")
epi_rodent  <- add_module_scores(epi_rodent,  subtype_markers, "Rodent")

# =============================================================================
# STEP 6 — FindAllMarkers (save CSVs for inspection)
# =============================================================================
cat("\n── Step 6: FindAllMarkers ──\n")

Idents(epi_primate) <- "seurat_clusters"
markers_primate <- FindAllMarkers(epi_primate, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.3,
                                  test.use = "wilcox", verbose = FALSE)
write.csv(markers_primate,
          file.path(outdir, "markers_epi_primate_clusters.csv"),
          row.names = FALSE)

Idents(epi_rodent) <- "seurat_clusters"
markers_rodent <- FindAllMarkers(epi_rodent, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.3,
                                 test.use = "wilcox", verbose = FALSE)
write.csv(markers_rodent,
          file.path(outdir, "markers_epi_rodent_clusters.csv"),
          row.names = FALSE)

cat("  Marker CSVs saved.\n")

# =============================================================================
# STEP 7 — Exploratory plots (pre-annotation)
# =============================================================================
cat("\n── Step 7: Exploratory plots ──\n")

# Key feature genes for FeaturePlot
key_features <- c(
  "KRT5", "KRT14", "TP63",          # basal
  "KRT4", "KRT13", "IVL",           # suprabasal
  "KRT1", "KRT10", "FLG",           # cornified
  "KRT8", "KRT18", "EPCAM",         # simple
  "MKI67", "TOP2A",                 # proliferating
  "SOX2", "LGR5",                   # progenitor
  "TRPM5", "GNAT3",                 # taste II
  "MUC5B", "TFF1"                   # secretory
)

# DotPlot genes (representative subset)
dotplot_genes <- c(
  "EPCAM", "CDH1", "KRT8", "KRT18",
  "KRT5", "KRT14", "KRT15", "TP63", "ITGA6",
  "SOX2", "LGR5", "LRIG1", "ALDH1A1",
  "MKI67", "TOP2A", "PCNA", "STMN1",
  "KRT1", "KRT10", "KRT4", "KRT13", "IVL", "SPRR1A",
  "KRT6A", "KRT16", "S100A2",
  "KRT2", "FLG", "LOR", "KRTDAP",
  "MUC5B", "MUC5AC", "LTF", "BPIFB1",
  "FOXI1", "CFTR",
  "TRPM5", "PLCB2", "GNAT3",
  "SYP", "SNAP25"
)

make_plots <- function(seu, group_name, spc_colors, dot_genes, feat_genes) {
  
  # UMAP cluster + species
  p_clust <- DimPlot(seu, group.by = "seurat_clusters",
                     label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle(paste(group_name, "— clusters")) + theme_pub() +
    theme(legend.position = "none")
  p_spc <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  
  cairo_pdf(file.path(outdir,
                      paste0("01_umap_", tolower(gsub(" ", "_", group_name)),
                             "_exploratory.pdf")),
            width = 14, height = 6)
  print(p_clust | p_spc)
  dev.off()
  
  # DotPlot
  genes_ok <- intersect(dot_genes, rownames(seu))
  Idents(seu) <- "seurat_clusters"
  p_dot <- suppressWarnings(
    DotPlot(seu, features = genes_ok, dot.scale = 5) +
      RotatedAxis() +
      scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                            midpoint = 0) +
      ggtitle(paste(group_name, "— canonical epithelial markers")) +
      theme_pub(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7,
                                       face = "italic"),
            legend.key.size = unit(0.35, "cm"))
  )
  
  cairo_pdf(file.path(outdir,
                      paste0("02_dotplot_", tolower(gsub(" ", "_", group_name)),
                             ".pdf")),
            width = 18, height = 7)
  print(p_dot)
  dev.off()
  
  # FeaturePlot
  feats_ok <- intersect(feat_genes, rownames(seu))
  p_feat <- FeaturePlot(seu, features = feats_ok,
                        order = TRUE, ncol = 4,
                        cols = c("lightgrey", "#a50026"), pt.size = 0.3) &
    theme_pub(base_size = 9)
  
  cairo_pdf(file.path(outdir,
                      paste0("03_featureplot_", tolower(gsub(" ", "_", group_name)),
                             ".pdf")),
            width = 18, height = 12)
  print(p_feat)
  dev.off()
  
  cat(sprintf("  [%s] Exploratory plots saved.\n", group_name))
}

make_plots(epi_primate, "Primate", species_colors, dotplot_genes, key_features)
make_plots(epi_rodent,  "Rodent",  species_colors, dotplot_genes, key_features)

# =============================================================================
# STEP 8 — ANNOTATION
# ─────────────────────────────────────────────────────────────────────────────
# Inspect CSV markers + PDFs from Steps 6-7 then fill in below
# =============================================================================
cat("\n── Step 8: Annotation ──\n")
cat("Primate clusters:", levels(epi_primate$seurat_clusters), "\n")
cat("Rodent clusters: ", levels(epi_rodent$seurat_clusters),  "\n")

# ─── Fill in after marker inspection ─────────────────────────────────────────
primate_annotations <- c(
  "0"  = "Differentiating",
  "1"  = "Superficial/Cornified",
  "2"  = "Differentiating",
  "3"  = "Suprabasal/Spinous",
  "4"  = "Suprabasal/Spinous",
  "5"  = "FibroImmune signatures?",        # Fibroblast/Immune — VIM, DCN, IL1B
  "6"  = "T cell?",        # T cells — LTB, IL10, ICOS, IL22
  "7"  = "Proliferating",
  "8"  = "Taste receptor (Type II)",  # CDH2, TRPM3 — confirm below
  "9"  = "??",        # Neuronal — SLC6A1, CCSER1
  "10" = "??",        # Neuronal — DLGAP2, ADARB2
  "11" = "Basal",
  "12" = "Lung alveolar???",        # Lung alveolar — SFTPC, NAPSA
  "13" = "stromal?",        # Stromal — EBF1, LAMA2
  "14" = "Superficial/Cornified",
  "15" = "Secretory/Mucous",
  "16" = "Progenitor/Stem",
  "17" = "??",        # Neuronal — GRIA3, PCDH11X
  "18" = "Simple Epithelial",
  "19" = "Mesenchymal",        # Mesenchymal — CCN2, IGF2
  "20" = "???"                 # need markers
)

primate_annotations <- c(
  "0"  = "Differentiating",
  "1"  = "Superficial/Cornified",
  "2"  = "Differentiating",
  "3"  = "Suprabasal/Spinous",
  "4"  = "Suprabasal/Spinous",
  "5"  = "Contaminant",
  "6"  = "Contaminant",
  "7"  = "Proliferating",
  "8"  = "Taste receptor (Type II)",
  "9"  = "Contaminant",
  "10" = "Contaminant",
  "11" = "Basal",
  "12" = "Contaminant",
  "13" = "Contaminant",
  "14" = "Superficial/Cornified",
  "15" = "Secretory/Mucous",
  "16" = "Progenitor/Stem",
  "17" = "Contaminant",
  "18" = "Simple Epithelial",
  "19" = "Contaminant",
  "20" = "Contaminant"   # placeholder — still need those markers
)

rodent_annotations <- c(
  "0"  = "Contaminant",
  "1"  = "Contaminant",
  "2"  = "Contaminant",
  "3"  = "Taste receptor (Type II)",
  "4"  = "Proliferating",
  "5"  = "Suprabasal/Spinous",
  "6"  = "Progenitor/Stem",
  "7"  = "Proliferating",
  "8"  = "Basal",
  "9"  = "Differentiating",
  "10" = "Contaminant",
  "11" = "Secretory/Mucous",
  "12" = "Contaminant",
  "13" = "Contaminant"
)

# Apply annotations
if (length(primate_annotations) > 0) {
  epi_primate$epi_subtype <- unname(
    primate_annotations[as.character(epi_primate$seurat_clusters)]
  )
} else {
  epi_primate$epi_subtype <- paste0("Cluster_", epi_primate$seurat_clusters)
  cat("  Primate: placeholder labels set — fill annotations above.\n")
}

if (length(rodent_annotations) > 0) {
  epi_rodent$epi_subtype <- unname(
    rodent_annotations[as.character(epi_rodent$seurat_clusters)]
  )
} else {
  epi_rodent$epi_subtype <- paste0("Cluster_", epi_rodent$seurat_clusters)
  cat("  Rodent: placeholder labels set — fill annotations above.\n")
}

# ─── Broad category grouping ──────────────────────────────────────────────────
category_map <- c(
  "Basal"                     = "Stratified squamous",
  "Progenitor/Stem"           = "Stratified squamous",
  "Proliferating"             = "Cycling",
  "Suprabasal/Spinous"        = "Stratified squamous",
  "Differentiating"           = "Stratified squamous",
  "Superficial/Cornified"     = "Stratified squamous",
  "Simple Epithelial"         = "Simple/Glandular",
  "Secretory/Mucous"          = "Simple/Glandular",
  "Ionocyte"                  = "Simple/Glandular",
  "Filiform papillae"         = "Specialized",
  "Taste receptor (Type II)"  = "Specialized",
  "Taste receptor (Type III)" = "Specialized",
  "Taste bud (Type I/Glial)"  = "Specialized",
  "Pan-epithelial"            = "Pan-epithelial"
)

assign_category <- function(subtypes, cmap) {
  cat <- cmap[subtypes]
  cat[is.na(cat)] <- "Unannotated"
  unname(cat)
}

epi_primate$epi_category <- assign_category(epi_primate$epi_subtype, category_map)
epi_rodent$epi_category  <- assign_category(epi_rodent$epi_subtype,  category_map)
epi_primate$epi_category <- assign_category(epi_primate$epi_subtype, category_map)
epi_rodent$epi_category  <- assign_category(epi_rodent$epi_subtype,  category_map)

# ↓ ADD HERE — verify categories look right
cat("\nPrimate epi_category distribution:\n")
print(table(epi_primate$epi_category, useNA = "ifany"))
cat("\nRodent epi_category distribution:\n")
print(table(epi_rodent$epi_category, useNA = "ifany"))

# ── Check Rat subtype counts ──────────────────────────────────────────────────
cat("\nRat subtype counts (from rodent object):\n")
print(table(epi_rodent$seurat_clusters[epi_rodent$species == "Rat"], useNA = "ifany"))

# =============================================================================
# STEP 9 — Annotated UMAPs
# =============================================================================
cat("\n── Step 9: Annotated UMAPs ──\n")

plot_annotated <- function(seu, group_name, palette, spc_colors,
                           subtype_col = "epi_subtype") {
  p1 <- DimPlot(seu, group.by = subtype_col, cols = palette,
                label = TRUE, label.size = 3, repel = TRUE) +
    ggtitle(paste(group_name, "— epithelial subtypes")) + theme_pub() +
    theme(legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  p2 <- DimPlot(seu, group.by = "species", cols = spc_colors) +
    ggtitle("Species") + theme_pub()
  p3 <- DimPlot(seu, group.by = "epi_category",
                label = TRUE, label.size = 3.5, repel = TRUE) +
    ggtitle("Category") + theme_pub()
  p1 | p2 | p3
}

cairo_pdf(file.path(outdir, "04_umap_primate_annotated.pdf"),
          width = 20, height = 6)
print(plot_annotated(epi_primate, "Primate (Human + Marmoset)",
                     subtype_palette, species_colors))
dev.off()

cairo_pdf(file.path(outdir, "04_umap_rodent_annotated.pdf"),
          width = 20, height = 6)
print(plot_annotated(epi_rodent, "Rodent (Mouse + Rat)",
                     subtype_palette, species_colors))
dev.off()

cat("  Annotated UMAPs saved.\n")

# =============================================================================
# STEP 10 — Merge primate + rodent → integrate all 4 species
# =============================================================================
cat("\n── Step 10: Merging primate + rodent for full integration ──\n")

epi_integrated <- merge(epi_primate, y = epi_rodent, merge.data = TRUE)
cat("Total cells after merge:", ncol(epi_integrated), "\n")

epi_integrated[["RNA"]] <- JoinLayers(epi_integrated[["RNA"]])
epi_integrated[["RNA"]] <- split(epi_integrated[["RNA"]],
                                 f = epi_integrated$species)

epi_integrated <- NormalizeData(epi_integrated, verbose = FALSE)
epi_integrated <- FindVariableFeatures(epi_integrated, nfeatures = 3000,
                                       verbose = FALSE)
epi_integrated <- ScaleData(epi_integrated, verbose = FALSE)
epi_integrated <- RunPCA(epi_integrated, npcs = N_PCS, verbose = FALSE)

cat("Running Harmony on merged object...\n")
epi_integrated <- IntegrateLayers(
  object         = epi_integrated,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  group.by.vars  = "species",
  verbose        = TRUE
)

epi_integrated[["RNA"]] <- JoinLayers(epi_integrated[["RNA"]])
epi_integrated <- FindNeighbors(epi_integrated, reduction = "harmony",
                                dims = 1:N_HARMONY_DIMS, verbose = FALSE)
epi_integrated <- FindClusters(epi_integrated, resolution = 0.5,
                               verbose = FALSE)
epi_integrated <- RunUMAP(epi_integrated, reduction = "harmony",
                          dims = 1:N_HARMONY_DIMS, verbose = FALSE)

cat("Integrated clusters:", length(levels(epi_integrated$seurat_clusters)), "\n")

# =============================================================================
# STEP 11 — Propagate annotations to integrated object
# =============================================================================
cat("\n── Step 11: Propagating annotations ──\n")

all_meta <- rbind(
  data.frame(cell         = colnames(epi_primate),
             epi_subtype  = epi_primate$epi_subtype,
             epi_category = epi_primate$epi_category,
             stringsAsFactors = FALSE),
  data.frame(cell         = colnames(epi_rodent),
             epi_subtype  = epi_rodent$epi_subtype,
             epi_category = epi_rodent$epi_category,
             stringsAsFactors = FALSE)
)

idx <- match(colnames(epi_integrated), all_meta$cell)
cat("Cells matched:", sum(!is.na(idx)), "/", ncol(epi_integrated), "\n")

epi_integrated$epi_subtype  <- all_meta$epi_subtype[idx]
epi_integrated$epi_category <- all_meta$epi_category[idx]

cat("Subtype distribution:\n")
print(table(epi_integrated$epi_subtype, useNA = "ifany"))

# =============================================================================
# STEP 12 — Annotated DotPlot per species (with gene-group annotation bar)
# =============================================================================
cat("\n── Step 12: Annotated DotPlots per species ──\n")

# Gene → subtype lookup table
gene_annotation <- stack(subtype_markers) %>%
  rename(gene = values, subtype = ind) %>%
  mutate(subtype = as.character(subtype)) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  filter(gene %in% genes_ordered)

gene_annotation$gene    <- factor(gene_annotation$gene,    levels = genes_ordered)
gene_annotation$subtype <- factor(gene_annotation$subtype, levels = names(subtype_markers))

species_list        <- c("Human", "Marmoset", "Mouse", "Rat")
subtype_order_present <- intersect(subtype_order,
                                   unique(epi_integrated$epi_subtype[
                                     !is.na(epi_integrated$epi_subtype)]))

# Check gene availability
all_markers    <- unique(unlist(subtype_markers))
gene_presence  <- sapply(species_list, function(sp) {
  cells_sp <- colnames(epi_integrated)[epi_integrated$species == sp]
  rn <- rownames(epi_integrated[, cells_sp])
  all_markers %in% rn
})
rownames(gene_presence) <- all_markers
genes_keep    <- all_markers[rowSums(gene_presence) >= 2]
genes_ordered_filt <- intersect(genes_ordered, genes_keep)

cat("\nGenes available for DotPlot:", length(genes_ordered_filt), "/",
    length(genes_ordered), "\n")

make_annotated_dotplot <- function(seu, species_name, genes,
                                   subtype_ord, ann_palette,
                                   gene_ann, subtype_col = "epi_subtype",
                                   base_size = 9) {
  cells_sp  <- colnames(seu)[seu$species == species_name]
  sub_obj   <- subset(seu, cells = cells_sp)
  genes_sp  <- intersect(genes, rownames(sub_obj))
  
  subtypes_sp <- intersect(subtype_ord,
                           unique(sub_obj[[subtype_col]][
                             !is.na(sub_obj[[subtype_col]])]))
  sub_obj[[subtype_col]] <- factor(sub_obj[[subtype_col]], levels = subtypes_sp)
  Idents(sub_obj) <- subtype_col
  
  cat(sprintf("  [%s] %d genes | %d subtypes\n",
              species_name, length(genes_sp), length(subtypes_sp)))
  
  dot_data <- suppressWarnings(
    DotPlot(sub_obj, features = genes_sp, scale = TRUE,
            col.min = -2, col.max = 2)$data
  ) %>%
    rename(gene = features.plot, subtype_y = id) %>%
    mutate(gene      = factor(gene, levels = genes_sp),
           subtype_y = factor(subtype_y, levels = subtypes_sp))
  
  # Main DotPlot
  p_main <- ggplot(dot_data,
                   aes(x = gene, y = subtype_y,
                       size = pct.exp, color = avg.exp.scaled)) +
    geom_point() +
    scale_color_gradient2(low = "#313695", mid = "white", high = "#a50026",
                          midpoint = 0, name = "Avg\nexp",
                          limits = c(-2, 2), oob = scales::squish) +
    scale_size_continuous(name = "% exp", range = c(0.2, 5),
                          breaks = c(10, 25, 50, 75, 100)) +
    labs(x = NULL, y = NULL, title = species_name) +
    theme_pub(base_size = base_size) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1,
                                          size = 6.5, face = "italic"),
          axis.text.y      = element_text(size = 8),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
          legend.key.size  = unit(0.3, "cm"),
          legend.text      = element_text(size = 6.5),
          legend.title     = element_text(size = 7),
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 10),
          plot.margin      = margin(t = 2, r = 5, b = 2, l = 5))
  
  # Annotation color bar
  ann_data <- gene_ann %>%
    filter(gene %in% genes_sp) %>%
    mutate(gene = factor(gene, levels = genes_sp), y = 1)
  
  p_ann <- ggplot(ann_data, aes(x = gene, y = y, fill = subtype)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = ann_palette, drop = FALSE) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(t = 2, r = 5, b = 0, l = 5))
  
  # Subtype label strip
  label_data <- ann_data %>%
    mutate(gene_num = as.integer(gene)) %>%
    group_by(subtype) %>%
    summarise(xmid = mean(gene_num),
              xmin = min(gene_num) - 0.5,
              xmax = max(gene_num) + 0.5,
              .groups = "drop") %>%
    mutate(y = 1)
  
  p_label <- ggplot(label_data,
                    aes(x = xmid, y = y, label = subtype, color = subtype)) +
    geom_text(size = 2.2, fontface = "bold", hjust = 0.5, vjust = 0.5,
              show.legend = FALSE) +
    geom_segment(aes(x = xmin, xend = xmax, y = 0.45, yend = 0.45,
                     color = subtype), linewidth = 0.5, show.legend = FALSE) +
    scale_color_manual(values = ann_palette) +
    scale_x_continuous(limits = c(0.5, length(genes_sp) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(t = 0, r = 5, b = 1, l = 5))
  
  p_ann / p_label / p_main + plot_layout(heights = c(0.06, 0.1, 1))
}

# Generate and save
plots_ann <- lapply(species_list, function(sp) {
  make_annotated_dotplot(
    seu          = epi_integrated,
    species_name = sp,
    genes        = genes_ordered_filt,
    subtype_ord  = subtype_order_present,
    ann_palette  = subtype_palette,
    gene_ann     = gene_annotation
  )
})
names(plots_ann) <- species_list
plots_ann <- Filter(Negate(is.null), plots_ann)

for (sp in names(plots_ann)) {
  cairo_pdf(file.path(outdir, paste0("05_dotplot_annotated_", tolower(sp), ".pdf")),
            width = 18, height = 8)
  print(plots_ann[[sp]])
  dev.off()
  cat(sprintf("  Saved: 05_dotplot_annotated_%s.pdf\n", tolower(sp)))
}

cairo_pdf(file.path(outdir, "05_dotplot_annotated_all_species.pdf"),
          width = 20, height = 30)
print(wrap_plots(plots_ann, ncol = 1) +
        plot_annotation(
          title = "Epithelial subtype marker expression per species",
          theme = theme(plot.title = element_text(face = "bold",
                                                  hjust = 0.5, size = 13))
        ))
dev.off()

# =============================================================================
# STEP 13 — Cell proportion barplot
# =============================================================================
cat("\n── Step 13: Cell proportion barplot ──\n")

prop_df <- epi_integrated@meta.data %>%
  filter(!is.na(epi_subtype)) %>%
  group_by(species, epi_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(total      = sum(n),
         proportion = n / total * 100) %>%
  ungroup() %>%
  mutate(epi_subtype = factor(epi_subtype, levels = subtype_order),
         species     = factor(species, levels = c("Human", "Marmoset", "Mouse", "Rat")))

totals_df <- prop_df %>% distinct(species, total)
label_df  <- prop_df %>%
  filter(proportion >= 3) %>%
  group_by(species) %>%
  arrange(epi_subtype) %>%
  mutate(label_y = cumsum(proportion) - proportion / 2) %>%
  ungroup()

p_bar <- ggplot(prop_df, aes(x = species, y = proportion, fill = epi_subtype)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(data = label_df,
            aes(y = label_y, label = n),
            size = 2.8, color = "white", fontface = "bold") +
  geom_text(data = totals_df,
            aes(x = species, y = 103,
                label = paste0("n = ", comma(total)), fill = NULL),
            size = 3.2, fontface = "bold", color = "black") +
  scale_fill_manual(values = subtype_palette, name = "Epithelial subtype") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.07)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Cell proportion (%)",
       title = "Epithelial subtype composition per species") +
  theme_pub(base_size = 11) +
  theme(legend.position    = "right",
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(0.45, "cm"),
        axis.text.x        = element_text(size = 11, face = "bold"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)) +
  guides(fill = guide_legend(ncol = 1))

cairo_pdf(file.path(outdir, "06_barplot_epi_proportion_by_species.pdf"),
          width = 10, height = 7)
print(p_bar)
dev.off()

write.csv(prop_df, file.path(outdir, "epi_subtype_counts_by_species.csv"),
          row.names = FALSE)
cat("  Barplot and count table saved.\n")

# =============================================================================
# STEP 14 — UMAP split by species
# =============================================================================
cat("\n── Step 14: UMAP split by species ──\n")

epi_integrated$epi_subtype <- factor(epi_integrated$epi_subtype,
                                     levels = subtype_order_present)

panels <- lapply(c("Human", "Marmoset", "Mouse", "Rat"), function(sp) {
  n_cells <- sum(epi_integrated$species == sp, na.rm = TRUE)
  DimPlot(epi_integrated,
          group.by = "epi_subtype",
          cells    = colnames(epi_integrated)[epi_integrated$species == sp],
          cols     = subtype_palette,
          label    = FALSE, pt.size = 0.4) +
    ggtitle(paste0(sp, "\n(n = ", comma(n_cells), ")")) +
    theme_pub(base_size = 9) +
    theme(legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
})

legend_plot <- DimPlot(epi_integrated, group.by = "epi_subtype",
                       cols = subtype_palette) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme(legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))
shared_legend <- cowplot::get_legend(legend_plot)

cairo_pdf(file.path(outdir, "07_umap_split_species_epi_subtype.pdf"),
          width = 20, height = 5)
print((panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) |
        wrap_elements(shared_legend))
dev.off()


#######################
## Best umap EVER #####
#######################
p_all <- DimPlot(epi_clean,
                 group.by  = "epi_subtype",
                 cols      = subtype_palette,
                 label     = TRUE, label.size = 3.2,
                 repel     = TRUE, pt.size = 0.3,
                 label.box = TRUE) +   # adds white box behind labels
  ggtitle("Epithelial subtypes — all species") +
  theme_pub() +
  theme(legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))



# =============================================================================
# STEP 15 — Save all objects
# =============================================================================
cat("\n── Step 15: Saving objects ──\n")
saveRDS(epi_primate,    file.path(outdir, "epi_primate_annotated.rds"))
saveRDS(epi_rodent,     file.path(outdir, "epi_rodent_annotated.rds"))
saveRDS(epi_integrated, file.path(outdir, "epi_integrated_all4.rds"))
cat("  All objects saved to:", outdir, "\n")

cat("\n── Session info ──\n")
sessionInfo()