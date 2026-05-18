library(Seurat)
library(ggplot2)

DefaultAssay(seurat_obj) <- "RNA"

chosen_palette <- c("#D3D1C7", "#00BFFF", "#005AFF")

genes   <- c("KRT5", "CDH1", "PTPRC", "FLT1", "PECAM1", "PDGFRA", "LUM", "ACTA2", "PLP1")
species <- unique(seurat_obj$species)

out_dir <- "feature_plots_by_species"
dir.create(out_dir, showWarnings = FALSE)

# Expression threshold to split gray vs colored dots
expr_threshold <- 0.1

for (sp in species) {

  obj_sp <- subset(seurat_obj, subset = species == sp)
  umap_coords <- as.data.frame(Embeddings(obj_sp, reduction = "umap"))
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

  for (gene in genes) {

    # Pull expression values
    expr <- FetchData(obj_sp, vars = gene)
    umap_coords$expr <- expr[[1]]

    # Split into low and high expression
    df_gray <- umap_coords[umap_coords$expr <= expr_threshold, ]
    df_blue <- umap_coords[umap_coords$expr >  expr_threshold, ]
    df_blue <- df_blue[order(df_blue$expr), ]  # order so brightest on top

    p <- ggplot() +
      geom_point(data = df_gray,
                 aes(x = UMAP_1, y = UMAP_2),
                 color = "#D3D1C7", size = 0.2, alpha = 0.6) +
      geom_point(data = df_blue,
                 aes(x = UMAP_1, y = UMAP_2, color = expr),
                 size = 0.7, alpha = 0.9) +
      scale_color_gradientn(
        colors = c("#D3D1C7", "#00BFFF", "#005AFF"),
        limits = c(0, 2.5),
        name   = "Expression"
      ) +
      ggtitle(paste0(gene, " — ", sp)) +
      theme_void() +
      theme(
        plot.title        = element_text(size = 12, hjust = 0.5, face = "bold"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.25, "cm"),
        legend.text       = element_text(size = 8),
        plot.margin       = margin(10, 10, 10, 10)
      )

    filename <- file.path(out_dir, paste0(sp, "_", gene, ".pdf"))
    ggsave(filename, plot = p, width = 6, height = 5, dpi = 300)
    message("Saved: ", filename)
  }
}

message("Done! Files saved to: ", out_dir)