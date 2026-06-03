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

outdir <- getwd()

subtype_order_v3 <- c(
  "Naive/CM T cells", "NK/Cytotoxic T cells", "\u03b3\u03b4 T cells",
  "Proliferating T cells", "B cells", "Plasma cells",
  "Tissue-resident Macrophages", "Scavenger Macrophages",
  "Classical Neutrophils", "AXL+ Dendritic cells",
  "Mature Mast Cells (MCTC)"
)

subtype_palette_final <- c(
  "Naive/CM T cells"            = "#F8766D",
  "NK/Cytotoxic T cells"        = "#E58700",
  "\u03b3\u03b4 T cells"      = "#C49A00",
  "Proliferating T cells"       = "#53B400",
  "B cells"                     = "#00C094",
  "Plasma cells"                = "#00B6EB",
  "Tissue-resident Macrophages" = "#06A4FF",
  "Scavenger Macrophages"       = "#9590FF",
  "Classical Neutrophils"       = "#D575FE",
  "AXL+ Dendritic cells"        = "#F962DD",
  "Mature Mast Cells (MCTC)"    = "#FF61C3"
)

species_colors <- c(
  Human    = "#1D9E75",
  Marmoset = "#377EB8",
  Mouse    = "#E41A1C",
  Rat      = "#984EA3"
)

