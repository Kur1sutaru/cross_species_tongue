###############################################################################
## Interactive cell selection: drag-select a region directly on the UMAP
## instead of guessing coordinate bounds, then rename those cells.
##
## IMPORTANT: this must be run interactively (RStudio, with a plot window/
## Plots pane active) -- it will NOT work in a non-interactive Rscript batch
## run, since it opens a live plot and waits for you to click/drag on it.
###############################################################################

library(Seurat)
library(ggplot2)
library(dplyr)

group_col <- "epi_subtype_v4"
stopifnot(group_col %in% colnames(epi_clean@meta.data))
stopifnot("umap" %in% Reductions(epi_clean))

## ---- Option A: Seurat's built-in CellSelector() -- rectangular drag-select ----
## 1. This opens the UMAP in your Plots pane (or a plot window).
## 2. Click-and-drag a box over the region you want (e.g. the Stress Responsive
##    patch sitting inside Basal Keratinocytes).
## 3. Double-click inside the box (or follow the on-screen prompt) to confirm.
## 4. The selected cell barcodes are returned into `selected_cells`.

p_umap <- DimPlot(epi_clean, reduction = "umap", group.by = group_col, label = FALSE)

selected_cells <- CellSelector(p_umap)   ## <-- interactive step happens here

message("Selected ", length(selected_cells), " cells")

## ---- Check what labels the selected cells actually carry before renaming ----
message("\n--- Label breakdown of your selection ---")
print(table(as.character(epi_clean[[group_col]][selected_cells, ])))

## ---- Rename only cells matching your confirmed target label ----
## Adjust TARGET_LABEL to whatever the table above shows is actually in there.
TARGET_LABEL <- "Stress Responsive"
NEW_LABEL    <- "Basal Keratinocytes"

RUN_RENAME <- FALSE   ## <-- set TRUE once you've confirmed the selection/label above

if (RUN_RENAME) {
  confirmed_cells <- selected_cells[as.character(epi_clean[[group_col]][selected_cells, ]) == TARGET_LABEL]
  message("\nReassigning ", length(confirmed_cells), " cells from '", TARGET_LABEL,
          "' to '", NEW_LABEL, "'")
  
  new_labels <- as.character(epi_clean[[group_col]][, 1])
  names(new_labels) <- colnames(epi_clean)
  new_labels[confirmed_cells] <- NEW_LABEL
  epi_clean[[group_col]] <- unname(new_labels)
  
  message("\n--- Updated ", group_col, " counts ---")
  print(table(epi_clean[[group_col]][, 1]))
  
  # saveRDS(epi_clean, "epi_clean_FINAL_v6.rds")
}

## =========================================================
## Option B: true freeform/lasso selection (not just a rectangle), via the
## `gatepoints` package -- closer to actually "drawing a circle."
## Install once if needed: install.packages("gatepoints")
## =========================================================

# library(gatepoints)
# umap_df <- as.data.frame(Embeddings(epi_clean, "umap"))
# colnames(umap_df) <- c("umap_1", "umap_2")
# plot(umap_df$umap_1, umap_df$umap_2, pch = 20, cex = 0.3,
#      col = factor(epi_clean[[group_col]][, 1]))
# selected_idx <- fhs(umap_df, mark = TRUE)   ## click points to draw a freehand polygon,
#                                               ## double-click/right-click to close it
# selected_cells_lasso <- rownames(umap_df)[selected_idx]