##convert rat, mouse and marmoset genes in human symbol to facilitate the integration
library(Seurat)
library(SeuratObject)
library(Matrix)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Rn.eg.db)

# --- sparse-safe row aggregation: sums rows with same group label ---
aggregate_rows_sparse <- function(mat, groups) {
  stopifnot(length(groups) == nrow(mat))
  f <- factor(groups)
  # Incidence matrix: (n_groups x n_rows) * (n_rows x n_cells) = (n_groups x n_cells)
  S <- Matrix::sparseMatrix(
    i = as.integer(f),
    j = seq_len(nrow(mat)),
    x = 1,
    dims = c(nlevels(f), nrow(mat))
  )
  out <- S %*% mat
  rownames(out) <- levels(f)
  out
}

convert_mouse_to_human_using_mouse_genes <- function(seu, mouse_genes,
                                                     assay = "RNA", layer = "counts",
                                                     mouse_ens_col = "Gene.stable.ID",
                                                     human_sym_col = "Symbol") {
  message("Mouse symbols → ENSMUSG (org.Mm.eg.db) → human Symbol (mouse_genes table)")
  
  DefaultAssay(seu) <- assay
  mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
  
  # 1) SYMBOL -> ENSMUSG (local)
  sym <- rownames(mat)
  map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = sym,
    keytype = "SYMBOL",
    columns = "ENSEMBL"
  )
  map <- map[!is.na(map$ENSEMBL) & nzchar(map$ENSEMBL), ]
  # Resolve 1:many by keeping the first ENSEMBL per SYMBOL
  map <- map[!duplicated(map$SYMBOL), c("SYMBOL", "ENSEMBL")]
  
  keep_sym <- intersect(sym, map$SYMBOL)
  if (length(keep_sym) == 0) stop("No SYMBOL→ENSEMBL matches found for mouse genes.")
  
  mat <- mat[keep_sym, , drop = FALSE]
  sym2ens <- setNames(map$ENSEMBL, map$SYMBOL)[rownames(mat)]
  rownames(mat) <- sym2ens
  
  # 2) ENSMUSG -> Human Symbol (YOUR mouse_genes table)
  idx <- match(rownames(mat), mouse_genes[[mouse_ens_col]])
  hum <- mouse_genes[[human_sym_col]][idx]
  
  keep <- !is.na(idx) & !is.na(hum) & nzchar(hum)
  if (!any(keep)) {
    stop("No overlap between ENSMUSG from Seurat and mouse_genes$", mouse_ens_col)
  }
  
  mat2 <- mat[keep, , drop = FALSE]
  hum2 <- hum[keep]
  rownames(mat2) <- hum2
  
  # 3) Aggregate duplicated human symbols (sparse-safe)
  mat2 <- aggregate_rows_sparse(mat2, rownames(mat2))
  
  seu_new <- Seurat::CreateSeuratObject(counts = mat2, meta.data = seu@meta.data)
  seu_new$species_original <- "mouse"
  return(seu_new)
}

# Run
seu_mouse_h <- convert_mouse_to_human_using_mouse_genes(seu_mouse, mouse_genes)

#Rat gene symbols for some reason are mixed up so the convertion will be a lil different

library(Seurat)
library(SeuratObject)
library(Matrix)

# sparse-safe row aggregation (sum rows with same label)
aggregate_rows_sparse <- function(mat, groups) {
  f <- factor(groups)
  S <- Matrix::sparseMatrix(
    i = as.integer(f),
    j = seq_len(nrow(mat)),
    x = 1,
    dims = c(nlevels(f), nrow(mat))
  )
  out <- S %*% mat
  rownames(out) <- levels(f)
  out
}

convert_rat_to_human_using_table_only <- function(seu, rat_genes, species_name = "rat",
                                                  assay = "RNA", layer = "counts",
                                                  rat_ens_col = "Gene.stable.ID",
                                                  rat_name_col = "Gene.name",
                                                  human_sym_col = "Symbol") {
  message("Converting ", species_name, " → human using rat_genes table only (handles SYMBOL + ENSRNOG)...")
  
  DefaultAssay(seu) <- assay
  mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
  
  g <- rownames(mat)
  
  # Build fast lookup maps from the table
  # (dedupe keys so we get a single mapping per key)
  tab_ens  <- rat_genes[!is.na(rat_genes[[rat_ens_col]])  & nzchar(rat_genes[[rat_ens_col]]), ]
  tab_name <- rat_genes[!is.na(rat_genes[[rat_name_col]]) & nzchar(rat_genes[[rat_name_col]]), ]
  
  tab_ens  <- tab_ens[!duplicated(tab_ens[[rat_ens_col]]), ]
  tab_name <- tab_name[!duplicated(tab_name[[rat_name_col]]), ]
  
  ens2human  <- setNames(tab_ens[[human_sym_col]],  tab_ens[[rat_ens_col]])
  name2human <- setNames(tab_name[[human_sym_col]], tab_name[[rat_name_col]])
  
  # First try Ensembl IDs
  hum_from_ens <- unname(ens2human[g])
  # For those that failed, try gene names
  missing <- is.na(hum_from_ens) | !nzchar(hum_from_ens)
  hum_from_name <- rep(NA_character_, length(g))
  hum_from_name[missing] <- unname(name2human[g[missing]])
  
  hum <- hum_from_ens
  hum[missing] <- hum_from_name[missing]
  
  keep <- !is.na(hum) & nzchar(hum)
  if (!any(keep)) {
    stop("No features could be mapped to human symbols using rat_genes table.\n",
         "Example Seurat genes: ", paste(head(g), collapse = ", "))
  }
  
  mat2 <- mat[keep, , drop = FALSE]
  rownames(mat2) <- hum[keep]
  
  # Sum duplicates (many-to-one orthology, etc.)
  mat2 <- aggregate_rows_sparse(mat2, rownames(mat2))
  
  seu_new <- Seurat::CreateSeuratObject(counts = mat2, meta.data = seu@meta.data)
  seu_new$species_original <- species_name
  return(seu_new)
}

# Run
seu_rat_h <- convert_rat_to_human_using_table_only(seu_rat, rat_genes, layer = "counts")

# Sanity check

dim(seu_rat[["RNA"]])
dim(seu_rat_h[["RNA"]])
head(rownames(seu_rat_h[["RNA"]]))
any(duplicated(rownames(seu_rat_h[["RNA"]])))  # FALSE expected


# Now for marmoset!
# Install if needed
## ============================================================
## NCBI: Build Callithrix jacchus (taxid 9483) -> Human (9606)
##       ortholog symbol map using gene_orthologs.gz + gene_info.gz
## ============================================================
download.file("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
              destfile = "gene_orthologs.gz", mode = "wb")

library(data.table)

go <- fread("gene_orthologs.gz", sep = "\t", header = TRUE)
setnames(go, names(go),
         c("tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"))

tax_marmoset <- 9483
tax_human    <- 9606

pairs <- rbind(
  go[tax_id == tax_marmoset & Other_tax_id == tax_human,
     .(marmoset_geneid = GeneID, human_geneid = Other_GeneID)],
  go[tax_id == tax_human & Other_tax_id == tax_marmoset,
     .(marmoset_geneid = Other_GeneID, human_geneid = GeneID)]
)

pairs <- unique(pairs)

library(httr)
library(jsonlite)

fetch_symbols <- function(gene_ids) {
  ids <- paste(gene_ids, collapse = ",")
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
    "db=gene&id=", ids, "&retmode=json"
  )
  res <- GET(url)
  js  <- fromJSON(content(res, "text", encoding = "UTF-8"))
  out <- js$result
  data.frame(
    GeneID = names(out),
    Symbol = sapply(out, function(x) x$name),
    stringsAsFactors = FALSE
  )
}

# fetch only needed IDs
marmoset_syms <- fetch_symbols(unique(pairs$marmoset_geneid))
human_syms    <- fetch_symbols(unique(pairs$human_geneid))


#Re-built seurat marmoset with human gene symbol
# pairs has integer IDs; marmoset_syms/human_syms have character GeneID
pairs_dt <- as.data.table(pairs)
pairs_dt[, marmoset_geneid := as.character(marmoset_geneid)]
pairs_dt[, human_geneid    := as.character(human_geneid)]

marmoset_syms_dt <- as.data.table(marmoset_syms)
human_syms_dt    <- as.data.table(human_syms)

setnames(marmoset_syms_dt, c("GeneID","Symbol"), c("marmoset_geneid","marmoset_symbol"))
setnames(human_syms_dt,    c("GeneID","Symbol"), c("human_geneid","human_symbol"))

orth <- merge(pairs_dt, marmoset_syms_dt, by = "marmoset_geneid", all.x = TRUE)
orth <- merge(orth, human_syms_dt,    by = "human_geneid",    all.x = TRUE)

marmoset_to_human <- orth[!is.na(marmoset_symbol) & nzchar(marmoset_symbol) &
                            !is.na(human_symbol) & nzchar(human_symbol),
                          .(gene_from = marmoset_symbol, Symbol = human_symbol)]

# optional: keep one human symbol per marmoset symbol
marmoset_to_human_1to1 <- marmoset_to_human[!duplicated(gene_from)]

dim(marmoset_to_human)
head(marmoset_to_human)

library(Seurat)
library(SeuratObject)
library(Matrix)

# sparse-safe row aggregation (sum rows with same label)
aggregate_rows_sparse <- function(mat, groups) {
  f <- factor(groups)
  S <- Matrix::sparseMatrix(
    i = as.integer(f),
    j = seq_len(nrow(mat)),
    x = 1,
    dims = c(nlevels(f), nrow(mat))
  )
  out <- S %*% mat
  rownames(out) <- levels(f)
  out
}

convert_seurat_symbols_to_human <- function(seu, map_table,
                                            species_name = "marmoset",
                                            assay = "RNA", layer = "counts",
                                            from_col = "gene_from",
                                            to_col   = "Symbol") {
  message("Converting ", species_name, " → human symbols using mapping table...")
  
  DefaultAssay(seu) <- assay
  mat <- SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
  
  # match Seurat gene symbols to mapping table
  idx <- match(rownames(mat), map_table[[from_col]])
  hum <- map_table[[to_col]][idx]
  
  keep <- !is.na(idx) & !is.na(hum) & nzchar(hum)
  if (!any(keep)) stop("No overlap between Seurat genes and mapping table for ", species_name)
  
  mat2 <- mat[keep, , drop = FALSE]
  rownames(mat2) <- hum[keep]
  
  # sum duplicates after mapping
  mat2 <- aggregate_rows_sparse(mat2, rownames(mat2))
  
  seu_new <- Seurat::CreateSeuratObject(counts = mat2, meta.data = seu@meta.data)
  seu_new$species_original <- species_name
  return(seu_new)
}

# Run for marmoset
seu_marmoset_h <- convert_seurat_symbols_to_human(seu_marmoset, marmoset_to_human, "marmoset")

dim(seu_marmoset[["RNA"]])
dim(seu_marmoset_h[["RNA"]])
head(rownames(seu_marmoset_h[["RNA"]]))
any(duplicated(rownames(seu_marmoset_h[["RNA"]])))  # should be FALSE

data.table::fwrite(marmoset_to_human, "marmoset_to_human_symbol.tsv", sep = "\t")
