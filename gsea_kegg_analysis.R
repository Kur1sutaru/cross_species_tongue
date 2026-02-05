#!/usr/bin/env Rscript

# ============================================================================
# GSEA and KEGG Pathway Analysis for Conserved Markers
# Analysis across all cell subtypes from 4-species comparison
# ============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(openxlsx)

# Set working directory
setwd("/home/claude")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading conserved markers data...\n")
# Read the RDS file - update this path with your actual filename
conserved_markers <- readRDS("/mnt/user-data/uploads/conserved_markers_all_4_species.rds")

cat("Subtypes found:", length(names(conserved_markers)), "\n")
cat("Subtype names:\n")
print(names(conserved_markers))

# ============================================================================
# 2. PREPARE GENE LISTS FOR EACH SUBTYPE
# ============================================================================

prepare_gene_list <- function(marker_data, species = "human") {
  """
  Prepare ranked gene list for GSEA
  Uses avg_log2FC from specified species (default: human)
  """
  
  # Select the appropriate log2FC column based on species
  log2fc_col <- paste0(species, "_avg_log2FC")
  pval_col <- paste0(species, "_p_val")
  
  if (!log2fc_col %in% colnames(marker_data)) {
    cat("Warning: Column", log2fc_col, "not found. Using human_avg_log2FC\n")
    log2fc_col <- "human_avg_log2FC"
  }
  
  # Create ranked gene list
  gene_list <- marker_data[[log2fc_col]]
  names(gene_list) <- rownames(marker_data)
  
  # Remove NA values
  gene_list <- gene_list[!is.na(gene_list)]
  
  # Sort by log2FC (descending)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# ============================================================================
# 3. CONVERT GENE SYMBOLS TO ENTREZ IDs
# ============================================================================

convert_to_entrez <- function(gene_list) {
  """
  Convert gene symbols to Entrez IDs for KEGG analysis
  """
  
  # Convert gene symbols to Entrez IDs
  gene_symbols <- names(gene_list)
  
  tryCatch({
    gene_list_entrez <- bitr(gene_symbols, 
                             fromType = "SYMBOL",
                             toType = "ENTREZID", 
                             OrgDb = org.Hs.eg.db)
    
    # Match and create ranked list with Entrez IDs
    gene_list_matched <- gene_list[gene_list_entrez$SYMBOL]
    names(gene_list_matched) <- gene_list_entrez$ENTREZID
    
    # Remove duplicates (keep highest fold change)
    gene_list_matched <- gene_list_matched[!duplicated(names(gene_list_matched))]
    
    cat("  - Converted", length(gene_list_matched), "genes to Entrez IDs\n")
    
    return(gene_list_matched)
    
  }, error = function(e) {
    cat("  Error in gene conversion:", e$message, "\n")
    return(NULL)
  })
}

# ============================================================================
# 4. PERFORM ENRICHMENT ANALYSES
# ============================================================================

perform_enrichment <- function(marker_data, subtype_name, species = "human") {
  """
  Perform comprehensive enrichment analysis for one subtype
  Includes: GSEA-GO, GSEA-KEGG, ORA-GO, ORA-KEGG
  """
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Processing:", subtype_name, "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Prepare gene list
  gene_list <- prepare_gene_list(marker_data, species)
  cat("  - Total genes:", length(gene_list), "\n")
  
  # Convert to Entrez IDs
  gene_list_entrez <- convert_to_entrez(gene_list)
  
  if (is.null(gene_list_entrez) || length(gene_list_entrez) < 10) {
    cat("  ! Insufficient genes for analysis\n")
    return(NULL)
  }
  
  results <- list()
  results$gene_list <- gene_list_entrez
  
  # -------------------------------------------------------------------------
  # GSEA - Gene Ontology (Biological Process, Molecular Function, Cellular Component)
  # -------------------------------------------------------------------------
  cat("\n  Running GSEA - Gene Ontology...\n")
  
  tryCatch({
    results$gsea_go <- gseGO(
      geneList = gene_list_entrez,
      OrgDb = org.Hs.eg.db,
      ont = "ALL",  # All ontologies
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    if (!is.null(results$gsea_go) && nrow(results$gsea_go) > 0) {
      cat("    ✓ Found", nrow(results$gsea_go), "significant GO terms\n")
    } else {
      cat("    - No significant GO terms found\n")
    }
  }, error = function(e) {
    cat("    ! Error in GSEA-GO:", e$message, "\n")
    results$gsea_go <- NULL
  })
  
  # -------------------------------------------------------------------------
  # GSEA - KEGG Pathways
  # -------------------------------------------------------------------------
  cat("  Running GSEA - KEGG...\n")
  
  tryCatch({
    results$gsea_kegg <- gseKEGG(
      geneList = gene_list_entrez,
      organism = 'hsa',  # Homo sapiens
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    if (!is.null(results$gsea_kegg) && nrow(results$gsea_kegg) > 0) {
      cat("    ✓ Found", nrow(results$gsea_kegg), "significant KEGG pathways\n")
    } else {
      cat("    - No significant KEGG pathways found\n")
    }
  }, error = function(e) {
    cat("    ! Error in GSEA-KEGG:", e$message, "\n")
    results$gsea_kegg <- NULL
  })
  
  # -------------------------------------------------------------------------
  # Over-Representation Analysis - GO
  # -------------------------------------------------------------------------
  cat("  Running ORA - Gene Ontology...\n")
  
  # Get top upregulated genes for ORA
  top_genes_entrez <- names(gene_list_entrez)[gene_list_entrez > 0]
  top_genes_entrez <- head(top_genes_entrez, min(500, length(top_genes_entrez)))
  
  if (length(top_genes_entrez) >= 10) {
    tryCatch({
      results$ora_go <- enrichGO(
        gene = top_genes_entrez,
        OrgDb = org.Hs.eg.db,
        ont = "ALL",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        readable = TRUE
      )
      
      if (!is.null(results$ora_go) && nrow(results$ora_go) > 0) {
        cat("    ✓ Found", nrow(results$ora_go), "significant GO terms\n")
      } else {
        cat("    - No significant GO terms found\n")
      }
    }, error = function(e) {
      cat("    ! Error in ORA-GO:", e$message, "\n")
      results$ora_go <- NULL
    })
  }
  
  # -------------------------------------------------------------------------
  # Over-Representation Analysis - KEGG
  # -------------------------------------------------------------------------
  cat("  Running ORA - KEGG...\n")
  
  if (length(top_genes_entrez) >= 10) {
    tryCatch({
      results$ora_kegg <- enrichKEGG(
        gene = top_genes_entrez,
        organism = 'hsa',
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      )
      
      if (!is.null(results$ora_kegg) && nrow(results$ora_kegg) > 0) {
        cat("    ✓ Found", nrow(results$ora_kegg), "significant KEGG pathways\n")
      } else {
        cat("    - No significant KEGG pathways found\n")
      }
    }, error = function(e) {
      cat("    ! Error in ORA-KEGG:", e$message, "\n")
      results$ora_kegg <- NULL
    })
  }
  
  return(results)
}

# ============================================================================
# 5. RUN ANALYSIS FOR ALL SUBTYPES
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("STARTING ENRICHMENT ANALYSIS FOR ALL SUBTYPES\n")
cat(rep("=", 70), "\n", sep = "")

enrichment_results <- list()

for (subtype in names(conserved_markers)) {
  enrichment_results[[subtype]] <- perform_enrichment(
    conserved_markers[[subtype]], 
    subtype,
    species = "human"  # Change to "mouse", "rat", or "marm" if needed
  )
}

# ============================================================================
# 6. CREATE VISUALIZATIONS
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("CREATING VISUALIZATIONS\n")
cat(rep("=", 70), "\n", sep = "")

dir.create("enrichment_plots", showWarnings = FALSE)

for (subtype in names(enrichment_results)) {
  
  if (is.null(enrichment_results[[subtype]])) {
    next
  }
  
  cat("\nCreating plots for:", subtype, "\n")
  
  # Clean subtype name for filenames
  subtype_clean <- gsub("[^A-Za-z0-9_]", "_", subtype)
  
  # -------------------------------------------------------------------------
  # GSEA GO Plots
  # -------------------------------------------------------------------------
  if (!is.null(enrichment_results[[subtype]]$gsea_go) && 
      nrow(enrichment_results[[subtype]]$gsea_go) > 0) {
    
    tryCatch({
      pdf(paste0("enrichment_plots/", subtype_clean, "_GSEA_GO.pdf"), 
          width = 14, height = 10)
      
      # Dotplot
      p1 <- dotplot(enrichment_results[[subtype]]$gsea_go, 
                    showCategory = 30,
                    split = "ONTOLOGY",
                    font.size = 10) +
        ggtitle(paste(subtype, "- GSEA GO Terms"))
      print(p1)
      
      # Ridge plot (top 20)
      if (nrow(enrichment_results[[subtype]]$gsea_go) >= 5) {
        p2 <- ridgeplot(enrichment_results[[subtype]]$gsea_go, 
                        showCategory = 20) + 
          ggtitle(paste(subtype, "- GSEA GO Ridge Plot"))
        print(p2)
      }
      
      dev.off()
      cat("  ✓ GSEA GO plots saved\n")
    }, error = function(e) {
      cat("  ! Error creating GSEA GO plots:", e$message, "\n")
      dev.off()
    })
  }
  
  # -------------------------------------------------------------------------
  # GSEA KEGG Plots
  # -------------------------------------------------------------------------
  if (!is.null(enrichment_results[[subtype]]$gsea_kegg) && 
      nrow(enrichment_results[[subtype]]$gsea_kegg) > 0) {
    
    tryCatch({
      pdf(paste0("enrichment_plots/", subtype_clean, "_GSEA_KEGG.pdf"), 
          width = 14, height = 10)
      
      # Dotplot
      p1 <- dotplot(enrichment_results[[subtype]]$gsea_kegg, 
                    showCategory = 25,
                    font.size = 10) +
        ggtitle(paste(subtype, "- GSEA KEGG Pathways"))
      print(p1)
      
      # Ridge plot
      if (nrow(enrichment_results[[subtype]]$gsea_kegg) >= 5) {
        p2 <- ridgeplot(enrichment_results[[subtype]]$gsea_kegg, 
                        showCategory = 20) + 
          ggtitle(paste(subtype, "- GSEA KEGG Ridge Plot"))
        print(p2)
      }
      
      # Network plot (if enough pathways)
      if (nrow(enrichment_results[[subtype]]$gsea_kegg) >= 5) {
        p3 <- tryCatch({
          cnetplot(enrichment_results[[subtype]]$gsea_kegg,
                   categorySize = "pvalue",
                   showCategory = 10,
                   foldChange = enrichment_results[[subtype]]$gene_list,
                   colorEdge = TRUE) +
            ggtitle(paste(subtype, "- KEGG Pathway Network"))
        }, error = function(e) NULL)
        
        if (!is.null(p3)) print(p3)
      }
      
      dev.off()
      cat("  ✓ GSEA KEGG plots saved\n")
    }, error = function(e) {
      cat("  ! Error creating GSEA KEGG plots:", e$message, "\n")
      dev.off()
    })
  }
  
  # -------------------------------------------------------------------------
  # ORA GO Plots
  # -------------------------------------------------------------------------
  if (!is.null(enrichment_results[[subtype]]$ora_go) && 
      nrow(enrichment_results[[subtype]]$ora_go) > 0) {
    
    tryCatch({
      pdf(paste0("enrichment_plots/", subtype_clean, "_ORA_GO.pdf"), 
          width = 12, height = 10)
      
      # Dotplot
      p1 <- dotplot(enrichment_results[[subtype]]$ora_go, 
                    showCategory = 30,
                    split = "ONTOLOGY",
                    font.size = 10) +
        ggtitle(paste(subtype, "- ORA GO Terms"))
      print(p1)
      
      # Barplot
      p2 <- barplot(enrichment_results[[subtype]]$ora_go, 
                    showCategory = 20) +
        ggtitle(paste(subtype, "- ORA GO Barplot"))
      print(p2)
      
      dev.off()
      cat("  ✓ ORA GO plots saved\n")
    }, error = function(e) {
      cat("  ! Error creating ORA GO plots:", e$message, "\n")
      dev.off()
    })
  }
  
  # -------------------------------------------------------------------------
  # ORA KEGG Plots
  # -------------------------------------------------------------------------
  if (!is.null(enrichment_results[[subtype]]$ora_kegg) && 
      nrow(enrichment_results[[subtype]]$ora_kegg) > 0) {
    
    tryCatch({
      pdf(paste0("enrichment_plots/", subtype_clean, "_ORA_KEGG.pdf"), 
          width = 12, height = 10)
      
      # Dotplot
      p1 <- dotplot(enrichment_results[[subtype]]$ora_kegg, 
                    showCategory = 25,
                    font.size = 10) +
        ggtitle(paste(subtype, "- ORA KEGG Pathways"))
      print(p1)
      
      # Barplot
      p2 <- barplot(enrichment_results[[subtype]]$ora_kegg, 
                    showCategory = 20) +
        ggtitle(paste(subtype, "- ORA KEGG Barplot"))
      print(p2)
      
      dev.off()
      cat("  ✓ ORA KEGG plots saved\n")
    }, error = function(e) {
      cat("  ! Error creating ORA KEGG plots:", e$message, "\n")
      dev.off()
    })
  }
}

# ============================================================================
# 7. EXPORT RESULTS TO EXCEL
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("EXPORTING RESULTS TO EXCEL\n")
cat(rep("=", 70), "\n", sep = "")

wb <- createWorkbook()

for (subtype in names(enrichment_results)) {
  
  if (is.null(enrichment_results[[subtype]])) {
    next
  }
  
  # Clean subtype name for sheet names (Excel limits)
  sheet_name_base <- substr(gsub("[^A-Za-z0-9]", "_", subtype), 1, 20)
  
  # GSEA GO results
  if (!is.null(enrichment_results[[subtype]]$gsea_go) && 
      nrow(enrichment_results[[subtype]]$gsea_go) > 0) {
    
    sheet_name <- paste0(sheet_name_base, "_GSEA_GO")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, 
              as.data.frame(enrichment_results[[subtype]]$gsea_go))
    cat("  Added sheet:", sheet_name, "\n")
  }
  
  # GSEA KEGG results
  if (!is.null(enrichment_results[[subtype]]$gsea_kegg) && 
      nrow(enrichment_results[[subtype]]$gsea_kegg) > 0) {
    
    sheet_name <- paste0(sheet_name_base, "_GSEA_KEGG")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, 
              as.data.frame(enrichment_results[[subtype]]$gsea_kegg))
    cat("  Added sheet:", sheet_name, "\n")
  }
  
  # ORA GO results
  if (!is.null(enrichment_results[[subtype]]$ora_go) && 
      nrow(enrichment_results[[subtype]]$ora_go) > 0) {
    
    sheet_name <- paste0(sheet_name_base, "_ORA_GO")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, 
              as.data.frame(enrichment_results[[subtype]]$ora_go))
    cat("  Added sheet:", sheet_name, "\n")
  }
  
  # ORA KEGG results
  if (!is.null(enrichment_results[[subtype]]$ora_kegg) && 
      nrow(enrichment_results[[subtype]]$ora_kegg) > 0) {
    
    sheet_name <- paste0(sheet_name_base, "_ORA_KEGG")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, 
              as.data.frame(enrichment_results[[subtype]]$ora_kegg))
    cat("  Added sheet:", sheet_name, "\n")
  }
}

saveWorkbook(wb, "enrichment_results_all_subtypes.xlsx", overwrite = TRUE)
cat("\n✓ Excel file saved: enrichment_results_all_subtypes.xlsx\n")

# ============================================================================
# 8. SAVE RDS FOR LATER USE
# ============================================================================

saveRDS(enrichment_results, "enrichment_results.rds")
cat("✓ RDS file saved: enrichment_results.rds\n")

# ============================================================================
# 9. CREATE SUMMARY REPORT
# ============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY REPORT\n")
cat(rep("=", 70), "\n", sep = "")

summary_df <- data.frame(
  Subtype = character(),
  Total_Genes = integer(),
  GSEA_GO_Terms = integer(),
  GSEA_KEGG_Pathways = integer(),
  ORA_GO_Terms = integer(),
  ORA_KEGG_Pathways = integer(),
  stringsAsFactors = FALSE
)

for (subtype in names(enrichment_results)) {
  if (is.null(enrichment_results[[subtype]])) {
    next
  }
  
  summary_df <- rbind(summary_df, data.frame(
    Subtype = subtype,
    Total_Genes = length(enrichment_results[[subtype]]$gene_list),
    GSEA_GO_Terms = ifelse(is.null(enrichment_results[[subtype]]$gsea_go), 0, 
                           nrow(enrichment_results[[subtype]]$gsea_go)),
    GSEA_KEGG_Pathways = ifelse(is.null(enrichment_results[[subtype]]$gsea_kegg), 0,
                                nrow(enrichment_results[[subtype]]$gsea_kegg)),
    ORA_GO_Terms = ifelse(is.null(enrichment_results[[subtype]]$ora_go), 0,
                         nrow(enrichment_results[[subtype]]$ora_go)),
    ORA_KEGG_Pathways = ifelse(is.null(enrichment_results[[subtype]]$ora_kegg), 0,
                               nrow(enrichment_results[[subtype]]$ora_kegg))
  ))
}

print(summary_df)

write.csv(summary_df, "enrichment_summary.csv", row.names = FALSE)
cat("\n✓ Summary saved: enrichment_summary.csv\n")

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nOutput files:\n")
cat("  1. enrichment_results_all_subtypes.xlsx - All enrichment results\n")
cat("  2. enrichment_results.rds - R object for further analysis\n")
cat("  3. enrichment_summary.csv - Summary table\n")
cat("  4. enrichment_plots/ - Directory with all plots (PDFs)\n")
cat("\n")
