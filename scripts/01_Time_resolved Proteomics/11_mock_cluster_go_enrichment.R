# ==============================================================================
# Script: 11_mock_cluster_go_enrichment.R
# Description: Dynamic GO (BP) enrichment for MOCK timeline clusters with index 
#              protection, path independence, and fixed coordinate axis mapping.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tibble)
library(dplyr)

# 1. IO Matrix Topography Setup
cluster_data_dir <- "./results/Mock_output/5_Cluster_output"
go_output_dir    <- file.path(cluster_data_dir, "5_cluster_GO")
go_plot_dir      <- file.path(go_output_dir, "5_cluster_GO_plots")

dir.create(go_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(go_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Ingest Cluster Profile Matrix
master_table_file <- file.path(cluster_data_dir, "5_cluster_merged_data_ST_annotations.tsv")
if (!file.exists(master_table_file)) {
  stop("Prerequisite cluster alignment file missing: ", master_table_file)
}
merge_data <- read_tsv(master_table_file, show_col_types = FALSE)
all_genes  <- unique(merge_data$Protein.Ids)

# 2. High-Throughput Accession Mapping (UNIPROT -> ENTREZID)
gene_conversion <- bitr(
  all_genes, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db
)
cat("✅ ID Translation complete. Registered common database intersections:", nrow(gene_conversion), "rows.\n")

# Storage infrastructure initialized for multi-step pipe handover
go_results_list <- list()
cluster_num     <- max(merge_data$Cluster)

# 3. Layer A: Automated Hypergeometric Distribution Enrichment Testing
for (i in 1:cluster_num) {
  cat("🚀 Running biological network enrichment on Cluster", i, "...\n")
  
  cluster_genes <- merge_data %>% filter(Cluster == i) %>% pull(Protein.Ids)
  
  cluster_genes_entrez <- gene_conversion %>%
    filter(UNIPROT %in% cluster_genes) %>%
    pull(ENTREZID) %>% unique()
  
  if (length(cluster_genes_entrez) < 5) {
    cat("⚠️ Notice: Insufficient input dimensions for Cluster", i, " (N < 5). Skipping.\n")
    next
  }
  
  # Execute targeted enrichment strategy (Biological Process focused)
  ego <- enrichGO(
    gene = cluster_genes_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("⚠️ Notice: Zero statically significant enrichments returned for Cluster", i, "\n")
    next
  }
  
  # Lock results into reference list context
  go_results_list[[paste0("Cluster_", i)]] <- ego
  
  # Write tabular data backup
  write.table(
    as.data.frame(ego),
    file = file.path(go_output_dir, paste0("Cluster_", i, "_GO_enrichment.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
cat("🎉 Tabular enrichment matrices saved cleanly.\n")

# 4. Layer B: Publication-Grade Visual Layout Rendering (Fixed Axis Orientation)
cat("Initiating high-resolution graphic generation...\n")

for (i in 1:cluster_num) {
  cluster_key <- paste0("Cluster_", i)
  
  # CRITICAL SECURITY FIX: Guarding against index breaches caused by NULL outputs
  if (!cluster_key %in% names(go_results_list)) {
    cat("⚠️ Skipping visualization for empty/non-significant node:", cluster_key, "\n")
    next
  }
  
  ego_object <- go_results_list[[cluster_key]]
  
  # Corrected layout map: terms occupy Y-axis natively, counts map to X-axis
  p <- barplot(
    ego_object, showCategory = 20, colorBy = "p.adjust",
    title = paste("GO Enrichment: Cluster", i)
  ) + 
    labs(x = "Gene Count", y = "Biological Process Terms") +
    theme_bw() + 
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      legend.title = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey92")
    )
  
  # File save blocks
  pdf_filename <- file.path(go_plot_dir, paste0("Cluster_", i, "_GO_enrichment_plot.pdf"))
  png_filename <- file.path(go_plot_dir, paste0("Cluster_", i, "_GO_enrichment_plot.png"))
  
  pdf(pdf_filename, width = 9, height = 7)
  print(p)
  dev.off()
  
  png(png_filename, width = 2700, height = 2100, res = 300)
  print(p)
  dev.off()
  
  cat("✅ Visualization successfully minted for Cluster", i, "\n")
}

cat("🎉 Pipeline 12 executed flawlessly. All viable outputs are located under your project hierarchy.\n")