# ==============================================================================
# Script: 03_functional_go.R
# Description: Downstream GO Biological Process (BP) enrichment analysis and 
#              visualization for individual Mfuzz temporal clusters.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 1. Define Directories
input_file <- "./results/time_proteomics_output/8_Cluster_output/8_cluster_merged_data_with_annotations.tsv"
output_go_dir <- "./results/time_proteomics_output/cluster_GO"

dir.create(output_go_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Required cluster output file missing: ", input_file)
}

# 2. Load Clustered Dataset
cluster_data <- read_tsv(input_file)
all_proteins <- unique(cluster_data$Protein.Ids)

# 3. Global ID Conversion: Uniprot to EntrezID
gene_conversion <- bitr(
  all_proteins,
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
cat("ID conversion completed. Valid mapped identifiers:", nrow(gene_conversion), "\n")

# 4. Graphical Parameters and Custom Color Palette
color_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

custom_theme <- theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
  )

# 5. Batch Run GO Enrichment Analysis
cluster_num <- max(cluster_data$Cluster)

for (i in 1:cluster_num) {
  cat("Processing GO enrichment for Cluster", i, "...\n")
  
  # Extract protein set for current cluster
  cluster_proteins <- cluster_data %>%
    filter(Cluster == i) %>%
    pull(Protein.Ids)
  
  # Map to mapped Entrez IDs
  cluster_entrez <- gene_conversion %>%
    filter(UNIPROT %in% cluster_proteins) %>%
    pull(ENTREZID) %>%
    unique()
  
  # Quality filter on gene count minimum thresholds
  if (length(cluster_entrez) < 5) {
    cat("Skipped Cluster", i, ": Insufficient gene identifiers.\n")
    next
  }
  
  # GO Biological Process (BP) Analysis
  ego <- enrichGO(
    gene = cluster_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  if (nrow(as.data.frame(ego)) > 0) {
    # Export GO Table
    write_tsv(as.data.frame(ego), file.path(output_go_dir, paste0("Cluster_", i, "_GO_enrichment.tsv")))
    
    # Plot and save GO figures
    p_go <- barplot(ego, showCategory = 20, colorBy = "p.adjust", 
                    palette = color_palette, title = paste("GO Enrichment for Cluster", i)) + custom_theme
    
    pdf(file.path(output_go_dir, paste0("Cluster_", i, "_GO_plot.pdf")), width = 10, height = 8)
    print(p_go)
    dev.off()
    
    png(file.path(output_go_dir, paste0("Cluster_", i, "_GO_plot.png")), width = 800, height = 600)
    print(p_go)
    dev.off()
    
    cat("Successfully exported GO results for Cluster", i, "\n")
  } else {
    cat("No significant GO terms identified for Cluster", i, "\n")
  }
}

cat("GO enrichment pipeline execution completed.\n")