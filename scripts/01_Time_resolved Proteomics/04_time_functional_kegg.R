# ==============================================================================
# Script: 04_functional_kegg.R
# Description: Downstream KEGG pathway enrichment analysis and visualization 
#              for individual Mfuzz temporal clusters.
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
output_kegg_dir <- "./results/time_proteomics_output/cluster_KEGG"

dir.create(output_kegg_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Required cluster output file missing: ", input_file)
}

# 2. Load Clustered Dataset
cluster_data <- read_tsv(input_file)
all_proteins <- unique(cluster_data$Protein.Ids)

# 3. Global ID Conversion: Uniprot to EntrezID (Required for KEGG hsa organism matching)
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

# 5. Batch Run KEGG Enrichment Analysis
cluster_num <- max(cluster_data$Cluster)

for (i in 1:cluster_num) {
  cat("Processing KEGG enrichment for Cluster", i, "...\n")
  
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
  
  # KEGG Pathway Enrichment Analysis
  ekegg <- enrichKEGG(
    gene = cluster_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (nrow(as.data.frame(ekegg)) > 0) {
    # Export KEGG Table
    write_tsv(as.data.frame(ekegg), file.path(output_kegg_dir, paste0("Cluster_", i, "_KEGG_enrichment.tsv")))
    
    # Plot and save KEGG figures
    p_kegg <- barplot(ekegg, showCategory = 20, colorBy = "p.adjust", 
                      palette = color_palette, title = paste("KEGG Enrichment for Cluster", i)) + custom_theme
    
    pdf(file.path(output_kegg_dir, paste0("Cluster_", i, "_KEGG_plot.pdf")), width = 10, height = 8)
    print(p_kegg)
    dev.off()
    
    png(file.path(output_kegg_dir, paste0("Cluster_", i, "_KEGG_plot.png")), width = 800, height = 600)
    print(p_kegg)
    dev.off()
    
    cat("Successfully exported KEGG results for Cluster", i, "\n")
  } else {
    cat("No significant KEGG pathways identified for Cluster", i, "\n")
  }
}

cat("KEGG enrichment pipeline execution completed.\n")