# ==============================================================================
# Script: 05_dep_kegg_down_plot.R
# Description: Batch generation of publication-quality barplots for significantly 
#              downregulated KEGG pathways from DEP analysis.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readxl)
library(ggplot2)
library(dplyr)

# 1. Define Input and Output Directories (Relative Path Paradigm)
input_dir <- "./data/Time_DEP_GO/KEGG_DOWN"
output_dir <- "./results/Time_DEP_Enrichment_plot/KEGGDOWN_plot"

if (!dir.exists(input_dir)) {
  stop("Input directory missing: ", input_dir)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Identify All Target Excel Workbooks
files <- list.files(input_dir, pattern = "\\.xlsx$", full.names = TRUE)

if (length(files) == 0) {
  cat("No Excel files found in the specified directory.\n")
}

# 3. Iterative Visual Generation Pipeline
for (file in files) {
  # Load current dataset
  kegg_df <- read_excel(file)
  
  # Structural validation check
  if (!all(c("pvalue", "Description") %in% colnames(kegg_df))) {
    cat("Skipped non-compliant file structure:", basename(file), "\n")
    next
  }
  
  # Filter, rank, and calculate negative log p-values for Top 10 pathways
  top10_kegg <- kegg_df %>%
    arrange(pvalue) %>%
    slice(1:10) %>% # Synced strictly to top 10 as noted in labels
    mutate(
      logp = -log10(pvalue),
      Description = factor(Description, levels = rev(Description))
    )
  
  # Construct ggplot graphics
  p <- ggplot(top10_kegg, aes(x = Description, y = logp)) +
    geom_bar(stat = "identity", fill = "#F8B606", width = 0.7) +
    coord_flip() +
    labs(
      title = "Top 10 Downregulated KEGG Pathways", 
      x = "Pathway", 
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(color = "black", face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Upgraded size -> linewidth
      axis.line = element_line(color = "black", linewidth = 0.5),              # Upgraded size -> linewidth
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),     # Optimized grid color visibility
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.25),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
    )
  
  # Extract core base name
  base_name <- tools::file_path_sans_ext(basename(file))
  
  # Export multi-format visualization plots
  ggsave(
    filename = file.path(output_dir, paste0(base_name, "_DOWN_KEGG_Barplot.png")),
    plot = p, width = 8, height = 6, dpi = 300
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(base_name, "_DOWN_KEGG_Barplot.pdf")),
    plot = p, width = 8, height = 6
  )
  
  cat("Successfully generated down-regulated KEGG barplot for:", base_name, "\n")
}

cat("DEP down-regulated KEGG visualization pipeline completed.\n")