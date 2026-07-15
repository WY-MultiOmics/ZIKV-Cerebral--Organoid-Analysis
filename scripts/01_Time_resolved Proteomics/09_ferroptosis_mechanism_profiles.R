# ==============================================================================
# Script: 09_ferroptosis_mechanism_profiles.R
# Description: Targeted extraction of ferroptosis-related biomarkers at 10dpi, 
#              automated t-test boxplotting, and high-resolution patchwork assembly.D(DPI),M(MOCK)
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(dplyr)
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggsignif)
library(patchwork)

# 1. Path Topology Setup
ferro_data_dir <- "./data/Ferroptosis_Ex"
box_output_dir <- "./results/Ferroptosis_DEP_10dpi_ggbox"

if (!dir.exists(ferro_data_dir)) {
  dir.create(ferro_data_dir, showWarnings = FALSE, recursive = TRUE)
}
dir.create(box_output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Part A: Target Protein Extraction from 10dpi Proteomics Landscape
diff_10dpi_file <- file.path(ferro_data_dir, "Diff_D10_M10_all.csv")

if (file.exists(diff_10dpi_file)) {
  df_10dpi <- read_csv(diff_10dpi_file, show_col_types = FALSE)
  
  # Crucial Ferroptosis Core Regulators (UniProt Accessions)
  target_proteins <- c(
    "P33121", "O95573", "P02794", "P02792", "P48507", "P36969",
    "P09601", "Q13772", "Q15365", "P04156", "Q96F10", "P49281",
    "Q15043", "P08195", "Q658P3", "P02786", "P04637", "P45880"
  )
  
  filtered_ferro_df <- df_10dpi %>% filter(Accession %in% target_proteins)
  write_tsv(filtered_ferro_df, file.path(ferro_data_dir, "filtered_DEP_Ferro_proteins.tsv"))
  cat("Successfully extracted and cached", nrow(filtered_ferro_df), "ferroptosis-related candidate rows.\n")
} else {
  warning("Base file 'Diff_D10_M10_all.csv' not detected. Skipping extraction phase.")
}

# 3. Part B: Individual Boxplot Rendering Engine
intensity_matrix_file <- file.path(ferro_data_dir, "Ferroptosis_DEP_10dpi.csv")
if (!file.exists(intensity_matrix_file)) {
  stop("Required intensity matrix file missing: ", intensity_matrix_file)
}

ferro_intensity_df <- read_csv(intensity_matrix_file, show_col_types = FALSE) %>% 
  mutate(group = as.character(group))

# Statistical hypothesis testing setup (ZIKV vs MOCK)
comparison_groups <- list(c("ZIKV", "MOCK"))

protein_list <- c(
  "ACSL3", "TFRC", "FTL", "FTH1", "PRNP", "TP53", "SLC3A2", "HMOX1",
  "ACSL1", "GPX4", "VDAC2", "GCLM", "SLC11A2", "NCOA4", "SLC39A14",
  "PCBP1", "STEAP3", "SAT2"
)

# Filter out targets not existing in the matrix headers to avoid crash
valid_proteins <- intersect(protein_list, colnames(ferro_intensity_df))

plot_storage <- list()
cat("Commencing single-gene boxplot rendering...\n")

# Elegant SCI Palette Mapping
sci_colors <- c("ZIKV" = "#C6595A", "MOCK" = "#00A289")

for (protein in valid_proteins) {
  p <- ggboxplot(
    ferro_intensity_df, x = "group", y = protein,
    color = "group", palette = sci_colors,
    add = "jitter", add.params = list(size = 0.8, alpha = 0.7),
    x.text.angle = 0
  ) + 
    labs(x = "Group", y = "Normalized Intensity", title = protein) +
    theme_test() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black", size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    ) +
    geom_signif(
      comparisons = comparison_groups, 
      step_increase = 0.1, 
      map_signif_level = FALSE, 
      test = "t.test", 
      test.args = list(var.equal = TRUE) # Standard Student's t-test matching standard pipelines
    )
  
  plot_storage[[protein]] <- p
  
  # Save individual high-res panels for separate figure building
  ggsave(filename = file.path(box_output_dir, paste0(protein, "_boxplot.png")), plot = p, width = 2.5, height = 4, dpi = 300)
  ggsave(filename = file.path(box_output_dir, paste0(protein, "_boxplot.pdf")), plot = p, width = 2.5, height = 4)
}

# 4. Part C: High-Throughput Matrix Panel Stitching via Patchwork
cat("Assembling panels into high-resolution multi-plot mosaic...\n")
combined_ferro_mosaic <- wrap_plots(plot_storage, ncol = 5) + 
  plot_annotation(
    title = "Expression Profiles of Ferroptosis-Related Markers (10 dpi)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

# Output final composite figures
ggsave(file.path(box_output_dir, "Ferroptosis_combined_plot.png"), plot = combined_ferro_mosaic, width = 13, height = 15, dpi = 300)
ggsave(file.path(box_output_dir, "Ferroptosis_combined_plot.pdf"), plot = combined_ferro_mosaic, width = 13, height = 15)

cat("Ferroptosis extension visualization pipeline completed flawlessly.\n")