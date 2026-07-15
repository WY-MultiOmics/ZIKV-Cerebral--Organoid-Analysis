# ==============================================================================
# Script:08 _ferroptosis_zscore_heatmaps.R
# Description: Mass matrix computation of Ferroptosis Z-score heatmaps across 
#              timepoints with dynamic group labeling and fallback image assembly.D(dpi),M(Mock)
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Verification and alignment of dependencies
required_pkgs <- c("pheatmap", "RColorBrewer", "magick", "gtools", "tidyverse")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(pheatmap)
library(RColorBrewer)
library(magick)
library(gtools)
library(tidyverse)

# 1. Directory Topology Configuration
input_dir  <- "./data/Ferroptosis/Z_Score"
output_dir <- file.path(input_dir, "Ferro_heatmap")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

txt_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(txt_files) == 0) {
  stop("Zero visualization .txt files detected in path: ", input_dir)
}

# 2. Extract Absolute Reference Gene Order (Locking Row Structure)
ref_df <- read.table(txt_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ref_gene_order <- ref_df$PG.Genes

# Define Signature Chromatographic Palette (Custom deep purple to soft gold)
heat_colors <- colorRampPalette(c("#3e0f73", "#50296e", "#2a9999", "#ccdf39", "#ffff33"))(100)

cat("Commencing automated dynamic heatmap pipeline...\n")

# 3. Dynamic Matrix Loop Processing
for (file_path in txt_files) {
  file_base <- tools::file_path_sans_ext(basename(file_path))
  
  # CRITICAL CORRECTION: Dynamically parse group labels from filename (e.g., "D1_M1" -> "D1", "M1")
  group_tags <- unlist(strsplit(file_base, "_"))
  if (length(group_tags) < 2) {
    group_tags <- c("Infected", "Mock") # Fallback naming convention
  }
  
  # Ingest expr matrix
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(df) <- df$PG.Genes
  df_matrix <- df[, 1:6]
  
  # Lock row sequence across different timepoint datasets
  df_matrix <- df_matrix[ref_gene_order[ref_gene_order %in% rownames(df_matrix)], , drop = FALSE]
  
  # Generate accurate metadata frames matching the exact samples
  group_info <- data.frame(Group = c(rep(group_tags[1], 3), rep(group_tags[2], 3)))
  rownames(group_info) <- colnames(df_matrix)
  
  # Establish clear mapping colors
  group_colors <- list(Group = setNames(c("#ffff33", "#74b346"), c(group_tags[1], group_tags[2])))
  
  # Define filenames
  png_out <- file.path(output_dir, paste0(file_base, "_Heatmap.png"))
  pdf_out <- file.path(output_dir, paste0(file_base, "_Heatmap.pdf"))
  
  # Output Render Layer A: Publication PNG
  png(filename = png_out, width = 2400, height = 2571, res = 300)
  pheatmap(df_matrix,
           scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = TRUE,
           annotation_col = group_info, annotation_colors = group_colors,
           fontsize_row = 8, fontsize_col = 8, border_color = NA,
           color = heat_colors, main = paste("Z-score Heatmap:", file_base))
  dev.off()
  
  # Output Render Layer B: Vector PDF
  pdf(file = pdf_out, width = 7, height = 7.5)
  pheatmap(df_matrix,
           scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
           show_rownames = TRUE, show_colnames = TRUE,
           annotation_col = group_info, annotation_colors = group_colors,
           fontsize_row = 8, fontsize_col = 8, border_color = NA,
           color = heat_colors, main = paste("Z-score Heatmap:", file_base))
  dev.off()
}
cat("Batch execution completed. All single matrices saved under precise group parameters.\n")

# 4. Multi-Panel Image Assembly Panel (Stitching Stage)
cat("Initiating multi-panel high-res stitching flow...\n")
pdf_list <- mixedsort(list.files(output_dir, pattern = "_Heatmap\\.pdf$", full.names = TRUE))
png_list <- mixedsort(list.files(output_dir, pattern = "_Heatmap\\.png$", full.names = TRUE))

# Use tryCatch to prevent crash if Ghostscript environment path is unconfigured
tryCatch({
  cat("Attempting vector PDF aggregation via Ghostscript device wrapper...\n")
  pdf_images <- lapply(pdf_list, function(f) image_read_pdf(f, density = 300))
  combined_mosaic <- image_append(image_join(pdf_images), stack = FALSE)
  combined_scaled <- image_scale(combined_mosaic, "3508") # Target standard A4 print-width limit
  
  image_write(combined_scaled, path = file.path(output_dir, "Combined_Heatmap_Row_A4_600dpi.pdf"), format = "pdf")
  image_write(combined_scaled, path = file.path(output_dir, "Combined_Heatmap_Row_A4_600dpi.png"), format = "png")
  cat("Success: PDF stitching complete.\n")
}, error = function(e) {
  cat("Notice: Ghostscript PDF library unavailable. Switching automatically to native PNG array stitching...\n")
  png_images <- lapply(png_list, function(f) image_read(f))
  combined_mosaic <- image_append(image_join(png_images), stack = FALSE)
  combined_scaled <- image_scale(combined_mosaic, "3508")
  
  image_write(combined_scaled, path = file.path(output_dir, "Combined_Heatmap_Row_A4_600dpi.png"), format = "png")
  cat("Success: Raster PNG stitching complete. Final multi-panel saved cleanly.\n")
})

cat("Pipeline 10 completed flawlessly.\n")