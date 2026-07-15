# ==============================================================================
# Script: 12_mock_cluster_time_heatmaps.R
# Description: Multi-scheme ComplexHeatmap visualization for MOCK dynamic clusters
#              featuring row scaling, automatic file handling, and %v% layout.
# ==============================================================================

rm(list = ls())  
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

# 1. Directory and Path Initializations
input_dir   <- "./results/Mock_output/5_Cluster_output"
heatmap_dir <- file.path(input_dir, "heatmap_color_set")
color_test_dir <- file.path(heatmap_dir, "color_tests")

if (!dir.exists(input_dir)) {
  stop("Input mapping directory not found: ", input_dir)
}
dir.create(color_test_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Dependency Loading
library(openxlsx)
library(readr)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(grid)
library(viridis)  
library(RColorBrewer)

# 3. Encapsulated Vector-Regulated Heatmap Engine
generate_cluster_heatmap <- function(cluster_data, cluster_id, output_dir, color_mapping) {
  expr_cols <- c("M1_Mean", "M2_Mean", "M3_Mean", "M4_Mean", "M6_Mean", "M8_Mean", "M10_Mean")
  expr_mat  <- as.matrix(cluster_data[, expr_cols])
  rownames(expr_mat) <- cluster_data$ProteinID
  
  # Row-wise standard normalization (Z-score transformation)
  expr_mat <- t(scale(t(expr_mat)))
  
  ht <- Heatmap(expr_mat,
                name = paste0("Cluster ", cluster_id),
                show_row_names = FALSE,
                cluster_rows = TRUE,
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D2",
                cluster_columns = FALSE,
                column_title = paste("Cluster", cluster_id),
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10),
                col = color_mapping,
                heatmap_legend_param = list(
                  title = "Z-score",
                  at = c(-2, -1, 0, 1, 2),
                  labels = c("-2", "-1", "0", "1", "2")
                ))
  
  pdf(file = file.path(output_dir, paste0("Cluster_", cluster_id, "_heatmap.pdf")), width = 6, height = 8)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  png(file = file.path(output_dir, paste0("Cluster_", cluster_id, "_heatmap.png")), width = 1200, height = 1600, res = 200)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  return(ht)
}

# 4. Data Ingestion & Target Infrastructure Auditing
master_file <- file.path(input_dir, "5_cluster_merged_data_ST_annotations.tsv")
if (!file.exists(master_file)) {
  stop("Prerequisite baseline cluster file missing: ", master_file)
}

df <- read_tsv(master_file, show_col_types = FALSE)
expr_cols <- c("M1_Mean", "M2_Mean", "M3_Mean", "M4_Mean", "M6_Mean", "M8_Mean", "M10_Mean")

if(!all(expr_cols %in% colnames(df))) {
  missing_cols <- expr_cols[!expr_cols %in% colnames(df)]
  stop(paste("Target intensity array intersection breach. Missing dimensions:", paste(missing_cols, collapse=", ")))
}

df$Cluster <- factor(df$Cluster)

# 5. Multipolar Analytical Color Schemes Dictionary
color_schemes <- list(
  classic_enhanced = colorRamp2(c(-2, -1, 0, 1, 2), c("#053061", "#6bacd0", "#f7f7f7", "#ef8a62", "#67001f")),
  bioinfo          = colorRamp2(c(-2, -1, 0, 1, 2), c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#b2182b")),
  nature           = colorRamp2(c(-2, -1, 0, 1, 2), c("#1a9850", "#91cf60", "#ffffbf", "#d73027", "#7f3b08")),
  high_contrast    = colorRamp2(c(-2, -1, 0, 1, 2), c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c")),
  viridis_colors   = colorRamp2(c(-2, -1, 0, 1, 2), viridis(5)),
  plasma_colors    = colorRamp2(c(-2, -1, 0, 1, 2), plasma(5)),
  temperature      = colorRamp2(c(-2, -1, 0, 1, 2), c("#313695", "#74add1", "#f7f7f7", "#f46d43", "#a50026")),
  gradient         = colorRamp2(c(-2, -1, 0, 1, 2), c("#3b4cc0", "#99c7ec", "#f7f7f7", "#f8a07e", "#b24745"))
)

# 6. Single Cluster Evaluation Benchmarking Runs
sample_cluster_id <- levels(df$Cluster)[1]
sample_data <- df %>% filter(Cluster == sample_cluster_id)
sample_mat  <- as.matrix(sample_data[, expr_cols])
rownames(sample_mat) <- sample_data$ProteinID
sample_mat  <- t(scale(t(sample_mat)))

for(scheme_name in names(color_schemes)) {
  pdf(file = file.path(color_test_dir, paste0("sample_", scheme_name, "_heatmap.pdf")), width = 7, height = 8)
  ht_test <- Heatmap(sample_mat,
                     name = "Z-score",
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE,
                     column_title = paste("Color Scheme:", scheme_name),
                     col = color_schemes[[scheme_name]])
  draw(ht_test, merge_legend = TRUE)
  dev.off()
}

# 7. Production Target Processing (Selected: high_contrast)
selected_scheme    <- "high_contrast"
scheme_folder_name <- paste0("selected_colors_", selected_scheme)
selected_colors_dir <- file.path(heatmap_dir, scheme_folder_name)
dir.create(selected_colors_dir, showWarnings = FALSE, recursive = TRUE)

ht_list <- list()
for (k in levels(df$Cluster)) {
  cluster_data <- df %>% filter(Cluster == k)
  tryCatch({
    ht <- generate_cluster_heatmap(cluster_data, k, selected_colors_dir, color_schemes[[selected_scheme]])
    ht_list[[k]] <- ht
  }, error = function(e) {
    print(paste("Fault event captured inside cluster mapping matrix loops on key [", k, "]: ", e$message))
  })
}

# 8. Composite Layout Construction Pipeline
combined_dir <- file.path(selected_colors_dir, "combined")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

tryCatch({
  if(length(ht_list) == length(levels(df$Cluster))) {
    pdf(file = file.path(combined_dir, "all_clusters_heatmap.pdf"), width = 8, height = 18)
    
    combined_ht <- NULL
    for(i in 1:length(ht_list)) {
      if(is.null(combined_ht)) {
        combined_ht <- ht_list[[i]]
      } else {
        # Vertical stacking of profiles matching the identical column coordinates
        combined_ht <- combined_ht %v% ht_list[[i]]
      }
    }
    
    draw(combined_ht, ht_gap = unit(0.5, "cm"), merge_legend = TRUE)
    dev.off()
    print("Integrated layout generated and written out safely.")
  }
}, error = function(e) {
  print(paste("Primary canvas allocation failed: ", e$message, ". Reverting to partitioned multi-page system..."))
  tryCatch({
    pdf(file = file.path(combined_dir, "selected_all_clusters_separate_pages.pdf"), width = 7, height = 8)
    for(k in names(ht_list)) {
      draw(ht_list[[k]], merge_legend = TRUE)
    }
    dev.off()
  }, error = function(ex) {
    print(paste("Critical error in secondary multi-page structural writing: ", ex$message))
  })
})