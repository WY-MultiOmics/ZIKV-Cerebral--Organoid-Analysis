# ==============================================================================
# Script: 13_zikv_cluster_time_heatmaps.R
# Description: Multi-scheme ComplexHeatmap visualization for ZIKV dynamic clusters
#              featuring row scaling, automatic file handling, and %v% layout.
# ==============================================================================

rm(list = ls())  
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

# 1. Directory and Path Initializations
input_dir   <- "D:/1.APM博士课题/1.1毕业论文/2.Time_MS_data_analysis/R_Time_zikv/R_project_Time_ZIKV_Infect/data"
output_dir  <- "D:/1.APM博士课题/1.1毕业论文/2.Time_MS_data_analysis/R_Time_zikv/R_project_Time_ZIKV_Infect/output/8_Cluster_output/8_cluster_heatmap"
color_test_dir <- file.path(output_dir, "color_tests")

if (!dir.exists(input_dir)) {
  stop("Input mapping directory not found: ", input_dir)
}
dir.create(color_test_dir, showWarnings = FALSE, recursive = TRUE)
setwd(input_dir)

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
generate_cluster_heatmap <- function(cluster_data, cluster_id, target_dir, color_mapping) {
  expr_cols <- c("D1", "D2", "D3", "D4", "D6", "D8", "D10")
  expr_mat  <- as.matrix(cluster_data[, expr_cols])
  rownames(expr_mat) <- cluster_data$ProteinID
  
  # Check if data appears to be Z-scored already by examining range
  data_range <- range(expr_mat, na.rm = TRUE)
  if(data_range[1] < -10 || data_range[2] > 10) {
    expr_mat <- t(scale(t(expr_mat)))
  }
  
  n_proteins <- nrow(expr_mat)
  dynamic_height <- min(max(6, n_proteins / 50), 12)
  
  ht <- Heatmap(expr_mat,
                name = "Expression\nZ-score",
                show_row_names = n_proteins <= 50,
                row_names_gp = gpar(fontsize = 8),
                cluster_rows = TRUE,
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D2",
                cluster_columns = FALSE,
                column_title = paste("Cluster", cluster_id, "-", n_proteins, "proteins"),
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10),
                col = color_mapping,
                heatmap_legend_param = list(
                  title = "Z-score",
                  at = c(-2, -1, 0, 1, 2),
                  labels = c("-2", "-1", "0", "1", "2")
                ),
                use_raster = n_proteins > 100)
  
  pdf(file = file.path(target_dir, paste0("Cluster_", cluster_id, "_heatmap.pdf")), width = 7, height = dynamic_height)
  print(draw(ht, merge_legend = TRUE))
  dev.off()
  
  png(file = file.path(target_dir, paste0("Cluster_", cluster_id, "_heatmap.png")), width = 1800, height = dynamic_height * 300, res = 300)
  print(draw(ht, merge_legend = TRUE))
  dev.off()
  
  return(ht)
}

# 4. Data Ingestion & Target Infrastructure Auditing
master_file <- "8_cluster_merged_data_ST_annotations.tsv"
if (!file.exists(master_file)) {
  stop("Prerequisite baseline cluster file missing: ", master_file)
}

df <- read_tsv(master_file, show_col_types = FALSE)
expr_cols <- c("D1", "D2", "D3", "D4", "D6", "D8", "D10")

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

data_range <- range(sample_mat, na.rm = TRUE)
if(data_range[1] < -10 || data_range[2] > 10) {
  sample_mat <- t(scale(t(sample_mat)))
}

for(scheme_name in names(color_schemes)) {
  pdf(file = file.path(color_test_dir, paste0("sample_", scheme_name, "_heatmap.pdf")), width = 7, height = 8)
  ht_test <- Heatmap(sample_mat,
                     name = "Expression\nZ-score",
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "ward.D2",
                     cluster_columns = FALSE,
                     column_title = paste("Color Scheme:", scheme_name),
                     column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                     column_names_gp = gpar(fontsize = 10),
                     col = color_schemes[[scheme_name]],
                     heatmap_legend_param = list(
                       title = "Z-score",
                       at = c(-2, -1, 0, 1, 2),
                       labels = c("-2", "-1", "0", "1", "2")
                     ))
  draw(ht_test, merge_legend = TRUE)
  dev.off()
}

# 7. Production Target Processing (Selected: viridis_colors)
selected_scheme    <- "viridis_colors"
scheme_folder_name <- paste0("selected_colors_", selected_scheme)
selected_colors_dir <- file.path(output_dir, scheme_folder_name)
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
    pdf(file = file.path(combined_dir, "all_clusters_heatmap.pdf"), width = 8, height = 14)
    
    combined_ht <- NULL
    for(i in 1:length(ht_list)) {
      if(is.null(combined_ht)) {
        combined_ht <- ht_list[[i]]
      } else {
        # Vertical stacking of profiles matching the identical column coordinates
        combined_ht <- combined_ht %v% ht_list[[i]]
      }
    }
    
    draw(combined_ht, ht_gap = unit(0.5, "cm"))
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