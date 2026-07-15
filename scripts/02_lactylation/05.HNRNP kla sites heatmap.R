
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(openxlsx)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(grid)

# ============================================================================
input_file <- "data/HNRNP_Kla_site_Heatmap/HNRNP_Kla.xlsx"
output_dir <- "output/HNRNP_Heatmap"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}


# ============================================================================
df <- read.xlsx(input_file, sheet = 1)

rownames(df) <- make.unique(df$PG.Genes)

expr_df <- df[, c("kla_ZIKV_1", "kla_ZIKV_2", "kla_ZIKV_3",
                  "kla_MOCK_1", "kla_MOCK_2", "kla_MOCK_3")]
expr_matrix <- as.matrix(expr_df)

# Horizontal Z-score Transformation
expr_zscore <- t(scale(t(expr_matrix)))

# Ensure continuous metadata matches perfectly
df$log2FC <- as.numeric(df$log2FC)
df$q_value <- as.numeric(df$q_value)
df <- df[match(rownames(expr_zscore), rownames(df)), ]

# ============================================================================
# ============================================================================
group_info <- data.frame(Group = c(rep("ZIKV", 3), rep("MOCK", 3)))
rownames(group_info) <- colnames(expr_zscore)

group_colors <- list(Group = c(ZIKV = "#ffff33", MOCK = "#74b346"))
heat_colors <- colorRampPalette(c("#3e0f73", "#50296e", "#2a9999", "#ccdf39", "#ffff33"))(100)


ha <- HeatmapAnnotation(df = group_info,
                        col = group_colors,
                        which = "column")

# Right Side Row Annotation (Barplot for FC + Points for Q-value)
row_anno <- rowAnnotation(
  log2FC = anno_barplot(df$log2FC,
                        gp = gpar(fill = "#74b346", col = NA),
                        border = FALSE,
                        width = unit(2, "cm")),
  q_value = anno_points(df$q_value,
                        gp = gpar(col = "#3e0f73"),
                        width = unit(2, "cm")),
  annotation_name_side = "top"
)

# ============================================================================
# ============================================================================
ht <- Heatmap(expr_zscore,
              name = "Z-score",
              top_annotation = ha,
              col = heat_colors,
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(title = "Z-score"),
              right_annotation = row_anno,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y,
                          width = width, height = height,
                          gp = gpar(col = "black", fill = NA, lwd = 0.5))
              }
)


draw(ht)

# ============================================================================

# ============================================================================
# Export PDF
pdf(file.path(output_dir, "HNRNP_Kla_Heatmap.pdf"), width = 8, height = 10)
draw(ht)
dev.off()

png(file.path(output_dir, "HNRNP_Kla_Heatmap.png"), width = 2400, height = 3000, res = 300)
draw(ht)
dev.off()