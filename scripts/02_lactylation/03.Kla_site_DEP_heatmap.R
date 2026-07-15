
library(openxlsx)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)

# ============================================================================
# ============================================================================
input_file <- "data/DEP_Kla_Site_Heatmap/2_1981_DE_Site_kla_peptide.xlsx"
output_dir <- "output/DE_kla_site_Heatmap"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================================
df <- read.xlsx(input_file, sheet = 1)

rownames(df) <- make.unique(df$PG.Genes_kla_site)

expr_df <- df[, c("kla_ZIKV_1", "kla_ZIKV_2", "kla_ZIKV_3",
                  "kla_MOCK_1", "kla_MOCK_2", "kla_MOCK_3")]
expr_matrix <- as.matrix(expr_df)

# ============================================================================
# Row Z-score transformation to ensure correct clustering and representation
scaled_matrix <- t(scale(t(expr_matrix)))

group_info <- data.frame(Group = c(rep("ZIKV", 3), rep("MOCK", 3)))
rownames(group_info) <- colnames(scaled_matrix)

group_colors <- list(Group = c(ZIKV = "#ffff33", MOCK = "#74b346"))

ha <- HeatmapAnnotation(df = group_info,
                        col = group_colors,
                        which = "column")

# Custom continuous color mapping from purple to yellow
heat_colors <- colorRampPalette(c("#3e0f73", "#50296e", "#2a9999", "#ccdf39", "#ffff33"))(100)

# ============================================================================
ht <- Heatmap(scaled_matrix,
              name = "Z-score",
              top_annotation = ha,
              col = heat_colors,
              show_row_names = FALSE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(title = "Z-score"))

# Display plot in current R session
draw(ht)

# ============================================================================
# Export PDF format
pdf(file = file.path(output_dir, "DE_kla_site_heatmap.pdf"), width = 6, height = 8)
draw(ht)
dev.off()

# Export high-resolution PNG format
png(file = file.path(output_dir, "DE_kla_site_heatmap.png"),
    width = 1800, height = 2400, res = 300)
draw(ht)
dev.off()