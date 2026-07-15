# ==============================================================================
# Script: 10_mock_timepoint_mfuzz.R
# Description: Time-resolved dynamic fuzzy c-means clustering for MOCK baseline 
#              controls with automated indexing, error-free joins, and export.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readxl)
library(readr)
library(Mfuzz)
library(Biobase)
library(tibble)
library(dplyr)

# 1. Standard Topology Configuration
input_dir  <- "./data/MOCK"
output_dir <- "./results/Mock_output"
cluster_out_dir <- file.path(output_dir, "5_Cluster_output")

if (!dir.exists(input_dir)) {
  stop("Input channel mapping failed, path invalid: ", input_dir)
}
dir.create(cluster_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Ingest Reference Baseline Profile
data_file <- file.path(input_dir, "M_Average.tsv")
if (!file.exists(data_file)) {
  stop("Missing M_Average.tsv file inside data repository.")
}

mock_data <- read_tsv(data_file, na = c("NA", "NaN", ""), show_col_types = FALSE)

# Dissect Metadata and Raw Expression Grids Cleanly
annotation_cols <- mock_data[, 1:5]
expression_data <- mock_data[, -(1:5)]
expression_matrix <- as.matrix(expression_data)
rownames(expression_matrix) <- annotation_cols$Protein.Ids

# 3. Create ExpressionSet Layer & Implement Preprocessing
mfuzz_class <- new('ExpressionSet', exprs = expression_matrix)

# Robust Filtering and Value Imputation Pipeline
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
mfuzz_class <- standardise(mfuzz_class)

# Summary Verification Layer
cat("Dimension of validated baseline expression grid:", dim(exprs(mfuzz_class)), "\n")
if (any(is.na(exprs(mfuzz_class)))) stop("Critical Warning: Remaining NA values detected.")

# 4. Cluster Partition Optimization (Dmin Metric Estimator)
m_val <- mestimate(mfuzz_class)
cat("Optimized fuzzifier value m estimated as:", m_val, "\n")

# Save Selection Dmin Curvature for Publication Quality Verification
pdf_curve_path <- file.path(output_dir, "Dmin_curve.pdf")
pdf(file = pdf_curve_path, width = 10, height = 6)
Dmin(mfuzz_class, m = m_val, crange = seq(2, 10, 1), repeats = 3, visu = TRUE)
dev.off()
cat("Dmin optimization curve safely written out to:", pdf_curve_path, "\n")

# 5. Execute Non-Overlapping Fuzzy C-Means Strategy
set.seed(123)
cluster_num <- 8
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = m_val)

# 6. High-Fidelity Trend Vis Rendering Stage
pdf_trend_path <- file.path(output_dir, "Mock_Cluster_Trends.pdf")
pdf(file = pdf_trend_path, width = 12, height = 8)
mfuzz.plot2(
  mfuzz_class, cl = mfuzz_cluster,
  time.labels = colnames(expression_matrix),
  new.window = FALSE, xlab = "Time point", ylab = "Z-score expression",
  cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.8, lwd = 2.5, plot = TRUE
)
dev.off()

# 7. Safe Joining Layer (CRITICAL FIX: Overcoming Tibble Rowname Loss)
cluster_assignments <- tibble(
  Protein.Ids = rownames(exprs(mfuzz_class)),
  Cluster = mfuzz_cluster$cluster
) %>% left_join(annotation_cols, by = "Protein.Ids")

# Export Profile Matrix Type A: Original Intensity Values Joined
expression_df <- as.data.frame(expression_matrix) %>%
  rownames_to_column(var = "Protein.Ids")

cluster_full_raw <- left_join(cluster_assignments, expression_df, by = "Protein.Ids")

# Multi-file dynamic dumping
for (i in 1:cluster_num) {
  sub_cluster <- cluster_full_raw %>% filter(Cluster == i)
  write.table(
    sub_cluster, file = file.path(cluster_out_dir, paste0("Cluster_", i, "_protein_list_raw.txt")),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
}

# Export Profile Matrix Type B: Normalized Matrix Output for Downstream Heatmap Input
exprs_norm_df <- as.data.frame(exprs(mfuzz_class)) %>%
  rownames_to_column(var = "Protein.Ids")

cluster_full_norm <- left_join(cluster_assignments, exprs_norm_df, by = "Protein.Ids")

# Dump master table configuration
write.table(
  cluster_full_norm, 
  file = file.path(cluster_out_dir, "5_cluster_merged_data_ST_annotations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Dump sub-files matching normalized distribution profiles
for (i in unique(cluster_full_norm$Cluster)) {
  sub_cluster_norm <- cluster_full_norm %>% filter(Cluster == i)
  write.table(
    sub_cluster_norm, file = file.path(cluster_out_dir, paste0("Cluster_", i, "_protein_list_norm.tsv")),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
}

cat("✅ Success: Mfuzz matrix mapping pipelines terminated without errors.\n")