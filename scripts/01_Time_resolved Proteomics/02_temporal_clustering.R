# ==============================================================================
# Script: 02_temporal_clustering.R
# Description: Soft clustering using Mfuzz to identify temporal expression profiles.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

library(Mfuzz)
library(Biobase)
library(tidyverse)

# Define IO directories
input_matrix <- "./data/time_proteomics/D_Average.tsv" 
output_dir <- "./results/time_proteomics_output/8_Cluster_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_matrix)) {
  stop("Required input matrix missing: ", input_matrix)
}

# Load intensity matrix and extract annotations
data <- read_tsv(input_matrix, na = c("NA", "NaN", ""))

annotation_cols <- data[, 1:5]
expression_data <- data[, -(1:5)]
expression_matrix <- as.matrix(expression_data)
rownames(expression_matrix) <- annotation_cols$Protein.Ids

# Initialize ExpressionSet object and execute standard Mfuzz preprocessing
mfuzz_class <- new('ExpressionSet', exprs = expression_matrix)
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
mfuzz_class <- standardise(mfuzz_class)

# Parameter estimation
m <- mestimate(mfuzz_class)
cat("Estimated parameter m:", m, "\n")

# Save Dmin curve optimization plot
pdf_path_dmin <- file.path(output_dir, "Dmin_curve.pdf")
pdf(pdf_path_dmin, width = 8, height = 6)
Dmin(mfuzz_class, m = m, crange = seq(2, 10, 1), repeats = 3, visu = TRUE)
dev.off()
cat("Dmin optimization plot saved to:", pdf_path_dmin, "\n")

# Execute fuzzy c-means clustering
set.seed(123)
cluster_num <- 8
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = m)

# Save cluster trend profiles
pdf_path_trends <- file.path(output_dir, "Mfuzz_Temporal_Clusters.pdf")
pdf(pdf_path_trends, width = 12, height = 6)
mfuzz.plot2(
  mfuzz_class,
  cl = mfuzz_cluster,
  mfrow = c(2, 4),
  time.labels = colnames(expression_matrix),
  xlab = "Time point",
  ylab = "Z-score expression",
  cex.lab = 1.5,
  cex.axis = 1.2,
  cex.main = 1.6,
  lwd = 2.5
)
dev.off()
cat("Temporal cluster profile plots saved to:", pdf_path_trends, "\n")

# Process cluster assignments and metadata merge
cluster_assignments <- tibble(
  Protein.Ids = rownames(exprs(mfuzz_class)),
  Cluster = mfuzz_cluster$cluster
) %>% left_join(annotation_cols, by = "Protein.Ids")

# Extract standardized expression data
exprs_zscore <- as.data.frame(exprs(mfuzz_class))
exprs_zscore$Protein.Ids <- rownames(exprs_zscore)

merged_cluster_data <- left_join(cluster_assignments, exprs_zscore, by = "Protein.Ids")

# Export comprehensive summary and partitioned files
write_tsv(merged_cluster_data, file.path(output_dir, "8_cluster_merged_data_with_annotations.tsv"))

for (i in 1:cluster_num) {
  cluster_subset <- merged_cluster_data %>% filter(Cluster == i)
  write_tsv(cluster_subset, file.path(output_dir, paste0("Cluster_", i, "_protein_list.tsv")))
}

cat("Clustering execution and file export completed under:", output_dir, "\n")