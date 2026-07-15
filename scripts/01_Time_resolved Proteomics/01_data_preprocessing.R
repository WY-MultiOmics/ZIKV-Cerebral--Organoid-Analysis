# ==============================================================================
# Script: 01_data_preprocessing.R
# Description: Data cleaning, normalization, and Z-score calculation for 
#              time-resolved proteomics matrix.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# Define input and output paths
input_file <- "./data/time_proteomics/ZIKV_MOCK_unique_PG.tsv"
output_file <- "./data/time_proteomics/Final_Matrix_with_Mean_and_Zscore.tsv"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

data <- read_tsv(input_file)

# Automatically remove all unnamed placeholder columns
data <- data %>% select(-starts_with("..."))

# Separate annotation and expression data
annotation_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description")
annotation <- data %>% select(all_of(annotation_cols))
expr_data <- data %>% select(-all_of(annotation_cols))

# Transform data into long format
expr_long <- expr_data %>%
  mutate(Row_ID = row_number()) %>%
  pivot_longer(-Row_ID, names_to = "Sample", values_to = "Expression")

# Set factor levels for chronological groups
group_levels <- c("D1", "D2", "D3", "D4", "D6", "D8", "D10", 
                  "M1", "M2", "M3", "M4", "M6", "M8", "M10")

expr_long <- expr_long %>%
  mutate(Group = gsub("_[0-9]+$", "", Sample),
         Group = factor(Group, levels = group_levels))

# Compute group means and apply row-wise filter
expr_avg_long <- expr_long %>%
  group_by(Row_ID, Group) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

expr_quality <- expr_avg_long %>%
  group_by(Row_ID) %>%
  filter(!all(is.na(Mean_Expression))) %>%
  ungroup()

expr_avg_quality <- expr_quality %>%
  pivot_wider(names_from = Group, values_from = Mean_Expression) %>%
  arrange(Row_ID)

# Perform row-wise Z-score transformation
expr_scaled <- expr_avg_quality %>%
  select(-Row_ID) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame() %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()

# Merge metrics and apply uniform column naming
expr_combined <- bind_cols(expr_avg_quality %>% select(-Row_ID), expr_scaled)
colnames(expr_combined) <- c(
  paste0(colnames(expr_avg_quality %>% select(-Row_ID)), "_Mean"),
  paste0(colnames(expr_avg_quality %>% select(-Row_ID)), "_Zscore")
)

# Combine expressions with original metadata
final_data <- expr_quality %>%
  distinct(Row_ID) %>%
  left_join(annotation %>% mutate(Row_ID = row_number()), by = "Row_ID") %>%
  select(-Row_ID) %>%
  bind_cols(expr_combined)

# Export preprocessed matrix
dir.create("./data/time_proteomics/", recursive = TRUE, showWarnings = FALSE)
write_tsv(final_data, output_file)
cat("Data preprocessing completed. Output saved to:", output_file, "\n")