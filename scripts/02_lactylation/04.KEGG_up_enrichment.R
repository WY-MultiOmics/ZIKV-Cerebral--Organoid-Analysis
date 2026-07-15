library(readxl)
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)


input_file <- "data/kla_enrichment/KEGG_up_enrichment_results.csv"
output_dir <- "output/GO_KEGG_enrichment"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}


kegg_df <- read_csv(input_file)

if (anyNA(kegg_df)) {
  stop("Dataset contains missing values. Please check the input file!")
}

expected_columns <- c("Description", "pvalue")
if (!all(expected_columns %in% colnames(kegg_df))) {
  stop("Missing required columns! Columns 'Description' and 'pvalue' are mandatory.")
}


top10_kegg <- kegg_df %>%
  arrange(pvalue) %>%
  slice(1:10) %>%
  mutate(logp = -log10(pvalue),
         Description = factor(Description, levels = rev(Description)))


p <- ggplot(top10_kegg, aes(x = Description, y = logp)) +
  geom_bar(stat = "identity", fill = "#0072B2", show.legend = FALSE) +
  coord_flip() +
  labs(title = "Top 10 KEGG Pathways", x = "Pathway", y = "-log10(p-value)") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "gray", linewidth = 0.25),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.25)
  )

print(p)


# ============================================================================
ggsave(file.path(output_dir, "Kla_up_KEGG.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Kla_up_KEGG.pdf"), plot = p, width = 10, height = 8, dpi = 300)