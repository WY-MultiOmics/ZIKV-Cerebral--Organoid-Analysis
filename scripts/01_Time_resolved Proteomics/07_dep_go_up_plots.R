# ==============================================================================
# Script: 07_dep_go_up_plot.R
# Description: Fine-tuned single-file visualization for D1_M1(ZIKV/Mock) critical group 
#              upregulated GO categories with precise scaling and rendering.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)

# 1. Routing Path Architecture
input_file <- "./data/Time_DEP_GO/GO_UP/Diff_D1_M1_up_ALL_GOEnrichmenRsult.xlsx"
output_dir <- "./results/Time_DEP_Enrichment_plot/upGo_plot"

if (!file.exists(input_file)) {
  stop("Critical verification dataset not found: ", input_file)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Dataset Ingestion & Preprocessing
goup_data <- read_excel(input_file)
goup_data$log_p <- -log10(goup_data$pvalue)

# Safeguarded slice mapping for BP, CC, and MF terms
extract_top_terms <- function(df, ontology, top_n = 4) {
  sub_df <- df %>% filter(ONTOLOGY == ontology) %>% arrange(pvalue)
  if (nrow(sub_df) == 0) return(data.frame())
  slice(sub_df, 1:min(nrow(sub_df), top_n))
}

go_BP <- extract_top_terms(goup_data, "BP")
go_CC <- extract_top_terms(goup_data, "CC")
go_MF <- extract_top_terms(goup_data, "MF")

# Merge elements into dedicated layout matrices
combined_df <- bind_rows(go_BP, go_CC, go_MF)

go_enrich_df <- data.frame(
  ID = combined_df$ID,
  Description = combined_df$Description,
  log10pvalue = combined_df$log_p,
  type = factor(combined_df$ONTOLOGY, levels = c("MF", "CC", "BP"))
) %>% arrange(desc(type), log10pvalue)

# Assign inverse factor strings to match coordinate flip orientations
go_enrich_df$number <- factor(1:nrow(go_enrich_df))

# 3. Robust Description Shortening Mechanics
shorten_names <- function(x, n_word = 4, n_char = 40) {
  if (is.na(x)) return("")
  if (nchar(x) > n_char) {
    x <- substr(x, 1, n_char)
  }
  words <- unlist(strsplit(x, " "))
  if (length(words) > n_word) {
    x <- paste(paste(words[1:n_word], collapse = " "), "...", sep = "")
  }
  return(x)
}

plot_labels <- sapply(go_enrich_df$Description, shorten_names)
names(plot_labels) <- go_enrich_df$number

# 4. ggplot Rendering Blueprint
CPCOLS <- c("MF" = "#74b346", "CC" = "#009E73", "BP" = "#0072B2")

p_single <- ggplot(data = go_enrich_df, aes(x = number, y = log10pvalue, fill = type)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS, name = "Ontology") +
  scale_x_discrete(labels = plot_labels) +
  labs(
    title = "Upregulated GO Enrichment (D1_M1)", 
    x = "GO Terms", 
    y = "-log10(p-value)"
  ) +
  theme_test() +
  theme(
    axis.text = element_text(face = "bold", color = "black", size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # size upgraded to linewidth
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
    panel.grid.minor = element_line(color = "gray96", linewidth = 0.25),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# Display graphic inside console environment
print(p_single)

# 5. Targeted Dimensional Image Archiving
ggsave(
  filename = file.path(output_dir, "D1_M1_up_ALL_GOEnrichplot.png"), 
  plot = p_single, width = 10, height = 8, dpi = 300
)

ggsave(
  filename = file.path(output_dir, "D1_M1_up_ALL_GOEnrichplot.pdf"), 
  plot = p_single, width = 10, height = 8
)

cat("Fine-tuned single GO plot generation successfully saved.\n")