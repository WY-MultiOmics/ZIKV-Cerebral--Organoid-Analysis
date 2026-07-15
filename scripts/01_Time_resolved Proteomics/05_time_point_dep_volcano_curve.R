# ==============================================================================
# Script: 06_dep_volcano_curve.R
# Description: Advanced Volcano Plot featuring custom hyperbolic thresholds 
#              and ggrepel gene labeling for publication-ready figures.
# ==============================================================================

rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 1. Environment Infrastructure Setup
pkg <- c("openxlsx", "readr", "tidyverse", "ggplot2", "ggrepel", "patchwork", "cowplot", "dplyr")
newpkg <- pkg[!pkg %in% installed.packages()[, "Package"]]
if (length(newpkg)) {
  sapply(newpkg, function(x) install.packages(x, dependencies = TRUE))
}
sapply(pkg, require, character.only = TRUE)

# 2. File Path Routing (Relative Network Paradigm)
input_file <- "./data/Volcano/Diff_D1_M1_all.csv"
output_dir <- "./results/Time_volcano"

if (!file.exists(input_file)) {
  stop("Input differentiation data file missing: ", input_file)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Load and Curate Raw Dataset
df <- read.csv(input_file, header = TRUE)

# Standardize column naming topologies
df$log2FC <- df$log2.foldchange.
df$pvalue_adj <- df$P_value_adjust
df$Symbol <- ifelse(is.na(df$SYMBOL), df$Gene, df$SYMBOL)

# Drop missing critical mapping elements
df <- df %>% filter(!is.na(log2FC) & !is.na(pvalue))

# 4. Statistical Cut-off Definitions
pvalue_cut <- 0.05
log2FC_cut <- 0.58

# Add descriptive standard -log10 metric
df$`-log10(pvalue)` <- -log10(df$pvalue)

# 5. Hyperbolic Threshold Modeling
# Function to generate custom hyperbolic arrays based on y = c / (x - x0) + y0
generate_hyperbola <- function(max_x) {
  input_x <- seq(0.0001, max_x, by = 0.0001)
  y_pos <- 0.05 / input_x + (-log10(pvalue_cut))
  
  rbind(
    data.frame(x = input_x + log2FC_cut, y = y_pos),
    data.frame(x = -(input_x + log2FC_cut), y = y_pos)
  )
}

dff_curve <- generate_hyperbola(max_x = 5)

# Calculate expected y-coordinate for each real data point along the hyperbola
df$curve_y <- case_when(
  df$log2FC > 0 ~ 0.05 / (df$log2FC - log2FC_cut) + (-log10(pvalue_cut)),
  df$log2FC <= 0 ~ 0.05 / (-df$log2FC - log2FC_cut) + (-log10(pvalue_cut))
)

# Allocate multi-layered classification categories using the hyperbola
df$group2 <- case_when(
  df$`-log10(pvalue)` > df$curve_y & df$log2FC >= log2FC_cut ~ 'up',
  df$`-log10(pvalue)` > df$curve_y & df$log2FC <= -log2FC_cut ~ 'down',
  TRUE ~ 'none'
)
df$group2 <- factor(df$group2, levels = c("up", "down", "none"))

# 6. Top 10 Significant Gene Selection for Labeling
top10_genes <- df %>%
  filter(group2 != "none") %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  top_n(10, abs(log2FC))

# 7. ggplot Technical Drafting
mycol2 <- c("#F8B606", "#4A1985", "#d8d8d8")

mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# Dynamically parse global limits to avoid clipping points
max_y_val <- max(df$`-log10(pvalue)`, na.rm = TRUE)
y_limit_upper <- ceiling(max_y_val + 0.5)

p_volcano <- ggplot(data = df, aes(x = log2FC, y = `-log10(pvalue)`, color = group2)) +
  geom_point(size = 2.2, alpha = 0.7) +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
  scale_y_continuous(
    expand = expansion(add = c(0.1, 0.5)),
    limits = c(0, y_limit_upper)
  ) +
  scale_colour_manual(name = "", values = mycol2) +
  geom_line(
    data = dff_curve, aes(x = x, y = y), 
    color = "black", lty = "dashed", linewidth = 0.7
  ) +
  geom_text_repel(
    data = top10_genes, aes(x = log2FC, y = `-log10(pvalue)`, label = Symbol),
    force = 60, color = 'black', size = 3.5, fontface = "bold",
    point.padding = 0.5, box.padding = 0.5,
    arrow = arrow(length = unit(0.015, "npc"), type = "open", ends = "last"),
    segment.color = "gray30", segment.linewidth = 0.4,
    max.overlaps = Inf
  ) +
  labs(title = "Volcano Distribution Profile", x = "log2(Fold Change)", y = "-log10(p-value)") +
  mytheme

# 8. High-Resolution Multi-Format Graphics Export
ggsave(
  filename = file.path(output_dir, "V_D1_M1_SCIENCE.pdf"),
  plot = p_volcano, width = 7.5, height = 6.0, device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "V_D1_M1_SCIENCE.png"),
  plot = p_volcano, width = 7.5, height = 6.0, dpi = 600
)

# 9. Structured Sub-Dataset Extraction & Excel Sheet Generation
write.xlsx(
  df %>% filter(group2 == "up") %>% arrange(desc(log2FC)),
  file = file.path(output_dir, "Upregulated_Genes_List.xlsx")
)

write.xlsx(
  df %>% filter(group2 == "down") %>% arrange(log2FC),
  file = file.path(output_dir, "Downregulated_Genes_List.xlsx")
)

cat("Advanced hyperbolic volcano modeling and sublist generation completed.\n")