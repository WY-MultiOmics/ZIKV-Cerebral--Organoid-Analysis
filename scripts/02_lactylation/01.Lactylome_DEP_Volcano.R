# ============================================================================
# 1. Environment Infrastructure & Libraries Setup
# ============================================================================
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")

library(openxlsx)
library(readr)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(cowplot)
library(conflicted)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
input_file <- "data/volcano/1.1981_DE_Site_kla_peptide_2.xlsx"
output_dir <- "output/plot/volcano"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

df <- read.xlsx(input_file, sheet = 1, colNames = TRUE)

df$pvalue <- df$q_value
df$Symbol <- ifelse(is.na(df$PG.Genes), df$PG.ProteinGroups, df$PG.Genes)
df <- df[!is.na(df$log2FC) & !is.na(df$pvalue), ]

# Threshold Definitions
pvalue_cutoff <- 0.05
log2FC_cutoff <- 0.58

df$'-log10(pvalue)' <- df$NG_lg_Q 

df$curve_y <- case_when(
  df$log2FC > 0 ~ 0.5 / (df$log2FC - log2FC_cutoff) + (-log10(pvalue_cutoff)),
  df$log2FC <= 0 ~ 0.5 / (-df$log2FC - log2FC_cutoff) + (-log10(pvalue_cutoff))
)

df$group2 <- case_when(
  df$`-log10(pvalue)` > df$curve_y & df$log2FC >= log2FC_cutoff ~ 'up',
  df$`-log10(pvalue)` > df$curve_y & df$log2FC <= -log2FC_cutoff ~ 'down',
  TRUE ~ 'none'
)
df$group2 <- factor(df$group2, levels = c("up", "down", "none"))

f_curve <- function(x_max = 10) {
  inputx <- seq(0.01, x_max, by = 0.01)
  y <- 0.5 / (inputx) + (-log10(pvalue_cutoff))
  rbind(
    data.frame(x = inputx + log2FC_cutoff, y = y),
    data.frame(x = -(inputx + log2FC_cutoff), y = y)
  )
}
dff_curve <- f_curve(8)


top10_proteins <- df %>%
  filter(group2 != "none") %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  slice_max(order_by = abs(log2FC), n = 10) %>%
  pull(Symbol)


interest_proteins <- c(
  "NCL(K116)", "NCL(K403)", "NCL(K63)", "NCL(K223)",
  "HNRNPA2B1(K113)", "HNRNPD(K183)", "HNRNPK(K52)",
  "HNRNPK(K34)", "HNRNPU(K516)"
)

# 3. Combine both lists and remove duplicates seamlessly
final_mark_list <- unique(c(top10_proteins, interest_proteins))
df_label <- df %>% filter(Symbol %in% final_mark_list)

# 6. Volcano Plot Construction (ggplot2)
mycol2 <- c("#C6595A", "#878586", "#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

up_count <- sum(df$group2 == "up")
down_count <- sum(df$group2 == "down")

p_final <- ggplot(data = df, aes(x = log2FC, y = `-log10(pvalue)`, color = group2)) +
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-15, 15), breaks = seq(-15, 15, by = 3)) +
  scale_y_continuous(expand = expansion(add = c(0.1, 0)), limits = c(0, 8), breaks = seq(0, 10, by = 2)) +
  scale_colour_manual(name = "", values = alpha(mycol2, 0.7)) +
  geom_line(data = dff_curve, aes(x = x, y = y), color = "black", lty = "dashed", size = 0.7) +
  geom_text_repel(data = df_label, aes(label = Symbol),
                  size = 2.8, fontface = "italic", color = "#333333",
                  max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
                  segment.color = "grey40", segment.size = 0.2, segment.linetype = "dashed",
                  force = 80, nudge_y = 0.4) +
  # Add Quantified Count Annotation
  annotate("text", x = 14, y = 7.5, label = paste0("Up: ", up_count), color = "#C6595A", size = 6, fontface = "bold", hjust = 1) +
  annotate("text", x = -14, y = 7.5, label = paste0("Down: ", down_count), color = "#878586", size = 6, fontface = "bold", hjust = 0) +
  mytheme

# ============================================================================
# ============================================================================
ggsave(filename = file.path(output_dir, "L_Lactylome_ZIKV_MOCK.pdf"),
       plot = p_final, width = 8, height = 7, dpi = 500, units = "in", device = cairo_pdf)

ggsave(filename = file.path(output_dir, "L_Lactylome_ZIKV_MOCK.png"),
       plot = p_final, width = 8, height = 7, dpi = 500, units = "in")