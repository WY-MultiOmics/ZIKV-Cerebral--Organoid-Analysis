
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(readxl)
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)


input_file <- "data/kla_enrichment/GO_up_enrichment_results.csv"
output_dir <- "output/GO_KEGG_enrichment"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}


goup_data <- read_csv(input_file)
goup_data$log10p <- -log10(goup_data$pvalue)


go_MF <- goup_data %>% filter(ONTOLOGY == "MF") %>% slice(1:10)
go_CC <- goup_data %>% filter(ONTOLOGY == "CC") %>% slice(1:10)
go_BP <- goup_data %>% filter(ONTOLOGY == "BP") %>% slice(1:10)


go_enrich_df <- data.frame(
  ID = c(go_BP$ID, go_CC$ID, go_MF$ID),
  Description = c(go_BP$Description, go_CC$Description, go_MF$Description),
  log10p = c(go_BP$log10p, go_CC$log10p, go_MF$log10p),
  type = factor(c(rep("BP", nrow(go_BP)), rep("CC", nrow(go_CC)), rep("MF", nrow(go_MF))),
                levels = c("MF", "CC", "BP"))
)


go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))


shorten_names <- function(x, n_word = 4, n_char = 40) {
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char)) {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x, " ")[[1]]), n_word)],
                     collapse = " "), "...", sep = "")
    return(x)
  } else {
    return(x)
  }
}

go_enrich_df$short_Description <- sapply(go_enrich_df$Description, shorten_names)


labels <- sapply(go_enrich_df$Description, shorten_names)
names(labels) <- rev(1:nrow(go_enrich_df))

# ============================================================================
CPCOLS <- c("#1f77b4", "#00A289", "#C6595A")

p <- ggplot(data = go_enrich_df, aes(x = number, y = log10p, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = CPCOLS) +
  theme_test() +  
  scale_x_discrete(labels = labels) +  
  labs(title = "Up GO Enrichment", x = "GO term", y = "-log10(p-value)") +
  theme(
    axis.text = element_text(face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  
    axis.line = element_line(color = "black", linewidth = 0.5),  
    panel.grid.major = element_line(color = "gray", linewidth = 0.25),  
    panel.grid.minor = element_line(color = "gray", linewidth = 0.25)   
  )

print(p)

# ============================================================================
# 6. High-Resolution Output Export
# ============================================================================
ggsave(file.path(output_dir, "Kla_up_site_GOEnrichplot.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "Kla_up_site_GOEnrichplot.pdf"), plot = p, width = 10, height = 8, dpi = 300)