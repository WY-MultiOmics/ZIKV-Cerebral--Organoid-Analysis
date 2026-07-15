
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggsignif)
library(patchwork)


# ============================================================================
input_file <- "data/HNRNP_GGBOX/HNRNP_GGBOX.xlsx"
output_dir <- "output/HNRNP_ggbox"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}


# ============================================================================
df <- read.xlsx(input_file, sheet = 1)

df$group <- factor(df$group, levels = c("ZIKV", "MOCK"))
compaired <- list(c("ZIKV", "MOCK"))

protein_list <- c(
  "HNRNPA0_K99", "HNRNPA1_K3", "HNRNPA2B1_K113", "HNRNPC_K157", 
  "HNRNPD_K183", "HNRNPDL_K311", "HNRNPF_K224", "HNRNPH1_K185", 
  "HNRNPK_K52", "HNRNPK_K34", "HNRNPL_K302", "HNRNPM_K698", 
  "HNRNPM_K651", "HNRNPR_K366", "HNRNPU_K516", "HNRNPU_K265", 
  "HNRNPUL2_K58"
)

plot_list <- list()

# ============================================================================
# ============================================================================
for (protein in protein_list) {
  if (protein %in% colnames(df)) {
    p <- ggboxplot(df, x = "group", y = protein,
                   color = "group", add = "jitter", add.params = list(size = 0.5),
                   x.text.angle = 0) + 
      labs(x = "Group", y = "Intensity", title = protein) +
      theme(legend.position = "none") +
      geom_signif(comparisons = compaired, step_increase = 0.1, 
                  map_signif_level = FALSE, test = t.test) + 
      scale_color_manual(values = c("ZIKV" = "#3C5487", "MOCK" = "#00A289"))
    
    print(p)
    plot_list[[protein]] <- p
    
    # Export individual high-resolution figures
    ggsave(filename = file.path(output_dir, paste0(protein, "_boxplot.pdf")), 
           plot = p, width = 2.5, height = 4)
    
    ggsave(filename = file.path(output_dir, paste0(protein, "_boxplot.png")), 
           plot = p, width = 2.5, height = 4, dpi = 300)
  } else {
    cat("Protein", protein, "not found in the dataset\n")
  }
}


# ============================================================================
combined_plot <- wrap_plots(plot_list, ncol = 5) +
  plot_annotation(title = "HNRNP Family Protein Expression")

print(combined_plot)

# Export publication-ready summary figures
ggsave(file.path(output_dir, "HNRNP_Family_combined_plot.pdf"),
       plot = combined_plot, width = 12.5, height = 16)

ggsave(file.path(output_dir, "HNRNP_Family_combined_plot.png"),
       plot = combined_plot, width = 12.5, height = 16, dpi = 300)