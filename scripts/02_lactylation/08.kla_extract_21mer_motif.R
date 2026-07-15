# ============================================================================
rm(list = ls())
Sys.setenv(LANGUAGE = "en")
options(encoding = "UTF-8")
options('repos' = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", update = FALSE)
}

library(openxlsx)
library(Biostrings)
library(stringr)
library(tidyverse)

# ============================================================================
input_xlsx  <- "data/motif/1.1981_DE_Site_kla_peptide_2.xlsx"
input_fasta <- "data/motif/uniprotkb_human_AND_reviewed_true_AND_m_2025_06_09.fasta"
output_dir  <- "output/data/motif"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}


# ============================================================================
df <- read.xlsx(input_xlsx, sheet = 1)

df <- df %>%
  filter(!is.na(PG.Genes)) %>%
  mutate(
    Protein = str_extract(PG.Genes, "^[^\\(]+"),
    RawSite = str_extract(PG.Genes, "(?<=\\()[^\\)]+"),
    K_site = str_extract(RawSite, "K\\d+"),            
    Position = as.numeric(str_extract(K_site, "\\d+")) 
  ) %>%
  drop_na(Protein, Position)


# ============================================================================
fasta <- readAAStringSet(input_fasta)
names(fasta) <- sub(".*GN=([^ ]+).*", "\\1", names(fasta)) 


get_21mer <- function(protein, pos, seq_db, window = 10) {
  if (!protein %in% names(seq_db)) return(NA)
  seq <- as.character(seq_db[[protein]])
  n <- nchar(seq)
  start <- max(1, pos - window)
  end <- min(n, pos + window)
  substr(seq, start, end)
}

df$motif_21mer <- mapply(get_21mer, df$Protein, df$Position, MoreArgs = list(seq_db = fasta))


# ============================================================================
df_motif <- df %>%
  filter(!is.na(motif_21mer)) %>%
  filter(nchar(motif_21mer) == 21) %>%
  mutate(
    left  = substr(motif_21mer, 1, 10),
    mid   = substr(motif_21mer, 11, 11),
    right = substr(motif_21mer, 12, 21)
  ) %>%
  # Filter sequences containing only standard 20 amino acids
  filter(str_detect(motif_21mer, "^[ACDEFGHIKLMNPQRSTVWY]+$")) %>%
  mutate(Unique_ID = paste0(PG.Genes, "_", row_number())) %>%
  select(-left, -mid, -right) 


write.fasta <- function(sequences, names, file) {
  con <- file(file, open = "w")
  for (i in seq_along(sequences)) {
    writeLines(paste0(">", names[i]), con)
    writeLines(sequences[i], con)
  }
  close(con)
}

# ============================================================================

write.fasta(df_motif$motif_21mer, df_motif$Unique_ID, 
            file.path(output_dir, "Kla_filtered_21mer.fasta"))


df_export <- df_motif %>%
  select(PG.Genes, motif_21mer) %>%
  rename(Identifier = PG.Genes, Peptide_21mer = motif_21mer)

write.xlsx(df_export,
           file = file.path(output_dir, "Kla_21mer_motifs.xlsx"),
           sheetName = "Motif_21mer",
           rowNames = FALSE)

# ============================================================================
# ============================================================================
df_motif <- df_motif %>%
  mutate(motif_21mer_noK = paste0(substr(motif_21mer, 1, 10), "-", substr(motif_21mer, 12, 21)))

write.fasta(df_motif$motif_21mer_noK, df_motif$Unique_ID,
            file = file.path(output_dir, "Kla_filtered_21mer_K_to_dash.fasta"))

df_export_noK <- df_motif %>%
  select(PG.Genes, motif_21mer_noK) %>%
  rename(Identifier = PG.Genes, Peptide_21mer_noK = motif_21mer_noK)

write.xlsx(df_export_noK,
           file = file.path(output_dir, "Kla_21mer_motifs_K_to_dash.xlsx"),
           sheetName = "Motif_21mer_noK",
           rowNames = FALSE)