library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)
# ---------------------------------------------------------------------------------- #

# Data implementation - T2D Lipodystrophy

lipodys <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/t2dlipodys_mahajan_summary_stat.txt")

colnames(lipodys)

lipodys <- lipodys %>%
  dplyr::select("chromosome(b37)","position(b37)","rsID","effect_allele","other_allele","effect_allele_frequency","Fixed-effects_beta","Fixed-effects_SE","Fixed-effects_p-value")

colnames(lipodys) <- c("chromosome","base_pair_location","rsid","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")

# ---------------------------------------------------------------------------------- #

# P-value filtering 

lipodys <- lipodys %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ---------------------------------------------------------------------------------- #

# Function application

lipodys <- pos_aligner_df(lipodys)    
lipodys <- alleloi(lipodys)            
lipodys <- eaf_chooser(lipodys)       
lipodys <- mhc_cleaner(lipodys)

# ---------------------------------------------------------------------------------- #

# chr_pos column creation

lipodys <- lipodys %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

lipodys <- lipodys %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# Data exporting

fwrite(lipodys, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/clean_gwas_sum_stats/filtered_chrposAdded_gwas/T2Dlipodys_mahajan_gwas_cleaned.txt")

# ---------------------------------------------------------------------------------- #