library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

# Data implementation

vatgfat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/vatGfatRatio_agrawal_summary_stat")

colnames(vatgfat)

vatgfat <- vatgfat %>%
  dplyr::select("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")

colnames(vatgfat) <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","info","beta","standard_error","p_value")

# ---------------------------------------------------------------------------------- #

# P-value filtering 

vatgfat <- vatgfat %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ---------------------------------------------------------------------------------- #

# Function application

vatgfat <- pos_aligner_df(vatgfat)    
vatgfat <- alleloi(vatgfat)            
vatgfat <- eaf_chooser(vatgfat)       
vatgfat <- mhc_cleaner(vatgfat)

# ---------------------------------------------------------------------------------- #

# chr_pos column creation

vatgfat <- vatgfat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

vatgfat <- vatgfat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# Data exporting

fwrite(vatgfat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/clean_gwas_sum_stats/filtered_chrposAdded_gwas/vatgfat_agrawal_gwas_cleaned.txt")

# ---------------------------------------------------------------------------------- #