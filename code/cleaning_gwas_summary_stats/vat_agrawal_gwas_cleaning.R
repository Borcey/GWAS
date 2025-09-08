library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

# Data implementation

vat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/vat_agrawal_summary_stat")

colnames(vat)

vat <- vat %>%
  dplyr::select("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")

colnames(vat) <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","info","beta","standard_error","p_value")

# ---------------------------------------------------------------------------------- #

# P-value filtering 

vat <- vat %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ---------------------------------------------------------------------------------- #

# Function application

vat <- pos_aligner_df(vat)    
vat <- alleloi(vat)            
vat <- eaf_chooser(vat)       
vat <- mhc_cleaner(vat)

# ---------------------------------------------------------------------------------- #

# chr_pos column creation

vat <- vat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

vat <- vat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# Data exporting

fwrite(vat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/clean_gwas_sum_stats/filtered_chrposAdded_gwas/vat_agrawal_gwas_cleaned.txt")

# ---------------------------------------------------------------------------------- #