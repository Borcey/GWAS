library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

# Data implementation

vatasat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/vatAsatRatio_agrawal_summary_stat")

colnames(vatasat)

vatasat <- vatasat %>%
  dplyr::select("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")

colnames(vatasat) <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","info","beta","standard_error","p_value")

# ---------------------------------------------------------------------------------- #

# P-value filtering 

vatasat <- vatasat %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ---------------------------------------------------------------------------------- #

# Function application

vatasat <- pos_aligner_df(vatasat)    
vatasat <- alleloi(vatasat)            
vatasat <- eaf_chooser(vatasat)       
vatasat <- mhc_cleaner(vatasat)

# ---------------------------------------------------------------------------------- #

# chr_pos column creation

vatasat <- vatasat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

vatasat <- vatasat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# Data exporting

fwrite(vatasat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/clean_gwas_sum_stats/filtered_chrposAdded_gwas/vatasat_agrawal_gwas_cleaned.txt")

# ---------------------------------------------------------------------------------- #