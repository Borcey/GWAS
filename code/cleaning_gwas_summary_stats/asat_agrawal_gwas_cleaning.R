library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

# Data implementation

asat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/asat_agrawal_summary_stat")

colnames(asat)

asat <- asat %>%
  dplyr::select("SNP","CHR","BP","ALLELE1","ALLELE0","A1FREQ","INFO","BETA","SE","P_BOLT_LMM_INF")

colnames(asat) <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","info","beta","standard_error","p_value")

# ---------------------------------------------------------------------------------- #

# P-value filtering 

asat <- asat %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ---------------------------------------------------------------------------------- #

# Function application

asat <- pos_aligner_df(asat)    
asat <- alleloi(asat)            
asat <- eaf_chooser(asat)       
asat <- mhc_cleaner(asat)

# ---------------------------------------------------------------------------------- #

# chr_pos column creation

asat <- asat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

asat <- asat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# Data exporting

fwrite(asat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/clean_gwas_sum_stats/filtered_chrposAdded_gwas/asat_agrawal_gwas_cleaned.txt")

# ---------------------------------------------------------------------------------- #
