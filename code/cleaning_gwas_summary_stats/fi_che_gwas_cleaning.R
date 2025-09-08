library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)


fi_che <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/fi_che/fi_che_raw_gwas.tsv")

# ----------------------------------------------------------------------- #

# P-value 

fi_che <- fi_che %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)  

# ----------------------------------------------------------------------- #

# chr_pos column creation

fi_che <- fi_che %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

fi_che <- fi_che %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# Function application

fi_che <- pos_aligner_df(fi_che)     
fi_che <- alleloi(fi_che)            
fi_che <- eaf_chooser(fi_che)        
fi_che <- mhc_cleaner(fi_che)

# ----------------------------------------------------------------------- #

fwrite(fi_che, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/chr_pos_gwas/fi_che_chrpos.txt")
check<- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/chr_pos_gwas/fi_che_chrpos.txt")

# ----------------------------------------------------------------------- #