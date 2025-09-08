library(data.table)
library(tidyverse)
library(dplyr)

# leptin_bmi


#leptin_bmi <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/leptin/merged_leptinBMI.tsv")



leptin_bmi$p_value <- as.numeric(leptin_bmi$LOG10P)
leptin_bmi$p_value <- 10^(-leptin_bmi$p_value)


leptin_bmi <- leptin_bmi %>%
  dplyr::select("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","CHISQ","p_value")

colnames(leptin_bmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","info","sample_size","beta","standard_error","chisq","p_value")


leptin_bmi <- leptin_bmi %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) 

fwrite(leptin_bmi, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/leptin/output/leptinBMI_putty_chrpos_filtered.txt")

