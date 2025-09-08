library(data.table)
library(tidyverse)
library(dplyr)


# resistin_bmi


#resistin_bmi <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/resistin/input/merged_retn_bmi.tsv")



resistin_bmi$p_value <- as.numeric(resistin_bmi$LOG10P)
resistin_bmi$p_value <- 10^(-resistin_bmi$p_value)


resistin_bmi <- resistin_bmi %>%
  dplyr::select("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","CHISQ","p_value")

colnames(resistin_bmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","info","sample_size","beta","standard_error","chisq","p_value")


resistin_bmi <- resistin_bmi %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) 

# Function application
resistin_bmi <- pos_aligner_df(resistin_bmi)
resistin_bmi <- alleloi(resistin_bmi)
resistin_bmi <- eaf_chooser(resistin_bmi)
resistin_bmi <- mhc_cleaner(resistin_bmi)

fwrite(resistin_bmi, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/resistin/output/resistinBMI_putty_chrpos_filtered.txt")



check <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/resistinBMI_putty_chrpos_filtered.txt")
ret