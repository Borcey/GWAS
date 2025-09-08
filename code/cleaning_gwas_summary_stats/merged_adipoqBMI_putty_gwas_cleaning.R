library(data.table)
library(tidyverse)
library(dplyr)

# adipoq_bmi


#adipoq_bmi <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/merged_all_chr.regenie.gz")



adipoq_bmi$p_value <- as.numeric(adipoq_bmi$LOG10P)
adipoq_bmi$p_value <- 10^(-adipoq_bmi$p_value)


adipoq_bmi <- adipoq_bmi %>%
  dplyr::select("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","CHISQ","p_value")

colnames(adipoq_bmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","info","sample_size","beta","standard_error","chisq","p_value")


adipoq_bmi <- adipoq_bmi %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) 

# Function application
adipoq_bmi <- pos_aligner_df(adipoq_bmi)
adipoq_bmi <- alleloi(adipoq_bmi)
adipoq_bmi <- eaf_chooser(adipoq_bmi)
adipoq_bmi <- mhc_cleaner(adipoq_bmi)

fwrite(adipoq_bmi, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/adiponectin/output/adipoqBMI_putty_chrpos_filtered.txt")

