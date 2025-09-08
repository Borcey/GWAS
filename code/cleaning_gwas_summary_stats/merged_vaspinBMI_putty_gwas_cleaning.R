library(data.table)
library(tidyverse)
library(dplyr)


# vaspin_bmi


#vaspin_bmi <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/merged_all_chr.regenie.gz")



vaspin_bmi$p_value <- as.numeric(vaspin_bmi$LOG10P)
vaspin_bmi$p_value <- 10^(-vaspin_bmi$p_value)


vaspin_bmi <- vaspin_bmi %>%
  dplyr::select("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","CHISQ","p_value")

colnames(vaspin_bmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","info","sample_size","beta","standard_error","chisq","p_value")


vaspin_bmi <- vaspin_bmi %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) 

# Function application
vaspin_bmi <- pos_aligner_df(vaspin_bmi)
vaspin_bmi <- alleloi(vaspin_bmi)
vaspin_bmi <- eaf_chooser(vaspin_bmi)
vaspin_bmi <- mhc_cleaner(vaspin_bmi)

fwrite(vaspin_bmi, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/vaspin/output/vaspinBMI_putty_chrpos_filtered.txt")



check <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/resistinBMI_putty_chrpos_filtered.txt")
ret