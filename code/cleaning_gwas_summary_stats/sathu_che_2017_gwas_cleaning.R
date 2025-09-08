library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

# ----------------------------------------------------------------------- #

sathu_sum <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/CHE_2017/SATHU_che_2017/SATHU_che_2017_chrpos.txt")

colnames(sathu_sum)

sathu_sum <- sathu_sum %>%
  dplyr::select("rsID","allele1","allele2","N","zscore","pval","chr_b37","pos_b37","chr_b38","pos_b38","ref","alt")

colnames(sathu_sum) <- c("variant","effect_allele","other_allele","sample_size","zscore","p_value","chromosome","base_pair_location","chr_b38","pos_b38","ref","alt")

# ----------------------------------------------------------------------- #

# P-value 

sathu_sum <- sathu_sum %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)  # NONE OF THE P VALUE IS UNDER 5x10^-8

# ----------------------------------------------------------------------- #

# chr_pos column creation

sathu_sum <- sathu_sum %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

sathu_sum <- sathu_sum %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# Upper the alleles
sathu_sum <- sathu_sum %>%
  mutate(effect_allele = toupper(effect_allele))

sathu_sum <- sathu_sum %>%
  mutate(other_allele = toupper(other_allele))

# ----------------------------------------------------------------------- #

# Function application

#sathu_sum <- pos_aligner_df(sathu_sum)     NO BETA VALUE
sathu_sum <- alleloi(sathu_sum)            
#sathu_sum <- eaf_chooser(sathu_sum)        NO EAF
sathu_sum <- mhc_cleaner(sathu_sum)

# ----------------------------------------------------------------------- #

fwrite(sathu_sum, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/sathu_che_2017_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/sathu_che_2017_chrpos.txt")
