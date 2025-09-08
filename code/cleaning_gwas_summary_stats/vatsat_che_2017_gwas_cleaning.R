library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

# ----------------------------------------------------------------------- #

vatsat_sum <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/CHE_2017/VATSAT_che_2017/VATSAT_che_2017_chrpos.txt")

colnames(vatsat_sum)

vatsat_sum <- vatsat_sum %>%
  dplyr::select("rsID","allele1","allele2","N","zscore","pval","chr_b37","pos_b37","chr_b38","pos_b38","ref","alt")

colnames(vatsat_sum) <- c("variant","effect_allele","other_allele","sample_size","zscore","p_value","chromosome","base_pair_location","chr_b38","pos_b38","ref","alt")

# ----------------------------------------------------------------------- #

# chr_pos column creation

vatsat_sum <- vatsat_sum %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

vatsat_sum <- vatsat_sum %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# P-value 

vatsat_sum <- vatsat_sum %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)    

# ----------------------------------------------------------------------- #

# Upper the alleles
vatsat_sum <- vatsat_sum %>%
  mutate(effect_allele = toupper(effect_allele))

vatsat_sum <- vatsat_sum %>%
  mutate(other_allele = toupper(other_allele))

# ----------------------------------------------------------------------- #

# Function application

#vatsat_sum <- pos_aligner_df(vatsat_sum)     NO BETA VALUE
vatsat_sum <- alleloi(vatsat_sum)            
#vatsat_sum <- eaf_chooser(vatsat_sum)        NO EAF
vatsat_sum <- mhc_cleaner(vatsat_sum)

# ----------------------------------------------------------------------- #
fwrite(vatsat_sum, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/vatsat_che_2017_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/vatsat_che_2017_chrpos.txt")
