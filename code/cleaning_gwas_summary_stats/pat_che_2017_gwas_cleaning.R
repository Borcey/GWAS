library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

# ----------------------------------------------------------------------- #

pat_sum <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/CHE_2017/PAT_che_2017/PAT_che_2017_chrpos.txt")

colnames(pat_sum)

pat_sum <- pat_sum %>%
  dplyr::select("rsID","allele1","allele2","N","zscore","pval","chr_b37","pos_b37","chr_b38","pos_b38","ref","alt")

colnames(pat_sum) <- c("variant","effect_allele","other_allele","sample_size","zscore","p_value","chromosome","base_pair_location","chr_b38","pos_b38","ref","alt")

# ----------------------------------------------------------------------- #

# chr_pos column creation

pat_sum <- pat_sum %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

pat_sum <- pat_sum %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# P-value 

pat_sum <- pat_sum %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)    

# ----------------------------------------------------------------------- #

# Upper the alleles
pat_sum <- pat_sum %>%
  mutate(effect_allele = toupper(effect_allele))

pat_sum <- pat_sum %>%
  mutate(other_allele = toupper(other_allele))

# ----------------------------------------------------------------------- #

# Function application

#pat_sum <- pos_aligner_df(pat_sum)     NO BETA VALUE
pat_sum <- alleloi(pat_sum)            
#pat_sum <- eaf_chooser(pat_sum)        NO EAF
pat_sum <- mhc_cleaner(pat_sum)

# ----------------------------------------------------------------------- #
fwrite(pat_sum, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/pat_che_2017_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/pat_che_2017_chrpos.txt")
