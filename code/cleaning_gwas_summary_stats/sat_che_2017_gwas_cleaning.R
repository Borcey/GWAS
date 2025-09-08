library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

# ----------------------------------------------------------------------- #

sat_sum <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/CHE_2017/SAT_che_2017/SAT_che_2017_chrpos.txt")

colnames(sat_sum)

sat_sum <- sat_sum %>%
  dplyr::select("rsID","allele1","allele2","N","zscore","pval","chr_b37","pos_b37","chr_b38","pos_b38","ref","alt")

colnames(sat_sum) <- c("variant","effect_allele","other_allele","sample_size","zscore","p_value","chromosome","base_pair_location","chr_b38","pos_b38","ref","alt")

# ----------------------------------------------------------------------- #

# chr_pos column creation

sat_sum <- sat_sum %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

sat_sum <- sat_sum %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# P-value 

sat_sum <- sat_sum %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)    

# ----------------------------------------------------------------------- #

# Upper the alleles
sat_sum <- sat_sum %>%
  mutate(effect_allele = toupper(effect_allele))

sat_sum <- sat_sum %>%
  mutate(other_allele = toupper(other_allele))

# ----------------------------------------------------------------------- #

# Function application

#sat_sum <- pos_aligner_df(sat_sum)     NO BETA VALUE
sat_sum <- alleloi(sat_sum)            
#sat_sum <- eaf_chooser(sat_sum)        NO EAF
sat_sum <- mhc_cleaner(sat_sum)

# ----------------------------------------------------------------------- #
fwrite(sat_sum, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/sat_che_2017_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/sat_che_2017_chrpos.txt")
