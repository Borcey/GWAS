library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)

# ----------------------------------------------------------------------- #

vathu_sum <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/CHE_2017/VATHU_che_2017/VATHU_che_2017_chrpos.txt") 

colnames(vathu_sum)

vathu_sum <- vathu_sum %>%
  dplyr::select("rsID","allele1","allele2","N","zscore","pval","chr_b37","pos_b37","chr_b38","pos_b38","ref","alt")

colnames(vathu_sum) <- c("variant","effect_allele","other_allele","sample_size","zscore","p_value","chromosome","base_pair_location","chr_b38","pos_b38","ref","alt")

# ----------------------------------------------------------------------- #

# P-value 

vathu_sum <- vathu_sum %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)    # NONE OF THE P VALUE IS UNDER 5x10^-8


# ----------------------------------------------------------------------- #



# chr_pos column creation

vathu_sum <- vathu_sum %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

vathu_sum <- vathu_sum %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# Upper the alleles
vathu_sum <- vathu_sum %>%
  mutate(effect_allele = toupper(effect_allele))

vathu_sum <- vathu_sum %>%
  mutate(other_allele = toupper(other_allele))

# ----------------------------------------------------------------------- #

# Function application

#vathu_sum <- pos_aligner_df(vathu_sum)     NO BETA VALUE
vathu_sum <- alleloi(vathu_sum)            
#vathu_sum <- eaf_chooser(vathu_sum)        NO EAF
vathu_sum <- mhc_cleaner(vathu_sum)

# ----------------------------------------------------------------------- #
fwrite(vathu_sum, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/vathu_che_2017_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/vathu_che_2017_chrpos.txt")
