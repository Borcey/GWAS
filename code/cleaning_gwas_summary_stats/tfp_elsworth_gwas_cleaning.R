# GOZLEMLEME
library(data.table)
library(tidyverse)
library(dplyr)

tfr <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/raw_tfr_filtered_gwas_elsworth.txt")


# rename col
colnames(tfr)

tfr <- tfr %>%
  dplyr::select("chromosome","position","variant","effect_allele","other_allele","beta","standard_error","p_value","effect_allele_frequency")
colnames(tfr) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","beta","standard_error","p_value","effect_allele_frequency")

tfr$sample_size <- NA


# create chr_pos column

tfr <- tfr %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

tfr <- tfr %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  


# application of 4 functions
tfr <- pos_aligner_df(tfr)
tfr <- alleloi(tfr)
tfr <- eaf_chooser(tfr)
tfr <- mhc_cleaner(tfr)

fwrite(tfr, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/tfr_elsworth_chrpos.txt")
