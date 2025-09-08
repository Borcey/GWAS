library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)
library(readr)

post_thigh <- readr::read_tsv("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/post_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.tsv")
colnames(post_thigh)[colnames(post_thigh) == "n"] <- "sample_size"

# ----------------------------------------------------------------------- #

# P-value 

post_thigh <- post_thigh %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)  

# ----------------------------------------------------------------------- #

# chr_pos column creation

post_thigh <- post_thigh %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

post_thigh <- post_thigh %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))


# ----------------------------------------------------------------------- #

# Function application
post_thigh$effect_allele_frequency <- NA

post_thigh <- pos_aligner_df(post_thigh)     
post_thigh <- alleloi(post_thigh)            
#post_thigh <- eaf_chooser(post_thigh) there is no eaf        
post_thigh <- mhc_cleaner(post_thigh)

# ----------------------------------------------------------------------- #

fwrite(post_thigh, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/post_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/post_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.txt")

# ----------------------------------------------------------------------- #