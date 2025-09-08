library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)
library(readr)

ant_thigh <- readr::read_tsv("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/ant_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.tsv")
colnames(ant_thigh)[colnames(ant_thigh) == "n"] <- "sample_size"

# ----------------------------------------------------------------------- #

# P-value 

ant_thigh <- ant_thigh %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)  

# ----------------------------------------------------------------------- #

# chr_pos column creation

ant_thigh <- ant_thigh %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

ant_thigh <- ant_thigh %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))


# ----------------------------------------------------------------------- #

# Function application
ant_thigh$effect_allele_frequency <- NA

ant_thigh <- pos_aligner_df(ant_thigh)     
ant_thigh <- alleloi(ant_thigh)            
#ant_thigh <- eaf_chooser(ant_thigh) there is no eaf        
ant_thigh <- mhc_cleaner(ant_thigh)

# ----------------------------------------------------------------------- #

fwrite(ant_thigh, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/ant_thigh_muscle_fat_vanDerMeer_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/ant_thigh_muscle_fat_vanDerMeer_chrpos.txt")

# ----------------------------------------------------------------------- #