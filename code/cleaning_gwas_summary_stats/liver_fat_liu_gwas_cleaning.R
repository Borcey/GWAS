library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)
library(readr)

liver_fat <- readr::read_tsv("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/liver_fat_liu_gwas_summary_stat.tsv")

# ----------------------------------------------------------------------- #

# chr_pos column creation

liver_fat <- liver_fat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

liver_fat <- liver_fat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ----------------------------------------------------------------------- #

# P-value 

liver_fat <- liver_fat %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)  

# ----------------------------------------------------------------------- #

# Function application

liver_fat <- pos_aligner_df(liver_fat)     
liver_fat <- alleloi(liver_fat)            
liver_fat <- eaf_chooser(liver_fat)        
liver_fat <- mhc_cleaner(liver_fat)

# ----------------------------------------------------------------------- #

fwrite(liver_fat, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/liver_fat_liu_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/liver_fat_liu_chrpos.txt")

# ----------------------------------------------------------------------- #