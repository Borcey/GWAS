library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(dplyr)
library(readr)

whradjbmi <- readr::read_tsv("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/whradjbmi_pulit_summary_stat.txt")
#colnames(whradjbmi)[colnames(ant_thigh) == "n"] <- "sample_size"

# ----------------------------------------------------------------------- #

# P-value 

whradjbmi <- whradjbmi %>%
  dplyr::mutate(P = as.numeric(P)) %>%
  dplyr::filter(!is.na(P), P < 5e-8)  

# ----------------------------------------------------------------------- #

whradjbmi <- whradjbmi %>%
  dplyr::select("CHR","POS","rsid","Tested_Allele","Other_Allele","Freq_Tested_Allele","INFO","N","BETA","SE","P")

colnames(whradjbmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","info","sample_size","beta","standard_error","p_value")

# ----------------------------------------------------------------------- #


# chr_pos column creation

whradjbmi <- whradjbmi %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

whradjbmi <- whradjbmi %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))


# ----------------------------------------------------------------------- #

# Function application

whradjbmi <- pos_aligner_df(whradjbmi)     
whradjbmi <- alleloi(whradjbmi)            
whradjbmi <- eaf_chooser(whradjbmi) 
whradjbmi <- mhc_cleaner(whradjbmi)

# ----------------------------------------------------------------------- #

fwrite(whradjbmi, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/ant_thigh_muscle_fat_vanDerMeer_chrpos.txt")
check<- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/ant_thigh_muscle_fat_vanDerMeer_chrpos.txt")

# ----------------------------------------------------------------------- #