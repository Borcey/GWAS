library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

main_table <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_afr_lfr_tfr_andersen_suppl_1.xlsx")
colnames(main_table)

# There are no standard error, other allele, effect allele frequency # 
# A1 is effect or other allele? Not sure, so implemented as A1 #

# ******************************************************************************************************************************************************* #
# ******************************************************************************************************************************************************* #

afr <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_afr_lfr_tfr_andersen_suppl_1.xlsx")
afr <- afr %>%
  dplyr::select("CHR","size (kb)","RSID","BP","A1_AFR","N_AFR","BETA_AFR","p_corr_AFR") # size (kb) is size of what? There is also N. 
colnames(afr) <- c("chromosome","size_kb","rsid","base_pair_location","A1","sample_size","beta","p_value")  

afr <- afr%>%
  rename(effect_allele = A1)    # ASSUMED A1 AS EFFECT ALLELE

afr$other_allele <- NA
afr$effect_allele_frequency <- NA
afr$standard_error <- NA

afr$p_value <- as.numeric(afr$p_value)
afr <- afr %>%
  filter(p_value < 5e-8)

afr <- afr %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

afr <- afr %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  

afr <- pos_aligner_df(afr)
#afr <- alleloi(afr)
#afr <- eaf_chooser(afr)
afr <- mhc_cleaner(afr)



fwrite(afr, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/afr_andersen_curated.txt")
check <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/afr_andersen_curated.txt")

# ******************************************************************************************************************************************************* #
# ******************************************************************************************************************************************************* #

lfr <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_afr_lfr_tfr_andersen_suppl_1.xlsx")
lfr <- lfr %>%
  dplyr::select("CHR","Locus_range","size (kb)","RSID","BP","A1_LFR","N_LFR","BETA_LFR","p_corr_LFR")
colnames(lfr) <- c("chromosome","locus_range","size_kb","rsid","base_pair_location","A1","sample_size","beta","p_value")

lfr <- lfr%>%
  rename(effect_allele = A1)    # ASSUMED A1 AS EFFECT ALLELE

lfr$other_allele <- NA
lfr$effect_allele_frequency <- NA
lfr$standard_error <- NA


lfr$p_value <- as.numeric(lfr$p_value)
lfr <- lfr %>%
  filter(p_value < 5e-8)

lfr <- lfr %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

lfr <- lfr %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  


lfr <- pos_aligner_df(lfr)
#lfr <- alleloi(lfr)
#afr <- eaf_chooser(lfr)
lfr <- mhc_cleaner(lfr)


fwrite(lfr, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/lfr_andersen_curated.txt")
check <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/lfr_andersen_curated.txt")

# ******************************************************************************************************************************************************* #
# ******************************************************************************************************************************************************* #

tfr <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_afr_lfr_tfr_andersen_suppl_1.xlsx")
tfr <- tfr %>%
  dplyr::select("CHR","Locus_range","size (kb)","RSID","BP","A1_TFR","N_TFR","BETA_TFR","p_corr_TFR")
colnames(tfr) <- c("chromosome","locus_range","size_kb","rsid","base_pair_location","A1","sample_size","beta","p_value")

tfr <- tfr%>%
  rename(effect_allele = A1)    # ASSUMED A1 AS EFFECT ALLELE

tfr$other_allele <- NA
tfr$effect_allele_frequency <- NA
tfr$standard_error <- NA


lfr$p_value <- as.numeric(lfr$p_value)
tfr <- tfr %>%
  filter(p_value < 5e-8)

tfr <- tfr %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

tfr <- tfr %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  

tfr <- pos_aligner_df(tfr)
#tfr <- alleloi(tfr)
#afr <- eaf_chooser(tfr)
tfr <- mhc_cleaner(tfr)



fwrite(tfr, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/tfr_andersen_curated.txt")
check <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/tfr_andersen_curated.txt")

# ******************************************************************************************************************************************************* #
# ******************************************************************************************************************************************************* #

main_table$A1_AFR <- as.character(main_table$A1_AFR)
main_table$A1_LFR <- as.character(main_table$A1_LFR)
main_table$A1_TFR <- as.character(main_table$A1_TFR)

identical(main_table$A1_TFR, main_table$A1_LFR)
