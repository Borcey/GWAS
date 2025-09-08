library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# Supplementary Table 2. Fourthin variants associated with "favorable adiposity". The imputation information is for the UK Biobank study. 
# BETA, SE and P-value are reported for the GWAS of body fat % in UK Biobank. SE: standard error; HWE: Hardy–Weinberg equilibrium.

# EAF is not included in the suppl. table

adipocardio <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_favoradip_ji_suppl_2.txt")

# ********************************************************************************************************************************* #

# Column renaming # 

colnames(adipocardio)
adipocardio_clean <- adipocardio %>%
  dplyr::select("RSID","Chr_pos","Effect allele","Other allele","BETA (SD)","SE","p-value","MetaCCA p-value")

colnames(adipocardio_clean) <- c("variant","chr_pos","effect_allele","other_allele","beta","standard_error","p_value","p_value_meta")

# ********************************************************************************************************************************* #

# Chr:bp seperation # 

adipocardio_clean <- adipocardio_clean %>%
  separate(col = chr_pos, into = c("chromosome", "base_pair_location"), sep = ":", remove = FALSE)

adipocardio_clean <- adipocardio_clean %>%
  mutate(chr_pos = paste0("chr", chr_pos))


# Functions #
adipo_pos <- pos_aligner_df(adipocardio_clean)
adipo_pos <- alleloi(adipo_pos)
#adipo_pos <- eaf_chooser(adipo_pos)
adipo_pos <- mhc_cleaner(adipo_pos)


adipo_pos$chromosome <- as.numeric(adipo_pos$chromosome)
adipo_pos$base_pair_location <- as.numeric(adipo_pos$base_pair_location)


# ********************************************************************************************************************************* #
fwrite(adipo_pos, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/favoradip_ji_suppl_2_curated.txt")


