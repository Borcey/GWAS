library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)
library(ieugwasr)
library(TwoSampleMR)

agrawal <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_vat_sat_gfat_agrawal_suppl_3.xlsx")

# Column renaming #
colnames(agrawal)
agrawal <- agrawal %>%
  dplyr::select("Trait","CHR","BP","SNP","Effect Allele","Other Allele","EAF","BETA","SE","P-value")

colnames(agrawal) <- c("trait","chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")

#  > Adjusted ones and sex specific traits were not included. < #

# ********************************************************************************************************************************* #

                                                              # For ASAT volume #

# -> Trait selection <-
asat <- agrawal[agrawal$trait =="ASAT",]
  
asat <- asat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

asat <- asat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  

asat_pos <- asat


# Function application
asat_pos <- pos_aligner_df(asat)
asat_pos <- alleloi(asat_pos)
asat_pos <- eaf_chooser(asat_pos)
asat_pos <- mhc_cleaner(asat_pos)

# Sample size column addition
asat_pos$sample_size <- NA

# ********************************************************************************************************************************* #

# Export and check the data
fwrite(asat_pos,"/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/asat_agrawal_curated.txt")
check <- fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/asat_agrawal_curated.txt")

