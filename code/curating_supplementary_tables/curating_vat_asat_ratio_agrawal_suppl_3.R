library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

agrawal <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_vat_sat_gfat_agrawal_suppl_3.xlsx")

# Column renaming # 
colnames(agrawal)
agrawal <- agrawal %>%
  dplyr::select("Trait","CHR","BP","SNP","Effect Allele","Other Allele","EAF","BETA","SE","P-value")

colnames(agrawal) <- c("trait","chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")

#  > Adjusted ones and sex specific traits were not included. < #

# ********************************************************************************************************************************* #

                                                      # For VAT/ASAT ratio #

# -> Trait selection <-
vat_asat_ratio <- agrawal[agrawal$trait =="VAT/ASAT",]

vat_asat_ratio <- vat_asat_ratio %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

vat_asat_ratio <- vat_asat_ratio %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

vat_asat_pos <- vat_asat_ratio

# Function application
vat_asat_pos <- pos_aligner_df(vat_asat_pos)
vat_asat_pos <- alleloi(vat_asat_pos)
vat_asat_pos <- eaf_chooser(vat_asat_pos)
vat_asat_pos <- mhc_cleaner(vat_asat_pos)
                            
# Sample size column addition
vat_asat_pos$sample_size <- NA

# Export and check the data
fwrite(vat_asat_pos,"/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/vat_asat_ratio_agrawal_curated.txt")
check <- fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/vat_asat_ratio_agrawal_curated.txt")

# ********************************************************************************************************************************* #