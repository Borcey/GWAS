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

                                                              # For GFAT volume #

# -> Trait selection <-
gfat <- agrawal[agrawal$trait == "GFAT",]

# Function application
gfat_pos <- pos_aligner_df(gfat)
gfat_pos <- alleloi(gfat_pos)
gfat_pos <- eaf_chooser(gfat_pos)
gfat_pos <- mhc_cleaner(gfat_pos)

# Sample size column addition
gfat_pos$sample_size <- NA

gfat_pos <- gfat_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

gfat_pos <- gfat_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))



# ********************************************************************************************************************************* #

# Export and check the data
fwrite(gfat_pos,"/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/gfat_agrawal_curated.txt")
check <- fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/gfat_agrawal_curated.txt")