library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

agrawal <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_vat_sat_gfat_agrawal_suppl_3.xlsx")
colnames(agrawal)

agrawal <- agrawal %>%
  dplyr::select("Trait","CHR","BP","SNP","Effect Allele","Other Allele","EAF","BETA","SE","P-value")

colnames(agrawal) <- c("trait","chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")


#  adjusted ones and sex specific traits were not included. 

# ************************************************************************************************************************************* #
                                    
                                                               # For VAT volume
# -> Trait selection <-
vat <- agrawal[agrawal$trait =="VAT",]

vat <- vat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

vat <- vat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

vat_pos <- vat

# Function application
vat_pos <- pos_aligner_df(vat)
vat_pos <- alleloi(vat_pos)
vat_pos <- eaf_chooser(vat_pos)
vat_pos <- mhc_cleaner(vat_pos)


# Sample size column addition
vat_pos$sample_size <- NA

# ************************************************************************************************************************************* #

fwrite(vat_pos,"/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/vat_agrawal_curated.txt")
