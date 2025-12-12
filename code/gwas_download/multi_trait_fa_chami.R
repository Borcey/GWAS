# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# Refference allele = other allele ===> should be asked later

# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

options(timeout = 600) # Surpass of timeout risk during download
url_lead <- url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-025-03931-0/MediaObjects/41591_2025_3931_MOESM2_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

multi_trait_fa_chami_lead <- readxl::read_excel(tmp_lead, sheet = 6, skip = 1)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

multi_trait_fa_chami_lead_clean <- multi_trait_fa_chami_lead %>%
  dplyr::select("SNP","Chr","Position","Reference Allele","Effect Allele","Effect Weight")

colnames(multi_trait_fa_chami_lead_clean) <- c("variant","chromosome","base_pair_location","other_allele","effect_allele","effect_weight")


# Function applications
#multi_trait_fa_chami_lead_pos <- pos_aligner_df(multi_trait_fa_chami_lead_clean)
multi_trait_fa_chami_lead_pos <- alleloi(multi_trait_fa_chami_lead_clean)
#multi_trait_fa_chami_lead_pos <- eaf_chooser(multi_trait_fa_chami_lead_pos)
multi_trait_fa_chami_lead_pos <- mhc_cleaner(multi_trait_fa_chami_lead_pos)



multi_trait_fa_chami_lead_pos <- multi_trait_fa_chami_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

multi_trait_fa_chami_lead_pos <- multi_trait_fa_chami_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(multi_trait_fa_chami_lead_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_chami_leads.txt")

# ---------------------------------------------------------------------------------- #