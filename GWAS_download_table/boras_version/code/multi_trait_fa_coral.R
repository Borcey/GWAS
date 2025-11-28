# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                  # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

options(timeout = 600) # Surpass of timeout risk during download
url_lead <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-022-00731-5/MediaObjects/42255_2022_731_MOESM3_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

multi_trait_fa_coral_lead <- readxl::read_excel(tmp_lead, sheet = 2)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation 

multi_trait_fa_coral_lead_clean <- multi_trait_fa_coral_lead %>%
  dplyr::select("CHR", "POSITION", "RSID", "EA", "NEA", "EAF_BMI", "BETA_BMI", "SE_BMI" , "PVAL_BMI", "N_BMI", "PROFILE")

colnames(multi_trait_fa_coral_lead_clean) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size","profile")

multi_trait_fa_coral_lead_clean <- multi_trait_fa_coral_lead_clean[multi_trait_fa_coral_lead_clean$profile == "Discordant",]


multi_trait_fa_coral_lead_pos <- pos_aligner_df(multi_trait_fa_coral_lead_clean)     
multi_trait_fa_coral_lead_pos <- alleloi(multi_trait_fa_coral_lead_pos)            
multi_trait_fa_coral_lead_pos <- eaf_chooser(multi_trait_fa_coral_lead_pos)        
multi_trait_fa_coral_lead_pos <- mhc_cleaner(multi_trait_fa_coral_lead_pos) # Function applications


multi_trait_fa_coral_lead_pos <- multi_trait_fa_coral_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

multi_trait_fa_coral_lead_pos <- multi_trait_fa_coral_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) # chr_pos column addition

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(multi_trait_fa_coral_lead_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_coral_lead.txt")

# ---------------------------------------------------------------------------------- #