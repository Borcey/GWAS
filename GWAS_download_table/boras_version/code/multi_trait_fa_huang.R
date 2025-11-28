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
url_lead <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-021-00346-2/MediaObjects/42255_2021_346_MOESM2_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

multi_trait_huang_lead <- readxl::read_excel(tmp_lead, sheet = 2, skip=3)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2.1 - Data Manipulation (Singe Trait Association Results)

multi_trait_huang_lead_clean_bmi <- multi_trait_huang_lead %>%
  dplyr::select("...1", "...2","...3","...4","...5","...6","beta BMI", "SE.BMI", "P.BMI", "N.BMI") # Target column selection

colnames(multi_trait_huang_lead_clean_bmi) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size") # Standardization of target column names

multi_trait_huang_lead_clean_bmi <- multi_trait_huang_lead_clean_bmi %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)


multi_trait_huang_lead_bmi_pos <- pos_aligner_df(multi_trait_huang_lead_clean_bmi)
multi_trait_huang_lead_bmi_pos <- alleloi(multi_trait_huang_lead_bmi_pos)
multi_trait_huang_lead_bmi_pos <- eaf_chooser(multi_trait_huang_lead_bmi_pos)
multi_trait_huang_lead_bmi_pos <- mhc_cleaner(multi_trait_huang_lead_bmi_pos)



multi_trait_huang_lead_bmi_pos <- multi_trait_huang_lead_bmi_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

multi_trait_huang_lead_bmi_pos <- multi_trait_huang_lead_bmi_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3.1 - Data Exportation

fwrite(multi_trait_huang_lead_bmi_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_huang_leads_single_trait.txt")

# ---------------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------------- #

# 2.2 - Data Manipulation (Pairwise association results)

multi_trait_huang_lead_clean_pairwise <- multi_trait_huang_lead %>%
  dplyr::select("...1", "...2","...3","...4","...5","...6","P.BMI_HDL", "P.BMI_LDL", "P.BMI_TG", "P.BMI_FI", "P.BMI_FG", "P.BMI_SBP", "P.BMI_CAD", "P.BMI_T2D") # Target column selection

colnames(multi_trait_huang_lead_clean_pairwise) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","pval_BMI_HDL", "pval_BMI_LDL", "pval_BMI_TG", "pval_BMI_FI", "pval_BMI_FG", "pval_BMI_SBP", "pval_BMI_CAD", "pval_BMI_T2D") # Standardization of target column names


#multi_trait_huang_lead_pairwise_pos <- pos_aligner_df(multi_trait_huang_lead_clean_pairwise)
multi_trait_huang_lead_pairwise_pos <- alleloi(multi_trait_huang_lead_clean_pairwise)
multi_trait_huang_lead_pairwise_pos <- eaf_chooser(multi_trait_huang_lead_pairwise_pos)
multi_trait_huang_lead_pairwise_pos <- mhc_cleaner(multi_trait_huang_lead_pairwise_pos)



multi_trait_huang_lead_pairwise_pos <- multi_trait_huang_lead_pairwise_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

multi_trait_huang_lead_pairwise_pos <- multi_trait_huang_lead_pairwise_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3.2 - Data Exportation 

fwrite(multi_trait_huang_lead_pairwise_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_huang_leads_pairwise.txt")

# ---------------------------------------------------------------------------------- #