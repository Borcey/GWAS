# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------------- #

                # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

options(timeout = 600) # Surpass of timeout risk during download
url_lead <- url <- "https://static-content.springer.com/esm/art%3A10.1007%2Fs00125-022-05848-6/MediaObjects/125_2022_5848_MOESM2_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

t2d_kim_lead <- readxl::read_excel(tmp_lead, sheet = 5, skip = 2)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

t2d_kim_lead_clean <- t2d_kim_lead %>%
  dplyr::select("rsID", "VAR_ID_hg19", "T2D risk increasing allele")

colnames(t2d_kim_lead_clean) <- c("variant","VAR_ID_hg19","T2D risk increasing allele")




t2d_kim_lead_clean <- t2d_kim_lead_clean %>%
  separate(
    VAR_ID_hg19,
    into   = c("chromosome", "base_pair_location", "allele1", "allele2"),
    sep    = "_",
    remove = FALSE
  ) %>%
  mutate(
    effect_allele = case_when(
      `T2D risk increasing allele` == allele1 ~ allele1,
      `T2D risk increasing allele` == allele2 ~ allele2,
      TRUE ~ NA_character_
    ),
    other_allele = case_when(
      `T2D risk increasing allele` == allele1 ~ allele2,
      `T2D risk increasing allele` == allele2 ~ allele1,
      TRUE ~ NA_character_
    )
  )

t2d_kim_lead_clean <- t2d_kim_lead_clean %>%
  dplyr::select("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele")

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(t2d_kim_lead_clean, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/t2d_kim_leads.txt")

# ---------------------------------------------------------------------------------- #