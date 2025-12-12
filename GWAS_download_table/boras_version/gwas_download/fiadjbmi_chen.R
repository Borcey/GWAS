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
url_lead <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00852-9/MediaObjects/41588_2021_852_MOESM4_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

fiadjbmi_lead <- readxl::read_excel(tmp_lead, sheet = 2, skip = 3)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation 

fiadjbmi_lead <- fiadjbmi_lead[, 1:20] # Lead variants columns

fiadjbmi_lead <- fiadjbmi_lead[fiadjbmi_lead$Trait == "FI",] # Fasting insulin (FI) trait selection
fiadjbmi_lead <- fiadjbmi_lead[fiadjbmi_lead$`Analysis in which variant association discovered` %in% c("EUR","TA, EUR"),] # European ancestry selection


fiadjbmi_lead$sample_size <- NA # sample_size column addition


fiadjbmi_lead_clean <- fiadjbmi_lead %>%
  dplyr::select("rsID","Chr...3","Pos (bp)","Effect Allele","Other Allele", "EAF...16","Effect...17","SE...18","GCTA Pvalue...20","sample_size") # Target column selection

colnames(fiadjbmi_lead_clean) <- c("variant","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size") # Standardization of target column names


# Function applications
fiadjbmi_lead_pos <- pos_aligner_df(fiadjbmi_lead_clean)
fiadjbmi_lead_pos <- alleloi(fiadjbmi_lead_pos)
fiadjbmi_lead_pos <- eaf_chooser(fiadjbmi_lead_pos)
fiadjbmi_lead_pos <- mhc_cleaner(fiadjbmi_lead_pos)



# Chromosome position (chr_pos) column addition

fiadjbmi_lead_pos <- fiadjbmi_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

fiadjbmi_lead_pos <- fiadjbmi_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(fiadjbmi_lead_pos,"C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/fiadjbmi_chen_leads.txt")

# ---------------------------------------------------------------------------------- #




######################################################################################




# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF GWAS SUMMARY STATS #

# ---------------------------------------------------------------------------------- #

# 1 - DATA FETCH 

url_gwas <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz"


tmp_gwas <- tempfile(fileext = ".tsv.gz")
download.file(url_gwas, tmp_gwas, mode = "wb")

fiadjbmi_gwas <- fread(tmp_gwas)

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation 

fiadjbmi_gwas_p <- fiadjbmi_gwas %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) # P-value filtration



fiadjbmi_gwas_clean <- pos_aligner_df(fiadjbmi_gwas_p)     
fiadjbmi_gwas_clean <- alleloi(fiadjbmi_gwas_clean)            
fiadjbmi_gwas_clean <- eaf_chooser(fiadjbmi_gwas_clean)        
fiadjbmi_gwas_clean <- mhc_cleaner(fiadjbmi_gwas_clean) # Function applications



fiadjbmi_gwas_clean <- fiadjbmi_gwas_clean %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

fiadjbmi_gwas_clean <- fiadjbmi_gwas_clean %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) # chr_pos addition

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(fiadjbmi_lead_clean,"C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/full_sumstats/fiadjbmi_chen_full_sumstats.txt")

# ---------------------------------------------------------------------------------- #