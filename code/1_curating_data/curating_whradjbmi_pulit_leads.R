# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)
options(timeout = 100000) # Surpass of timeout risk during download

# ---------------------------------------------------------------------------------- #

                  # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

url_lead <- "https://raw.githubusercontent.com/lindgrengroup/fatdistnGWAS/master/SuppTable1/whradjbmi.giant-ukbb.meta.1.merged.indexSnps.combined.parsed.txt"
whradjbmi_pulit_lead <- readr::read_table(url_lead)

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

# A1 = effect allele, A2 = other allele

whradjbmi_pulit_lead_clean <- whradjbmi_pulit_lead %>%
  dplyr::select("SNP", "Chr.ref.males","Pos.ref.males", "A1.combined", "A2.combined", "frqA1.combined", "beta.combined", "se.combined", "pval.combined", "info.combined", "nmeta.combined")

colnames(whradjbmi_pulit_lead_clean) <- c("variant", "chromosome", "base_pair_location","effect_allele","other_allele", "effect_allele_frequency","beta","standard_error","p_value","info", "sample_size")

# Function applications

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R")

whradjbmi_pulit_lead_pos <- alleloi(whradjbmi_pulit_lead_clean)
whradjbmi_pulit_lead_pos <- eaf_chooser(whradjbmi_pulit_lead_pos)
whradjbmi_pulit_lead_pos <- pos_aligner_df(whradjbmi_pulit_lead_pos)
whradjbmi_pulit_lead_pos <- mhc_cleaner(whradjbmi_pulit_lead_pos)
whradjbmi_pulit_lead_pos$variant <- as.character(unlist(sapply(whradjbmi_pulit_lead_pos$variant, rsid_parser)))


whradjbmi_pulit_lead_pos <- whradjbmi_pulit_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

whradjbmi_pulit_lead_pos <- whradjbmi_pulit_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) # chr_pos addition


# ---------------------------------------------------------------------------------- #

# 3 - saving data

fwrite(whradjbmi_pulit_lead_pos, "output/1_curated_data/whradjbmi_pulit_leads.txt")

# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF GWAS SUM STATS #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

url_gwas <- "https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"

tmp <- tempfile(fileext = ".txt.gz")    
download.file(url_gwas, tmp, mode = "wb")  

whradjbmi_pulit_gwas <- data.table::fread(tmp) 

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

# Standardization of target column names
whradjbmi_pulit_gwas_clean <- whradjbmi_pulit_gwas %>%
  dplyr::select("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele", "BETA", "SE", "P", "N")

colnames(whradjbmi_pulit_gwas_clean) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size") 

# Function applications
whradjbmi_pulit_gwas_clean <- eaf_chooser(whradjbmi_pulit_gwas_clean) #remove rare variants
whradjbmi_pulit_gwas_clean <- mhc_cleaner(whradjbmi_pulit_gwas_clean) #remove MHC region

whradjbmi_pulit_gwas_clean$variant <- sub(":.*", "", whradjbmi_pulit_gwas_clean$variant)

# chr_pos addition
whradjbmi_pulit_gwas_clean <- whradjbmi_pulit_gwas_clean %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

whradjbmi_pulit_gwas_clean <- whradjbmi_pulit_gwas_clean %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) 

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(whradjbmi_pulit_gwas_clean, "output/1_curated_data/whradjbmi_pulit_full_sumstats.txt")

# ---------------------------------------------------------------------------------- #