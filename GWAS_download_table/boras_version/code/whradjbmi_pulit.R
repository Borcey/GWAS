# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

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
  dplyr::select("SNP", "Chr.ref.males","Pos.ref.males", "A1.combined", "A2.combined", "frqA1.combined", "beta.combined", "se.combined", "pval.combined", "info.combined")

colnames(whradjbmi_pulit_lead_clean) <- c("variant", "chromosome", "base_pair_location","effect_allele","other_allele", "effect_allele_frequency","beta","standard_error","p_value","info")


# Function applications
whradjbmi_pulit_lead_pos <- alleloi(whradjbmi_pulit_lead_clean)
whradjbmi_pulit_lead_pos <- eaf_chooser(whradjbmi_pulit_lead_pos)
whradjbmi_pulit_lead_pos <- pos_aligner_df(whradjbmi_pulit_lead_pos)
whradjbmi_pulit_lead_pos <- mhc_cleaner(whradjbmi_pulit_lead_pos)



whradjbmi_pulit_lead_pos <- whradjbmi_pulit_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

whradjbmi_pulit_lead_pos <- whradjbmi_pulit_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) # chr_pos addition

# ---------------------------------------------------------------------------------- #

# 3 - Data Extraction

fwrite(whradjbmi_pulit_lead_pos, "C:/Users/borac/Desktop/CBMR/GWAS_download_table/output/leadswhradjbmi_pulit_leads.txt")

# ---------------------------------------------------------------------------------- #

                    # STANDARDIZATION OF GWAS SUM STATS #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

url_gwas <- "https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1"

tmp <- tempfile(fileext = ".txt.gz")    
download.file(url_gwas, tmp, mode = "wb")  

whradjbmi_pulit_gwas <- read_table(tmp)   

glimpse(whradjbmi_pulit_gwas)


# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

whradjbmi_pulit_gwas_clean <- whradjbmi_pulit_gwas %>%
  dplyr::select("rsID","Chr...3","Pos (bp)","Effect Allele","Other Allele", "EAF...16","Effect...17","SE...18","GCTA Pvalue...20","sample_size") # Target column selection

colnames(fiadjbmi_lead_clean) <- c("variant","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size") # Standardization of target column names


# P-value filtration
whradjbmi_pulit_gwas_p <- whradjbmi_pulit_gwas %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) 


# Function applications
whradjbmi_pulit_gwas_pos <- pos_aligner_df(whradjbmi_pulit_gwas_p)     
whradjbmi_pulit_gwas_pos <- alleloi(whradjbmi_pulit_gwas_pos)            
whradjbmi_pulit_gwas_pos <- eaf_chooser(whradjbmi_pulit_gwas_pos)        
whradjbmi_pulit_gwas_pos <- mhc_cleaner(whradjbmi_pulit_gwas_pos) 


# chr_pos addition
whradjbmi_pulit_gwas_pos <- whradjbmi_pulit_gwas_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

whradjbmi_pulit_gwas_pos <- whradjbmi_pulit_gwas_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) 

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(whradjbmi_pulit_gwas_pos, "C:/Users/borac/Desktop/CBMR/GWAS_download_table/output/full_sumstats/whradjbmi_pulit_full_sumstats.txt")

# ---------------------------------------------------------------------------------- #