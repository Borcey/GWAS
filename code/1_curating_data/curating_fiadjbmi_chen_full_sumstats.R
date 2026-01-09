# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF GWAS SUMMARY STATS #

# ---------------------------------------------------------------------------------- #

# 1 - DATA FETCH 

url_gwas <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002238/GCST90002238_buildGRCh37.tsv.gz"


tmp_gwas <- tempfile(fileext = ".tsv.gz")
download.file(url_gwas, tmp_gwas, mode = "wb")

fiadjbmi_gwas <- fread(tmp_gwas)

# ---------------------------------------------------------------------------------- #

# 2 - Data formatting

fiadjbmi_gwas_clean <- alleloi(fiadjbmi_gwas) #removing indels           
fiadjbmi_gwas_clean <- eaf_chooser(fiadjbmi_gwas_clean)        
fiadjbmi_gwas_clean <- mhc_cleaner(fiadjbmi_gwas_clean) # Function applications

fiadjbmi_gwas_clean <- fiadjbmi_gwas_clean %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

fiadjbmi_gwas_clean <- fiadjbmi_gwas_clean %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) # chr_pos addition

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(fiadjbmi_gwas_clean,"output/1_curated_data/fiadjbmi_chen_full_sumstats.txt")
