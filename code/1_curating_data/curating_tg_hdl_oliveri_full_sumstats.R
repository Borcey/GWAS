# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

######################################################################################

# ---------------------------------------------------------------------------------- #

                    # STANDARDIZATION OF GWAS SUMMARY STATS #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch

url_gwas <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90295001-GCST90296000/GCST90295949/GCST90295949.tsv" #correct data

tmp_gwas <- tempfile(fileext = ".tsv") # Temporary file creation
download.file(url_gwas, tmp_gwas, mode = "wb") # Downloading excel file from URL 

tg_hdl_oliveri_gwas <- data.table::fread(tmp_gwas)

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

tg_hdl_oliveri_gwas_clean <- tg_hdl_oliveri_gwas %>%
  dplyr::select("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value", "rs_id", "n")

colnames(tg_hdl_oliveri_gwas_clean) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value", "variant", "sample_size")


# Function applications
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here
tg_hdl_oliveri_gwas_pos <- alleloi(tg_hdl_oliveri_gwas_clean) #removing indels
tg_hdl_oliveri_gwas_pos <- eaf_chooser(tg_hdl_oliveri_gwas_pos)
tg_hdl_oliveri_gwas_pos <- mhc_cleaner(tg_hdl_oliveri_gwas_pos)

# p-value conversion
tg_hdl_oliveri_gwas_pos$p_value <- NA
tg_hdl_oliveri_gwas_pos$p_value <- 10^(-(tg_hdl_oliveri_gwas_pos$neg_log_10_p_value))

tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_pos %>%
  dplyr::select("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "variant", "sample_size", "p_value")

# chr_pos column addition
tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_p %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_p %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) 

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(tg_hdl_oliveri_gwas_p, "output/1_curated_data/tg_hdl_oliveri_full_sumstats.txt")
