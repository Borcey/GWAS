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
url_lead <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01625-2/MediaObjects/41588_2023_1625_MOESM3_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

tg_hdl_oliveri_lead <- readxl::read_excel(tmp_lead, sheet = 3, skip = 2)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

tg_hdl_oliveri_lead_clean <- tg_hdl_oliveri_lead %>%
  dplyr::select("rsIDa", "CHR", "POS (hg19/b37)", "EA", "OA", "EAF", "BETA", "SE","P")

colnames(tg_hdl_oliveri_lead_clean) <- c("variant","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")

tg_hdl_oliveri_lead_clean$p_value <- as.numeric(tg_hdl_oliveri_lead_clean$p_value) # p_value column has to be converted into "numeric" from "character"

# Function applications
tg_hdl_oliveri_lead_pos <- alleloi(tg_hdl_oliveri_lead_clean)
tg_hdl_oliveri_lead_pos <- pos_aligner_df(tg_hdl_oliveri_lead_pos)
tg_hdl_oliveri_lead_pos <- eaf_chooser(tg_hdl_oliveri_lead_pos)
tg_hdl_oliveri_lead_pos <- mhc_cleaner(tg_hdl_oliveri_lead_pos)


# chr_pos column addition
tg_hdl_oliveri_lead_pos <- tg_hdl_oliveri_lead_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

tg_hdl_oliveri_lead_pos <- tg_hdl_oliveri_lead_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation 

fwrite(tg_hdl_oliveri_lead_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/tg_hdl_oliveri_leads.txt")

# # ---------------------------------------------------------------------------------- #




######################################################################################




# ---------------------------------------------------------------------------------- #

                    # STANDARDIZATION OF GWAS SUMMARY STATS #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch

url_gwas <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90295001-GCST90296000/GCST90295949/GCST90295949.tsv"

tmp_gwas <- tempfile(fileext = ".tsv") # Temporary file creation
download.file(url_gwas, tmp_gwas, mode = "wb") # Downloading excel file from URL 

tg_hdl_oliveri_gwas <- data.table::fread(tmp_gwas)

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

tg_hdl_oliveri_gwas_clean <- tg_hdl_oliveri_gwas %>%
  dplyr::select("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value", "rs_id", "n")

colnames(tg_hdl_oliveri_gwas_clean) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value", "variant", "sample_size")


# Function applications
tg_hdl_oliveri_gwas_pos <- alleloi(tg_hdl_oliveri_gwas_clean)
tg_hdl_oliveri_gwas_pos <- eaf_chooser(tg_hdl_oliveri_gwas_pos)
tg_hdl_oliveri_gwas_pos <- mhc_cleaner(tg_hdl_oliveri_gwas_pos)
tg_hdl_oliveri_gwas_pos <- pos_aligner_df(tg_hdl_oliveri_gwas_pos)


# p-value conversion
tg_hdl_oliveri_gwas_pos$p_value <- NA
tg_hdl_oliveri_gwas_pos$p_value <- 10^(-(tg_hdl_oliveri_gwas_pos$neg_log_10_p_value))
                    
                  
# p-value filtration
tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_pos %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_p %>%
  dplyr::select("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "variant", "sample_size", "p_value")



# chr_pos column addition
tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_p %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

tg_hdl_oliveri_gwas_p <- tg_hdl_oliveri_gwas_p %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location)) 

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(tg_hdl_oliveri_gwas_p, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/full_sumstats/tg_hdl_oliveri_full_sumstats.txt")

# ---------------------------------------------------------------------------------- #