# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                    # CLEAN LEAD SNPS FROM SUPPLEMENTARY TABLE 2 #

# ---------------------------------------------------------------------------------- #

# 1 - the data is in PDF form. We manually extracted from https://ada.silverchair-cdn.com/ada/content_public/journal/diabetes/68/1/10.2337_db18-0708/5/db180708supplementarydata.pdf?Expires=1765876255&Signature=Ef7vEQXZS-59ORF3494FPfPF4dDvbeC6k04NMegNilP7abNnPsy-NReFEHKR4UK6k4Hwrd1w8Ehp9mfhJWowB-TZSkSix5cM1a1XzWnTEH7GNyLR0mbJowfIUxtnRaQFX9hYucTl6UkiMO-dOzXNYyKihZNS2Y~cyxFllyYDfcu1GiRlfWZdH-sOmYm7iTI0rQkycoPdJjembstEiZettGDFVJm6l4txh3oBcfatrOTXNDqgN~Tx1CKBIcNZzErd-5aXTZyDXKjDvbLD5cFGfVfaJHy7U10gbKW~5t9y1AA2OtkjiVpAJNFh5pPrPKglfEPyz1ycgo5IReakigNAUw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here

ji = fread("output/1_curated_data/multi_trait_fa_ji_leads.txt")

#Let's format it and update it

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

clean_multi_trait_fa_ji_leads_st2 <- ji %>%
  dplyr::select("RSID", "Chr_pos", "Effect allele", "Other allele", "BETA (SD)", "SE", "MetaCCA p-value", "Imputation info")

colnames(clean_multi_trait_fa_ji_leads_st2) <- c("variant","chr_pos","effect_allele","other_allele","beta","standard_error","p_value","info")

# Chromosome column addition
clean_multi_trait_fa_ji_leads_st2 <- clean_multi_trait_fa_ji_leads_st2 %>%
  mutate(
    chromosome = sub(":.*", "", chr_pos),        
    position   = sub(".*:", "", chr_pos)         
  )

clean_multi_trait_fa_ji_leads_st2$chr_pos <- paste0("chr", clean_multi_trait_fa_ji_leads_st2$chr_pos)

clean_multi_trait_fa_ji_leads_st2 <- clean_multi_trait_fa_ji_leads_st2 %>%
  dplyr::rename(base_pair_location = position)


# Function applications
clean_multi_trait_fa_ji_leads_st2 <- alleloi(clean_multi_trait_fa_ji_leads_st2)
#clean_multi_trait_fa_ji_leads_st2 <- eaf_chooser(clean_multi_trait_fa_ji_leads_st2)
clean_multi_trait_fa_ji_leads_st2 <- mhc_cleaner(clean_multi_trait_fa_ji_leads_st2)

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(clean_multi_trait_fa_ji_leads_st2, "output/1_curated_data/multi_trait_fa_ji_leads.txt")