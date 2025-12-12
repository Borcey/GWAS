# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                        # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

                              # For Supplementary Table 2

# 1 - Data Fetch 

multi_trait_fa_ji_leads_st2 <- readxl::read_xlsx("C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/multi_trait_fa_ji_st2.xlsx")

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

clean_multi_trait_fa_ji_leads_st2 <- multi_trait_fa_ji_leads_st2 %>%
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
clean_multi_trait_fa_ji_leads_st2 <- pos_aligner_df(clean_multi_trait_fa_ji_leads_st2)

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(clean_multi_trait_fa_ji_leads_st2, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_ji_leads.txt")

# ---------------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------------- #

                              # For Supplementary Table 4

# 1 - Data Fetch

multi_trait_fa_ji_leads_st4 <- readxl::read_xlsx("C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/multi_trait_fa_ji_st4.xlsx")

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

clean_multi_trait_fa_ji_leads_st4 <- multi_trait_fa_ji_leads_st4 %>%
  dplyr::select("RSID", "Chr:Position", "A1", "A0", "Beta Adipo", "P Adipo", "Beta ALT", "P ALT", "Beta HDL", "P HDL", "Beta SHBG", "P SHBG", "Beta TG", "P TG", "Beta FladjBMI", "P FladjBMI",)

colnames(clean_multi_trait_fa_ji_leads_st4) <- c("variant","base_pair_location","effect_allele","other_allele",)


clean_multi_trait_fa_ji_leads_st4 <- clean_multi_trait_fa_ji_leads_st4 %>%
  mutate(
    chromosome = sub(":.*", "", base_pair_location),        
    position   = sub(".*:", "", base_pair_location)         
  )

# Function applications
clean_multi_trait_fa_ji_leads_st4 <- alleloi(clean_multi_trait_fa_ji_leads_st4)
clean_multi_trait_fa_ji_leads_st4 <- eaf_chooser(clean_multi_trait_fa_ji_leads_st4)
clean_multi_trait_fa_ji_leads_st4 <- mhc_cleaner(clean_multi_trait_fa_ji_leads_st4)
clean_multi_trait_fa_ji_leads_st4 <- pos_aligner_df(clean_multi_trait_fa_ji_leads_st4)

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(clean_multi_trait_fa_ji_leads_st4, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_ji_leads_single_trait_association.txt")

# ---------------------------------------------------------------------------------- #