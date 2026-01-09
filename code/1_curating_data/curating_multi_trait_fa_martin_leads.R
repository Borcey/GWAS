# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                      # CLEANING OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - The original data is in a PDF!!!! We extracted the data manually, making sure all was good. We will use this code to standardize the data columns and do some check-ups.

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here
multi_trait_fa_martin <- readxl::read_xlsx("output/1_curated_data/multi_trait_fa_martin_leads.xlsx")

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

multi_trait_fa_martin_clean <- multi_trait_fa_martin %>%
  dplyr::select("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

multi_trait_fa_martin_clean$sample_size=442278 #BFP sample size! Then they computed multi-trait associations with other traits. Since they are special cases of BFP - we can use the BFP p-value. This is only relevant for LD-clumping and works well for us

multi_trait_fa_martin_clean$chr_pos <- paste0("chr", multi_trait_fa_martin_clean$chromosome, ":", multi_trait_fa_martin_clean$base_pair_location, sep = "")

# Function applications
multi_trait_fa_martin_clean <- alleloi(multi_trait_fa_martin_clean)
multi_trait_fa_martin_clean <- mhc_cleaner(multi_trait_fa_martin_clean)
multi_trait_fa_martin_clean$effect_allele_frequency=NA
# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(multi_trait_fa_martin_clean, "output/1_curated_data/multi_trait_fa_martin_leads.txt")

# ---------------------------------------------------------------------------------- #