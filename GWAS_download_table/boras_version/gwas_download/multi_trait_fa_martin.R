# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

multi_trait_fa_martin <- readxl::read_xlsx("C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/input/multi_trait_fa_martin_leads.xlsx")

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

multi_trait_fa_martin_clean <- multi_trait_fa_martin %>%
  dplyr::select("RSID", "Chr:Position", "EA", "EAF", "metaCCA P", "BFP beta", "BFP P")

colnames(multi_trait_fa_martin_clean) <- c("variant","chr_pos","effect_allele","effect_allele_frequency","p_value","beta_BFP","p_value_BFP")

multi_trait_fa_martin_clean <- multi_trait_fa_martin_clean %>%
  mutate(
    chromosome = sub(":.*", "", chr_pos),        
    base_pair_location   = sub(".*:", "", chr_pos)         
  )

multi_trait_fa_martin_clean$chr_pos <- paste0("chr", multi_trait_fa_martin_clean$chr_pos)

# Function applications
multi_trait_fa_martin_pos <- eaf_chooser(multi_trait_fa_martin_clean)
multi_trait_fa_martin_pos$other_allele <- "A" # Temporary assignment to eliminate the multi alleles in effect allele. Function would not be worked because there is no other allele.

multi_trait_fa_martin_pos <- alleloi(multi_trait_fa_martin_pos)
multi_trait_fa_martin_pos <- mhc_cleaner(multi_trait_fa_martin_pos)
#multi_trait_fa_martin_pos <- pos_aligner_df(multi_trait_fa_martin_pos) # All of the beta values of BFP are positive. There is no just "beta"

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(multi_trait_fa_martin_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_martin_leads.txt")

# ---------------------------------------------------------------------------------- #