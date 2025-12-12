# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                        # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 -  Data Fetch 

ir_lotta <- readxl::read_xlsx("C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/input/multi_trait_ir_lotta_leads.xlsx")

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation

ir_lotta_clean <- ir_lotta %>%
  dplyr::select("SNP", "Genomic coordinate", "Alleles (effect / other)", "Beta FIadjBMI per allele", "FIadjBMI p-value")

colnames(ir_lotta_clean) <- c("variant", "chr_pos","effect_allele_other_allele","beta","p_value")


ir_lotta_clean <- ir_lotta_clean %>%
  tidyr::separate(chr_pos, into = c("chromosome", "base_pair_location"), sep = ":", remove = FALSE)

ir_lotta_clean$chromosome <- as.numeric(gsub("Chr", "", ir_lotta_clean$chromosome))
ir_lotta_clean$base_pair_location <- as.numeric(ir_lotta_clean$base_pair_location)

ir_lotta_clean$chr_pos <- tolower(ir_lotta_clean$chr_pos)

ir_lotta_clean <- ir_lotta_clean %>%
  tidyr::separate(effect_allele_other_allele, into = c("effect_allele", "other_allele"), sep = " / ", remove = FALSE)

ir_lotta_clean <- ir_lotta_clean %>%
  dplyr::select("variant", "chr_pos", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "p_value")

ir_lotta_clean$variant <- gsub("\\*", "", ir_lotta_clean$variant)
ir_lotta_clean$effect_allele_frequency <- 0.5


# Function applications
ir_lotta_pos <- alleloi(ir_lotta_clean)
#ir_lotta_pos <- eaf_chooser(ir_lotta_pos)
ir_lotta_pos <- mhc_cleaner(ir_lotta_pos)
ir_lotta_pos <- pos_aligner_df(ir_lotta_pos)

# P-value filtration
ir_lotta_p <- ir_lotta_pos %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8) # P-value filtration

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(ir_lotta_p, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_ir_lotta_leads.txt")

# ---------------------------------------------------------------------------------- #