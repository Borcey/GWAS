library(data.table)
library(tidyverse)
library(dplyr)

tfp <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/raw_gwas_sum_stats/ready_to_use/trunk_fat_percentage/ukb_b_16407_raw/ukb_b_16407_raw.tsv")


tfp <- tfp %>%
  separate(
    col = `UKB-b-16407`,
    into = c("beta","standard_error","p_value","effect_allele_frequency","ID"),
    sep = ":",
    remove = FALSE
  )

colnames(tfp)

tfp <- tfp%>%
  dplyr::select("#CHROM","POS","REF","ALT","FORMAT","UKB-b-16407","beta","standard_error","p_value","effect_allele_frequency","ID")

colnames(tfp) <- c("chromosome","base_pair_location","effect_allele","other_allele","FORMAT","UKB-b-16407","beta","standard_error","p_value","effect_allele_frequency","variant")

tfp$sample_size <- NA
nrow(tfp)


#LOG10P TRANSFORMING
tfp$p_value <- as.numeric(tfp$p_value)
tfp$p_value <- 10^(-tfp$p_value)


# p_value filtering
tfp <- tfp %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# application of 4 functions

tfp <- pos_aligner_df(tfp)
nrow(tfp)
tfp <- alleloi(tfp)
nrow(tfp)
tfp <- eaf_chooser(tfp)
nrow(tfp)
tfp <- mhc_cleaner(tfp)
nrow(tfp)

# create chr_pos column

tfp <- tfp %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

tfp <- tfp %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  
nrow(tfp)



fwrite(tfp, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/cleaned_gwas_sum_stats/tfp_elsworth_gwas_cleaned.txt")
