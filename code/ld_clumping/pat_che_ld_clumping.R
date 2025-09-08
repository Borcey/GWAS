library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)

# DATA IMPLEMENTATION - PAT_che_2017

pat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/cleaned_gwas_sum_stats/pat_che_gwas_cleaned.txt")

pat <- pat %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(pat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_pat <- ld_clump_local(
  dat = pat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/plink_files/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_pat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/ld_clumped_gwas_sum_stats/lead_snp_gwas/pat_ld_clumped.txt")
