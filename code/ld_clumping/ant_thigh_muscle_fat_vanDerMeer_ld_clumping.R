library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)

# DATA IMPLEMENTATION -   Ant Thigh Muscle Fat Van Der Meer

atm <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/cleaned_gwas_sum_stats/ant_thigh_muscle_fat_vanDerMeer_gwas_cleaned.txt")

atm <- atm %>%
  dplyr::select("variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(atm) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_atm <- ld_clump_local(
  dat = atm,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/plink_files/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_atm, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/ld_clumped_gwas_sum_stats/lead_snp_gwas/anterior_thigh_muscle_fat_ld_clumped.txt")
