library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)

# DATA IMPLEMENTATION -   Liver Fat Liu

livf <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/cleaned_gwas_sum_stats/liver_fat_liu_gwas_cleaned.txt")

livf <- livf %>%
  dplyr::select("variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(livf) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_livf <- ld_clump_local(
  dat = livf,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/plink_files/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_livf, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/ld_clumped_gwas_sum_stats/lead_snp_gwas/liver_fat_ld_clumped.txt")
