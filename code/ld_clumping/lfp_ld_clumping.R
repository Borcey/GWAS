library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)


lfp <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/cleaned_gwas_sum_stats/leg_fat_percentage/lfp_merged_elsworth_gwas_cleaned.txt")


lfp <- lfp %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(lfp) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script




clumped_lfp <- ld_clump_local(
  dat = lfp,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/plink_files/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_lfp, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/ld_clumped_gwas_sum_stats/lead_snp_gwas/lfp_ld_clumped.txt")