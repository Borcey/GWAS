library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)


# DATA IMPLEMENTATION - asat

asat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/chr_pos_gwas/asat_agrawal_chrpos.txt")

asat <- asat %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(asat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_asat <- ld_clump_local(
  dat = asat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_asat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/asat_ld_clumped.txt")
