library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)


# DATA IMPLEMENTATION - gfat

gfat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/chr_pos_gwas/gfat_agrawal_chrpos.txt")

gfat <- gfat %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(gfat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_gfat <- ld_clump_local(
  dat = gfat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_gfat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/gfat_ld_clumped.txt")



