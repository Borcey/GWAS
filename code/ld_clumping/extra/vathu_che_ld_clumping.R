library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)

# DATA IMPLEMENTATION - VATHU_che_2017

vathu <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/VATHU_che_2017/VATHU_che_2017_chrpos.txt")

vathu <- vathu %>%
  dplyr::select("rsID", "pval", "chr_b37", "pos_b37", "allele1", "allele2")

colnames(vathu) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_vathu <- ld_clump_local(
  dat = vathu,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 1e-5,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_vathu, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/vathu_ld_clumped.txt")