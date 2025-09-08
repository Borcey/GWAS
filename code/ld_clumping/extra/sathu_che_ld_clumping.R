library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)

# DATA IMPLEMENTATION - SATHU_che_2017

sathu <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/SATHU_che_2017/SATHU_che_2017_chrpos.txt")

sathu <- sathu %>%
  dplyr::select("rsID", "pval", "chr_b37", "pos_b37", "allele1", "allele2")

colnames(sathu) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_sathu <- ld_clump_local(
  dat = sathu,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 1e-5,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_sathu, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/sathu_ld_clumped.txt")