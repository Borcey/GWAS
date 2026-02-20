# LD-clumping of 34 merged traits

library(data.table)
library(tidyverse)
library(dplyr)
library(ieugwasr)

merged <- data.table::fread("/home/boraceylan/Desktop/34_traits_ld/merged_34_traits_with_lowest_pval.txt")

merged <- merged %>%
  dplyr::select("variant","p_value","chromosome","base_pair_location","effect_allele","other_allele")

colnames(merged) <- c("rsid","pval","chr","pos","A1","A2")

clumped_merged <- ld_clump_local(
  dat = merged,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p = 5e-8,
  bfile = "/home/boraceylan/Desktop/ld_clumping/EUR",
  plink_bin = "/usr/local/bin/plink"
)

fwrite(clumped_merged, "/home/boraceylan/Desktop/output/2nd_ld_clumping_of_merged_34_traits.txt")
