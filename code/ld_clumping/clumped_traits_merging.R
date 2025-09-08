# Harmonising ld_clumped data


library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)


adipoq_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/adipoq_ld_clumped.txt")
nrow(adipoq_ld)

afr_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/afr_ld_clumped.txt")
nrow(afr_ld)

anterior_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/ant_thigh_muscle_fat_vanDerMeer_ld_clumped.txt")
nrow(anterior_ld)

asat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/asat_ld_clumped.txt")
nrow(asat_ld)

fi_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/fi_ld_clumped.txt")
nrow(fi_ld)

gfat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/gfat_ld_clumped.txt")
nrow(gfat_ld)

leptin_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/leptin_ld_clumped.txt")
nrow(leptin_ld)

lfr_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/lfr_ld_clumped.txt")
nrow(lfr_ld)

liverfat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/liver_fat_liu_ld_clumped.txt")
nrow(liverfat_ld)

pat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/pat_ld_clumped.txt")
nrow(pat_ld)

posterior_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/post_thigh_muscle_fat_vanDerMeer_ld_clumped.txt")
nrow(posterior_ld)

resistin_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/resistin_ld_clumped.txt")
nrow(resistin_ld)

sat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/sat_ld_clumped.txt")
nrow(sat_ld)

t2dlipo_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/t2dlipodys_mahajan_ld_clumped.txt")
nrow(t2dlipo_ld)

tfr_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/tfr_ld_clumped.txt")
nrow(tfr_ld)

vaspin_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/vaspin_ld_clumped.txt")
nrow(vaspin_ld)

vat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/vat_ld_clumped.txt")
nrow(vat_ld)

vatasat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/vatAsatRatio_ld_clumped.txt")
nrow(vatasat_ld)

vatgfat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/vatGfatRatio_ld_clumped.txt")
nrow(vatgfat_ld)

vatsat_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/vatsat_ld_clumped.txt")
nrow(vatsat_ld)

whradjbmi_ld <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/output/whradjbmi_ld_clumped.txt")
nrow(whradjbmi_ld)



final_df <- list(adipoq_ld,afr_ld,anterior_ld,asat_ld,fi_ld,gfat_ld,leptin_ld,liverfat_ld,pat_ld,posterior_ld,resistin_ld,sat_ld,t2dlipo_ld,tfr_ld,vaspin_ld,
                 vat_ld, vatasat_ld,vatgfat_ld,vatsat_ld,whradjbmi_ld)

final_df <- rbindlist(final_df)



# LD - CLUMPING

clumped_all <- ld_clump_local(
  dat = final_df,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_all, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/ld_clumping/all_traits_merged_and_second_clumped.txt")



gormek <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/input/leg_fat_percentage/all_traits_merged_and_second_clumped.txt")

