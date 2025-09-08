library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(dplyr)
library(ieugwasr)
library(ggplot2)
library(dplyr)

# p = 5e-8

# ************************************************************************************************************************************************ # 
# ------------------------------------------------------------------------------------------------------------------------------------------------ # 

#                                                                   A D I P O K I N E S                                                            #



# DATA IMPLEMENTATION - ADIPONECTIN


#adipoq <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
adipoq <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/adiponectin/adipoqBMI_putty_chrpos_filtered.txt")

adipoq <- adipoq %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(adipoq) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_adipoq <- ld_clump_local(
  dat = adipoq,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_adipoq, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/adipoq_ld_clumped.txt")


# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - LEPTIN 


#leptin <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/leptinBMI_putty_chrpos_filtered.txt")
leptin <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/leptin/leptinBMI_putty_chrpos_filtered.txt")

leptin <- leptin %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(leptin) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_leptin <- ld_clump_local(
  dat = leptin,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_leptin, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/leptin_ld_clumped.txt")


# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - OMENTIN

# MISSING DATA !!!!!! #

#omentin <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
omentin <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/omentin/omentinBMI_putty_chrpos_filtered.txt")

omentin <- omentin %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(omentin) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_omentin <- ld_clump_local(
  dat = omentin,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_omentin, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/omentin_ld_clumped.txt")


# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - RESISTIN


#resistin <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
resistin <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/resistin/output/resistinBMI_putty_chrpos_filtered.txt")

resistin <- resistin %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(resistin) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_resistin <- ld_clump_local(
  dat = resistin,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_resistin, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/resistin_ld_clumped.txt")


# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - VASPIN


#vaspin <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
vaspin <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/adipokines/vaspin/output/vaspinBMI_putty_chrpos_filtered.txt")

vaspin <- vaspin %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(vaspin) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_vaspin <- ld_clump_local(
  dat = vaspin,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_vaspin, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/vaspin_ld_clumped.txt")


# ------------------------------------------------------------------------------------------------------------------------------------------------ # 
# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - PAT_che_2017


#pat <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
pat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/PAT_che_2017/PAT_che_2017_chrpos.txt")

pat <- pat %>%
  dplyr::select("rsID", "pval", "chr_b37", "pos_b37", "allele1", "allele2")

colnames(pat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_pat <- ld_clump_local(
  dat = pat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_pat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/pat_ld_clumped.txt")


# ************************************************************************************************************************************************ # 
# DATA IMPLEMENTATION - SAT_che_2017


#sat <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
sat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/SAT_che_2017/SAT_che_2017_chrpos.txt")

sat <- sat %>%
  dplyr::select("rsID", "pval", "chr_b37", "pos_b37", "allele1", "allele2")

colnames(sat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_sat <- ld_clump_local(
  dat = sat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_sat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/sat_ld_clumped.txt")

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - SATHU_che_2017


#sathu <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
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
# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - VATHU_che_2017


#vathu <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
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

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - VATSAT_che_2017


#vatsat <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
vatsat <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/VATSAT_che_2017/VATSAT_che_2017_chrpos.txt")

vatsat <- vatsat %>%
  dplyr::select("rsID", "pval", "chr_b37", "pos_b37", "allele1", "allele2")

colnames(vatsat) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_vatsat <- ld_clump_local(
  dat = vatsat,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_vatsat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/vatsat_ld_clumped.txt")

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION - Trunk Fat Ratio


#tfr <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoqBMI_putty_chrpos_filtered.txt")
tfr <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/trunk_fat_percentage/tfr_elsworth_filtered_chrpos.txt")

tfr <- tfr %>%
  dplyr::select("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(tfr) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_tfr <- ld_clump_local(
  dat = tfr,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_tfr, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/tfr_ld_clumped.txt")

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION -   Ant Thigh Muscle Fat Van Der Meer

atm <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ant_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.tsv")

atm <- atm %>%
  dplyr::select("variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(atm) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_atm <- ld_clump_local(
  dat = atm,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_atm, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/ant_thigh_muscle_fat_vanDerMeer_ld_clumped.txt")

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION -   Liver Fat Liu

livf <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/liver_fat_liu_gwas_summary_stat.tsv")

livf <- livf %>%
  dplyr::select("variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(livf) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_livf <- ld_clump_local(
  dat = livf,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_livf, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/liver_fat_liu_ld_clumped.txt")

# ************************************************************************************************************************************************ # 

# DATA IMPLEMENTATION -   Post Thigh Muscle Fat Van Der Meer

ptm <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/post_thigh_muscle_fat_vanDerMeer_gwas_summary_stat.tsv")

ptm <- ptm %>%
  dplyr::select("variant_id", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele")

colnames(ptm) <- c("rsid", "pval", "chr", "pos", "A1", "A2") # Required form for LD-clumping script

clumped_ptm <- ld_clump_local(
  dat = ptm,
  clump_kb = 1000,
  clump_r2 = 0.001,
  clump_p = 5e-8,
  bfile = "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/EUR",
  plink_bin = "/opt/software/plink/1.9.0/bin/plink"
)

fwrite(clumped_ptm, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/post_thigh_muscle_fat_vanDerMeer_ld_clumped.txt")

# ************************************************************************************************************************************************ # 
# ************************************************************************************************************************************************ # 
# ************************************************************************************************************************************************ # 
# ************************************************************************************************************************************************ # 
# ************************************************************************************************************************************************ # 
# ************************************************************************************************************************************************ # 


                                                            # VISUALIZATION SCRIPTS 


clumped <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/adipoq_ld_clumped.txt")

# MANHATTAN PLOT 


# orijinal veri: adipoq
# clumped SNP listesi
indep_snps <- clumped$rsid

adipoq <- adipoq %>%
  mutate(is_clumped = ifelse(rsid %in% indep_snps, "Clumped", "Other"))

p <- ggplot(adipoq, aes(x = pos, y = -log10(pval), color = is_clumped)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~chr, scales = "free_x") +
  scale_color_manual(values = c("Clumped" = "red", "Other" = "grey")) +
  theme_bw() +
  labs(title = "LD-clumped SNPs", y = "-log10(p-value)", x = "Position")

ggsave("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/adipoq_manhattan.png",
       plot = p, width = 12, height = 6, dpi = 300)

# ************************************************************************************************************************************************ # 

# P_VAL DISTRIBUTION 

p_hist <- ggplot(adipoq, aes(x = -log10(pval))) +
  geom_histogram(bins = 50, fill = "grey", alpha = 0.5) +
  geom_vline(data = adipoq %>% filter(rsid %in% indep_snps),
             aes(xintercept = -log10(pval)), color = "red", alpha = 0.7) +
  labs(title = "P-value distribution with clumped SNPs",
       x = "-log10(p)", y = "Count") +
  theme_minimal()

# PNG kaydet
ggsave("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/ld_clumping_studies/adipoq_pval_hist.png",
       plot = p_hist, width = 10, height = 6, dpi = 300)

# ************************************************************************************************************************************************ # 

# CLUMP BUYUKLUKLERI (BARPLOT)

clump_sizes <- clumped %>%
  count(clump)   # 'clump' kolonu varsa (her LD grubunu numaralar)

ggplot(clump_sizes, aes(x = factor(clump), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Number of SNPs per clump", x = "Clump ID", y = "SNP count") +
  theme_minimal()

# ************************************************************************************************************************************************ # 

# CHROMOSOME BASED SUMMARY 

clumped %>%
  group_by(chr) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = chr, y = n)) +
  geom_col(fill = "darkgreen") +
  labs(title = "Independent SNPs per chromosome", y = "Count") +
  theme_bw()

# ************************************************************************************************************************************************ # 
