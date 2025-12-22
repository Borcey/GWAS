##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for ASATadjBMI from Agrawal et al, available in https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-30931-2/MediaObjects/41467_2022_30931_MOESM4_ESM.xlsx

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/0321_asatadjbmi3_bgen_stats.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#SNP CHR    BP GENPOS                ALLELE1 ALLELE0     A1FREQ     INFO CHISQ_LINREG P_LINREG         BETA        SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM P_BOLT_LMM
#1:                      rs367896724   1 10177      0                      A      AC 0.59808500 0.467935    0.2940910     0.59  0.006999410 0.0106559        4.31462e-01           0.51      0.3837400       0.54
#2:                      rs201106462   1 10352      0                      T      TA 0.60547500 0.447895    0.4267010     0.51 -0.003384030 0.0109373        9.57297e-02           0.76      0.1711800       0.68
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C   1 10616      0 CCGCCGTTGCAAAGGCGCGCCG       C 0.00591215 0.468098    0.0131984     0.91  0.026051900 0.0685064        1.44616e-01           0.70      0.0320796       0.86
#4:                      rs575272151   1 11008      0                      C       G 0.91378400 0.495023    0.1749780     0.68 -0.001726710 0.0181481        9.05269e-03           0.92      0.0102343       0.92
#5:                      rs544419019   1 11012      0                      C       G 0.91378400 0.495023    0.1749780     0.68 -0.001726710 0.0181481        9.05269e-03           0.92      0.0102343       0.92
#6:                      rs540538026   1 13110      0                      G       A 0.94266500 0.391804    0.0389807     0.84  0.000153695 0.0246552        3.88602e-05           1.00      0.0536600       0.82

#No need for Z-scores here, let's filter and clean:

raw_sumstat <- raw_sumstat %>%
  dplyr::select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value") 

raw_sumstat$sample_size=37641

#Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

head(raw_sumstat)

#variant chromosome base_pair_location          effect_allele effect_allele_frequency other_allele         beta standard_error p_value sample_size    chr_pos
#1:                      rs367896724          1              10177                      A                      AC   0.59808500  0.006999410      0.0106559    0.51       37641 chr1:10177
#2:                      rs201106462          1              10352                      T                      TA   0.60547500 -0.003384030      0.0109373    0.76       37641 chr1:10352
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C          1              10616 CCGCCGTTGCAAAGGCGCGCCG                       C   0.00591215  0.026051900      0.0685064    0.70       37641 chr1:10616
#4:                      rs575272151          1              11008                      C                       G   0.91378400 -0.001726710      0.0181481    0.92       37641 chr1:11008
#5:                      rs544419019          1              11012                      C                       G   0.91378400 -0.001726710      0.0181481    0.92       37641 chr1:11012
#6:                      rs540538026          1              13110                      G                       A   0.94266500  0.000153695      0.0246552    1.00       37641 chr1:13110

#######################################################################
#STEP 2: remove insertions and deletions, rare variants and MHC region#
#######################################################################

#First let's remove indels:

yes_vect <- c("A", "G", "T", "C")

raw_sumstat_no_indel <- raw_sumstat[which(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect),] 
#raw_sumstat_indel <- raw_sumstat[which(!(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect)),] #it seems we don't have any

#Second let's remove rare variant:

summary(raw_sumstat_no_indel$effect_allele_frequency)

raw_sumstat_no_indel_common = raw_sumstat_no_indel[which(as.numeric(raw_sumstat_no_indel$effect_allele_frequency) > 0.01 & as.numeric(raw_sumstat_no_indel$effect_allele_frequency) < 0.99),]

summary(raw_sumstat_no_indel_common$effect_allele_frequency)

#Third, let's remove extended MHC region

raw_sumstat_no_indel_common_mhc = raw_sumstat_no_indel_common[which(raw_sumstat_no_indel_common$chromosome == 6 & raw_sumstat_no_indel_common$base_pair_location >= 26000000 & raw_sumstat_no_indel_common$base_pair_location <= 34000000),]

summary(raw_sumstat_no_indel_common_mhc$chromosome) #perfect
summary(raw_sumstat_no_indel_common_mhc$base_pair_location) #perfect

raw_sumstat_no_indel_common_mhc_removed = raw_sumstat_no_indel_common[which(!(raw_sumstat_no_indel_common$chr_pos%in%raw_sumstat_no_indel_common_mhc$chr_pos)),]

#Let's be one sure

length(raw_sumstat_no_indel_common$chr_pos) - length(raw_sumstat_no_indel_common_mhc_removed$chr_pos) #perfect

####################################
#STEP 3: save the data! We are done#
####################################

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/asatadjbmi_agrawal_full_sumstats.txt")
