##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for gsatadjBMI from Agrawal et al, available in https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-30931-2/MediaObjects/41467_2022_30931_MOESM4_ESM.xlsx

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/0321_gfatadjbmi3_bgen_stats.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#SNP CHR    BP GENPOS                ALLELE1 ALLELE0     A1FREQ     INFO CHISQ_LINREG P_LINREG        BETA        SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM P_BOLT_LMM
#1:                      rs367896724   1 10177      0                      A      AC 0.59808500 0.467935   0.02917800     0.86 -0.00343743 0.0105930          0.1053000           0.75     0.05614190       0.81
#2:                      rs201106462   1 10352      0                      T      TA 0.60547500 0.447895   1.58531000     0.21  0.01484500 0.0108728          1.8641400           0.17     1.33053000       0.25
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C   1 10616      0 CCGCCGTTGCAAAGGCGCGCCG       C 0.00591215 0.468098   0.00221426     0.96  0.00975707 0.0681023          0.0205265           0.89     0.00435309       0.95
#4:                      rs575272151   1 11008      0                      C       G 0.91378400 0.495023   0.23540000     0.63 -0.01189670 0.0180410          0.4348430           0.51     0.66165700       0.42
#5:                      rs544419019   1 11012      0                      C       G 0.91378400 0.495023   0.23540000     0.63 -0.01189670 0.0180410          0.4348430           0.51     0.66165700       0.42
#6:                      rs540538026   1 13110      0                      G       A 0.94266500 0.391804   0.83152400     0.36  0.02239780 0.0245097          0.8350890           0.36     1.25641000       0.26

#No need for Z-scores here, let's filter and clean:

raw_sumstat <- raw_sumstat %>%
  dplyr::select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value") 

raw_sumstat$sample_size=37641

#Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

head(raw_sumstat)

#variant chromosome base_pair_location          effect_allele other_allele effect_allele_frequency        beta standard_error p_value sample_size    chr_pos
#1:                      rs367896724          1              10177                      A           AC              0.59808500 -0.00343743      0.0105930    0.75       37641 chr1:10177
#2:                      rs201106462          1              10352                      T           TA              0.60547500  0.01484500      0.0108728    0.17       37641 chr1:10352
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C          1              10616 CCGCCGTTGCAAAGGCGCGCCG            C              0.00591215  0.00975707      0.0681023    0.89       37641 chr1:10616
#4:                      rs575272151          1              11008                      C            G              0.91378400 -0.01189670      0.0180410    0.51       37641 chr1:11008
#5:                      rs544419019          1              11012                      C            G              0.91378400 -0.01189670      0.0180410    0.51       37641 chr1:11012
#6:                      rs540538026          1              13110                      G            A              0.94266500  0.02239780      0.0245097    0.36       37641 chr1:13110

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/gsatadjbmi_agrawal_full_sumstats.txt")
