##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for VAT/ASAT ratio from Agrawal et al, available in https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-30931-2/MediaObjects/41467_2022_30931_MOESM4_ESM.xlsx

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/0321_vatAsatRatio_bgen_stats.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#SNP CHR    BP GENPOS                ALLELE1 ALLELE0     A1FREQ     INFO CHISQ_LINREG P_LINREG         BETA        SE CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM P_BOLT_LMM
#1:                      rs367896724   1 10177      0                      A      AC 0.59767100 0.467935    0.0044993     0.95 -0.000789656 0.0101979         0.00599587           0.94      0.1128510       0.74
#2:                      rs201106462   1 10352      0                      T      TA 0.60568400 0.447895    0.0304250     0.86  0.000896248 0.0104593         0.00734260           0.93      0.0397492       0.84
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C   1 10616      0 CCGCCGTTGCAAAGGCGCGCCG       C 0.00595562 0.468098    1.1061100     0.29 -0.083886000 0.0653143         1.64954000           0.20      2.0965500       0.15
#4:                      rs575272151   1 11008      0                      C       G 0.91377100 0.495023    0.1461500     0.70 -0.004852580 0.0173645         0.07809460           0.78      0.0489216       0.82
#5:                      rs544419019   1 11012      0                      C       G 0.91377100 0.495023    0.1461500     0.70 -0.004852580 0.0173645         0.07809460           0.78      0.0489216       0.82
#6:                      rs540538026   1 13110      0                      G       A 0.94285300 0.391804    0.8781100     0.35 -0.024138800 0.0236645         1.04049000           0.31      1.2549900       0.26

#No need for Z-scores here, let's filter and clean:

raw_sumstat <- raw_sumstat %>%
  dplyr::select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value") 

raw_sumstat$sample_size=38965

#Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

head(raw_sumstat)

#variant chromosome base_pair_location          effect_allele other_allele effect_allele_frequency        beta standard_error p_value sample_size    chr_pos
#1:                      rs367896724          1              10177                      A           AC              0.59767100  0.00720764      0.0104009   0.490       38965 chr1:10177
#2:                      rs201106462          1              10352                      T           TA              0.60568400 -0.00122601      0.0106675   0.910       38965 chr1:10352
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C          1              10616 CCGCCGTTGCAAAGGCGCGCCG            C              0.00595562 -0.00450296      0.0666141   0.950       38965 chr1:10616
#4:                      rs575272151          1              11008                      C            G              0.91377100  0.00847971      0.0177101   0.630       38965 chr1:11008
#5:                      rs544419019          1              11012                      C            G              0.91377100  0.00847971      0.0177101   0.630       38965 chr1:11012
#6:                      rs540538026          1              13110                      G            A              0.94285300 -0.04448090      0.0241354   0.065       38965 chr1:13110

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/vat_asat_ratio_agrawal_full_sumstats.txt")
