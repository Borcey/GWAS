##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for Epicardial+Pericardial adipose tissue area from Ramo et al, available in https://api.kpndataregistry.org/api/d/DoNNBJ

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/EPAT_area_regenie_GWAS_summary_statistics.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#CHROM GENPOS          ID ALLELE0 ALLELE1    A1FREQ     INFO     N TEST       BETA       SE     CHISQ    LOG10P EXTRA        P
#1:     1  10177 rs367896724       A      AC 0.4021310 0.467607 41494  ADD  0.0222519 0.106077 0.0440036 0.0789137    NA 0.833847
#2:     1  10352 rs201106462       T      TA 0.3947070 0.445780 41494  ADD -0.1423240 0.109114 1.7013600 0.7164480    NA 0.192111
#3:     1  11008 rs575272151       C       G 0.0858994 0.488697 41494  ADD  0.2147180 0.181543 1.3988800 0.6254130    NA 0.236912
#4:     1  11012 rs544419019       C       G 0.0858994 0.488697 41494  ADD  0.2147180 0.181543 1.3988800 0.6254130    NA 0.236912
#5:     1  13110 rs540538026       G       A 0.0568116 0.383327 41494  ADD  0.2557710 0.248403 1.0602000 0.5183160    NA 0.303168
#6:     1  13116  rs62635286       T       G 0.1863530 0.402023 41494  ADD  0.1590480 0.144324 1.2144600 0.5679100    NA 0.270452

#No need for Z-scores here, let's filter and clean:

raw_sumstat <- raw_sumstat %>%
  dplyr::select(ID, CHROM, GENPOS, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size") 

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

tail(raw_sumstat) ##Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

#variant chromosome base_pair_location          effect_allele other_allele effect_allele_frequency        beta standard_error p_value sample_size    chr_pos
#1:                      rs367896724          1              10177                      A           AC              0.59767100  0.00169645      0.0103367    0.87       38965 chr1:10177
#2:                      rs201106462          1              10352                      T           TA              0.60568400 -0.00119692      0.0106017    0.91       38965 chr1:10352
#3: 1:10616_CCGCCGTTGCAAAGGCGCGCCG_C          1              10616 CCGCCGTTGCAAAGGCGCGCCG            C              0.00595562 -0.06775070      0.0662033    0.31       38965 chr1:10616
#4:                      rs575272151          1              11008                      C            G              0.91377100  0.00348227      0.0176009    0.84       38965 chr1:11008
#5:                      rs544419019          1              11012                      C            G              0.91377100  0.00348227      0.0176009    0.84       38965 chr1:11012
#6:                      rs540538026          1              13110                      G            A              0.94285300 -0.03892340      0.0239866    0.10       38965 chr1:13110

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/epat_area_ramo_full_sumstats.txt")
