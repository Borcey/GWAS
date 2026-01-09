##############
#INTRODUCTION#
##############

#Let's clean the full GWAS summary statistics for T2D from Suzuki.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/Suzuki.Nature2024.T2DGGI.EUR.sumstats/EUR_Metal_LDSC-CORR_Neff.v2.txt")

#Let's clean data:

colnames(raw_sumstat)=c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "sample_cases", "sample_controls", "sample_size")

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

#Is this build 37:

check=raw_sumstat[which(as.numeric(raw_sumstat$p_value) < 5e-08),]
check=check[order(as.numeric(check$p_value)),] #YESSS they are

#######################################################################
#STEP 2: remove insertions and deletions, rare variants and MHC region#
#######################################################################

#First let's remove indels:

yes_vect <- c("A", "G", "T", "C")

raw_sumstat_no_indel <- raw_sumstat[which(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect),] #we go from 9M to 8M
#raw_sumstat_indel <- raw_sumstat[which(!(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect)),] #we did it correctly

#Second let's remove rare variant:

summary(raw_sumstat_no_indel$effect_allele_frequency)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0006  0.2228  0.4341  0.9883  1.0000 

raw_sumstat_no_indel_common = raw_sumstat_no_indel[which(as.numeric(raw_sumstat_no_indel$effect_allele_frequency) > 0.01 & as.numeric(raw_sumstat_no_indel$effect_allele_frequency) < 0.99),]

summary(raw_sumstat_no_indel_common$effect_allele_frequency)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0101  0.0865  0.4027  0.4565  0.8356  0.9899 

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/t2d_suzuki_full_sumstats.txt")
