##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for Anterior Thigh fat infiltration from Van der meer et al, available in https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267351/GCST90267351.tsv.gz

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/GCST90267351.tsv.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

# variant_id   p_value chromosome base_pair_location effect_allele other_allele     z_value         beta standard_error     n
# 1: rs533090414 0.7950760          1              18849             G            C -0.25972492 -0.013728669     0.05285850 30951
# 2:  rs28619159 0.7620750          1              91421             T            C -0.30275710 -0.016391647     0.05414125 31755
# 3: rs539032812 0.9217795          1             665266             T            C -0.09819242 -0.005126519     0.05220891 31943
# 4:  rs12238997 0.7580255          1             693731             A            G -0.30807470 -0.004255064     0.01381179 30771
# 5: rs200531508 0.7233943          1             711310             G            A -0.35392600 -0.008859894     0.02503318 30394
# 6: rs373285923 0.7973986          1             712547             G            C  0.25671530  0.013606740     0.05300323 31327

#No need for Z-scores here, let's filter and clean:

raw_sumstat <- raw_sumstat %>%
  dplyr::select(variant_id, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value, n)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "sample_size") 

#We do not have effect_allele_frequencies, we will have to add them afterwards

#Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

head(raw_sumstat)

# variant chromosome base_pair_location effect_allele other_allele         beta standard_error   p_value sample_size     chr_pos
# 1: rs533090414          1              18849             G            C -0.013728669     0.05285850 0.7950760       30951  chr1:18849
# 2:  rs28619159          1              91421             T            C -0.016391647     0.05414125 0.7620750       31755  chr1:91421
# 3: rs539032812          1             665266             T            C -0.005126519     0.05220891 0.9217795       31943 chr1:665266
# 4:  rs12238997          1             693731             A            G -0.004255064     0.01381179 0.7580255       30771 chr1:693731
# 5: rs200531508          1             711310             G            A -0.008859894     0.02503318 0.7233943       30394 chr1:711310
# 6: rs373285923          1             712547             G            C  0.013606740     0.05300323 0.7973986       31327 chr1:712547

#######################################################################
#STEP 2: remove insertions and deletions, rare variants and MHC region#
#######################################################################

#First let's remove indels:

yes_vect <- c("A", "G", "T", "C")

raw_sumstat_no_indel <- raw_sumstat[which(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect),] #it seems we don't have any
#raw_sumstat_indel <- raw_sumstat[which(!(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect)),] #it seems we don't have any

#Second let's remove rare variant:

raw_sumstat_no_indel$effect_allele_frequency = NA
summary(raw_sumstat_no_indel$effect_allele_frequency)

#raw_sumstat_no_indel_common = raw_sumstat_no_indel[which(as.numeric(raw_sumstat_no_indel$effect_allele_frequency) > 0.01 & as.numeric(raw_sumstat_no_indel$effect_allele_frequency) < 0.99),]

#summary(raw_sumstat_no_indel_common$effect_allele_frequency)

#Third, let's remove extended MHC region

raw_sumstat_no_indel_common_mhc = raw_sumstat_no_indel[which(raw_sumstat_no_indel$chromosome == 6 & raw_sumstat_no_indel$base_pair_location >= 26000000 & raw_sumstat_no_indel$base_pair_location <= 34000000),]

summary(raw_sumstat_no_indel_common_mhc$chromosome) #perfect
summary(raw_sumstat_no_indel_common_mhc$base_pair_location) #perfect

raw_sumstat_no_indel_common_mhc_removed = raw_sumstat_no_indel[which(!(raw_sumstat_no_indel$chr_pos%in%raw_sumstat_no_indel_common_mhc$chr_pos)),]

#Let's be one sure

length(raw_sumstat_no_indel$chr_pos) - length(raw_sumstat_no_indel_common_mhc_removed$chr_pos) #perfect

####################################
#STEP 3: save the data! We are done#
####################################

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/anterior_thigh_fat_infiltration_van_der_meer_full_sumstats.txt")
