##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for ADIPOQadjBMI from Sanz-Martinez et al, available in in-house

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

for(chr in seq(1, 22)){
  
  tmp_df= fread(paste("input/raw_full_sumstats/Step2_combined_retn_BMI_Chr", chr, "_pheno.regenie.gz", sep = ""))
  
  if(!(exists("raw_sumstat"))){
    
    raw_sumstat=tmp_df
    
  } else {
    
    raw_sumstat <- rbind(raw_sumstat, tmp_df)
    
  }
  
  
}

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#CHROM GENPOS          ID ALLELE0 ALLELE1   A1FREQ     INFO     N TEST         BETA         SE       CHISQ     LOG10P EXTRA
#1:     1  10177 rs367896724      AC       A 0.602445 0.464697 43576  ADD -0.007453660 0.00926067 0.647820000 0.37582800    NA
#2:     1  10352 rs201106462      TA       T 0.608107 0.443765 43576  ADD -0.000175001 0.00946458 0.000341884 0.00645449    NA
#3:     1  11008 rs575272151       G       C 0.914038 0.490110 43576  ADD  0.005209770 0.01583270 0.108274000 0.12952700    NA
#4:     1  11012 rs544419019       G       C 0.914038 0.490110 43576  ADD  0.005209770 0.01583270 0.108274000 0.12952700    NA
#5:     1  13110 rs540538026       A       G 0.940524 0.398849 43576  ADD  0.016415700 0.02072260 0.627522000 0.36828700    NA
#6:     1  13116  rs62635286       G       T 0.810001 0.404093 43576  ADD  0.008556440 0.01242030 0.474590000 0.30902200    NA

#Let's convert to p-value:

raw_sumstat$P=10^(as.numeric(raw_sumstat$LOG10P)*(-1)) #NAs are added because in the merging, someone added the column names!

raw_sumstat <- raw_sumstat[which(!(is.na(raw_sumstat$P))),]

raw_sumstat <- raw_sumstat %>%
  dplyr::select(ID, CHROM, GENPOS, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size") 

#raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "") #the data is too big, let's do it later ;)

tail(raw_sumstat) ##Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

#variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error   p_value sample_size
#1: rs200507571         22           51236013             A           AT                0.748568  0.007043750     0.00791453 0.3734780       43576
#2:   rs3896457         22           51237063             T            C                0.703952  0.005357460     0.00730165 0.4631120       43576
#3: rs200607599         22           51237364             A            G                0.984923  0.021623900     0.03493690 0.5359546       43576
#4: rs370652263         22           51237712             G            A                0.944619 -0.003917110     0.01459480 0.7883989       43576
#5: rs202228854         22           51240820             C            T                0.973479  0.000699635     0.02574650 0.9783209       43576
#6: rs575160859         22           51244237             C            T                0.986963  0.024074100     0.03890410 0.5360435       43576

#######################################################################
#STEP 2: remove insertions and deletions, rare variants and MHC region#
#######################################################################

#First let's remove indels:

yes_vect <- c("A", "G", "T", "C")

raw_sumstat_no_indel <- raw_sumstat[which(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect),] 
#raw_sumstat_indel <- raw_sumstat[which(!(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect)),] #it seems we don't have any

#Second let's remove rare variant:

summary(as.numeric(raw_sumstat_no_indel$effect_allele_frequency))

raw_sumstat_no_indel_common = raw_sumstat_no_indel[which(as.numeric(raw_sumstat_no_indel$effect_allele_frequency) > 0.01 & as.numeric(raw_sumstat_no_indel$effect_allele_frequency) < 0.99),]

summary(as.numeric(raw_sumstat_no_indel_common$effect_allele_frequency))

raw_sumstat_no_indel_common$chr_pos = paste("chr", raw_sumstat_no_indel_common$chromosome, ":", raw_sumstat_no_indel_common$base_pair_location, sep = "")

#Third, let's remove extended MHC region

raw_sumstat_no_indel_common_mhc = raw_sumstat_no_indel_common[which(raw_sumstat_no_indel_common$chromosome == 6 & as.numeric(raw_sumstat_no_indel_common$base_pair_location) >= 26000000 & as.numeric(raw_sumstat_no_indel_common$base_pair_location) <= 34000000),]

summary(as.numeric(raw_sumstat_no_indel_common_mhc$chromosome)) #perfect
summary(as.numeric(raw_sumstat_no_indel_common_mhc$base_pair_location)) #perfect

raw_sumstat_no_indel_common_mhc_removed = raw_sumstat_no_indel_common[which(!(raw_sumstat_no_indel_common$chr_pos%in%raw_sumstat_no_indel_common_mhc$chr_pos)),]

#Let's be one sure

length(raw_sumstat_no_indel_common$chr_pos) - length(raw_sumstat_no_indel_common_mhc_removed$chr_pos) #perfect

####################################
#STEP 3: save the data! We are done#
####################################

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/retnadjbmi_sanz_martinez_full_sumstats.txt")
