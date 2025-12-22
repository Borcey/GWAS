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
  
  tmp_df= fread(paste("input/raw_full_sumstats/Step2_combined_lep_BMI_Chr", chr, "_pheno.regenie.gz", sep = ""))
  
  if(!(exists("raw_sumstat"))){
    
    raw_sumstat=tmp_df
    
  } else {
    
    raw_sumstat <- rbind(raw_sumstat, tmp_df)
    
  }
  
  
}

unique(raw_sumstat$CHROM)

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

#CHROM GENPOS          ID ALLELE0 ALLELE1   A1FREQ     INFO     N TEST        BETA         SE     CHISQ    LOG10P EXTRA
#1:     1  10177 rs367896724      AC       A 0.602595 0.464567 43673  ADD  0.00647817 0.00642152 1.0177200 0.5043720    NA
#2:     1  10352 rs201106462      TA       T 0.608084 0.443499 43673  ADD -0.01020670 0.00656607 2.4163600 0.9205510    NA
#3:     1  11008 rs575272151       G       C 0.914360 0.489533 43673  ADD -0.01744970 0.01099950 2.5166800 0.9482780    NA
#4:     1  11012 rs544419019       G       C 0.914360 0.489533 43673  ADD -0.01744970 0.01099950 2.5166800 0.9482780    NA
#5:     1  13110 rs540538026       A       G 0.940493 0.399897 43673  ADD  0.00217175 0.01435060 0.0229023 0.0556598    NA
#6:     1  13116  rs62635286       G       T 0.810185 0.404106 43673  ADD  0.00445858 0.00861665 0.2677410 0.2183520    NA

#Let's convert to p-value:

raw_sumstat$P=10^(as.numeric(raw_sumstat$LOG10P)*(-1)) #NAs are added because in the merging, someone added the column names!

raw_sumstat <- raw_sumstat[which(!(is.na(raw_sumstat$P))),]

raw_sumstat <- raw_sumstat %>%
  dplyr::select(ID, CHROM, GENPOS, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size") 

#raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "") #the data is too big, let's do it later ;)

tail(raw_sumstat) ##Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

#variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error   p_value sample_size
#1: rs200507571         22           51236013             A           AT                0.748539  0.004503530     0.00548493 0.4116044       43673
#2:   rs3896457         22           51237063             T            C                0.704109 -0.000870815     0.00506144 0.8633991       43673
#3: rs200607599         22           51237364             A            G                0.984918  0.019901300     0.02417580 0.4103997       43673
#4: rs370652263         22           51237712             G            A                0.944552  0.000253895     0.01011310 0.9799708       43673
#5: rs202228854         22           51240820             C            T                0.973493 -0.010674000     0.01785890 0.5500523       43673
#6: rs575160859         22           51244237             C            T                0.986969  0.041259400     0.02694420 0.1256985       43673

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/lepadjbmi_sanz_martinez_full_sumstats.txt")
