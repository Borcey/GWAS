##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for VASPINadjBMI from Sanz-Martinez et al, available in in-house

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
  
  tmp_df= fread(paste("input/raw_full_sumstats/Step2_combined_serpina12_BMI_Chr", chr, "_pheno.regenie.gz", sep = ""))
  
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
#1:     1  10177 rs367896724      AC       A 0.602481 0.464428 42412  ADD  0.01185090 0.00825525 2.0608400 0.8206580    NA
#2:     1  10352 rs201106462      TA       T 0.608194 0.443430 42412  ADD -0.00158112 0.00844675 0.0350387 0.0698079    NA
#3:     1  11008 rs575272151       G       C 0.914448 0.489245 42412  ADD -0.00401040 0.01415520 0.0802679 0.1096150    NA
#4:     1  11012 rs544419019       G       C 0.914448 0.489245 42412  ADD -0.00401040 0.01415520 0.0802679 0.1096150    NA
#5:     1  13110 rs540538026       A       G 0.940471 0.399932 42412  ADD  0.00915804 0.01843610 0.2467540 0.2080500    NA
#6:     1  13116  rs62635286       G       T 0.809993 0.404114 42412  ADD  0.01774130 0.01106250 2.5719400 0.9634740    NA

#Let's convert to p-value:

raw_sumstat$P=10^(as.numeric(raw_sumstat$LOG10P)*(-1)) #NAs are added because in the merging, someone added the column names!

raw_sumstat <- raw_sumstat[which(!(is.na(raw_sumstat$P))),]

raw_sumstat <- raw_sumstat %>%
  dplyr::select(ID, CHROM, GENPOS, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size") 

#raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "") #the data is too big, let's do it later ;)

tail(raw_sumstat) ##Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

#variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency        beta standard_error   p_value sample_size
#1: rs200507571         22           51236013             A           AT                0.748715  0.00460783     0.00706252 0.5141218       42412
#2:   rs3896457         22           51237063             T            C                0.704017  0.00513554     0.00651458 0.4305127       42412
#3: rs200607599         22           51237364             A            G                0.984941 -0.00446578     0.03122390 0.8862707       42412
#4: rs370652263         22           51237712             G            A                0.944576 -0.00914042     0.01303500 0.4831656       42412
#5: rs202228854         22           51240820             C            T                0.973526  0.03590040     0.02299890 0.1185340       42412
#6: rs575160859         22           51244237             C            T                0.986975 -0.01274940     0.03470110 0.7133144       42412

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/vaspinadjbmi_sanz_martinez_full_sumstats.txt")
