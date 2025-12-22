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
  
  tmp_df= fread(paste("input/raw_full_sumstats/Step2_combined_adipoq_BMI_Chr", chr, "_pheno.regenie.gz", sep = ""))
  
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

#CHROM GENPOS          ID ALLELE0 ALLELE1   A1FREQ     INFO     N TEST        BETA         SE    CHISQ   LOG10P EXTRA
#1:     1  10177 rs367896724      AC       A 0.602619 0.463354 37635  ADD  0.00417136 0.00961220 0.188326 0.177627    NA
#2:     1  10352 rs201106462      TA       T 0.608492 0.442322 37635  ADD  0.01356450 0.00982751 1.905110 0.775965    NA
#3:     1  11008 rs575272151       G       C 0.914363 0.489477 37635  ADD -0.01726660 0.01643140 1.104240 0.532631    NA
#4:     1  11012 rs544419019       G       C 0.914363 0.489477 37635  ADD -0.01726660 0.01643140 1.104240 0.532631    NA
#5:     1  13110 rs540538026       A       G 0.940567 0.397955 37635  ADD -0.01392120 0.02152100 0.418434 0.285906    NA
#6:     1  13116  rs62635286       G       T 0.810550 0.402595 37635  ADD  0.00482195 0.01290390 0.139638 0.149574    NA

#Let's convert to p-value:

raw_sumstat$P=10^(as.numeric(raw_sumstat$LOG10P)*(-1)) #NAs are added because in the merging, someone added the column names!

raw_sumstat <- raw_sumstat[which(!(is.na(raw_sumstat$P))),]

raw_sumstat <- raw_sumstat %>%
  dplyr::select(ID, CHROM, GENPOS, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, N)

colnames(raw_sumstat) = c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size") 

#raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "") #the data is too big, let's do it later ;)

tail(raw_sumstat) ##Let's also add chr_pos. Just checked in chr22 that variants are in build 37!

#variant chromosome base_pair_location effect_allele other_allele effect_allele_frequency        beta standard_error    p_value sample_size
#1: rs200507571         22           51236013             A           AT                0.748160  0.00107655     0.00819837 0.89552827       37635
#2:   rs3896457         22           51237063             T            C                0.704955  0.00147895     0.00756054 0.84491204       37635
#3: rs200607599         22           51237364             A            G                0.984869 -0.06661970     0.03616250 0.06544101       37635
#4: rs370652263         22           51237712             G            A                0.945033  0.02137080     0.01517330 0.15899812       37635
#5: rs202228854         22           51240820             C            T                0.973445 -0.00558273     0.02671850 0.83449051       37635
#6: rs575160859         22           51244237             C            T                0.986920  0.02427390     0.04050680 0.54900336       37635

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

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/adipoqadjbmi_sanz_martinez_full_sumstats.txt")
