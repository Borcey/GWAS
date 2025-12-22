##############
#INTRODUCTION#
##############

#This code performs ld_clump_local to UKBB data to get the leads for our analyses on trunk fat %

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(ieugwasr)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

clean_data<- fread("output/1_curated_data/leg_fat_perc_right_elseworth_full_sumstats.txt")

#Let's make a df that has rsID and pval as columns:

data_4_clumping = clean_data
data_4_clumping$rsid = data_4_clumping$variant
data_4_clumping$pval = data_4_clumping$p_value

#For triallelic SNPs, let's just get the combo of alleles with the best p-value:

data_4_clumping <- data_4_clumping[order(as.numeric(data_4_clumping$pval)),]
data_4_clumping_nondupl  <- data_4_clumping[which(duplicated(data_4_clumping$chr_pos) == FALSE),]

#We just need the genome-wide significant hits as leads:

data_4_clumping_nondupl_gw <- data_4_clumping_nondupl[which(as.numeric(data_4_clumping_nondupl$pval) < 5e-08),] #48610

###########################
#Let's perform LD-clumping#
###########################

clumped_data = ieugwasr::ld_clump_local(dat=data_4_clumping_nondupl_gw, clump_kb = 1000, clump_r2 = 0.01, bfile = "C:/Users/zlc436/Desktop/1kg.v3/EUR", plink_bin = "C:/Users/zlc436/Documents/R/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", clump_p = 5e-08)

#Let's get the data in the format we want:

leads_df = data_4_clumping_nondupl[which(data_4_clumping_nondupl$variant%in%clumped_data$variant),] #perfect match

leads_df <- leads_df %>%
  dplyr::select(-c("rsid", "pval"))

#####################
#Let's save the data#
#####################

fwrite(leads_df, "output/1_curated_data/leg_fat_perc_right_elseworth_leads.txt")
