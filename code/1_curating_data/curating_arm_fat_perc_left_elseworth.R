##############
#INTRODUCTION#
##############

#This code is to curate trunk fat % GWAS from Elseworth available in gwas.mrcieu: Links are only available for a certain amount of time. Download VCF here:https://opengwas.io/datasets/ukb-b-20188#

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)

###################
#Loading functions#
###################

beta_retriever <- function(header_info){
  
  beta <- strsplit(header_info, ":")[[1]][1]
  
  return(beta)
  
}

se_retriever <- function(header_info){
  
  se <- strsplit(header_info, ":")[[1]][2]
  
  return(se)
  
}

logp_retriever <- function(header_info){
  
  logp <- strsplit(header_info, ":")[[1]][3]
  
  return(logp)
  
}

EAF_retriever <- function(header_info){
  
  EAF <- strsplit(header_info, ":")[[1]][4]
  
  return(EAF)
  
}


rsid_retriever <- function(header_info){
  
  rsid <- strsplit(header_info, ":")[[1]][5]
  
  return(rsid)
  
}

##############
#Loading data#
##############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/"

setwd(path_2_input)

raw_data<-read.table("input/raw_full_sumstats/ukb-b-20188.vcf.gz", stringsAsFactors = FALSE)

#We are also gonna load the genome-wide significant SNPs to be able to see if the data
#matches and check whether we have weird results or not with this format.

#######################
#REFORMATTING THE DATA#
#######################

#Now that we know the basics of the data.
#Let's start by making some functions to make it work:

raw_data_copy <- raw_data

colnames(raw_data_copy) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure", "dot", "PASS", "AF", "Headers", "Headers_info")

#Let's see if we can just use the functions from BFP...

check <- which(duplicated(raw_data_copy$Headers) == FALSE) #YES they are the same

#Now we are gonna make functions to retrieve the info on headers info:

raw_data_copy$beta.exposure <- as.numeric(as.character(unlist(sapply(raw_data_copy$Headers_info, beta_retriever))))
raw_data_copy$se.exposure <- as.numeric(as.character(unlist(sapply(raw_data_copy$Headers_info, se_retriever))))
raw_data_copy$logp <- as.numeric(as.character(unlist(sapply(raw_data_copy$Headers_info, logp_retriever))))
raw_data_copy$eaf.exposure <- as.numeric(as.character(unlist(sapply(raw_data_copy$Headers_info, EAF_retriever))))
raw_data_copy$N <- 	454724 #no need to filter by N cuz all of them have the same
raw_data_copy$rsid <- as.character(unlist(sapply(raw_data_copy$Headers_info, rsid_retriever)))

#To convert to p-value I need to transform them again:

raw_data_copy$logp_neg <- raw_data_copy$logp*(-1)

#And thus...

raw_data_copy$pval.exposure <- 10^(raw_data_copy$logp_neg)

############
#FINISHED!!#
############

##########
#CURATION#
##########

#We don't have INFO, so we cannot rely on that. 
#1. But we can remove those variants that are MAF < 0.01

raw_data_maf <- raw_data_copy[which(as.numeric(raw_data_copy$eaf.exposure) > 0.01),]
raw_data_maf <- raw_data_maf[which(as.numeric(raw_data_maf$eaf.exposure) < 0.99),]

summary(raw_data_maf$eaf.exposure)

#We go from 9.2M to 7.7M. Good enough!!

#2. Now we remove the MHC region:

raw_data_maf_mhc <- raw_data_maf[which(as.numeric(raw_data_maf$chr.exposure) == 6 & as.numeric(raw_data_maf$pos.exposure) >= 26000000 & as.numeric(raw_data_maf$pos.exposure) <= 34000000),]

summary(as.numeric(raw_data_maf_mhc$chr.exposure)) #perfect.
summary(as.numeric(raw_data_maf_mhc$pos.exposure)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(raw_data_maf_mhc$pval.exposure)) #we actually do. But we said we should remove it so...

raw_data_maf_no_mhc <- raw_data_maf[which(!(raw_data_maf$SNP%in%raw_data_maf_mhc$SNP)),]

#3. We should have done that before. But let's check whether we are in build 37. I already know that it is, hence why I just check it here.

head(raw_data_maf_no_mhc)

#Perfect.

#4. Let's check the RSIDs:

raw_data_maf_no_mhc <- raw_data_maf_no_mhc[order(raw_data_maf_no_mhc$SNP),]

head(raw_data_maf_no_mhc)

#Some of the SNPs do not have rsID, if we have mismatches in the future, we will try to correct them, but so far there is no need to do it ;)

#Let's get chr_pos:

raw_data_maf_no_mhc$chr_pos <- paste("chr", raw_data_maf_no_mhc$chr.exposure, ":", raw_data_maf_no_mhc$pos.exposure, sep = "")

data_end <- raw_data_maf_no_mhc %>%
  select(rsid, chr.exposure, pos.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, N, chr_pos)

colnames(data_end) <-c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(data_end, "output/1_curated_data/arm_fat_perc_left_elseworth_full_sumstats.txt")
