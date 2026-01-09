# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

#0.1 special function to get data:

effect_allele_parser = function(alleles_){
  
  #STEP 1: let's split:
  
  a1_vect=as.character(unlist(str_split(alleles_, "/")))
  
  a1_data=a1_vect[1]
  
  a1_data=as.character(unlist(str_remove(a1_data, " ")))
  
  return(a1_data)
  
}

other_allele_parser = function(alleles_){
  
  #STEP 1: let's split:
  
  a2_vect=as.character(unlist(str_split(alleles_, "/")))
  
  a2_data=a2_vect[2]
  
  a2_data=as.character(unlist(str_remove(a2_data, " ")))
  
  return(a2_data)
  
}

chr_parser = function(chr_pos){
  
  #STEP 1: let's split:
  
  chr_vect=as.character(unlist(str_split(chr_pos, ":")))
  
  chr_vect=chr_vect[1]
  
  chr_clean=as.character(unlist(str_split(chr_vect, "chr")))
  
  chr_clean=chr_clean[2]
  
  return(as.numeric(chr_clean))
  
}

pos_parser = function(chr_pos){
  
  #STEP 1: let's split:
  
  chr_vect=as.character(unlist(str_split(chr_pos, ":")))
  
  pos_clean=chr_vect[2]
  
  return(as.numeric(pos_clean))
  
}

rsid_parser = function(rsid){
  
  #STEP 1: let's split:
  
  rsid_vect=as.character(unlist(str_split(rsid, "[*]")))
  
  rsid_clean=rsid_vect[1]
  
  return(rsid_clean)
  
}

# ---------------------------------------------------------------------------------- #

                        # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 -  Data Fetch 

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here

ir_lotta <- readxl::read_xlsx("output/1_curated_data/multi_trait_ir_lotta_leads.xlsx")

# ---------------------------------------------------------------------------------- #

#Let's clean the data:

ir_lotta=ir_lotta[-c(1,12,13, 36,37,38),]

ir_lotta$chr_pos=tolower(ir_lotta$`Genomic coordinate`)

#Let's get effect and other allele:

ir_lotta$effect_allele=as.character(unlist(sapply(ir_lotta$`Alleles (effect / other)`, effect_allele_parser)))
ir_lotta$other_allele=as.character(unlist(sapply(ir_lotta$`Alleles (effect / other)`, other_allele_parser)))

ir_lotta$`Beta FIadjBMI per allele` <- as.numeric(ir_lotta$`Beta FIadjBMI per allele`)
ir_lotta$`Beta HDL cholesterol per allele` <- as.numeric(ir_lotta$`Beta HDL cholesterol per allele`)
ir_lotta$`Beta triglycerides per allele` <- as.numeric(ir_lotta$`Beta triglycerides per allele`)

ir_lotta$`FIadjBMI p-value` <- as.numeric(ir_lotta$`FIadjBMI p-value`)
ir_lotta$`HDL cholesterol p-value` <- as.numeric(ir_lotta$`HDL cholesterol p-value`)
ir_lotta$`Triglycerides p-value` <- as.numeric(ir_lotta$`Triglycerides p-value`)

ir_lotta$beta=NA
ir_lotta$p_value=NA
ir_lotta$best_univariate_trait=NA

for(index in seq(1, length(ir_lotta$`Genomic coordinate`))){
  
  #STEP 1: get data
  
  trait_vect <- c("FIadjBMI", "HDL","TG")
  
  beta_vect <- c(ir_lotta$`Beta FIadjBMI per allele`[index],
                 ir_lotta$`Beta HDL cholesterol per allele`[index],
                 ir_lotta$`Beta triglycerides per allele`[index])
  
  p_vect <- c(ir_lotta$`FIadjBMI p-value`[index],
              ir_lotta$`HDL cholesterol p-value`[index],
              ir_lotta$`Triglycerides p-value`[index])
  
  
  #Let's get the index:
  
  index_best=which.min(as.numeric(p_vect))
  
  trait_best=trait_vect[index_best]
  beta_best=beta_vect[index_best]
  p_best=p_vect[index_best]
  
  ir_lotta$beta[index]=beta_best
  ir_lotta$p_value[index]=p_best
  ir_lotta$best_univariate_trait[index]=trait_best
  
}

#This worked!

ir_lotta$chromosome=as.numeric(as.character(unlist(sapply(ir_lotta$chr_pos, chr_parser))))
ir_lotta$base_pair_location=as.numeric(as.character(unlist(sapply(ir_lotta$chr_pos, pos_parser))))

#Let's get all possible data:

ir_lotta_clean = ir_lotta %>%
  dplyr::select(SNP, chromosome, base_pair_location, effect_allele, other_allele, beta, p_value, best_univariate_trait, chr_pos)

ir_lotta_clean$standard_error=NA
ir_lotta_clean$sample_size=NA

# Function applications
ir_lotta_clean <- alleloi(ir_lotta_clean)
ir_lotta_clean <- mhc_cleaner(ir_lotta_clean)

#Let's clean the data:

ir_lotta_clean$variant=as.character(unlist(sapply(ir_lotta_clean$SNP, rsid_parser))) #this works only when I run everything and load it again???? Bizarre behaviour.

ir_lotta_clean = ir_lotta_clean %>%
  dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, beta, p_value, best_univariate_trait, chr_pos)

# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(ir_lotta_clean, "output/1_curated_data/multi_trait_ir_lotta_leads.txt")

# ---------------------------------------------------------------------------------- #