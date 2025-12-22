##############
#INTRODUCTION#
##############

#This code reads the supplementary table 3 and retrieves that lead posterior thigh fat infiltration from van der meer et al. Data is available here: https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-022-04237-4/MediaObjects/42003_2022_4237_MOESM4_ESM.xlsx

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(readxl)

###################
#Loading functions#
###################

chr_pos_parser = function(id){
  
  vect=as.character(unlist(str_split(id, ":")))
  
  chr_pos_ = paste("chr", vect[1], ":", vect[2], sep = "")
  
  return(chr_pos_)
  
}

a1_parser = function(id){
  
  vect=as.character(unlist(str_split(id, ":")))
  
  allele=vect[3]
  
  return(allele)
  
}

a2_parser = function(id){
  
  vect=as.character(unlist(str_split(id, ":")))
  
  allele=vect[4]
  
  return(allele)
  
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_data<- as.data.frame(readxl::read_excel("input/raw_supplementary_tables/42003_2022_4237_MOESM4_ESM.xlsx")) #dunno why it is sheet 10 according to this. This is Suppl. Table 1h

#The data here is quite difficult to parse, let's just filter for leads and get them through full sumstats:

clean_data = raw_data[which(raw_data$Trait == "PTMFI"),] #17 - no indels!

##################################################################################
#Let's try get the data for pancreas fat percentage: recover the indels via proxy#
##################################################################################

#We have several indels that we should take care of with LDLink and our newly curated data!

posterior_thigh = fread("output/1_curated_data/posterior_thigh_fat_infiltration_van_der_meer_full_sumstats.txt") #full sumstats here do not have INDELS

lead_df <- posterior_thigh[which(posterior_thigh$variant%in%clean_data$`Lead SNP`),] #17

#####################
#Let's save the data#
#####################

fwrite(lead_df, "output/1_curated_data/posterior_thigh_fat_infiltration_van_der_meer_leads.txt")
