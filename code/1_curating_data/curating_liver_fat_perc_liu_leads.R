##############
#INTRODUCTION#
##############

#This code reads the supplementary table 3 and retrieves that lead Liver fat % from Liu et al. Data is available here: https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjU1NTQvZWxpZmUtNjU1NTQtc3VwcDEtdjEueGxzeA--/elife-65554-supp1-v1.xlsx?_hash=HK79XqeA2thYaX%2F24w7j4obnmtM%2FXJEL7Wvp7uNxxug%3D

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

raw_data<- as.data.frame(readxl::read_excel("input/raw_supplementary_tables/elife-65554-supp1-v1.xlsx", sheet = 10)) #dunno why it is sheet 10 according to this. This is Suppl. Table 1h

#Let's remove the first 4 rows, which are just an explanation of the columns and traits included

colnames(raw_data)=c("trait", "organ", "var_id", "var_index", "variant", "var_conditional", "p_value", "pp", "beta", "closest_gene", "distance")

#Let's use the full summary statistics for this, because the format is horrible

clean_data = raw_data[which(raw_data$trait == "Fat" & raw_data$organ == "Liver"),] #12

#We have one indel, let's take chr_pos and run our pipeline:

clean_data$chr_pos = as.character(unlist(sapply(clean_data$var_id, chr_pos_parser)))
clean_data$effect_allele = as.character(unlist(sapply(clean_data$var_id, a1_parser)))
clean_data$other_allele = as.character(unlist(sapply(clean_data$var_id, a2_parser)))

###############################################################################
#Let's try get the data for liver fat percentage: recover the indels via proxy#
###############################################################################

#We have several indels that we should take care of with LDLink and our newly curated data!

liver_fat_perc = fread("output/1_curated_data/liver_fat_perc_liu_full_sumstats.txt") #full sumstats here do not have INDELS

#Let's loop over the indels and retrieve the data:

yes_vect = c("A", "G", "C", "T")

indel <- clean_data[which(!(clean_data$effect_allele%in%yes_vect) | !(clean_data$other_allele%in%yes_vect)),]

for(index in seq(1, length(indel$chr_pos))){
  
  #STEP 1: retrieve the chr_pos:
  
  snp=indel$chr_pos[index]
  
  #STEP 2: get proxies:
  
  proxies=tryCatch(LDlinkR::LDproxy(snp, pop = "EUR", token = "04cad4ca4374", r2d = "r2", genome_build = "grch37"),   error = function(e) return(NA))
  
  if(length(proxies) > 1){
    
    proxies_r2_08=proxies[which(as.numeric(proxies$R2) > 0.8),]
    
    #Now let's find it in the best dataset:
    
    liver_fat_perc_match <- liver_fat_perc[which(liver_fat_perc$chr_pos%in%proxies_r2_08$Coord),]
    
    #Great, let's find the one with the highest r2 in the original one:
    
    proxies_match <- proxies_r2_08[which(proxies_r2_08$Coord%in%liver_fat_perc_match$chr_pos),]
    
    proxies_best = proxies_match[which.max(as.numeric(proxies_match$R2)),]
    
    #And retrieve it:
    
    liver_fat_perc_best = liver_fat_perc_match[which(liver_fat_perc_match$chr_pos%in%proxies_best$Coord),]
    
    #Just in case we have several that have R2:
    
    liver_fat_perc_best <- liver_fat_perc_best[which.min(as.numeric(liver_fat_perc_best$p_value)),]
    
    #It is almost impossible that this happens, but there is chance that the p-values have been approximated and show the same. 
    #Let's take the first row, randomly
    
    liver_fat_perc_best <- liver_fat_perc_best[1,]
    
    liver_fat_perc_best$lead=snp #we will remove later, just to trace how many we ended up recovering
    
    #Awesome -  let's go and save this:
    
    if(!(exists("recovered_snps"))){
      
      recovered_snps=liver_fat_perc_best
      
    } else {
      
      recovered_snps <- rbind(recovered_snps, liver_fat_perc_best)
      
    }
    
  }
  
}

#We only found one other... shame!
#Okay, then let's add the recovered snps to the data:

non_indels= clean_data[which(clean_data$effect_allele%in%yes_vect & clean_data$other_allele%in%yes_vect),]
non_indels$sample_size = liver_fat_perc$sample_size[1] #we can do this since all SNPs have the same sample size for this case

#And now let's add the recovered_snps, which come from the full_sumstats

recovered_snps <- recovered_snps[which(is.na(recovered_snps$variant) == FALSE),]

#Final set of leads:

lead_df <- liver_fat_perc[which(liver_fat_perc$chr_pos%in%recovered_snps$chr_pos | liver_fat_perc$chr_pos%in%non_indels$chr_pos),] #11 - we lose one

#This is actually a very clean way to deal with this, since the full_sumstats are already clean!
#We are done!

#####################
#Let's save the data#
#####################

fwrite(lead_df, "output/1_curated_data/liver_fat_perc_agrawal_leads.txt")
