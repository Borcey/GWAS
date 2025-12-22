##############
#INTRODUCTION#
##############

#This code reads the supplementary table 3 and retrieves that lead gsatadjBMI from Agrawal et al. Data is available here: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-30931-2/MediaObjects/41467_2022_30931_MOESM4_ESM.xlsx

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(readxl)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_data<- readxl::read_excel("input/raw_supplementary_tables/41467_2022_30931_MOESM4_ESM.xlsx", sheet = 3)

#Let's remove the first 4 rows, which are just an explanation of the columns and traits included

col_names = raw_data[5,]
raw_data_ = as.data.frame(raw_data[6:length(rownames(raw_data)),])
colnames(raw_data_) <- col_names[1,]

raw_data_ <- raw_data_ %>%
  dplyr::select("Trait","CHR","BP","SNP","Effect Allele","Other Allele","EAF","BETA","SE","P-value")

colnames(raw_data_) <- c("trait","chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value")

raw_data_$chr_pos = paste("chr", raw_data_$chromosome, ":", raw_data_$base_pair_location, sep = "")

###################################################
#Let's filter the data for gsatadjbmi - here VATadj#
###################################################

clean_data <- raw_data_[raw_data_$trait =="GFATadj",]

#We have several indels that we should take care of with LDLink and our newly curated data!

gsatadjbmi = fread("output/1_curated_data/gsatadjbmi_agrawal_full_sumstats.txt") #full sumstats here do not have INDELS

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
    
    gsatadjbmi_match <- gsatadjbmi[which(gsatadjbmi$chr_pos%in%proxies_r2_08$Coord),]
    
    #Great, let's find the one with the highest r2 in the original one:
    
    proxies_match <- proxies_r2_08[which(proxies_r2_08$Coord%in%gsatadjbmi_match$chr_pos),]
    
    proxies_best = proxies_match[which.max(as.numeric(proxies_match$R2)),]
    
    #And retrieve it:
    
    gsatadjbmi_best = gsatadjbmi_match[which(gsatadjbmi_match$chr_pos%in%proxies_best$Coord),]
    
    #Just in case we have several that have R2:
    
    gsatadjbmi_best <- gsatadjbmi_best[which.min(as.numeric(gsatadjbmi_best$p_value)),]
    
    #It is almost impossible that this happens, but there is chance that the p-values have been approximated and show the same. 
    #Let's take the first row, randomly
    
    gsatadjbmi_best <- gsatadjbmi_best[1,]
    
    gsatadjbmi_best$lead=snp #we will remove later, just to trace how many we ended up recovering
    
    #Awesome -  let's go and save this:
    
    if(!(exists("recovered_snps"))){
      
      recovered_snps=gsatadjbmi_best
      
    } else {
      
      recovered_snps <- rbind(recovered_snps, gsatadjbmi_best)
      
    }
    
  }
  
}

#We only found one other... shame!
#Okay, then let's add the recovered snps to the data:

non_indels= clean_data[which(clean_data$effect_allele%in%yes_vect & clean_data$other_allele%in%yes_vect),]
non_indels$sample_size = gsatadjbmi$sample_size[1] #we can do this since all SNPs have the same sample size for this case

#And now let's add the recovered_snps, which come from the full_sumstats

recovered_snps <- recovered_snps[which(is.na(recovered_snps$variant) == FALSE),]

#Final set of leads:

lead_df <- gsatadjbmi[which(gsatadjbmi$chr_pos%in%recovered_snps$chr_pos | gsatadjbmi$chr_pos%in%non_indels$chr_pos),] #19 - also removed MHC region!

#This is actually a very clean way to deal with this, since the full_sumstats are already clean!
#We are done!

#####################
#Let's save the data#
#####################

fwrite(lead_df, "output/1_curated_data/gsatadjbmi_agrawal_leads.txt")
