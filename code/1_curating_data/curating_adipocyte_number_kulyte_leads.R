# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

#0.1 cleaner function:

cleaner_=function(pos){
  
  pos_clean = str_remove_all(pos, " ")
  
  return(pos_clean)
  
}

# ---------------------------------------------------------------------------------- #

                  # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here

#The data has mostly being cleaned manually since it has been extracted from a word document.

leads <- as.data.frame(readxl::read_excel("input/raw_supplementary_tables/19383830/adipocyte_number_kulyte_leads.xlsx"))  # data is in build 37! Just checked 

leads$base_pair_location <- as.numeric(as.character(unlist(sapply(leads$base_pair_location, cleaner_))))

leads$effect_allele_frequency <- as.numeric(leads$effect_allele_frequency) #this worked


# Function applications
clean_data <- alleloi(leads)
#clean_data <- eaf_chooser(clean_data) #Let's use all of them, which are borderline rare. Nonetheless they behave well under LD.
clean_data <- mhc_cleaner(leads)


# chr_pos addition
clean_data <- clean_data %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

clean_data <- clean_data %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

clean_data <- clean_data %>%
  dplyr::select(-c("...6")) #removing weird column

# ---------------------------------------------------------------------------------- #

# 3.2 - Data Exportation 

fwrite(clean_data, "output/1_curated_data/adipocyte_number_kulyte_leads.txt")

