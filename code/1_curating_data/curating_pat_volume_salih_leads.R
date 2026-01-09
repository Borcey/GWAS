# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                  # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here

#The data has mostly being cleaned manually since it has been extracted from a word document.

leads <- as.data.frame(readxl::read_excel("input/raw_supplementary_tables/pat_volume_salih_leads.xlsx"))  #it was easier to extract the data manually - here we will just do a cleaning.

#Just checked and the effect alleles are the MAF alleles, so we can switch the MAF column to eaf:

colnames(leads)[which(colnames(leads) == "minimum_allele_frequency")]="effect_allele_frequency"

leads$base_pair_location <- as.numeric(leads$base_pair_location)
leads$beta <- as.numeric(leads$beta)
leads$standard_error <- as.numeric(leads$standard_error)
leads$p_value <- as.numeric(leads$p_value)

# Function applications
clean_data <- alleloi(leads)
clean_data <- eaf_chooser(clean_data) #Let's use all of them, which are borderline rare. Nonetheless they behave well under LD.
clean_data <- mhc_cleaner(clean_data)

# chr_pos addition
clean_data <- clean_data %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

clean_data <- clean_data %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3.2 - Data Exportation 

fwrite(clean_data, "output/1_curated_data/pat_volume_salih_leads.txt")

