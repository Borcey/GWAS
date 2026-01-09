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

leads <- as.data.frame(readxl::read_excel("input/raw_supplementary_tables/jamacardiol-e240080-s002.xlsx", sheet = 10))  # first page is an index

leads <- leads %>%
  dplyr::select(rsID, Chr, Position, EA, NEA, EAF, Beta, SE, P)

colnames(leads)=c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")

leads$sample_size=41494

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

fwrite(clean_data, "output/1_curated_data/epat_area_ramo_leads.txt")

