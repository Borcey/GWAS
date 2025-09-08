library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

#Seeting working directory and loading data

fiadj <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_fiadjbmi_che_suppl_2.xlsx", skip = 3) #the data has weird formatting due to complex Excel formatting. 
fiadj <- fiadj[, 1:20]  # Take the first 20 columns which includes everything we need. #


# ********* FI Trait Selection ********* #
fiadj <- fiadj[fiadj$Trait == "FI",] # , is important

# ********* European Ancestry Selection ********* #
fiadj <- fiadj[fiadj$`Analysis in which variant association discovered` %in% c("EUR","TA, EUR"),] # %>% is essential for multiple selection. #

#Let's check the data:
View(fiadj)
table(fiadj$Trait)

#   2hGlu    FG     FI   HbA1c 
#   28       182    95   218 


#AA      EUR     HISP       TA       TA, EUR     TA, HISP 
#1       30      2          47       14          1 

#We are interested in EUR and TA, EUR -  a total of 44 signals


dim(fiadj)
#[1] 44 20 #PERFECT!

# ********************************************************************************************************************************* #

# Column renaming #
colnames(fiadj)
fiadj$sample_size <- NA

fiadj_clean <- fiadj %>%
  dplyr::select("rsID","Chr...3","Pos (bp)","Effect Allele","Other Allele", "EAF...16","Effect...17","SE...18","GCTA Pvalue...20","sample_size")

colnames(fiadj_clean)
colnames(fiadj_clean) <- c("variant","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size")

# ********************************************************************************************************************************* #

fiadj_clean$effect_allele_frequency <- as.numeric(fiadj_clean$effect_allele_frequency)
fiadj_clean$p_value <- as.numeric(fiadj_clean$p_value)
fiadj_clean$beta <- as.numeric(fiadj_clean$beta)
fiadj_clean$standard_error <- as.numeric(fiadj_clean$standard_error)

# ********************************************************************************************************************************* #

# Function application
fiadj_pos <- pos_aligner_df(fiadj_clean)
fiadj_pos <- alleloi(fiadj_pos)
fiadj_pos <- eaf_chooser(fiadj_pos)
fiadj_pos <- mhc_cleaner(fiadj_pos)

fiadj_pos <- fiadj_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

fiadj_pos <- fiadj_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ********************************************************************************************************************************* #

# Export the data as .txt #
fwrite(fiadj_pos,"N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/output/1_clean_suppl_tables/fiadj_che_curated.txt")
check <- fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/output/1_clean_suppl_tables/fiadj_che_curated.txt")




