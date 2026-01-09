# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# 0.2 functions:

rsid_special_cleaner <- function(rsid){
  
  if(str_detect(rsid, "Affx")){
    
    vect=as.character(unlist(str_split(rsid, "-")))
    
    rsid_clean=vect[3]
    
  } else {
    
    return(rsid)
    
  }
  
}

# ---------------------------------------------------------------------------------- #

                      # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

options(timeout = 600) # Surpass of timeout risk during download
url_lead <- url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-025-03931-0/MediaObjects/41591_2025_3931_MOESM2_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

multi_trait_fa_chami_lead <- as.data.frame(readxl::read_excel(tmp_lead, sheet = 2, skip = 1))  # Read the excel file, delete first 3 columns

#Let's clean the names first:

## 1. Promote first row to column names
clean_df <- multi_trait_fa_chami_lead
colnames(clean_df) <- as.character(unlist(clean_df[1, ]))
clean_df <- clean_df[-1, ]

## 2. Fix duplicated / empty column names
raw_names <- colnames(multi_trait_fa_chami_lead)
clean_names <- colnames(clean_df)

# Identify trait headers (e.g. BFP_HDL, BMI_SBP, WHR_TG, etc.)
is_trait <- str_detect(raw_names, "^[A-Z]+_[A-Z0-9]+$")

# Carry trait names forward
trait_context <- raw_names
trait_context[!is_trait] <- NA
trait_context <- zoo::na.locf(trait_context, na.rm = FALSE) #Woah this worked amazingly. I love this function

# Rename repeated stat columns
new_names <- map2_chr(clean_names, trait_context, function(nm, tr) {
  if (nm %in% c("Beta", "SE", "Pvalue", "N")) {
    paste(nm, tr, sep = ".")
  } else {
    nm
  }
})

colnames(clean_df) <- new_names

#We almost got it! But we have some format issues. Let's do a quick manual change.

colnames(clean_df) =  c("SNP",                         "Variant ID",                  "Chr",                         "Position",                    "Reference Allele",            "Effect Allele",               "Effect allele frequency (%)", "Traits Prioritized",          "Beta.BFP_HDL",                "SE.BFP_HDL",                  "Pvalue.BFP_HDL",              "N.BFP_HDL",                   "Beta.BFP_TG",                
                        "SE.BFP_TG",                   "Pvalue.BFP_TG",               "N.BFP_TG",                    "Beta.BFP_LDL",                "SE.BFP_LDL",                  "Pvalue.BFP_LDL",              "N.BFP_LDL",                   "Beta.BFP_TC",                 "SE.BFP_TC",                   "Pvalue.BFP_TC",               "N.BFP_TC",                    "Beta.BFP_DBP",                "SE.BFP_DBP",                 
                        "Pvalue.BFP_DBP",              "N.BFP_DBP",                   "Beta.BFP_SBP",                "SE.BFP_SBP",                  "Pvalue.BFP_SBP",              "N.BFP_SBP",                   "Beta.BFP_HbA1c",                "SE.BFP_HbA1c",                  "Pvalue.BFP_HbA1c",              "N.BFP_HbA1c",                   "Beta.BFP_Glucose",                "SE.BFP_Glucose",                  "Pvalue.BFP_Glucose",             
                        "N.BFP_Glucose",                   "Beta.BMI_HDL",                "SE.BMI_HDL",                  "Pvalue.BMI_HDL",              "N.BMI_HDL",                   "Beta.BMI_TG",                 "SE.BMI_TG",                   "Pvalue.BMI_TG",               "N.BMI_TG",                    "Beta.BMI_LDL",                "SE.BMI_LDL",                  "Pvalue.BMI_LDL",              "N.BMI_LDL",                  
                        "Beta.BMI_TC",                 "SE.BMI_TC",                   "Pvalue.BMI_TC",               "N.BMI_TC",                    "Beta.BMI_DBP",                "SE.BMI_DBP",                  "Pvalue.BMI_DBP",              "N.BMI_DBP",                   "Beta.BMI_SBP",                "SE.BMI_SBP",                  "Pvalue.BMI_SBP",              "N.BMI_SBP",                   "Beta.BMI_HbA1c",               
                        "SE.BMI_HbA1c",                  "Pvalue.BMI_HbA1c",              "N.BMI_HbA1c",                   "Beta.BMI_Glucose",                "SE.BMI_Glucose",                  "Pvalue.BMI_Glucose",              "N.BMI_Glucose",                   "Beta.WHR_HDL",                "SE.WHR_HDL",                  "Pvalue.WHR_HDL",              "N.WHR_HDL",                   "Beta.WHR_TG",               "SE.WHR_TG",                  
                        "Pvalue.WHR_TG",               "N.WHR_TG",                    "Beta.WHR_LDL",                "SE.WHR_LDL",                  "Pvalue.WHR_LDL",              "N.WHR_LDL",                   "Beta.WHR_TC",                 "SE.WHR_TC",                   "Pvalue.WHR_TC",               "N.WHR_TC",                    "Beta.WHR_DBP",                "SE.WHR_DBP",                  "Pvalue.WHR_DBP",             
                        "N.WHR_DBP",                   "Beta.WHR_SBP",                "SE.WHR_SBP",                  "Pvalue.WHR_SBP",              "N.WHR_SBP",                   "Beta.WHR_HbA1c",                "SE.WHR_HbA1c",                  "Pvalue.WHR_HbA1c",              "N.WHR_HbA1c",                  "Beta.WHR_Glucose",                "SE.WHR_Glucose",                  "Pvalue.WHR_Glucose",              "N.WHR_Glucose") 


#3. Let's get the best p-value then:

best_assoc <- clean_df %>%
  rowwise() %>%
  mutate(
    best_trait = {
      traits <- str_split(`Traits Prioritized`, ";")[[1]]
      
      pvals <- map_dbl(traits, ~ as.numeric(cur_data()[[paste0("Pvalue.", .x)]]))
      
      traits[which.min(pvals)]
    },
    
    best_beta = cur_data()[[paste0("Beta.", best_trait)]],
    best_se   = cur_data()[[paste0("SE.", best_trait)]],
    best_p    = cur_data()[[paste0("Pvalue.", best_trait)]],
    best_n    = cur_data()[[paste0("N.", best_trait)]]
  ) %>%
  ungroup()

#4. Let's get the best data:

best_data_clean <- best_assoc %>%
  dplyr::select("SNP",                            "Chr",                         "Position",                    "Reference Allele",            "Effect Allele",               "Effect allele frequency (%)",  best_trait, best_beta, best_se, best_p, best_n)

#Let's finish cleaning this:

colnames(best_data_clean) = c("variant", "chromosome", "base_pair_location", "other_allele", "effect_allele", "effect_allele_frequency", "best_trait", "beta", "standard_error", "p_value", "sample_size")

#We are almost ready to go, let's make all columns ready for analyses:

best_data_clean$effect_allele_frequency = as.numeric(best_data_clean$effect_allele_frequency)
best_data_clean$beta = as.numeric(best_data_clean$beta)
best_data_clean$standard_error = as.numeric(best_data_clean$standard_error)
best_data_clean$p_value = as.numeric(best_data_clean$p_value)
best_data_clean$sample_size = as.numeric(best_data_clean$sample_size)

#Let's clean the data (I checked - it is on build 37)

# Function applications

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")
source("code/1_curating_data/data_manipulation_functions.R") #all SNPs are positive already, no need to do anything else here
best_data_clean <- alleloi(best_data_clean)
best_data_clean <- eaf_chooser(best_data_clean)
best_data_clean <- mhc_cleaner(best_data_clean)

best_data_clean <- best_data_clean %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

best_data_clean <- best_data_clean %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

#We need to do a special cleaning of rsIDs:

best_data_clean$variant = as.character(sapply(best_data_clean$variant, rsid_special_cleaner))


# ---------------------------------------------------------------------------------- #

# 3 - Data Exportation

fwrite(best_data_clean, "output/1_curated_data/multi_trait_fa_chami_leads.txt")
