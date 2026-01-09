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

options(timeout = 100000) # Surpass of timeout risk during download
url_lead <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-021-00346-2/MediaObjects/42255_2021_346_MOESM2_ESM.xlsx"

tmp_lead <- tempfile(fileext = ".xlsx") # Temporary file creation
download.file(url_lead, tmp_lead, mode = "wb") # Downloading excel file from URL 

multi_trait_huang_lead <- readxl::read_excel(tmp_lead, sheet = 2, skip=3)  # Read the excel file, delete first 3 columns

# ---------------------------------------------------------------------------------- #

# 2 - Data Manipulation 

multi_trait_huang_lead_clean_pairwise <- multi_trait_huang_lead %>%
  dplyr::select("...1", "...2","...3","...4","...5","...6","beta BFP", "SE.BFP", "P.BFP", "N.BFP","beta BMI", "SE.BMI", "P.BMI", "N.BMI","beta WHR","SE.WHR", "P.WHR","N.WHR",
                "P.BFP_HDL" ,    
                "P.BFP_LDL"    ,  "P.BFP_TG"    ,   "P.BFP_FI"    ,   "P.BFP_FG"   ,    "P.BFP_SBP"   ,   "P.BFP_CAD"    ,  "P.BFP_T2D" ,    
                "P.BMI_HDL"   ,   "P.BMI_LDL"    ,  "P.BMI_TG"   ,    "P.BMI_FI"    ,   "P.BMI_FG"    ,   "P.BMI_SBP"   ,   "P.BMI_CAD"  ,   
                "P.BMI_T2D"   ,   "P.WHR_HDL"    ,  "P.WHR_LDL"  ,    "P.WHR_TG"     ,  "P.WHR_FI"     ,  "P.WHR_FG"    ,   "P.WHR_SBP"   ,  
                "P.WHR_CAD"   ,   "P.WHR_T2D") 

colnames(multi_trait_huang_lead_clean_pairwise) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency",
                                                     "beta_BFP","standard_error_BFP","p_value_BFP","sample_size_BFP",
                                                     "beta_BMI","standard_error_BMI","p_value_BMI","sample_size_BMI",
                                                     "beta_WHR","standard_error_WHR","p_value_WHR","sample_size_WHR",
                                                  
                                                     "P.BFP_HDL" ,   "P.BFP_LDL"  ,    "P.BFP_TG"   ,   "P.BFP_FI"    ,   "P.BFP_FG"   ,   "P.BFP_SBP"   ,   "P.BFP_CAD"  ,  "P.BFP_T2D",    
                                                     "P.BMI_HDL" ,   "P.BMI_LDL"  ,    "P.BMI_TG"   ,   "P.BMI_FI"    ,   "P.BMI_FG"   ,   "P.BMI_SBP"   ,   "P.BMI_CAD"  ,  "P.BMI_T2D",   
                                                     "P.WHR_HDL" ,   "P.WHR_LDL"  ,    "P.WHR_TG"   ,   "P.WHR_FI"    ,   "P.WHR_FG"   ,   "P.WHR_SBP"   ,   "P.WHR_CAD"  ,  "P.WHR_T2D") # Standardization of target column names




bfp_cols <- c("P.BFP_HDL", "P.BFP_LDL", "P.BFP_TG", "P.BFP_FI",
              "P.BFP_FG",  "P.BFP_SBP", "P.BFP_CAD", "P.BFP_T2D")

bmi_cols <- c("P.BMI_HDL", "P.BMI_LDL", "P.BMI_TG", "P.BMI_FI",
              "P.BMI_FG",  "P.BMI_SBP", "P.BMI_CAD", "P.BMI_T2D")

whr_cols <- c("P.WHR_HDL", "P.WHR_LDL", "P.WHR_TG", "P.WHR_FI",
              "P.WHR_FG",  "P.WHR_SBP", "P.WHR_CAD", "P.WHR_T2D")


###### F U N C T I O N ########

# Helper function that performs the same operation for each trait block
pick_min_p <- function(df, cols, prefix) {
  mat <- as.matrix(
    data.frame(lapply(df[cols], as.numeric))
  )
  
  # Row-wise extraction of the minimum p-value

  min_p <- apply(
    mat, 1,
    function(x) {
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    }
  )
  
  # Identify, for each row, the column index where the minimum p-value occurs
  min_idx <- apply(
    mat, 1,
    function(x) {
      if (all(is.na(x))) {
        NA_integer_
      } else {
        which.min(x)  # index of the minimum p-value
      }
    }
  )
  
  # Store both the minimum p-value and the corresponding trait name
  df[[paste0("multi_trait_p_value_", prefix)]] <- min_p
  df[[paste0("multi_trait_trait_", prefix)]]   <- cols[min_idx]
  
  df
}

# Apply it
multi_trait_huang_lead_clean_pairwise <- multi_trait_huang_lead_clean_pairwise |>
  pick_min_p(bfp_cols, "BFP") |>
  pick_min_p(bmi_cols, "BMI") |>
  pick_min_p(whr_cols, "WHR")

#####################################################################################

#Let's format the data to get the best data:

multi_trait_huang_lead_clean_pairwise$beta=NA
multi_trait_huang_lead_clean_pairwise$standard_error=NA
multi_trait_huang_lead_clean_pairwise$p_value.univariate=NA
multi_trait_huang_lead_clean_pairwise$p_value=NA
multi_trait_huang_lead_clean_pairwise$best_univariate_trait=NA
multi_trait_huang_lead_clean_pairwise$best_pairwise_trait=NA

#Let's loop first for the univariate data:

for(index in seq(1, length(multi_trait_huang_lead_clean_pairwise$chromosome))){
  
  #STEP 1: get data
  
  trait_vect <- c("BFP", "BMI","WHR")
  
  beta_vect <- c(multi_trait_huang_lead_clean_pairwise$beta_BFP[index],
              multi_trait_huang_lead_clean_pairwise$beta_BMI[index],
              multi_trait_huang_lead_clean_pairwise$beta_WHR[index])
  
  se_vect <- c(multi_trait_huang_lead_clean_pairwise$standard_error_BFP[index],
              multi_trait_huang_lead_clean_pairwise$standard_error_BMI[index],
              multi_trait_huang_lead_clean_pairwise$standard_error_WHR[index])
  
  p_vect <- c(multi_trait_huang_lead_clean_pairwise$p_value_BFP[index],
               multi_trait_huang_lead_clean_pairwise$p_value_BMI[index],
               multi_trait_huang_lead_clean_pairwise$p_value_WHR[index])
  
  
  #Let's get the index:
  
  index_best=which.min(as.numeric(p_vect))
  
  trait_best=trait_vect[index_best]
  beta_best=beta_vect[index_best]
  se_best=se_vect[index_best]
  p_best=p_vect[index_best]
  
  multi_trait_huang_lead_clean_pairwise$beta[index]=beta_best
  multi_trait_huang_lead_clean_pairwise$standard_error[index]=beta_best
  multi_trait_huang_lead_clean_pairwise$p_value.univariate[index]=p_best
  multi_trait_huang_lead_clean_pairwise$best_univariate_trait[index]=trait_best

}

#Seems like this work!! Let's change do the same on the multi-trait

for(index in seq(1, length(multi_trait_huang_lead_clean_pairwise$chromosome))){
  
  #STEP 1: get data
  
  trait_vect <- c("BFP", "BMI","WHR")
  
  multi_trait <- c(multi_trait_huang_lead_clean_pairwise$multi_trait_trait_BFP[index],
                 multi_trait_huang_lead_clean_pairwise$multi_trait_trait_BMI[index],
                 multi_trait_huang_lead_clean_pairwise$multi_trait_trait_WHR[index])
  
  p_vect <- c(multi_trait_huang_lead_clean_pairwise$multi_trait_p_value_BFP[index],
               multi_trait_huang_lead_clean_pairwise$multi_trait_p_value_BMI[index],
               multi_trait_huang_lead_clean_pairwise$multi_trait_p_value_WHR[index])
  
  
  #Let's get the index:
  
  index_best=which.min(as.numeric(p_vect))
  
  trait_best=trait_vect[index_best]
  multi_trait_best=multi_trait[index_best]
  p_best=p_vect[index_best]
  
  multi_trait_huang_lead_clean_pairwise$p_value[index]=p_best
  multi_trait_huang_lead_clean_pairwise$best_pairwise_trait[index]=multi_trait_best
  
}

#We've got this!!! Let's filter the data:

clean_data= multi_trait_huang_lead_clean_pairwise %>%
  dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value.univariate, best_univariate_trait, p_value, best_pairwise_trait)

clean_data$sample_size=NA


# Function applications
clean_data <- alleloi(clean_data)
clean_data <- eaf_chooser(clean_data)
clean_data <- mhc_cleaner(clean_data)


# chr_pos addition
clean_data <- clean_data %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

clean_data <- clean_data %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3.2 - Data Exportation 

fwrite(clean_data, "output/1_curated_data/multi_trait_fa_huang_leads.txt")

