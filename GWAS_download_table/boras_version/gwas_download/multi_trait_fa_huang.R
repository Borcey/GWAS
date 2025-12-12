# 0 - Libraries

library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

# ---------------------------------------------------------------------------------- #

                  # STANDARDIZATION OF LEAD SNPS TABLES #

# ---------------------------------------------------------------------------------- #

# 1 - Data Fetch 

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

# p-value filtration
multi_trait_huang_lead_clean_pairwise <- multi_trait_huang_lead_clean_pairwise %>%
  dplyr::mutate(
    p_value_BFP = as.numeric(p_value_BFP),
    p_value_BMI = as.numeric(p_value_BMI),
    p_value_WHR = as.numeric(p_value_WHR)
  ) %>%
  
  dplyr::filter(
    (p_value_BFP < 5e-8 & !is.na(p_value_BFP)) |
      (p_value_BMI < 5e-8 & !is.na(p_value_BMI)) |
      (p_value_WHR < 5e-8 & !is.na(p_value_WHR))
  )



multi_trait_huang_lead_pairwise_pos <- multi_trait_huang_lead_pairwise_pos %>%
  dplyr::rename()


# Function applications
#multi_trait_huang_lead_pairwise_pos <- pos_aligner_df(multi_trait_huang_lead_clean_pairwise)
multi_trait_huang_lead_pairwise_pos <- alleloi(multi_trait_huang_lead_clean_pairwise)
multi_trait_huang_lead_pairwise_pos <- eaf_chooser(multi_trait_huang_lead_pairwise_pos)
multi_trait_huang_lead_pairwise_pos <- mhc_cleaner(multi_trait_huang_lead_pairwise_pos)


# chr_pos addition
multi_trait_huang_lead_pairwise_pos <- multi_trait_huang_lead_pairwise_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

multi_trait_huang_lead_pairwise_pos <- multi_trait_huang_lead_pairwise_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# ---------------------------------------------------------------------------------- #

# 3.2 - Data Exportation 

fwrite(multi_trait_huang_lead_pairwise_pos, "C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_huang_leads_pairwise.txt")

# ---------------------------------------------------------------------------------- #

#check <- data.table::fread("C:/Users/Bora Ceylan/Desktop/CBMR/gwas_download/output/leads/multi_trait_fa_huang_leads_pairwise.txt")
