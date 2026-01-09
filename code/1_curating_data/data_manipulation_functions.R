# ---------------------------------------------------------------------------------- #

                          # DATA MANIPLUATION FUNCTIONS # 

# ---------------------------------------------------------------------------------- #

# obtaining rsID for Pulit data:

rsid_parser = function(rsid_){
  
  #Split data according to : and convert output to vector
  vect_ = as.character(unlist(stringr::str_split(rsid_, ":")))
  
  #Recover the clean rsid, which is the first object in the vector:
  
  rsid_clean=vect_[1]
  
  return(rsid_clean)
  
}

# Effect & Other Allele, Beta, EAF #

pos_aligner_df <- function(df) {
  # numeric beta & eaf
  beta <- as.numeric(df$beta)
  eaf  <- as.numeric(df$effect_allele_frequency)
  
  # Allelles
  new_a1 <- ifelse(beta < 0, df$other_allele, df$effect_allele)
  new_a2 <- ifelse(beta < 0, df$effect_allele, df$other_allele)
  
  # Beta & EAF
  new_beta <- abs(beta)
  new_eaf  <- ifelse(beta < 0, 1 - eaf, eaf)
  
  # Updating
  df$effect_allele <- new_a1
  df$other_allele  <- new_a2
  df$beta          <- new_beta
  df$effect_allele_frequency <- new_eaf
  
  return(df)
}

# ---------------------------------------------------------------------------------- #

# Indel Deletion #

alleloi <- function(data) {
  yes_vect <- c("A","G","C","T")
  data[   data$effect_allele %in% yes_vect 
          & data$other_allele  %in% yes_vect, 
  ]
}

# ---------------------------------------------------------------------------------- #

# EAF Filtering (<0.01 & >0.99) # 

eaf_chooser <- function(data) {
  low <- 0.01 
  up <- 0.99
  data %>%
    dplyr::filter(
      as.numeric(effect_allele_frequency) > low ,
      as.numeric(effect_allele_frequency) < up
    )               
}

# ---------------------------------------------------------------------------------- #

# MHC Cleaning

mhc_cleaner <- function(data) {
  data %>%
    dplyr::filter(
      !(as.numeric(chromosome) == 6 &
          as.numeric(base_pair_location) >= 26000000 &
          as.numeric(base_pair_location) <= 34000000)
    )
}

# ---------------------------------------------------------------------------------- #

# Trait comparison and min p_value selector

pick_min_p <- function(df, cols, prefix) {
  mat <- as.matrix(
    data.frame(lapply(df[cols], as.numeric))
  )
  
  min_p <- apply(
    mat, 1,
    function(x) {
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    }
  )
  
  min_idx <- apply(
    mat, 1,
    function(x) {
      if (all(is.na(x))) {
        NA_integer_
      } else {
        which.min(x)  # index of p value
      }
    }
  )
  
  df[[paste0("multi_trait_p_value_", prefix)]] <- min_p
  df[[paste0("multi_trait_trait_", prefix)]]   <- cols[min_idx]
  
  df
}

# ---------------------------------------------------------------------------------- #


