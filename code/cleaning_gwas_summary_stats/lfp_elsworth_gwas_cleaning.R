library(data.table)
library(tidyverse)
library(dplyr)

# ----------------------------------------------------------------------- #
# Data implementation
lfr_left <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/leg_fat_percentage/lfr_left.tsv")
lfr_right <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/leg_fat_percentage/lfr_right.tsv")

# ----------------------------------------------------------------------- #
# "FORMAT" column separation

lfr_left <- lfr_left %>%
  separate(
    col = `UKB-b-18377`,
    into = c("beta","standard_error","p_value","effect_allele_frequency","rsID"),
    sep = ":",
    remove = FALSE
  )


lfr_right <- lfr_right %>%
  separate(
    col = `UKB-b-20531`,
    into = c("beta","standard_error","p_value","effect_allele_frequency","rsID"),
    sep = ":",
    remove = FALSE
  )

# ----------------------------------------------------------------------- #
# Column names organizing

colnames(lfr_left)

lfr_left <- lfr_left %>%
  dplyr::select("#CHROM","POS","rsID","REF","ALT","FORMAT","UKB-b-18377","beta","standard_error","p_value","effect_allele_frequency")

colnames(lfr_left) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","FORMAT","UKB-b-18377","beta","standard_error","p_value","effect_allele_frequency")


colnames(lfr_right)

lfr_right <- lfr_right %>%
  dplyr::select("#CHROM","POS","rsID","REF","ALT","FORMAT","UKB-b-20531","beta","standard_error","p_value","effect_allele_frequency")

colnames(lfr_right) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","FORMAT","UKB-b-18377","beta","standard_error","p_value","effect_allele_frequency")

# ----------------------------------------------------------------------- #

lfr_left$sample_size <- NA
lfr_right$sample_size <- NA

# ----------------------------------------------------------------------- #
# LOG10P transform

lfr_left$p_value <- as.numeric(lfr_left$p_value)
lfr_right$p_value <- as.numeric(lfr_right$p_value)
lfr_left$p_value <- 10^(-lfr_left$p_value)
lfr_right$p_value <- 10^(-lfr_right$p_value)

# ----------------------------------------------------------------------- #
# p_value filtering

lfr_left <- lfr_left %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)


lfr_right <- lfr_right %>%
  dplyr::mutate(p_value = as.numeric(p_value)) %>%
  dplyr::filter(!is.na(p_value), p_value < 5e-8)

# ----------------------------------------------------------------------- #
# application of 4 functions

nrow(lfr_left)
lfr_left <- pos_aligner_df(lfr_left)
nrow(lfr_left)
lfr_left <- alleloi(lfr_left)
nrow(lfr_left)
lfr_left <- eaf_chooser(lfr_left)
nrow(lfr_left)
lfr_left <- mhc_cleaner(lfr_left)
nrow(lfr_left)


nrow(lfr_right)
lfr_right <- pos_aligner_df(lfr_right)
nrow(lfr_right)
lfr_right <- alleloi(lfr_right)
nrow(lfr_right)
lfr_right <- eaf_chooser(lfr_right)
nrow(lfr_right)
lfr_right <- mhc_cleaner(lfr_right)
nrow(lfr_right)

# ----------------------------------------------------------------------- #
# chr_pos column creation

lfr_left <- lfr_left %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

lfr_left <- lfr_left %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  



lfr_right <- lfr_right %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":" ))

lfr_right <- lfr_right %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))  

# ----------------------------------------------------------------------- #
# Export for separeted versions
fwrite(lfr_left, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/lfr_left_elsworth_chrpos.txt")
fwrite(lfr_right, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/raw_full_sumstats/output/chr_pos_gwas/lfr_right_elsworth_chrpos.txt")

# ----------------------------------------------------------------------- #

# These parts cover the merging of lfr_left and lfr_right steps

# A - FINDING THE MAX AND MIN STANDARD_ERROR DIFFERENCES

# 1. Find common variants for both 
common_variants <- intersect(lfr_left$variant, lfr_right$variant)

# 2. Just get the common ones
left_common  <- lfr_left  %>% filter(variant %in% common_variants) %>% select(variant, se_left = standard_error)
right_common <- lfr_right %>% filter(variant %in% common_variants) %>% select(variant, se_right = standard_error)

# 3. Merge them 
merged_se <- left_common %>%
  inner_join(right_common, by = "variant") %>%
  mutate(
    se_left = as.numeric(se_left),
    se_right = as.numeric(se_right),
    diff_se = abs(se_left - se_right)
  )

# 4. Find the Min ve Max differences 
min_diff_variant <- merged_se %>% filter(diff_se == min(diff_se, na.rm = TRUE))
max_diff_variant <- merged_se %>% filter(diff_se == max(diff_se, na.rm = TRUE))

min_diff_variant
max_diff_variant

# ----------------------------------------------------------------------- #
# FOR P VAL (May not be necessary)

library(dplyr)

# 1. Find common variants for both 
common_variants <- intersect(lfr_left$variant, lfr_right$variant)

# 2. Just get the common ones
left_common  <- lfr_left  %>% filter(variant %in% common_variants) %>% select(variant, p_left = p_value)
right_common <- lfr_right %>% filter(variant %in% common_variants) %>% select(variant, p_right = p_value)

# 3. Merge them 
merged_p <- left_common %>%
  inner_join(right_common, by = "variant") %>%
  mutate(
    p_left = as.numeric(p_left),
    p_right = as.numeric(p_right),
    diff_p = abs(p_left - p_right)
  )

# 4. Find the Min ve Max differences 
min_diff_variant <- merged_p %>% filter(diff_p == min(diff_p, na.rm = TRUE))
max_diff_variant <- merged_p %>% filter(diff_p == max(diff_p, na.rm = TRUE))

min_diff_variant
max_diff_variant

# ----------------------------------------------------------------------- #
# B - UNIQUE - DUPLICATED VARIANTS AND AVERAGE VALUES


# 1. Find unique and duplicated ones

common_variants <- intersect(lfr_left$variant, lfr_right$variant)

left_unique  <- lfr_left  %>% filter(!variant %in% common_variants)
right_unique <- lfr_right %>% filter(!variant %in% common_variants)

left_common  <- lfr_left  %>% filter(variant %in% common_variants)
right_common <- lfr_right %>% filter(variant %in% common_variants)


# 2. Merge commons

common_all <- bind_rows(left_common, right_common)


# 3. Numeric transforming

make_numeric <- function(df) {
  df %>%
    mutate(
      beta = as.numeric(beta),
      standard_error = as.numeric(standard_error),
      p_value = as.numeric(p_value),
      effect_allele_frequency = as.numeric(effect_allele_frequency)
    )
}

left_unique  <- make_numeric(left_unique)
right_unique <- make_numeric(right_unique)

common_all   <- make_numeric(common_all)


# 4. Calculate the averages for commons

common_avg <- common_all %>%
  group_by(chromosome, base_pair_location, variant, effect_allele, other_allele) %>%
  summarise(
    beta = mean(beta, na.rm = TRUE),
    standard_error = mean(standard_error, na.rm = TRUE),
    p_value = mean(p_value, na.rm = TRUE),
    effect_allele_frequency = mean(effect_allele_frequency, na.rm = TRUE),
    .groups = "drop"
  )


# 5. Final dataset (unique + averaged duplicates)

final_df <- bind_rows(left_unique, right_unique, common_avg)

# 6. CHECK

cat("Duplicate (common):", length(common_variants), "\n")
cat("Unique variant:", nrow(left_unique) + nrow(right_unique), "\n")
cat("Final dataset :", nrow(final_df), "\n")
cat("Real unique variant (final_df):", length(unique(final_df$variant)), "\n")

length(unique(final_df$variant))

# ----------------------------------------------------------------------- #
# Export the merged version of lfr
fwrite(final_df, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/raw_data/raw_full_sumstats/gwas_sums/leg_fat_percentage/merged_lfr.txt")

# ! Bon appetit !