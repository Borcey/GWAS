library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

whradjbmi <- data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/raw_data/clean_supplementary_tables/clean_whradjbmi_pulit_suppl_7.txt", fill = TRUE)
whradjbmi <- whradjbmi[, -c(7:26)]
colnames(whradjbmi)

# Beta values are already positive, just the columns are renamed and the data is ordered&cleaned (from NA values) #
# P < 5 × 10−9 #


whradjbmi_clean <- whradjbmi %>%
  dplyr::select("RSID","CHR","POS","A1","A0","WHRadjBMI.Combined.A1FREQ", "WHRadjBMI.Combined.BETA","WHRadjBMI.Combined.SE","WHRadjBMI.Combined.Pvalue","WHRadjBMI.Combined.N")

colnames(whradjbmi_clean) <- c("variant","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size")

whradjbmi_clean <- na.omit(whradjbmi_clean) # NA deletion


# Function application
whradjbmi_final <- pos_aligner_df(whradjbmi_clean)
whradjbmi_final <- alleloi(whradjbmi_final)
whradjbmi_final <- eaf_chooser(whradjbmi_final)
whradjbmi_final <- mhc_cleaner(whradjbmi_final)

whradjbmi_final <- whradjbmi_final %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

whradjbmi_final <- whradjbmi_final %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# Export the data and check
fwrite(whradjbmi_final, "N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/output/1_clean_suppl_tables/whradjbmi_pulit_curated.txt")
check <- fread("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora/output/1_clean_suppl_tables/whradjbmi_pulit_curated.txt")
