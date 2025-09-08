library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)


# ...associated with type 2 diabetes, we extracted variants reaching genome-wide significance (p<5×10−8) from large-scale type 2 diabetes studies. 
setwd("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Bora")
getwd()

# Import the data #
diamante_eur <- data.table::fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_t2dlipodystrophy_mahajan_suppl_1.txt")


# Column names rearrangement # 
colnames(diamante_eur)
diamante_eur$sample_size <- NA

diamante_clean <- diamante_eur %>%
  dplyr::select("chromosome(b37)","position(b37)","rsID","effect_allele","other_allele","effect_allele_frequency","Fixed-effects_beta","Fixed-effects_SE","Fixed-effects_p-value","sample_size")


colnames(diamante_clean) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value","sample_size") 
diamante_clean$variant <- as.character(diamante_clean$variant)

# Genome-wide significance
diamante_clean$p_value <- as.numeric(diamante_clean$p_value)
diamante_clean <- diamante_clean[diamante_clean$p_value <= 5e-8] # Manuel filtering #

diamante_clean$effect_allele <- toupper(diamante_clean$effect_allele)
diamante_clean$other_allele <- toupper(diamante_clean$other_allele)

# Function application
diamante_pos <- pos_aligner_df(diamante_clean)
diamante_pos <- alleloi(diamante_pos)
diamante_pos <- eaf_chooser(diamante_pos)
diamante_pos <- mhc_cleaner(diamante_pos)

diamante_pos <- diamante_pos %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

diamante_pos <- diamante_pos %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# Save and check the file #
fwrite(diamante_pos, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/t2dlipodystrophy_mahajan_curated.txt")
check <- fread("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/t2dlipodystrophy_mahajan_curated.txt")
