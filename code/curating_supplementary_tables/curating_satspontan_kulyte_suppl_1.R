library(ggplot2)
library(tidyverse)
library(data.table)
library(readxl)

spontan_sat <- readxl::read_excel("/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/input/clean_supplementary_tables/clean_satspontan_kulyte_suppl_1.xlsx")
colnames(spontan_sat)


spontan_sat <- spontan_sat %>%
  dplyr::select("Chr","Pos","SNP","Ref","Alt","Freq","N","Beta","Se","P")

colnames(spontan_sat) <- c("chromosome","base_pair_location","variant","effect_allele","other_allele","effect_allele_frequency","sample_size","beta","standard_error","p_value")

spontan_sat <- spontan_sat %>%
  filter(p_value<5e-8)


# Function application
spontan_sat <- pos_aligner_df(spontan_sat)
spontan_sat <- alleloi(spontan_sat)
spontan_sat <- eaf_chooser(spontan_sat)
spontan_sat <- mhc_cleaner(spontan_sat)

spontan_sat <- spontan_sat %>%
  mutate(chr_pos = paste(chromosome, base_pair_location, sep = ":"))

spontan_sat <- spontan_sat %>%
  mutate(chr_pos = paste0("chr",chromosome, ":",base_pair_location))

# Beta values are alreadily positive and base_pair_locations are really close each other.

fwrite(spontan_sat, "/maps/projects/kilpelainen-AUDIT/data/team_projects/Bora/output/curated_supplementary_tables/satspontan_kulyte_curated.txt")
