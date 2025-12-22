##############
#INTRODUCTION#
##############

#This code clean the full summary statistics for Liver fat % from Liu et al, available in https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016673/GCST90016673_buildGRCh37.tsv.gz

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

raw_sumstat <- fread("input/raw_full_sumstats/GCST90016673_buildGRCh37.tsv.gz")

###########################################
#STEP 1: let's check that we have all data#
###########################################

head(raw_sumstat)

# variant_id p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error
# 1: rs367896724    0.95          1              10177            AC            A                0.399875  0.000681779      0.0114655
# 2: rs201106462    0.52          1              10352            TA            T                0.394299  0.007605820      0.0118052
# 3: rs575272151    0.20          1              11008             G            C                0.085902 -0.025156200      0.0196611
# 4: rs544419019    0.20          1              11012             G            C                0.085902 -0.025156200      0.0196611
# 5:  rs62635286    0.61          1              13116             G            T                0.188596  0.007998440      0.0155328
# 6:  rs62028691    0.61          1              13118             G            A                0.188596  0.007998440      0.0155328

#Only lacking sample size - let's add it:

raw_sumstat$sample_size = 32858 

#Let's also add chr_pos. We know it is in build 37, so we do not have to worry about that at all. 

raw_sumstat$chr_pos = paste("chr", raw_sumstat$chromosome, ":", raw_sumstat$base_pair_location, sep = "")

#And we can change the name for variant_id to variant

colnames(raw_sumstat)[which(colnames(raw_sumstat) == "variant_id")] = "variant" #easy 

head(raw_sumstat)
# variant p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency         beta standard_error sample_size    chr_pos
# 1: rs367896724    0.95          1              10177            AC            A                0.399875  0.000681779      0.0114655       32858 chr1:10177
# 2: rs201106462    0.52          1              10352            TA            T                0.394299  0.007605820      0.0118052       32858 chr1:10352
# 3: rs575272151    0.20          1              11008             G            C                0.085902 -0.025156200      0.0196611       32858 chr1:11008
# 4: rs544419019    0.20          1              11012             G            C                0.085902 -0.025156200      0.0196611       32858 chr1:11012
# 5:  rs62635286    0.61          1              13116             G            T                0.188596  0.007998440      0.0155328       32858 chr1:13116
# 6:  rs62028691    0.61          1              13118             G            A                0.188596  0.007998440      0.0155328       32858 chr1:13118

#######################################################################
#STEP 2: remove insertions and deletions, rare variants and MHC region#
#######################################################################

#First let's remove indels:

yes_vect <- c("A", "G", "T", "C")

raw_sumstat_no_indel <- raw_sumstat[which(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect),] #we go from 9M to 8M
#raw_sumstat_indel <- raw_sumstat[which(!(raw_sumstat$effect_allele%in%yes_vect & raw_sumstat$other_allele%in%yes_vect)),] #we did it correctly

#Second let's remove rare variant:

summary(raw_sumstat_no_indel$effect_allele_frequency)

#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000166 0.039396 0.147175 0.253808 0.403627 0.999786 

raw_sumstat_no_indel_common = raw_sumstat_no_indel[which(as.numeric(raw_sumstat_no_indel$effect_allele_frequency) > 0.01 & as.numeric(raw_sumstat_no_indel$effect_allele_frequency) < 0.99),]

summary(raw_sumstat_no_indel_common$effect_allele_frequency)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01000 0.04219 0.15226 0.25675 0.40792 0.99000 

#Third, let's remove extended MHC region

raw_sumstat_no_indel_common_mhc = raw_sumstat_no_indel_common[which(raw_sumstat_no_indel_common$chromosome == 6 & raw_sumstat_no_indel_common$base_pair_location >= 26000000 & raw_sumstat_no_indel_common$base_pair_location <= 34000000),]

summary(raw_sumstat_no_indel_common_mhc$chromosome) #perfect
summary(raw_sumstat_no_indel_common_mhc$base_pair_location) #perfect

raw_sumstat_no_indel_common_mhc_removed = raw_sumstat_no_indel_common[which(!(raw_sumstat_no_indel_common$chr_pos%in%raw_sumstat_no_indel_common_mhc$chr_pos)),]

#Let's be one sure

length(raw_sumstat_no_indel_common$chr_pos) - length(raw_sumstat_no_indel_common_mhc_removed$chr_pos) #perfect

####################################
#STEP 3: save the data! We are done#
####################################

fwrite(raw_sumstat_no_indel_common_mhc_removed, "output/1_curated_data/liver_fat_perc_liu_full_sumstats.txt")
