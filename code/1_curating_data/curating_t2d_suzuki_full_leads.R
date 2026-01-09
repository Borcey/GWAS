##############
#INTRODUCTION#
##############

#Let's get the the lipodystrophy-like cluster for T2D in Suzuki et al.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/perturb_seq/GWAS-main/GWAS-main/")

lead_snps <- as.data.frame(readxl::read_xlsx("input/raw_supplementary_tables/41586_2024_7019_MOESM3_ESM.xlsx", sheet = 6)) #we read it properly!!

colnames(lead_snps)=c("Locus", "chromosome", "variant", "base_pair_location", "cluster", "distance")

#Let's take the rest:

lead_snps_clean = lead_snps[4:length(lead_snps$Locus),] #seems like all data.

#Let's get the lipodystrophy ones, we have some issues with the chromosomes:

lipodystrophy=lead_snps_clean[which(lead_snps_clean$cluster=="Lipodystrophy"),]

#We could do this with an algorhythm... but it is easier to just do it ourselves

chr_2= c("rs10184004", "rs9287852", "rs78058190", "rs12613239", "rs2203452") #double-checked all is good! Done
chr_3= c("rs308958",  "rs17036160", "rs4132228", "rs62271373")
chr_4= c("rs4450871", "rs7697644",   "rs4691378")
chr_5= c("rs9687846",  "rs4976033",  "rs2963470")
chr_6= c("rs10214450",  "rs6905288",   "rs668459",   "rs4709746")
chr_7= c("rs1534696",  "rs4731702",  "rs62492368")
chr_9= c("rs1044531")
chr_10= c("rs1907225")
chr_11= c("rs28599796")
chr_12= c("rs10842708", "rs34680764")
chr_17= c("rs56745821")
chr_19= c("rs2115107", "rs116843064")
chr_20= c("rs1800961", "rs6063046")


lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_2)] <- 2
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_3)] <- 3
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_4)] <- 4
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_5)] <- 5
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_6)] <- 6
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_7)] <- 7
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_9)] <- 9
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_10)] <- 10
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_11)] <- 11
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_12)] <- 12
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_17)] <- 17
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_19)] <- 19
lipodystrophy$chromosome[which(lipodystrophy$variant%in%chr_20)] <- 20

lipodystrophy$chr_pos = paste("chr", lipodystrophy$chromosome, ":", lipodystrophy$base_pair_location, sep = "")

#Let's get the data from the sumstats:

t2d_full_sumstats=fread("output/1_curated_data/t2d_suzuki_full_sumstats.txt")

t2d_match = t2d_full_sumstats[which(t2d_full_sumstats$chr_pos%in%lipodystrophy$chr_pos),]

#Let's organize the data:

t2d_match <- t2d_match[order(match(t2d_match$chr_pos, lipodystrophy$chr_pos)),]

length(which(t2d_match$chr_pos == lipodystrophy$chr_pos)) #45

t2d_match$variant=lipodystrophy$variant

fwrite(t2d_match, "output/1_curated_data/t2d_suzuki_leads.txt")
