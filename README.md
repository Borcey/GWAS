# GWAS
Bora Ceylan's GWAS studies and further

# "code" section
- ld_clumping = LD Clumping codes for LD-clumped GWAS'S. Also, includes rest of the GWAS's.
- cleaning_gwas_summary_stats = p-value filtering, EAF filtering, MHC cleaning, positive beta aligning and allele exchanging based on the beta signs, chr_pos addition.
- cleaning_supplementary_tables = EAF filtering, MHC cleaning, positive beta aligning and allele exchanging based on the beta signs, chr_pos addition


# "output" section for code section's output
- curated_supplementary_tables = All of the necessary curated versions in the curated_supplementary_tables, not the inner files.
- cleaned_gwas_sum_stats = Only "leg_fat_percentage" and "adipokines" are in separated files. Rest of them in .txt file
- ld_clumped_gwas_sum_stats ;
    - lead_snp_gwas = Whole necessary LD-clumping results
    - trait_by_trait_ld_clumped_ALL_GWASS_CLUMPED = Includes every LD-clumping results for each trait (Unnecessary in general)
    - merged_clumped_traits_second_ld_clumped.txt = Merged version of necessary LD-clumping results (Just includes LD-clumped GWAS'S results, not included supplementary tables' lead SNPs)
    - merged_clumped_traits.txt = merging of LD-clumping needed GWAS'S traits
 

# "input" section
- clean_supplementary_tables
- raw_supplementary_tables
- raw_gwas_sum_stats ;
    - ready_to_use = unziped versions of GWAS's (.txt, .tsv, or just file format)
    - zip_versions = the rawest version of GWAS'S


# "plink_files" 
The file where 2 of the LD-clumping's parameters come from.
