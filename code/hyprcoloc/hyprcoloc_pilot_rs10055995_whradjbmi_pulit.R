library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(hyprcoloc)

trait_aligner <- function(query_ss, other_ss) {

  check <- which(colnames(query_ss) == "effect_allele_frequency")
  if (is_empty(check)) query_ss$effect_allele_frequency <- NA

  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele,
           effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                          "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "samplesize.exposure", "rsid.exposure")
  exposure$exposure <- "fat_distr"
  exposure$id.exposure <- "fat_distr"

  check <- which(colnames(other_ss) == "effect_allele_frequency")
  if (is_empty(check)) other_ss$effect_allele_frequency <- NA

  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele,
           effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                         "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
                         "pval.outcome", "samplesize.outcome", "rsid.outcome")
  outcome$outcome <- "Outcome"
  outcome$id.outcome <- "Outcome"

  merged_df <- harmonise_data(exposure, outcome, action = 3)
  merged_df <- merged_df[which(merged_df$remove == FALSE), ]

  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome",
           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome",
           "pval.outcome", "samplesize.outcome", "SNP")
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location", "effect_allele",
                                  "other_allele", "effect_allele_frequency", "beta",
                                  "standard_error", "p_value", "sample_size", "chr_pos")
  return(other_ss_aligned)
}

# Step 0.1 -- chr/pos info of lead SNP
TARGET_LEAD <- "rs10055995"

PROJECT_ROOT <- "C:/Users/Bora Ceylan/Desktop/CBMR_Works_CLAUDE_works_in_there"  

lead_info_path <- file.path(PROJECT_ROOT, "merged_34_traits_ld_clumping/output/2nd_ld_clumping_of_merged_34_traits.txt")
lead_master <- fread(lead_info_path)
setnames(lead_master,
         c("rsid", "pval", "chr", "pos", "A1", "A2"),
         c("variant", "p_value", "chromosome", "base_pair_location", "effect_allele", "other_allele"))

lead_row <- lead_master[variant == TARGET_LEAD]
lead_chr <- as.numeric(lead_row$chromosome)
lead_pos <- as.numeric(lead_row$base_pair_location)   # build 37

# Step 0.2 -- genes and traits of target lead
overlap_report_path <- file.path(PROJECT_ROOT, "gene_day_distribution_with_overlap_report/output/qtl_sgbs_overlap/qtl_sgbs_overlap_report.xlsx")
overlap_report <- readxl::read_xlsx(overlap_report_path)

lead_rows_in_report <- overlap_report %>% filter(lead_snp == TARGET_LEAD)

genes_for_lead <- lead_rows_in_report %>% pull(ENSG) %>% unique()

# One lead may have multiple traits, this section automate that.
traits_for_lead <- lead_rows_in_report %>%
  pull(traits) %>%
  unique() %>%
  str_split(",\\s*") %>%
  unlist() %>%
  unique()

trait_sumstats_paths <- c(
  whradjbmi_pulit = file.path(PROJECT_ROOT, "whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt")
)
stopifnot(all(traits_for_lead %in% names(trait_sumstats_paths)))

# AdipoExpress marginal eQTL data -- ftp://mohlkeanon:anon@rc-ns-ftp.its.unc.edu/marginal_byChr_EURonly.tar.gz
qtl_chr_path <- file.path(
  PROJECT_ROOT,
  paste0("hyprcoloc/input/marginal_byChr_EURonly/EURonly_marginal_local_eQTL_meta_chr", lead_chr, ".txt")
)
qtl_tmp <- fread(qtl_chr_path)
qtl_tmp <- as.data.frame(qtl_tmp)

# Step 1+2 -- per trait cut/align/dedup exposure window, then run hyprcoloc per gene (no liftover, build 37)
start_ <- lead_pos - 500000
end_   <- lead_pos + 500000

coloc_df <- NULL

for (trait_name in traits_for_lead) {

  cat("\n=== TRAIT:", trait_name, "===\n")

  exposure_ss <- fread(trait_sumstats_paths[[trait_name]])
  # Rename raw GIANT/UKBB columns to the expected schema (Pulit format).
  setnames(exposure_ss,
           c("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele", "BETA", "SE", "P", "N"),
           c("chromosome", "base_pair_location", "variant", "effect_allele", "other_allele",
             "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size"),
           skip_absent = TRUE)

  exposure_tmp <- exposure_ss[which(as.numeric(exposure_ss$chromosome) == lead_chr &
                                      as.numeric(exposure_ss$base_pair_location) > start_ &
                                      as.numeric(exposure_ss$base_pair_location) < end_), ]

  # Align to positive beta + dedup by variant (avoids duplicate snp.id in hyprcoloc)
  exposure_pos <- exposure_tmp
  new_a1   <- ifelse(exposure_pos$beta < 0, exposure_pos$other_allele, exposure_pos$effect_allele)
  new_a2   <- ifelse(exposure_pos$beta < 0, exposure_pos$effect_allele, exposure_pos$other_allele)
  new_beta <- ifelse(exposure_pos$beta < 0, as.numeric(exposure_pos$beta) * (-1), as.numeric(exposure_pos$beta))
  exposure_pos$effect_allele <- new_a1
  exposure_pos$other_allele  <- new_a2
  exposure_pos$beta          <- new_beta

  exposure_pos <- exposure_pos[order(exposure_pos$p_value), ]
  exposure_pos <- exposure_pos[which(duplicated(exposure_pos$variant) == FALSE), ]
  exposure_pos <- exposure_pos[order(exposure_pos$chromosome, exposure_pos$base_pair_location), ]

  for (gene_ in genes_for_lead) {

    gene_tmp <- qtl_tmp[which(qtl_tmp$ENSG == gene_), ]
    if (is_empty(gene_tmp$ENSG)) next()

    gene_tmp$chr_pos <- paste0("chr", gene_tmp$chr, ":", gene_tmp$pos)

    exp_match <- exposure_pos %>%
      dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele,
                    effect_allele_frequency, beta, standard_error, p_value, sample_size)
    exp_match$chr_pos <- paste0("chr", exp_match$chromosome, ":", exp_match$base_pair_location)

    colnames(gene_tmp) <- c("chromosome", "base_pair_location", "other_allele", "effect_allele",
                            "variant_id", "ENSG", "gene", "studies", "beta", "standard_error",
                            "p_value", "sample_size", "chr_pos")
    gene_tmp$variant <- gene_tmp$chr_pos

    exp_match <- exp_match[which(exp_match$chr_pos %in% gene_tmp$chr_pos), ]

    if (nrow(exp_match) < 2) {
      cat("Gene:", gene_, "| fewer than 2 shared SNPs in region, hyprcoloc cannot run, skipping\n")
      next()
    }

    gene_aligned <- trait_aligner(exp_match, gene_tmp)

    exp_match <- exp_match[which(exp_match$chr_pos %in% gene_aligned$chr_pos), ]
    gene_aligned <- trait_aligner(exp_match, gene_aligned)  # drops leftover duplicates

    exp_ordered  <- exp_match[order(as.numeric(exp_match$chromosome), as.numeric(exp_match$base_pair_location)), ]
    gene_ordered <- gene_aligned[order(match(gene_aligned$chr_pos, exp_ordered$chr_pos)), ]

    cat("Gene:", gene_, "| matched positions:",
        length(which(exp_ordered$chr_pos == gene_ordered$chr_pos)),
        "| matched effect alleles:",
        length(which(exp_ordered$effect_allele == gene_ordered$effect_allele)), "\n")

    if (nrow(exp_ordered) < 2) next()

    betas <- as.matrix(as.data.frame(cbind(as.numeric(exp_ordered$beta), as.numeric(gene_ordered$beta))))
    ses   <- as.matrix(as.data.frame(cbind(as.numeric(exp_ordered$standard_error), as.numeric(gene_ordered$standard_error))))

    traits_ <- c(trait_name, gene_)
    colnames(betas) <- traits_
    colnames(ses)   <- traits_
    rownames(betas) <- exp_ordered$variant
    rownames(ses)   <- exp_ordered$variant

    res <- hyprcoloc::hyprcoloc(effect.est = betas, effect.se = ses,
                                snp.id = rownames(betas), trait.names = traits_)
    res_df <- res$results
    res_df$lead_snp <- TARGET_LEAD
    res_df$trait    <- trait_name
    res_df$gene     <- gene_

    print(res_df)

    coloc_df <- if (is.null(coloc_df)) res_df else rbind(coloc_df, res_df)
  }
}

# Step 3 -- save result
HYPRCOLOC_OUT_DIR <- file.path(PROJECT_ROOT, "hyprcoloc/output")
out_path <- file.path(HYPRCOLOC_OUT_DIR, paste0("hyprcoloc_pilot_", TARGET_LEAD, ".RDS"))
saveRDS(coloc_df, out_path)
cat("\nDone. Result:", out_path, "\n")