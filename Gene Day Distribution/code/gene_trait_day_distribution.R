#!/usr/bin/env Rscript
# =============================================================================
#  Lead SNP x Trait x ABC-Day Distribution
#
#  Mario's task 1: for each of the 1,603 lead SNPs, (a) flag which of the 34
#  traits it is genome-wide significant for (P<5e-8), (b) flag whether any of
#  its LD proxies overlap an ABC peak at day0 / day4 / day14, then use this to
#  answer: which day has more genes, and are WHRadjBMI SNPs enriched at day14.
#
#  STEP 1  Load the 1,603 lead SNPs (2nd LD-clumping output).
#  STEP 2  For each of the 34 traits, load its own lead-SNP list
#            (individual_traits/) and flag membership -> 34 binary columns.
#            NOTE: each trait's file is already that trait's own genome-wide
#            significant (P<5e-8), LD-clumped lead-SNP list, so membership
#            IS the P<5e-8 flag Mario asked for.
#  STEP 3  Sanity check against the pre-clump "winning trait" tag that
#            survived the 34-trait merge (merged_34_traits_with_lowest_pval.txt).
#  STEP 4  Add 3 binary columns: does this lead SNP have a proxy overlapping
#            an ABC peak at day0 / day4 / day14 (from ABC_ATAC_Analysis output)?
#  STEP 5  Summarize: which day has more genes? (raw ABC overlap vs eQTL-ABC
#            cross-validated gene counts, per day)
#  STEP 6  Enrichment test: are WHRadjBMI SNPs enriched for ABC overlap at
#            day14 vs day0/day4 (Fisher's exact test)?
#  STEP 7  Which traits have more genes? (distinct genes per trait, via raw
#            ABC overlap and via the eQTL-ABC cross-validated 74-gene set)
#  STEP 8  Write the summary report.
#
#  Required input files:
#    - 2nd_ld_clumping_of_merged_34_traits.txt          (34_merged_traits_LD_clumping/output)
#    - merged_34_traits_with_lowest_pval.txt             (34_merged_traits_LD_clumping/input)
#    - individual_traits/*_leads.txt                     (34 per-trait lead SNP lists)
#    - day0/4/14_proxy_abc_overlap.tsv                   (ABC_ATAC_Analysis/output)
#    - day0/4/14_eqtl_abc_cross_validation_summary.tsv   (ABC_ATAC_Analysis/output)
# =============================================================================

library(data.table)

# ---------------------------------------------------------------------------
# PATHS — adjust base_path if you run this on a different machine
# ---------------------------------------------------------------------------

base_path   <- "/home/boraceylan/Desktop"
ld_file     <- file.path(base_path, "output/2nd_ld_clumping_of_merged_34_traits.txt")
premerge_file <- file.path(base_path, "input/merged_34_traits_with_lowest_pval.txt")
trait_dir   <- file.path(base_path, "input/individual_traits")
abc_out_dir <- file.path(base_path, "ATAC_seq/sgbs_day_0_4_14/abc_analysis/output_EN")
out_dir     <- file.path(base_path, "gene_trait_day_distribution/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# STEP 1 — Lead SNPs
# =============================================================================

cat("\n── STEP 1: Loading lead SNPs ──\n")

lead <- fread(ld_file)
setnames(lead, "rsid", "lead_snp")
cat("  Lead SNPs:", nrow(lead), "\n")

matrix_dt <- lead[, .(lead_snp)]


# =============================================================================
# STEP 2 — Per-trait significance columns (P<5e-8)
# =============================================================================

cat("\n── STEP 2: Per-trait significance (from individual_traits/) ──\n")

trait_files <- list.files(trait_dir, pattern = "_leads\\.txt$", full.names = TRUE)
cat("  Trait files found:", length(trait_files), "\n")

trait_cols <- character(0)
for (f in trait_files) {
  trait_name <- sub("_leads\\.txt$", "", basename(f))
  df <- fread(f)
  stopifnot("variant" %chin% names(df))
  matrix_dt[, (trait_name) := as.integer(lead_snp %chin% df$variant)]
  trait_cols <- c(trait_cols, trait_name)
}
cat("  Trait columns added:", length(trait_cols), "\n")


# =============================================================================
# STEP 3 — Sanity check against the pre-clump "winning trait" tag
#
# Each lead SNP survived from merged_34_traits_with_lowest_pval.txt, where it
# was tagged with exactly ONE trait (the one with the lowest p-value across
# traits, at merge time). That trait's column in our new matrix must be 1.
# =============================================================================

cat("\n── STEP 3: Sanity check ──\n")

premerge <- fread(premerge_file)
premerge[, trait := sub("_leads$", "", trait)]
win_map <- premerge[variant %chin% matrix_dt$lead_snp, .(lead_snp = variant, trait)]

mismatches <- 0L
for (i in seq_len(nrow(win_map))) {
  snp <- win_map$lead_snp[i]
  tr  <- win_map$trait[i]
  if (tr %chin% trait_cols) {
    flag <- matrix_dt[lead_snp == snp][[tr]]
    if (flag != 1) mismatches <- mismatches + 1L
  }
}
cat("  Checked:", nrow(win_map), "lead SNPs | Mismatches (should be 0):", mismatches, "\n")


# =============================================================================
# STEP 4 — Proxy-ABC overlap flags per day
# =============================================================================

cat("\n── STEP 4: Proxy-ABC overlap flags ──\n")

for (day in c("day0", "day4", "day14")) {
  ov <- fread(file.path(abc_out_dir, paste0(day, "_proxy_abc_overlap.tsv")))
  overlap_snps <- unique(ov$lead_snp)
  matrix_dt[, (paste0("overlap_", day)) := as.integer(lead_snp %chin% overlap_snps)]
  cat("  ", day, ": ", length(overlap_snps), " lead SNPs with proxy-ABC overlap\n", sep = "")
}

matrix_dt[, n_traits := rowSums(.SD), .SDcols = trait_cols]

fwrite(matrix_dt, file.path(out_dir, "lead_snp_trait_day_matrix.tsv"), sep = "\t")
cat("\n  Saved:", file.path(out_dir, "lead_snp_trait_day_matrix.tsv"), "\n")


# =============================================================================
# STEP 5 — Which day has more genes?
# =============================================================================

cat("\n── STEP 5: Which day has more genes? ──\n")

cat("\n  Raw ABC overlap (any proxy in an ABC peak, regardless of eQTL evidence):\n")
raw_counts <- list()
for (day in c("day0", "day4", "day14")) {
  ov <- fread(file.path(abc_out_dir, paste0(day, "_proxy_abc_overlap.tsv")))
  raw_counts[[day]] <- uniqueN(ov$gene_name)
  cat("   ", day, ":", raw_counts[[day]], "genes\n")
}

cat("\n  Cross-validated (eQTL ∩ ABC) genes:\n")
cv_counts <- list()
for (day in c("day0", "day4", "day14")) {
  cv <- fread(file.path(abc_out_dir, paste0(day, "_eqtl_abc_cross_validation_summary.tsv")))
  cv_counts[[day]] <- uniqueN(cv$gene_name)
  cat("   ", day, ":", cv_counts[[day]], "genes\n")
}

day_gene_counts <- data.table(
  day               = c("day0", "day4", "day14"),
  raw_abc_genes     = unlist(raw_counts),
  cross_validated_genes = unlist(cv_counts)
)
fwrite(day_gene_counts, file.path(out_dir, "day_gene_counts.tsv"), sep = "\t")


# =============================================================================
# STEP 6 — WHRadjBMI enrichment by day (Fisher's exact test)
# =============================================================================

cat("\n── STEP 6: WHRadjBMI enrichment by day ──\n")

stopifnot("whradjbmi_pulit" %chin% names(matrix_dt))
enrichment_rows <- list()
for (day in c("day0", "day4", "day14")) {
  ov_col <- paste0("overlap_", day)
  tab <- table(matrix_dt$whradjbmi_pulit, matrix_dt[[ov_col]])
  ft <- fisher.test(tab)
  cat("  ", day, ": OR=", round(unname(ft$estimate), 3),
      " p=", format(ft$p.value, digits = 4), "\n", sep = "")
  enrichment_rows[[day]] <- data.table(
    day = day,
    OR = unname(ft$estimate),
    p_value = ft$p.value,
    n_whr_and_overlap = tab["1", "1"],
    n_whr_no_overlap = tab["1", "0"],
    n_not_whr_overlap = tab["0", "1"],
    n_neither = tab["0", "0"]
  )
}
enrichment_dt <- rbindlist(enrichment_rows)
fwrite(enrichment_dt, file.path(out_dir, "whradjbmi_day_enrichment.tsv"), sep = "\t")


# =============================================================================
# STEP 7 — Which traits have more genes?
#
# For each trait, take the lead SNPs flagged 1 for that trait, then count the
# distinct genes linked to them via (a) any proxy-ABC overlap (any day) and
# (b) the eQTL x ABC cross-validated (74-gene) set.
# =============================================================================

cat("\n── STEP 7: Which traits have more genes? ──\n")

raw_pairs <- rbindlist(lapply(c("day0", "day4", "day14"), function(day) {
  ov <- fread(file.path(abc_out_dir, paste0(day, "_proxy_abc_overlap.tsv")))
  unique(ov[, .(lead_snp, gene_name)])
}))
raw_pairs <- unique(raw_pairs)

cv_pairs <- unique(fread(file.path(abc_out_dir, "temporal_cross_validated_genes.tsv"))[, .(lead_snp, gene_name)])

trait_gene_counts <- rbindlist(lapply(trait_cols, function(trait) {
  snps <- matrix_dt[get(trait) == 1, lead_snp]
  data.table(
    trait = trait,
    n_lead_snps = length(snps),
    n_raw_abc_genes = uniqueN(raw_pairs[lead_snp %chin% snps, gene_name]),
    n_cross_validated_genes = uniqueN(cv_pairs[lead_snp %chin% snps, gene_name])
  )
}))
setorder(trait_gene_counts, -n_raw_abc_genes)

fwrite(trait_gene_counts, file.path(out_dir, "trait_gene_counts.tsv"), sep = "\t")
cat("  Top 5 traits by raw ABC-overlap gene count:\n")
print(head(trait_gene_counts, 5))


# =============================================================================
# STEP 8 — Summary report
# =============================================================================

lines <- c(
  "============================================================",
  "  Lead SNP x Trait x ABC-Day Distribution — Summary",
  paste("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "============================================================", "",
  "── INPUTS ──",
  paste("  Lead SNPs:", nrow(matrix_dt)),
  paste("  Traits:", length(trait_cols)),
  paste("  Sanity check mismatches (should be 0):", mismatches), "",
  "── N TRAITS PER LEAD SNP ──",
  capture.output(print(matrix_dt[, .N, by = n_traits][order(n_traits)])), "",
  "── WHICH DAY HAS MORE GENES ──",
  capture.output(print(day_gene_counts)), "",
  "── WHRadjBMI ENRICHMENT BY DAY ──",
  capture.output(print(enrichment_dt[, .(day, OR, p_value)])), "",
  "── TOP 10 TRAITS BY GENE COUNT ──",
  capture.output(print(head(trait_gene_counts, 10)))
)

writeLines(lines, file.path(out_dir, "pipeline_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

cat("\n============================================================\n")
cat("  Output folder:", out_dir, "\n")
cat("============================================================\n")
