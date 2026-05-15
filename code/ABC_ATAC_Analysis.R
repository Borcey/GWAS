#!/usr/bin/env Rscript
# =============================================================================
#  ABC / ATAC-seq Overlap + eQTL Cross-Validation
#
#  What this script does (step by step):
#
#  STEP 1  Load the proxy SNP table and eQTL gene map produced earlier.
#  STEP 2  For each ABC timepoint (Day 0 / 4 / 14):
#            - Read the ABC file (ATAC-seq peaks linked to genes by ABC score)
#            - Find every proxy whose genomic position falls inside a peak
#            - Save the overlapping rows
#  STEP 3  Cross-validation: keep only genes that appear in BOTH
#            the ABC overlap (proxy → peak → gene)  AND
#            the eQTL gene map (proxy → eQTL → gene)
#  STEP 4  Temporal comparison: classify genes by which timepoints they appear in.
#  STEP 5  Print a summary and write all outputs.
#
#  Required input files:
#    - all_proxies_r2ge0.8.tsv                  (from LDlinkR step)
#    - GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv (from eQTL matching step)
#    - sgbs_day_0_ABCpp_scoredInteractions_value_.txt
#    - sgbs_day_4_ABCpp_scoredInteractions_value_.txt
#    - sgbs_day_14_ABCpp_scoredInteractions_value_.txt
# =============================================================================

library(data.table)

# ---------------------------------------------------------------------------
# PATHS  — adjust base_path if you run this on a different machine
# ---------------------------------------------------------------------------

base_path <- "/home/boraceylan/Desktop/ATAC_seq"

proxy_file <- file.path(base_path, "all_proxies_r2ge0.8.tsv")
eqtl_file  <- file.path(base_path, "GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv")

abc_files <- list(
  day0  = file.path(base_path, "sgbs_day_0_4_14", "sgbs_day_0_ABCpp_scoredInteractions_value_.txt"),
  day4  = file.path(base_path, "sgbs_day_0_4_14", "sgbs_day_4_ABCpp_scoredInteractions_value_.txt"),
  day14 = file.path(base_path, "sgbs_day_0_4_14", "sgbs_day_14_ABCpp_scoredInteractions_value_.txt")
)

out_dir <- file.path(base_path, "sgbs_day_0_4_14", "abc_analysis", "output_EN")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# STEP 1 — Load input files
# =============================================================================

cat("\n── STEP 1: Loading input files ──\n")

# --- 1a. Proxy table ---------------------------------------------------------
# Columns we need:
#   query      = lead SNP (rs-number from LD clumping)
#   RS_Number  = proxy SNP rs-number
#   Coord      = proxy position in format "chr4:89752913"
#   R2         = LD r² between lead and proxy (already filtered ≥ 0.8)

proxies <- fread(proxy_file)
cat("  Proxies loaded:", nrow(proxies), "rows,", uniqueN(proxies$query), "lead SNPs\n")

# Parse Coord ("chr4:89752913") → two new columns: proxy_chr ("4") and proxy_pos (89752913)
# We strip the "chr" prefix so it matches the ABC file, which uses plain numbers (1, 2, 3…)
proxies[, c("proxy_chr", "proxy_pos") := {
  stripped <- sub("^chr", "", Coord)          # "chr4:89752913" → "4:89752913"
  parts    <- tstrsplit(stripped, ":", fixed = TRUE)
  list(parts[[1]], as.integer(parts[[2]]))    # chr = "4", pos = 89752913
}]

proxies <- proxies[!is.na(proxy_pos)]         # drop any rows where parsing failed
cat("  After coordinate parsing:", nrow(proxies), "rows remain\n")

# --- 1b. eQTL gene map -------------------------------------------------------
# This file was produced by matching proxies to adipose eQTL signals (pval_joint ≤ 1e-6).
# Columns we care about:
#   query      = lead SNP
#   gene       = gene name (e.g. "PPARG")
#   ENSG       = Ensembl gene ID
#   pval_joint = eQTL p-value

eqtl_map <- fread(eqtl_file)
cat("  eQTL gene map loaded:", nrow(eqtl_map), "rows,",
    uniqueN(eqtl_map$gene), "genes,", uniqueN(eqtl_map$query), "lead SNPs\n")

# --- 1c. Check which ABC files exist -----------------------------------------
cat("\n  ABC file check:\n")
for (label in names(abc_files)) {
  status <- if (file.exists(abc_files[[label]])) "OK" else "NOT FOUND"
  cat("    ", label, ":", status, "\n")
}
available_abc <- abc_files[sapply(abc_files, file.exists)]
if (length(available_abc) == 0) stop("No ABC files found — check the paths above.")


# =============================================================================
# STEP 2 — ABC overlap: find proxies that fall inside ATAC-seq peaks
# =============================================================================

cat("\n── STEP 2: ABC overlap ──\n")

# This function processes one ABC file (one timepoint).
#
# The ABC file has columns:
#   #chr   = chromosome (plain number: 1, 2, …)
#   start  = peak start position
#   end    = peak end position
#   Gene Name  = gene linked to this peak by ABC model
#   ABC-Score  = strength of the enhancer-gene link (higher = stronger)
#
# A "hit" is when:  proxy_chr == abc_chr  AND  start <= proxy_pos <= end
#
# We use a data.table non-equi join (start <= proxy_pos AND end >= proxy_pos)
# which is faster than a row-by-row loop.

run_overlap <- function(abc_path, day_label) {

  cat("  [", day_label, "] Reading ABC file… ")
  abc <- fread(abc_path)

  # Rename "#chr" → "chr" (the "#" is just a file-format artefact)
  setnames(abc, "#chr", "chr", skip_absent = TRUE)
  abc[, chr := as.character(chr)]   # make sure both sides are character for joining

  cat(format(nrow(abc), big.mark = ","), "peaks |")

  # Non-equi join: for each proxy, find ABC rows where chr matches AND
  # start <= proxy_pos <= end.
  # nomatch = 0L means: only keep rows that actually matched (inner join).
  overlap <- abc[
    proxies,
    on = .(chr = proxy_chr, start <= proxy_pos, end >= proxy_pos),
    nomatch = 0L,
    allow.cartesian = TRUE   # one proxy can overlap many peaks
  ]

  if (nrow(overlap) == 0) {
    cat(" No overlaps found.\n")
    return(data.table())
  }

  # Select and rename to meaningful column names
  result <- overlap[, .(
    lead_snp     = query,       # the original lead SNP from LD clumping
    proxy_rsid   = RS_Number,   # the proxy that actually overlaps the peak
    proxy_coord  = Coord,       # chr:pos of the proxy
    proxy_r2     = R2,          # LD r² with the lead
    peak_id      = PeakID,
    gene_name    = `Gene Name`,
    ensembl_id   = `Ensembl ID`,
    abc_score    = `ABC-Score`,
    tss_distance = `TSS-dist`,
    timepoint    = day_label
  )]

  fwrite(result, file.path(out_dir, paste0(day_label, "_proxy_abc_overlap.tsv")), sep = "\t")
  cat("", format(nrow(result), big.mark = ","), "overlaps |",
      uniqueN(result$gene_name), "unique genes |",
      uniqueN(result$lead_snp), "lead SNPs\n")

  return(result)
}

# Run for all available timepoints
abc_results <- lapply(names(available_abc), function(lbl) run_overlap(available_abc[[lbl]], lbl))
names(abc_results) <- names(available_abc)

# Merge all timepoints into one table and save
abc_all <- rbindlist(abc_results, fill = TRUE)
if (nrow(abc_all) > 0) {
  fwrite(abc_all, file.path(out_dir, "all_timepoints_proxy_abc_overlap.tsv"), sep = "\t")
  cat("  Combined file saved:", nrow(abc_all), "rows\n")
}


# =============================================================================
# STEP 3 — Cross-validation: eQTL ∩ ABC
#
# Key question: which genes are supported by BOTH lines of evidence?
#   Evidence A: proxy overlaps an ATAC-seq peak → ABC model says that peak
#               regulates gene X
#   Evidence B: proxy is in LD with an eQTL for gene X
#
# We match on (lead_snp, gene_name) — same lead SNP AND same gene must appear
# in both tables.
# =============================================================================

cat("\n── STEP 3: Cross-validation (eQTL ∩ ABC) ──\n")

cv_results <- list()

for (label in names(abc_results)) {

  abc_r <- abc_results[[label]]
  if (nrow(abc_r) == 0) { cv_results[[label]] <- data.table(); next }

  # Keep only the columns needed for merging
  abc_slim  <- abc_r[,  .(lead_snp, gene_name, proxy_rsid, proxy_r2,
                           abc_score, peak_id, tss_distance, timepoint)]
  eqtl_slim <- eqtl_map[, .(lead_snp = query, gene_name = gene,
                              eqtl_pval = pval_joint, ENSG, eqtl_variant)]

  # Inner join: only rows present in BOTH tables
  cross <- merge(abc_slim, eqtl_slim, by = c("lead_snp", "gene_name"),
                 allow.cartesian = TRUE)

  if (nrow(cross) == 0) {
    cat("  [", label, "] No cross-validated genes.\n")
    cv_results[[label]] <- data.table()
    next
  }

  # Per-gene summary (best ABC score + best eQTL p-value across all proxies/peaks)
  summary_dt <- cross[, .(
    max_abc_score = max(abc_score),
    min_eqtl_pval = min(as.numeric(eqtl_pval)),
    n_peaks       = uniqueN(peak_id),
    n_proxies     = uniqueN(proxy_rsid),
    timepoint     = label
  ), by = .(lead_snp, gene_name)][order(-max_abc_score)]

  fwrite(cross,      file.path(out_dir, paste0(label, "_eqtl_abc_cross_validation_full.tsv")),    sep = "\t")
  fwrite(summary_dt, file.path(out_dir, paste0(label, "_eqtl_abc_cross_validation_summary.tsv")), sep = "\t")

  cat("  [", label, "]", uniqueN(cross$gene_name), "cross-validated genes |",
      uniqueN(cross$lead_snp), "lead SNPs\n")
  cat("  Top 5 genes by ABC score:\n")
  print(head(summary_dt[, .(lead_snp, gene_name, max_abc_score, min_eqtl_pval)], 5))
  cat("\n")

  cv_results[[label]] <- summary_dt
}


# =============================================================================
# STEP 4 — Temporal comparison
#
# Among cross-validated genes, classify each (lead_snp, gene) pair by when it
# shows up:
#   all_stages  = present at Day 0, Day 4, and Day 14  → constitutively active
#   early_only  = Day 0 but not Day 14                 → preadipocyte-specific
#   late_only   = Day 14 but not Day 0                 → mature adipocyte-specific
#   mid_only    = Day 4 only                           → transitional
#   other       = any other combination
# =============================================================================

cat("── STEP 4: Temporal comparison ──\n")

cv_all <- rbindlist(cv_results, fill = TRUE)

if (nrow(cv_all) > 0) {

  # Reshape to wide format: one row per (lead_snp, gene_name), one column per timepoint
  temporal <- dcast(
    cv_all,
    lead_snp + gene_name ~ timepoint,
    value.var = "max_abc_score",
    fun.aggregate = max,
    fill = NA_real_
  )

  day_cols <- intersect(c("day0", "day4", "day14"), names(temporal))

  # Helper: TRUE if a value is present and positive
  present <- function(col) {
    if (col %in% names(temporal)) !is.na(temporal[[col]]) & temporal[[col]] > 0
    else rep(FALSE, nrow(temporal))
  }

  temporal[, pattern := {
    d0 <- present("day0"); d4 <- present("day4"); d14 <- present("day14")
    fifelse(d0 & d4 & d14,      "all_stages",
    fifelse(d0 & !d14,          "early_only",
    fifelse(!d0 & d14,          "late_only",
    fifelse(d4 & !d0 & !d14,   "mid_only", "other"))))
  }]

  fwrite(temporal, file.path(out_dir, "temporal_cross_validated_genes.tsv"), sep = "\t")

  cat("\n  Temporal pattern counts:\n")
  print(temporal[, .N, by = pattern][order(-N)])

  cat("\n  Genes present at all three stages:\n")
  all_stages <- temporal[pattern == "all_stages"]
  if (nrow(all_stages) > 0) print(all_stages) else cat("  (none)\n")

} else {
  cat("  No cross-validated genes found across timepoints.\n")
}


# =============================================================================
# STEP 5 — Summary report
# =============================================================================

cat("\n── STEP 5: Summary ──\n\n")

lines <- c(
  "============================================================",
  "  ABC / ATAC-seq Overlap + eQTL Cross-Validation — Summary",
  paste("  Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "============================================================", "",
  "── INPUTS ──",
  paste("  Proxy SNPs (r² ≥ 0.8):", format(nrow(proxies), big.mark = ",")),
  paste("  Lead SNPs:", uniqueN(proxies$query)),
  paste("  eQTL gene map rows:", format(nrow(eqtl_map), big.mark = ",")),
  paste("  eQTL genes:", uniqueN(eqtl_map$gene)), "",
  "── ABC OVERLAP ──"
)

for (lbl in names(abc_results)) {
  r <- abc_results[[lbl]]
  lines <- c(lines, paste0("  ", lbl, ": ",
    if (nrow(r) > 0) paste(format(nrow(r), big.mark=","), "overlaps |",
                            uniqueN(r$gene_name), "genes") else "0 overlaps"))
}

lines <- c(lines, "", "── CROSS-VALIDATION ──")
for (lbl in names(cv_results)) {
  r <- cv_results[[lbl]]
  lines <- c(lines, paste0("  ", lbl, ": ",
    if (!is.null(r) && nrow(r) > 0) paste(uniqueN(r$gene_name), "genes") else "0 genes"))
}

if (exists("temporal") && nrow(temporal) > 0) {
  lines <- c(lines, "", "── TEMPORAL PATTERNS ──",
    capture.output(print(temporal[, .N, by = pattern][order(-N)])))
}

writeLines(lines, file.path(out_dir, "pipeline_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")

cat("\n============================================================\n")
cat("  Output folder:", out_dir, "\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
