# Merging 34 traits and 2nd LD-Clumping

library(data.table)

trait_dir <- "/home/boraceylan/Desktop/34_traits_ld/curated_data"
files <- list.files(trait_dir, pattern="\\.(txt|tsv|csv)$", full.names=TRUE, ignore.case=TRUE)

std_names <- function(x) tolower(gsub("[^a-z0-9]+", "_", x))

trait_from_file <- function(f){
  sub("\\.(txt|tsv|csv)$", "", basename(f), ignore.case = TRUE)
}

read_one_with_p <- function(f) {
  dt <- fread(f, showProgress = FALSE)
  
  old <- names(dt)
  new <- std_names(old)
  setnames(dt, old, new)
  
  p_candidates <- c("p_value", "pval", "pvalue", "p")
  p_col <- intersect(p_candidates, names(dt))[1]
  
  if (is.na(p_col)) {
    stop(sprintf("P column is not found. File: %s\nAvailable columns: %s",
                 f, paste(names(dt), collapse=", ")))
  }
  
  need <- c("chromosome","base_pair_location","variant","effect_allele","other_allele")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop(sprintf("Absent column in (%s): %s\nAvailable columns are in: %s",
                 paste(miss, collapse=", "),
                 f,
                 paste(names(dt), collapse=", ")))
  }
  
  out <- dt[, .(
    chromosome,
    base_pair_location,
    variant,
    effect_allele,
    other_allele,
    p_value = as.numeric(get(p_col))
  )]
  
  out[p_value == 0, p_value := 1e-300]
  out <- out[is.finite(p_value) & !is.na(p_value)]
  
  out[, trait := trait_from_file(f)]
  out
}


merged_p <- rbindlist(lapply(files, read_one_with_p), use.names=TRUE, fill=TRUE)


cat("OK files:", length(files), "\n")
cat("Merged rows:", nrow(merged_p), "\n")
cat("Unique trait:", uniqueN(merged_p$trait), "\n")
cat("p_value number of NA :", sum(is.na(merged_p$p_value)), "\n")

fwrite(merged_p, "/home/boraceylan/Desktop/34_traits_ld/ALL34_merged_with_p.tsv")
