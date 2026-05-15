library(data.table)
library(LDlinkR)

# Paths
gwas_path <- "/home/boraceylan/Desktop/output/2nd_ld_clumping_of_merged_34_traits.txt"
eqtl_path <- "/home/boraceylan/Desktop/input/SuppTable4_EURonly_eQTL_signals.txt"

# Token
Sys.setenv(LDLINK_TOKEN = "7e9ecefb4951")
token <- Sys.getenv("LDLINK_TOKEN")

# 1) Lead SNPs from clumped GWAS (use SNP column if present, otherwise use rsid/variant)
gwas <- fread(gwas_path, fill = TRUE)
lead_col <- if ("SNP" %chin% names(gwas)) "SNP" else {
  cand <- intersect(c("variant","rsid","RSID","snp","Snp"), names(gwas))
  if (length(cand) == 0) stop("No SNP/rsid/variant column found in the GWAS file.")
  cand[1]
}
lead_snps <- unique(gwas[!is.na(get(lead_col)) & get(lead_col) != "", get(lead_col)])

# 2) LD proxies (EUR, GRCh37), keep R2 >= 0.8
get_proxies <- function(s) {
  x <- tryCatch(
    LDproxy(snp = s, pop = "EUR", r2d = "r2", token = token, genome_build = "grch37"),
    error = function(e) NULL
  )
  if (is.null(x)) return(NULL)
  x <- as.data.table(x)
  x[, query := s]
  x
}

proxies <- rbindlist(lapply(lead_snps, get_proxies), fill = TRUE)
proxies[, R2 := as.numeric(R2)]
proxies <- proxies[is.finite(R2) & R2 >= 0.8]

# 3) eQTL (SuppTable4): filter and build coord as chr:pos to match LDproxy Coord
eqtl <- fread(eqtl_path)
eqtl <- eqtl[pval_joint <= 1e-6]

parts <- tstrsplit(eqtl$variant, "_", fixed = TRUE)
eqtl[, coord := paste0("chr", parts[[1]], ":", parts[[2]])]

eqtl_small <- unique(eqtl[, .(ENSG, gene, signal, eqtl_variant = variant, coord, pval_joint)])

# 4) Match: query SNP -> proxy (Coord) -> eQTL gene
hits <- merge(proxies, eqtl_small, by.x = "Coord", by.y = "coord", allow.cartesian = TRUE)
final <- hits[, .(query, RS_Number, Coord, R2, ENSG, gene, signal, pval_joint, eqtl_variant)]

fwrite(final, "/home/boraceylan/Desktop/output/GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv", sep = "\t")
fwrite(proxies, "/home/boraceylan/Desktop/output/all_proxies_r2ge_08.tsv")
