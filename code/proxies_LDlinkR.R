library(data.table)
library(LDlinkR)
library(data.table)

dt <- fread("/home/boraceylan/Desktop/34_traits_ld/2nd_ld_clumped_34_merged_traits.txt")
setnames(dt, c("rsid","pval","chr","pos","A1","A2"),
         c("variant","p_value","chromosome","base_pair_location","effect_allele","other_allele"))

dt[, p_value := as.numeric(p_value)]
dt <- dt[is.finite(p_value) & !is.na(p_value)]

setkey(dt, variant, p_value)
dedup_minp <- dt[, .SD[1], by=variant]

# PLINK clump input (SNP, P is obligatory)  # This part is negligible, I put the clumped dataset.
clump_in <- dedup_minp[, .(SNP=variant, P=p_value, CHR=chromosome, BP=base_pair_location)]
fwrite(clump_in, "/home/boraceylan/Desktop/34_traits_ld/ALL34_for_clump.tsv", sep="\t")


cl <- fread("/home/boraceylan/Desktop/34_traits_ld/ALL34_for_clump.tsv", fill=TRUE)
lead_snps <- unique(cl$SNP)
lead_snps <- lead_snps[!is.na(lead_snps) & lead_snps != ""]


# LDlinkR Section
Sys.setenv(LDLINK_TOKEN = "7e9ecefb4951")
token <- Sys.getenv("LDLINK_TOKEN")


get_proxies <- function(s) {
  x <- tryCatch(LDproxy(snp=s, pop="EUR", r2d="r2", token=token, genome_build="grch37"),
                error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x <- as.data.table(x)
  x[, query := s]
  x
}

proxies <- rbindlist(lapply(lead_snps, get_proxies), fill=TRUE)
proxies[, R2 := as.numeric(R2)]
proxies <- proxies[is.finite(R2) & R2 >= 0.8]



eqtl <- fread("/home/boraceylan/Desktop/ST4_eQTL/SuppTable4_EURonly_eQTL_signals.txt")

eqtl <- eqtl[pval_joint <= 1e-6]

# S4 variant -> coord (chr:pos)
parts <- tstrsplit(eqtl$variant, "_")
eqtl[, coord := paste0("chr", parts[[1]], ":", parts[[2]])]

hits <- merge(proxies, unique(eqtl[, .(ENSG, gene, signal, variant, coord, pval_joint)]),
              by.x="Coord", by.y="coord", allow.cartesian=TRUE)

# Final: query SNP -> proxy -> gene
final <- hits[, .(query, RS_Number, Coord, R2, ENSG, gene, signal, pval_joint, eqtl_variant=variant)]
fwrite(final, "GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv", sep="\t")


