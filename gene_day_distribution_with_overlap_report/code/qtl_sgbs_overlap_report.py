#!/usr/bin/env python
# =============================================================================
#  Lead SNP x Gene report -- QTL (eQTL) AND SGBS (ABC) overlap
#
#  Mario's request: "a report with the leads AND genes that overlap with QTL
#  data and SGBS data". That overlap was already computed in the ABC/ATAC
#  cross-validation step (temporal_cross_validated_genes.tsv = genes with
#  BOTH eQTL evidence AND an ABC-ATAC peak in SGBS adipocytes). This script
#  just reshapes that result into a lead-SNP-centric report and attaches
#  which of the 34 traits each lead SNP is genome-wide significant for.
#
#  Grain: one row per (lead_snp, gene) pair.
#  Scope: only the 57 lead SNPs / 74 genes that have QTL+SGBS overlap.
#
#  Inputs:
#    - ABC_ATAC_Analysis/output/temporal_cross_validated_genes.tsv
#    - ABC_ATAC_Analysis/output/day{0,4,14}_eqtl_abc_cross_validation_summary.tsv
#    - ABC_ATAC_Analysis/input/GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv
#    - Gene_Trait_Day_Distribution/output/trait_matrix.xlsx
# =============================================================================

import pandas as pd

CBMR = r"c:\Users\Bora Ceylan\Desktop\CBMR_pipeline"
abc_out = rf"{CBMR}\ABC_ATAC_Analysis\output"
abc_in = rf"{CBMR}\ABC_ATAC_Analysis\input"
trait_matrix_xlsx = r"c:\Users\Bora Ceylan\Desktop\task2_plus\output\trait_matrix.xlsx"
out_dir = r"c:\Users\Bora Ceylan\Desktop\task2_plus\output"

# -----------------------------------------------------------------------
# STEP 1 -- cross-validated (QTL + SGBS) lead_snp x gene pairs, wide by day
# -----------------------------------------------------------------------

temporal = pd.read_csv(rf"{abc_out}\temporal_cross_validated_genes.tsv", sep="\t")
temporal = temporal.rename(columns={"day0": "day0_abc_score", "day4": "day4_abc_score", "day14": "day14_abc_score"})
print(f"STEP 1: {temporal.shape[0]} lead_snp x gene pairs | "
      f"{temporal.lead_snp.nunique()} lead SNPs | {temporal.gene_name.nunique()} genes")

# -----------------------------------------------------------------------
# STEP 2 -- attach eQTL p-value (same for a given pair across all days it
# appears in -- verified earlier: 0/74 pairs disagree)
# -----------------------------------------------------------------------

eqtl_pval_rows = []
for day in ["day0", "day4", "day14"]:
    d = pd.read_csv(rf"{abc_out}\{day}_eqtl_abc_cross_validation_summary.tsv", sep="\t")
    eqtl_pval_rows.append(d[["lead_snp", "gene_name", "min_eqtl_pval"]])
eqtl_pval = pd.concat(eqtl_pval_rows).drop_duplicates(subset=["lead_snp", "gene_name"])
eqtl_pval = eqtl_pval.rename(columns={"min_eqtl_pval": "eqtl_pval"})

report = temporal.merge(eqtl_pval, on=["lead_snp", "gene_name"], how="left")
print(f"STEP 2: eqtl_pval attached | missing: {report.eqtl_pval.isna().sum()}")

# -----------------------------------------------------------------------
# STEP 3 -- attach ENSG (gene -> ENSG is 1:1, verified earlier)
# -----------------------------------------------------------------------

eqtl_map = pd.read_csv(rf"{abc_in}\GWAS_to_AdipoExpress_gene_map_r2ge0.8.tsv", sep="\t")
gene_to_ensg = eqtl_map[["gene", "ENSG"]].drop_duplicates().rename(columns={"gene": "gene_name"})
report = report.merge(gene_to_ensg, on="gene_name", how="left")
print(f"STEP 3: ENSG attached | missing: {report.ENSG.isna().sum()}")

# -----------------------------------------------------------------------
# STEP 4 -- attach which trait(s) each lead SNP is genome-wide significant
# for, from the 34 binary trait columns in the trait x day matrix
# -----------------------------------------------------------------------

matrix = pd.read_excel(trait_matrix_xlsx)
non_trait_cols = {"lead_snp", "overlap_day0", "overlap_day4", "overlap_day14", "n_traits"}
trait_cols = [c for c in matrix.columns if c not in non_trait_cols]

matrix_subset = matrix[matrix.lead_snp.isin(report.lead_snp.unique())].set_index("lead_snp")
traits_per_lead = matrix_subset[trait_cols].apply(
    lambda row: ", ".join(row.index[row == 1]), axis=1
)
report["traits"] = report["lead_snp"].map(traits_per_lead)
print(f"STEP 4: traits attached | missing: {report.traits.isna().sum()}")

# -----------------------------------------------------------------------
# STEP 5 -- finalize column order, sort, write outputs
# -----------------------------------------------------------------------

report = report[[
    "lead_snp", "traits", "gene_name", "ENSG", "eqtl_pval",
    "day0_abc_score", "day4_abc_score", "day14_abc_score", "pattern",
]].sort_values(["lead_snp", "gene_name"]).reset_index(drop=True)

report.to_csv(rf"{out_dir}\qtl_sgbs_overlap_report.tsv", sep="\t", index=False)
report.to_excel(rf"{out_dir}\qtl_sgbs_overlap_report.xlsx", index=False)

print(f"\nSTEP 5: saved {report.shape[0]} rows to "
      f"qtl_sgbs_overlap_report.tsv / .xlsx in {out_dir}")
print("\nPreview:")
print(report.head(10).to_string())
