# ----------------------------- 11) Run manifest and session info -------------
manifest <- tibble::tribble(
  ~section, ~output, ~description,
  "metadata", "01_metadata_group_counts", "Sample count by inferred group.",
  "metadata", "02_metadata_overview_tiles", "Visual overview of sample metadata and detected batch field.",
  "QC", "03_qc_library_sizes", "Raw library size per sample.",
  "QC", "04_qc_detected_genes", "Detected genes per sample using count threshold.",
  "GEO2R-like", "05_geo2r_style_boxplot", "Expression distribution boxplot.",
  "GEO2R-like", "06_geo2r_style_density", "Expression density curves.",
  "batch", "07_batch_qc_pca_by_condition", "PCA colored by condition.",
  "GEO2R-like", "10_geo2r_style_umap", "UMAP sample relationship plot.",
  "QC", "11_sample_distance_heatmap", "Sample-to-sample distance heatmap.",
  "RNA-seq QC", "12_deseq2_dispersion_estimates", "DESeq2 dispersion estimates.",
  "GEO2R-like", "13_geo2r_style_mean_variance_trend", "Mean-variance trend.",
  "DE", "14_deg_counts_summary", "Counts of up/down/not significant genes.",
  "GEO2R-like", "15_geo2r_style_volcano", "Volcano plot.",
  "GEO2R-like", "16_geo2r_style_MD_MA_plot", "Mean-difference / MA plot.",
  "DE", "17_raw_pvalue_histogram", "Raw P-value histogram.",
  "GEO2R-like", "18_geo2r_style_adjusted_pvalue_histogram", "Adjusted P-value histogram.",
  "GEO2R-like", "19_geo2r_style_qq_plot_pvalues", "Q-Q plot of P-values; RNA-seq diagnostic analogue.",
  "publication", "20_top_de_genes_heatmap", "Top DE gene heatmap.",
  "GEO2R-like", "21_geo2r_style_gene_profile_top_genes", "Per-sample expression profiles for top genes.",
  "GEO2R-like", "22_geo2r_style_venn_diagram", "Generated only for >=2 contrasts. Not meaningful for default GSE160299 PD vs NC."
)
readr::write_csv(manifest, file.path(RESDIR, "output_manifest.csv"))

get_count_or_zero <- function(label) {
  idx <- match(label, as.character(summary_tbl$regulation))
  if (is.na(idx)) 0L else summary_tbl$n_genes[idx]
}

run_summary <- tibble::tibble(
  gse_id = GSE_ID,
  test_level = TEST_LEVEL,
  reference_level = REF_LEVEL,
  alpha = ALPHA,
  lfc_cutoff = LFC_CUTOFF,
  n_samples = ncol(count_mat),
  n_genes_before_filter = nrow(count_mat),
  n_genes_after_filter = nrow(dds),
  design_formula = paste(deparse(design_formula), collapse = ""),
  batch_column_detected = ifelse(is.na(batch_col), "none", batch_col),
  batch_adjusted_in_model = batch_adjusted,
  significant_up = get_count_or_zero("Up"),
  significant_down = get_count_or_zero("Down")
)
readr::write_csv(run_summary, file.path(RESDIR, "run_summary.csv"))

sink(file.path(RESDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

message_step("Analysis complete.")
message("Output directory: ", OUTDIR)
message("Main DEG table: ", file.path(RESDIR, paste0(GSE_ID, "_DESeq2_", TEST_LEVEL, "_vs_", REF_LEVEL, "_all_genes.csv")))
message("Figures directory: ", FIGDIR)
