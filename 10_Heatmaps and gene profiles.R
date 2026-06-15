# ----------------------------- 10) Heatmaps and gene profiles ----------------
ranked_genes <- res_tbl %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

# dplyr::n() cannot be used inside the n = argument of slice_head().
# Compute slice sizes explicitly to remain compatible with current dplyr versions.
n_heat <- min(TOP_HEATMAP_GENES, nrow(ranked_genes))
heat_genes <- ranked_genes %>%
  dplyr::slice_head(n = n_heat) %>%
  dplyr::pull(gene_id)

if (length(heat_genes) >= 2) {
  heat_mat <- expr_vst[heat_genes, , drop = FALSE]
  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[is.na(heat_mat_z)] <- 0
  row_labels <- res_tbl$gene_label[match(rownames(heat_mat_z), res_tbl$gene_id)]
  rownames(heat_mat_z) <- make.unique(row_labels)
  
  save_pheatmap(
    heat_mat_z,
    "20_top_de_genes_heatmap",
    width = 8,
    height = max(6, min(12, 0.14 * length(heat_genes) + 3)),
    annotation_col = ann_col,
    show_colnames = TRUE,
    show_rownames = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    fontsize_row = 6,
    main = paste0("Top ", length(heat_genes), " DE genes, z-scored VST")
  )
}

sig_ranked_genes <- ranked_genes %>%
  dplyr::filter(significant_fdr_lfc)

n_profile_sig <- min(TOP_PROFILE_GENES, nrow(sig_ranked_genes))
profile_genes <- sig_ranked_genes %>%
  dplyr::slice_head(n = n_profile_sig) %>%
  dplyr::pull(gene_id)

if (length(profile_genes) == 0) {
  n_profile_all <- min(TOP_PROFILE_GENES, nrow(ranked_genes))
  profile_genes <- ranked_genes %>%
    dplyr::slice_head(n = n_profile_all) %>%
    dplyr::pull(gene_id)
}

if (length(profile_genes) >= 1) {
  profile_df <- as.data.frame(log2(norm_counts[profile_genes, , drop = FALSE] + 1)) %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "log2_norm_count") %>%
    dplyr::left_join(res_tbl %>% dplyr::select(gene_id, gene_label, padj, log2FoldChange), by = "gene_id") %>%
    dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition), by = "sample_id")
  
  p_profile <- ggplot(profile_df, aes(x = sample_id, y = log2_norm_count, fill = condition)) +
    geom_col(color = "black", linewidth = 0.15, width = 0.72) +
    facet_wrap(~ gene_label, scales = "free_y", ncol = 4) +
    manual_fill_or_null(condition_palette) +
    labs(
      title = "Gene expression profile plots",
      subtitle = "GEO2R-style per-sample profiles for top-ranked genes",
      x = NULL,
      y = "log2(DESeq2 normalized count + 1)"
    ) +
    theme_pub(base_size = 8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
  save_gg(p_profile, "21_geo2r_style_gene_profile_top_genes", width = 10, height = 7)
}
  # Optional Venn support. GSE160299 has one contrast, so a Venn diagram is not meaningful.
  deg_sets <- list()
  deg_sets[[paste0(TEST_LEVEL, "_vs_", REF_LEVEL)]] <- res_tbl %>%
    dplyr::filter(significant_fdr_lfc) %>%
    dplyr::pull(gene_id)
  
  if (length(deg_sets) >= 2 && has_ggvenn) {
    p_venn <- ggvenn::ggvenn(deg_sets, stroke_size = 0.4, set_name_size = 3.5) +
      labs(title = "Venn diagram of significant DE genes") +
      theme_pub()
    save_gg(p_venn, "22_geo2r_style_venn_diagram", width = 6, height = 5)
  } else {
    writeLines(
      c(
        "Venn diagram not generated.",
        "Reason: GSE160299 default analysis has only one contrast (PD vs NC).",
        "If you adapt this script to >1 contrast, store significant gene vectors in deg_sets and rerun the ggvenn block."
      ),
      con = file.path(RESDIR, "venn_diagram_note.txt")
    )
  }