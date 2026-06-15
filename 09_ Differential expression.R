# ----------------------------- 9) Differential expression --------------------
message_step("Extracting differential expression results...")

coef_name <- paste0("condition_", TEST_LEVEL, "_vs_", REF_LEVEL)
available_results <- DESeq2::resultsNames(dds)
if (!coef_name %in% available_results) {
  message("Available DESeq2 result names: ", paste(available_results, collapse = ", "))
  stop("Expected coefficient not found: ", coef_name,
       ". Check TEST_LEVEL/REF_LEVEL and condition labels.")
}

res <- DESeq2::results(dds, contrast = c("condition", TEST_LEVEL, REF_LEVEL), alpha = ALPHA)

res_shrunk <- tryCatch({
  if (has_apeglm) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "normal")
  }
}, error = function(e) {
  warning("lfcShrink failed; using unshrunk DESeq2 results. Reason: ", conditionMessage(e))
  res
})

# Gene symbol annotation.
annotate_gene_symbols <- function(gene_ids, gene_annot) {
  gene_no_version <- stringr::str_remove(gene_ids, "\\.\\d+$")
  symbol <- rep(NA_character_, length(gene_ids))
  
  # Use symbol-like columns from the GEO count matrix if available.
  possible_symbol_cols <- names(gene_annot)[stringr::str_detect(
    names(gene_annot),
    stringr::regex("symbol|gene.?name|external.?gene.?name|gene.?symbol", ignore_case = TRUE)
  )]
  possible_symbol_cols <- setdiff(possible_symbol_cols, c("gene_id"))
  
  if (length(possible_symbol_cols) > 0) {
    sym_df <- gene_annot %>%
      dplyr::transmute(gene_id, symbol_from_matrix = as.character(.data[[possible_symbol_cols[1]]]))
    symbol <- sym_df$symbol_from_matrix[match(gene_ids, sym_df$gene_id)]
    symbol <- ifelse(is.na(symbol) | symbol == "", NA_character_, as.character(symbol))
  }
  
  # Fallback: org.Hs.eg.db mapping.
  if (mean(stringr::str_detect(gene_no_version, "^ENSG"), na.rm = TRUE) > 0.3) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = gene_no_version,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    symbol <- dplyr::coalesce(symbol, unname(mapped[gene_no_version]))
  } else if (mean(stringr::str_detect(gene_no_version, "^[0-9]+$"), na.rm = TRUE) > 0.3) {
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = gene_no_version,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )
    symbol <- dplyr::coalesce(symbol, unname(mapped[gene_no_version]))
  }
  
  dplyr::coalesce(symbol, gene_ids)
}

res_tbl <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(
    symbol = annotate_gene_symbols(gene_id, gene_annot),
    gene_label = ifelse(is.na(symbol) | symbol == "", gene_id, symbol),
    regulation = dplyr::case_when(
      !is.na(padj) & padj < ALPHA & log2FoldChange >=  LFC_CUTOFF ~ "Up",
      !is.na(padj) & padj < ALPHA & log2FoldChange <= -LFC_CUTOFF ~ "Down",
      TRUE ~ "Not significant"
    ),
    significant_fdr = !is.na(padj) & padj < ALPHA,
    significant_fdr_lfc = !is.na(padj) & padj < ALPHA & abs(log2FoldChange) >= LFC_CUTOFF,
    neg_log10_pvalue = -log10(pmax(pvalue, .Machine$double.xmin)),
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin))
  ) %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

readr::write_csv(res_tbl, file.path(RESDIR, paste0(GSE_ID, "_DESeq2_", TEST_LEVEL, "_vs_", REF_LEVEL, "_all_genes.csv")))
readr::write_csv(res_tbl %>% dplyr::filter(significant_fdr_lfc),
                 file.path(RESDIR, paste0(GSE_ID, "_DESeq2_", TEST_LEVEL, "_vs_", REF_LEVEL, "_significant_FDR", ALPHA, "_absLFC", LFC_CUTOFF, ".csv")))

summary_tbl <- res_tbl %>%
  dplyr::count(regulation, name = "n_genes") %>%
  dplyr::mutate(regulation = factor(regulation, levels = c("Up", "Down", "Not significant"))) %>%
  dplyr::arrange(regulation)
readr::write_csv(summary_tbl, file.path(RESDIR, paste0(GSE_ID, "_DEG_summary.csv")))

p_deg_counts <- ggplot(summary_tbl, aes(x = regulation, y = n_genes, fill = regulation)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.7) +
  geom_text(aes(label = scales::comma(n_genes)), vjust = -0.35, size = 3.4) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not significant" = "grey75")) +
  labs(
    title = "Differential expression summary",
    subtitle = paste0(TEST_LEVEL, " vs ", REF_LEVEL, "; FDR < ", ALPHA, ", |log2FC| >= ", LFC_CUTOFF),
    x = NULL,
    y = "Genes"
  ) +
  theme_pub() +
  theme(legend.position = "none")
save_gg(p_deg_counts, "14_deg_counts_summary", width = 5.5, height = 4.5)

# GEO2R-like volcano plot.
top_labels <- res_tbl %>%
  dplyr::filter(significant_fdr_lfc) %>%
  dplyr::slice_min(order_by = padj, n = 15, with_ties = FALSE)

p_volcano <- ggplot(res_tbl, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_point(aes(color = regulation), alpha = 0.65, size = 0.85, na.rm = TRUE) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = gene_label),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not significant" = "grey72")) +
  labs(
    title = "Volcano plot",
    subtitle = paste0(TEST_LEVEL, " vs ", REF_LEVEL, "; colors use FDR < ", ALPHA, " and |log2FC| >= ", LFC_CUTOFF),
    x = paste0("log2 fold change (", TEST_LEVEL, " / ", REF_LEVEL, ")"),
    y = "-log10(raw P value)",
    color = "Class"
  ) +
  theme_pub()
save_gg(p_volcano, "15_geo2r_style_volcano", width = 7, height = 5.5)

# GEO2R-like mean-difference / MA plot.
p_ma <- ggplot(res_tbl, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = regulation), alpha = 0.65, size = 0.85, na.rm = TRUE) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = gene_label),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(values = c("Up" = "#D62728", "Down" = "#1F77B4", "Not significant" = "grey72")) +
  labs(
    title = "Mean-difference / MA plot",
    subtitle = "GEO2R-style log2FC versus average expression",
    x = "log2(baseMean + 1)",
    y = paste0("log2 fold change (", TEST_LEVEL, " / ", REF_LEVEL, ")"),
    color = "Class"
  ) +
  theme_pub()
save_gg(p_ma, "16_geo2r_style_MD_MA_plot", width = 7, height = 5.5)

# P-value and adjusted P-value histograms.
p_p_hist <- res_tbl %>%
  dplyr::filter(!is.na(pvalue)) %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(bins = 40, color = "black", linewidth = 0.2) +
  labs(
    title = "Raw P-value histogram",
    subtitle = "Useful for diagnosing global signal and model behavior",
    x = "Raw P value",
    y = "Genes"
  ) +
  theme_pub()
save_gg(p_p_hist, "17_raw_pvalue_histogram", width = 6.5, height = 4.5)

p_adj_hist <- res_tbl %>%
  dplyr::filter(!is.na(padj)) %>%
  ggplot(aes(x = padj)) +
  geom_histogram(bins = 40, color = "black", linewidth = 0.2) +
  labs(
    title = "Adjusted P-value histogram",
    subtitle = "GEO2R-style summary of multiple-testing-adjusted results",
    x = "Adjusted P value / FDR",
    y = "Genes"
  ) +
  theme_pub()
save_gg(p_adj_hist, "18_geo2r_style_adjusted_pvalue_histogram", width = 6.5, height = 4.5)

# Q-Q plot analogue for RNA-seq P values.
pvals <- res_tbl$pvalue[!is.na(res_tbl$pvalue) & res_tbl$pvalue > 0 & res_tbl$pvalue <= 1]
pvals <- sort(pvals)
qq_df <- tibble::tibble(
  expected = -log10(ppoints(length(pvals))),
  observed = -log10(pvals)
) %>% dplyr::arrange(expected)
lim_qq <- max(c(qq_df$expected, qq_df$observed), na.rm = TRUE)

p_qq <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.35) +
  coord_equal(xlim = c(0, lim_qq), ylim = c(0, lim_qq)) +
  labs(
    title = "Q-Q plot of differential-expression P values",
    subtitle = "RNA-seq analogue of GEO2R q-q diagnostics",
    x = "Expected -log10(P)",
    y = "Observed -log10(P)"
  ) +
  theme_pub()
save_gg(p_qq, "19_geo2r_style_qq_plot_pvalues", width = 5.5, height = 5.5)
