# ----------------------------- 8) GEO2R-like QC plots -------------------------
expr_vst <- assay(vsd)
expr_long <- as.data.frame(expr_vst) %>%
  tibble::rownames_to_column("gene_id") %>%
  tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "vst_expression") %>%
  dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition, batch), by = "sample_id")

# GEO2R-like boxplot of expression distributions.
p_box <- ggplot(expr_long, aes(x = sample_id, y = vst_expression, fill = condition)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.25) +
  manual_fill_or_null(condition_palette) +
  labs(
    title = "Expression distribution boxplot",
    subtitle = "Variance-stabilized expression; GEO2R-style sample comparability check",
    x = NULL,
    y = "VST expression"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_box, "05_geo2r_style_boxplot", width = 8, height = 4.8)

# GEO2R-like density plot.
p_density <- ggplot(expr_long, aes(x = vst_expression, color = sample_id, group = sample_id)) +
  geom_density(linewidth = 0.45, alpha = 0.85) +
  labs(
    title = "Expression density",
    subtitle = "GEO2R-style check for comparable expression distributions",
    x = "VST expression",
    y = "Density",
    color = "Sample"
  ) +
  theme_pub(base_size = 9)
save_gg(p_density, "06_geo2r_style_density", width = 8, height = 5)

# PCA for group separation and possible batch effects.
pca <- prcomp(t(expr_vst), scale. = FALSE)
pvar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition, batch, title), by = "sample_id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
  geom_point(size = 3.2, alpha = 0.95) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  manual_color_or_null(condition_palette) +
  labs(
    title = "PCA of samples",
    subtitle = "Primary unsupervised check for group structure and hidden batch effects",
    x = paste0("PC1 (", pvar[1], "%)"),
    y = paste0("PC2 (", pvar[2], "%)")
  ) +
  theme_pub()
save_gg(p_pca, "07_batch_qc_pca_by_condition", width = 6.5, height = 5.2)

if (has_real_batch) {
  p_pca_batch <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = batch, label = sample_id)) +
    geom_point(size = 3.2, alpha = 0.95) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
    manual_color_or_null(condition_palette) +
    labs(
      title = "PCA colored by condition and shaped by batch",
      subtitle = "Batch is included in the model only if non-confounded and full-rank",
      x = paste0("PC1 (", pvar[1], "%)"),
      y = paste0("PC2 (", pvar[2], "%)")
    ) +
    theme_pub()
  save_gg(p_pca_batch, "08_batch_qc_pca_by_detected_batch", width = 6.8, height = 5.2)
  
  expr_batch_removed <- limma::removeBatchEffect(expr_vst, batch = colData(dds)$batch)
  pca_adj <- prcomp(t(expr_batch_removed), scale. = FALSE)
  pvar_adj <- round(100 * pca_adj$sdev^2 / sum(pca_adj$sdev^2), 1)
  pca_adj_df <- as.data.frame(pca_adj$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition, batch), by = "sample_id")
  p_pca_adj <- ggplot(pca_adj_df, aes(x = PC1, y = PC2, color = condition, shape = batch, label = sample_id)) +
    geom_point(size = 3.2, alpha = 0.95) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
    manual_color_or_null(condition_palette) +
    labs(
      title = "PCA after batch removal for visualization",
      subtitle = "removeBatchEffect is used only for this QC plot, not for DESeq2 counts",
      x = paste0("PC1 (", pvar_adj[1], "%)"),
      y = paste0("PC2 (", pvar_adj[2], "%)")
    ) +
    theme_pub()
  save_gg(p_pca_adj, "09_batch_qc_pca_after_visual_batch_removal", width = 6.8, height = 5.2)
} else {
  writeLines("No explicit non-confounded batch-like metadata column was detected. Batch adjustment was not applied.",
             con = file.path(RESDIR, "batch_adjustment_note.txt"))
}

# UMAP, as in GEO2R visualizations.
n_neighbors <- max(2, min(5, nrow(coldata) - 1))
umap_mat <- uwot::umap(t(expr_vst), n_neighbors = n_neighbors, min_dist = 0.1, metric = "euclidean", verbose = FALSE)
umap_df <- as.data.frame(umap_mat)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- umap_df %>%
  dplyr::mutate(sample_id = rownames(t(expr_vst))) %>%
  dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition, batch), by = "sample_id")

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition, label = sample_id)) +
  geom_point(size = 3.3, alpha = 0.95) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  manual_color_or_null(condition_palette) +
  labs(
    title = "UMAP of samples",
    subtitle = paste0("GEO2R-style sample relationship plot; n_neighbors = ", n_neighbors),
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_pub()
save_gg(p_umap, "10_geo2r_style_umap", width = 6.5, height = 5.2)

# Sample distance heatmap.
sample_dist <- dist(t(expr_vst))
ann_col <- coldata %>% as.data.frame() %>% dplyr::select(condition, batch)
rownames(ann_col) <- coldata$sample_id
ann_colors <- list(condition = condition_palette)
batch_palette <- safe_named_palette(ann_col$batch)
if (!is.null(batch_palette) && length(batch_palette) > 0) ann_colors$batch <- batch_palette

save_pheatmap(
  as.matrix(sample_dist),
  "11_sample_distance_heatmap",
  width = 7,
  height = 6.5,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  main = "Sample-to-sample distances (VST)"
)

# DESeq2 dispersion estimates: RNA-seq equivalent of inspecting mean-variance behavior.
save_base_plot(
  "12_deseq2_dispersion_estimates",
  expr = DESeq2::plotDispEsts(dds, main = "DESeq2 dispersion estimates"),
  width = 6.5,
  height = 5
)

# Mean-variance trend plot, GEO2R-style quality check.
log_norm <- log2(norm_counts + 1)
mean_var_df <- tibble::tibble(
  mean_log_norm = rowMeans(log_norm),
  var_log_norm = matrixStats::rowVars(log_norm),
  mean_vst = rowMeans(expr_vst),
  sd_vst = matrixStats::rowSds(expr_vst)
)

p_meanvar <- ggplot(mean_var_df, aes(x = mean_log_norm, y = var_log_norm)) +
  geom_point(alpha = 0.18, size = 0.35) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  labs(
    title = "Mean-variance trend",
    subtitle = "Computed from log2 normalized counts; complements DESeq2 dispersion estimates",
    x = "Mean log2 normalized count",
    y = "Variance of log2 normalized count"
  ) +
  theme_pub()
save_gg(p_meanvar, "13_geo2r_style_mean_variance_trend", width = 6.5, height = 5)
