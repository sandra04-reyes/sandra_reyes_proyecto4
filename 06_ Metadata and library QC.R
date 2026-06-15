# ----------------------------- 6) Metadata and library QC --------------------
preferred_condition_palette <- c("NC" = "#4C78A8", "PD" = "#F58518")
condition_palette <- safe_named_palette(coldata$condition, preferred = preferred_condition_palette)
if (is.null(condition_palette) || length(condition_palette) == 0) {
  stop("Could not build a non-empty condition palette. Check coldata$condition.")
}

p_meta_counts <- coldata %>%
  dplyr::count(condition) %>%
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.35, size = 3.5) +
  manual_fill_or_null(condition_palette) +
  labs(
    title = paste0(GSE_ID, " sample composition"),
    subtitle = "Groups inferred from GEO sample titles / count-matrix labels",
    x = NULL,
    y = "Number of samples"
  ) +
  theme_pub() +
  theme(legend.position = "none")
save_gg(p_meta_counts, "01_metadata_group_counts", width = 5, height = 4)

meta_tile <- coldata %>%
  dplyr::select(sample_id, geo_accession, title, condition, batch) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
  tidyr::pivot_longer(-sample_id, names_to = "field", values_to = "value") %>%
  dplyr::mutate(value = ifelse(is.na(value) | value == "", "missing", value))

p_meta_tile <- ggplot(meta_tile, aes(x = field, y = sample_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.25) +
  labs(
    title = "Metadata overview",
    subtitle = "Useful for checking group assignment and available batch fields",
    x = NULL,
    y = NULL,
    fill = "Value"
  ) +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
save_gg(p_meta_tile, "02_metadata_overview_tiles", width = 8, height = 4.8)

lib_df <- tibble::tibble(
  sample_id = colnames(count_mat),
  library_size = colSums(count_mat),
  detected_genes_count_ge_10 = colSums(count_mat >= MIN_COUNT)
) %>%
  dplyr::left_join(coldata %>% tibble::rownames_to_column("sample_id2") %>% dplyr::select(sample_id = sample_id2, condition, batch), by = "sample_id") %>%
  dplyr::arrange(condition, sample_id) %>%
  dplyr::mutate(sample_id = factor(sample_id, levels = sample_id))

p_lib <- ggplot(lib_df, aes(x = sample_id, y = library_size, fill = condition)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.75) +
  scale_y_continuous(labels = label_number_compact_safe()) +
  manual_fill_or_null(condition_palette) +
  labs(
    title = "Library sizes",
    subtitle = "Total raw counts per sample",
    x = NULL,
    y = "Total counts"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_lib, "03_qc_library_sizes", width = 7, height = 4.5)

p_detected <- ggplot(lib_df, aes(x = sample_id, y = detected_genes_count_ge_10, fill = condition)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.75) +
  scale_y_continuous(labels = scales::comma) +
  manual_fill_or_null(condition_palette) +
  labs(
    title = "Detected genes per sample",
    subtitle = paste0("Genes with raw count >= ", MIN_COUNT),
    x = NULL,
    y = "Detected genes"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_gg(p_detected, "04_qc_detected_genes", width = 7, height = 4.5)
