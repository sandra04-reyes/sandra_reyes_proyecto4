# ----------------------------- 5) Count matrix parsing -----------------------
message_step("Reading and parsing count matrix...")
counts_raw <- readr::read_tsv(count_file, show_col_types = FALSE, progress = FALSE)
counts_raw <- as.data.frame(counts_raw, check.names = FALSE)

is_numeric_like <- function(x) {
  suppressWarnings(y <- as.numeric(as.character(x)))
  mean(!is.na(y)) > 0.95
}

numeric_cols <- names(counts_raw)[vapply(counts_raw, is_numeric_like, logical(1))]
if (length(numeric_cols) < 2) stop("Could not identify numeric count columns.")

# Prefer numeric columns whose names look like sample identifiers; otherwise take the last N numeric columns.
name_hits <- numeric_cols[stringr::str_detect(
  stringr::str_to_upper(numeric_cols), "GSM|(^|[^A-Z])(NC|PD)[ _.-]*[0-9]+|PLASMA"
)]

expected_n <- nrow(metadata)
if (length(name_hits) >= 2) {
  sample_cols <- name_hits
} else if (length(numeric_cols) >= expected_n) {
  sample_cols <- tail(numeric_cols, expected_n)
} else {
  sample_cols <- numeric_cols
}

annotation_cols <- setdiff(names(counts_raw), sample_cols)
if (length(annotation_cols) == 0) {
  counts_raw$gene_id_auto <- paste0("gene_", seq_len(nrow(counts_raw)))
  annotation_cols <- "gene_id_auto"
}

gene_id_col <- annotation_cols[1]
gene_ids <- as.character(counts_raw[[gene_id_col]])
gene_ids[is.na(gene_ids) | gene_ids == ""] <- paste0("gene_", which(is.na(gene_ids) | gene_ids == ""))

count_mat <- counts_raw[, sample_cols, drop = FALSE] %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ suppressWarnings(as.numeric(as.character(.x))))) %>%
  as.matrix()
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)

# Aggregate duplicated gene IDs by summing counts.
count_mat <- rowsum(count_mat, group = gene_ids, reorder = FALSE)
mode(count_mat) <- "integer"

# Keep the first annotation row per gene ID.
gene_annot <- counts_raw[, annotation_cols, drop = FALSE] %>%
  dplyr::mutate(gene_id = gene_ids) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)

# Build colData aligned to count matrix columns.
make_coldata_row <- function(count_col) {
  short <- infer_short_id(count_col)
  hit <- which(metadata$geo_accession == count_col)
  if (length(hit) == 0 && !is.na(short)) hit <- which(metadata$meta_short_id == short)
  
  if (length(hit) > 0) {
    row <- metadata[hit[1], , drop = FALSE]
    tibble::tibble(
      sample_id = count_col,
      geo_accession = row$geo_accession,
      title = row$title,
      source_name_ch1 = row$source_name_ch1,
      condition = dplyr::coalesce(row$condition, infer_condition_from_text(count_col)),
      inferred_short_id = dplyr::coalesce(row$meta_short_id, short),
      batch = row$batch_detected
    )
  } else {
    tibble::tibble(
      sample_id = count_col,
      geo_accession = NA_character_,
      title = count_col,
      source_name_ch1 = NA_character_,
      condition = infer_condition_from_text(count_col),
      inferred_short_id = short,
      batch = "not_available"
    )
  }
}

coldata <- purrr::map_dfr(colnames(count_mat), make_coldata_row)

if (any(is.na(coldata$condition))) {
  warning("Some samples could not be assigned to condition. They will be removed: ",
          paste(coldata$sample_id[is.na(coldata$condition)], collapse = ", "))
}
keep_samples <- !is.na(coldata$condition)
count_mat <- count_mat[, coldata$sample_id[keep_samples], drop = FALSE]
coldata <- coldata[keep_samples, , drop = FALSE]

coldata <- coldata %>%
  dplyr::mutate(
    condition = factor(condition, levels = c(REF_LEVEL, TEST_LEVEL)),
    batch = factor(ifelse(is.na(batch) | batch == "", "not_available", as.character(batch)))
  ) %>%
  droplevels()
rownames(coldata) <- coldata$sample_id

if (nrow(coldata) == 0 || ncol(count_mat) == 0) {
  stop("No samples remained after condition parsing. Check sample labels and the REF_LEVEL/TEST_LEVEL settings.")
}
if (nlevels(coldata$condition) < 2) {
  stop("Fewer than two condition levels were detected after parsing: ",
       paste(levels(coldata$condition), collapse = ", "),
       ". Differential expression requires at least two groups.")
}

# Drop genes with all-zero counts and save parsed inputs.
count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]
readr::write_csv(coldata %>% tibble::rownames_to_column("count_matrix_col"),
                 file.path(RESDIR, paste0(GSE_ID, "_sample_metadata_final.csv")))
readr::write_csv(as.data.frame(count_mat) %>% tibble::rownames_to_column("gene_id"),
                 file.path(RESDIR, paste0(GSE_ID, "_raw_counts_parsed.csv")))
readr::write_csv(gene_annot, file.path(RESDIR, paste0(GSE_ID, "_gene_annotation_from_matrix.csv")))

message_step("Samples retained: ", ncol(count_mat), "; genes retained before filtering: ", nrow(count_mat))
print(table(coldata$condition))
