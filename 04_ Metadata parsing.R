# ----------------------------- 4) Metadata parsing ---------------------------
infer_condition_from_text <- function(x) {
  y <- stringr::str_to_upper(ifelse(is.na(x), "", x))
  dplyr::case_when(
    stringr::str_detect(y, "(^|[^A-Z])NC[ _.-]*[0-9]*([^A-Z]|$)|CONTROL|HEALTHY|NON[ -]?PD") ~ REF_LEVEL,
    stringr::str_detect(y, "(^|[^A-Z])PD[ _.-]*[0-9]*([^A-Z]|$)|PARKINSON") ~ TEST_LEVEL,
    TRUE ~ NA_character_
  )
}

infer_short_id <- function(x) {
  y <- stringr::str_to_upper(ifelse(is.na(x), "", x))
  y <- stringr::str_replace_all(y, "[^A-Z0-9]", "")
  stringr::str_extract(y, "(NC|PD)[0-9]+")
}

metadata <- metadata_raw %>%
  dplyr::mutate(
    title = if ("title" %in% names(.)) title else geo_accession,
    source_name_ch1 = if ("source_name_ch1" %in% names(.)) source_name_ch1 else NA_character_,
    condition_from_title = infer_condition_from_text(title),
    condition_from_source = infer_condition_from_text(source_name_ch1),
    condition = dplyr::coalesce(condition_from_title, condition_from_source),
    meta_short_id = infer_short_id(title)
  )

# Detect likely batch/lane/run metadata, but use it only if statistically usable.
batch_candidates <- names(metadata)[stringr::str_detect(
  names(metadata),
  stringr::regex("batch|lane|flowcell|flow_cell|run|library|lib|center|site|date|instrument|chip", ignore_case = TRUE)
)]

batch_col <- NA_character_
if (length(batch_candidates) > 0) {
  for (bc in batch_candidates) {
    vals <- metadata[[bc]]
    vals <- ifelse(is.na(vals) | vals == "", NA, as.character(vals))
    nvals <- length(unique(stats::na.omit(vals)))
    # A useful batch variable should not be all identical or all unique.
    if (nvals > 1 && nvals < nrow(metadata)) {
      batch_col <- bc
      break
    }
  }
}

metadata <- metadata %>%
  dplyr::mutate(batch_detected = if (!is.na(batch_col)) as.character(.data[[batch_col]]) else "not_available")

readr::write_csv(metadata, file.path(RESDIR, paste0(GSE_ID, "_metadata_raw_parsed.csv")))
