#!/usr/bin/env Rscript
# =============================================================================
# Bulk RNA-seq differential expression from GEO: GSE160299
# GEO2R-style + publication-grade DESeq2 workflow
# Revision: namespace-safe + scales-compatible + palette-safe version; avoids dplyr masking, removed defunct label_number_si(), and prevents empty fill/color scales.
#
# Dataset default:
#   GSE160299: human plasma/serum RNA-seq, NC vs Parkinson's disease (PD)
#   Samples expected: Plasma NC1-NC4 and Plasma PD1-PD4
#
# Main outputs:
#   results/   DESeq2 tables, normalized counts, sample metadata, manifest
#   figures/   PDF + PNG plots for QC, metadata, batch effects and DE results
#
# Notes:
#   - Uses raw gene count matrix provided as GEO supplementary file.
#   - Adjusts for batch only if a usable, non-confounded batch field is detected.
#   - For GSE160299 no explicit batch column is expected, so batch assessment is QC-only.
#   - GEO2R-like plots included: volcano, MD/MA, UMAP, boxplot, density,
#     adjusted P-value histogram, q-q plot analogue, mean-variance trend,
#     gene profile plot, and Venn diagram support for multi-contrast analyses.
# =============================================================================

# ----------------------------- 0) Configuration ------------------------------
GSE_ID      <- "GSE160299"
TEST_LEVEL  <- "PD"
REF_LEVEL   <- "NC"
ALPHA       <- 0.05
LFC_CUTOFF  <- 1.0
MIN_COUNT   <- 10
MIN_SAMPLES <- 2
ADJUST_FOR_BATCH_IF_POSSIBLE <- TRUE
TOP_HEATMAP_GENES <- 50
TOP_PROFILE_GENES <- 12
SEED <- 160299

OUTDIR <- file.path(getwd(), paste0(GSE_ID, "_DESeq2_analysis"))
RAWDIR <- file.path(OUTDIR, "data_raw")
RESDIR <- file.path(OUTDIR, "results")
FIGDIR <- file.path(OUTDIR, "figures")
DIRS <- c(OUTDIR, RAWDIR, RESDIR, FIGDIR)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

set.seed(SEED)
options(stringsAsFactors = FALSE)

# ----------------------------- 1) Packages -----------------------------------
INSTALL_MISSING <- TRUE

cran_pkgs <- c(
  "tidyverse", "ggrepel", "pheatmap", "uwot", "matrixStats",
  "scales", "janitor", "ggvenn"
)
bioc_pkgs <- c(
  "GEOquery", "DESeq2", "limma", "AnnotationDbi", "org.Hs.eg.db",
  "apeglm"
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!INSTALL_MISSING) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         ". Set INSTALL_MISSING <- TRUE or install them manually.")
  }
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
  library(uwot)
  library(matrixStats)
  library(scales)
  library(janitor)
  library(GEOquery)
  library(DESeq2)
  library(limma)
  # AnnotationDbi is used with explicit AnnotationDbi:: calls to avoid dplyr::select() masking.
  library(org.Hs.eg.db)
})

has_apeglm <- requireNamespace("apeglm", quietly = TRUE)
has_ggvenn <- requireNamespace("ggvenn", quietly = TRUE)

# Use explicit namespaces for data manipulation because Bioconductor packages,
# plyr, Hmisc or other attached packages can mask dplyr::count() and dplyr::select().

# ----------------------------- 2) Plot helpers -------------------------------
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey88"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey25"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

save_gg <- function(p, name, width = 7, height = 5, dpi = 320) {
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf_file, png_file))
}


# scales::label_number_si() was deprecated in scales 1.2.0 and is defunct in
# recent versions. This wrapper uses the modern API when available and falls
# back to comma labels on older/incomplete scales installations.
label_number_compact_safe <- function(...) {
  if (exists("cut_si", where = asNamespace("scales"), inherits = FALSE)) {
    scales::label_number(scale_cut = scales::cut_si(""), ...)
  } else {
    scales::label_comma(...)
  }
}

# Robust discrete palette helper. The grid error
# "the 'gpar' element 'fill' must not have length 0" is commonly triggered
# when a manual fill/color scale receives an empty named vector. These helpers
# always return either a non-empty named palette for the values actually present
# in the plotted data, or NULL. ggplot can safely add NULL layers.
safe_named_palette <- function(x, preferred = NULL, fallback_palette = scales::hue_pal()) {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  vals <- sort(vals)
  if (length(vals) == 0) return(NULL)

  if (!is.null(preferred) && length(preferred) > 0 && all(vals %in% names(preferred))) {
    pal <- preferred[vals]
  } else {
    pal <- setNames(fallback_palette(length(vals)), vals)
  }

  if (length(pal) == 0 || any(is.na(pal)) || any(!nzchar(as.character(pal)))) {
    pal <- setNames(fallback_palette(length(vals)), vals)
  }
  pal
}

manual_fill_or_null <- function(pal, drop = FALSE) {
  if (is.null(pal) || length(pal) == 0) return(NULL)
  ggplot2::scale_fill_manual(values = pal, drop = drop)
}

manual_color_or_null <- function(pal, drop = FALSE) {
  if (is.null(pal) || length(pal) == 0) return(NULL)
  ggplot2::scale_color_manual(values = pal, drop = drop)
}

save_base_plot <- function(name, expr, width = 7, height = 5, dpi = 320) {
  expr_sub <- substitute(expr)
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  eval(expr_sub, envir = parent.frame())
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = dpi)
  eval(expr_sub, envir = parent.frame())
  dev.off()
  invisible(c(pdf_file, png_file))
}

save_pheatmap <- function(mat, name, width = 7, height = 7, ...) {
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  print(pheatmap::pheatmap(mat, silent = TRUE, ...))
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 320)
  print(pheatmap::pheatmap(mat, silent = TRUE, ...))
  dev.off()
  invisible(c(pdf_file, png_file))
}

message_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

# ----------------------------- 3) GEO download -------------------------------
message_step("Downloading GEO metadata and supplementary files for ", GSE_ID, "...")

gse_list <- GEOquery::getGEO(GSE_ID, GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gse_list) == 0) stop("No Series Matrix object retrieved from GEO.")
eset <- gse_list[[1]]
metadata_raw <- Biobase::pData(eset) %>% as.data.frame()

if (!"geo_accession" %in% names(metadata_raw)) {
  metadata_raw <- metadata_raw %>% tibble::rownames_to_column("geo_accession")
}
metadata_raw <- metadata_raw %>% janitor::clean_names()

# Download supplementary files. GEOquery may download both the raw tar and the count matrix.
supp <- GEOquery::getGEOSuppFiles(GSE_ID, makeDirectory = FALSE, baseDir = RAWDIR)

# Extract tar files if present.
tar_files <- list.files(RAWDIR, pattern = "\\.tar$", full.names = TRUE, ignore.case = TRUE)
if (length(tar_files) > 0) {
  for (tf in tar_files) {
    exdir <- file.path(RAWDIR, tools::file_path_sans_ext(basename(tf)))
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    try(utils::untar(tf, exdir = exdir), silent = TRUE)
  }
}

count_files <- list.files(
  RAWDIR,
  pattern = "(raw.*gene.*count|gene.*count|count.*matrix|counts).*\\.txt(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
if (length(count_files) == 0) {
  stop("No count matrix found in supplementary files. Inspect ", RAWDIR,
       " and update the count_files pattern manually.")
}

# Prefer the GEO-provided raw gene counts matrix if present.
count_file <- count_files[which.max(grepl("Raw_gene_counts_matrix", basename(count_files), ignore.case = TRUE))]
message_step("Using count matrix: ", count_file)

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

# ----------------------------- 7) DESeq2 model -------------------------------
message_step("Constructing DESeq2 object...")

dds0 <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ condition
)

keep_genes <- rowSums(DESeq2::counts(dds0) >= MIN_COUNT) >= MIN_SAMPLES
dds <- dds0[keep_genes, ]
message_step("Genes retained after expression filter: ", nrow(dds))

# Use batch in the design only if it exists and the model matrix is full-rank.
has_real_batch <- nlevels(colData(dds)$batch) > 1 && !all(colData(dds)$batch == "not_available")
design_formula <- ~ condition
batch_adjusted <- FALSE

if (ADJUST_FOR_BATCH_IF_POSSIBLE && has_real_batch) {
  mm_try <- try(model.matrix(~ batch + condition, data = as.data.frame(colData(dds))), silent = TRUE)
  if (!inherits(mm_try, "try-error") && qr(mm_try)$rank == ncol(mm_try)) {
    design_formula <- ~ batch + condition
    batch_adjusted <- TRUE
  } else {
    warning("A batch-like column was detected but is confounded or rank-deficient. Using ~ condition only.")
  }
}

design(dds) <- design_formula
message_step("DESeq2 design formula: ", paste(deparse(design_formula), collapse = ""))

dds <- DESeq2::DESeq(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)

norm_counts <- DESeq2::counts(dds, normalized = TRUE)
readr::write_csv(as.data.frame(norm_counts) %>% tibble::rownames_to_column("gene_id"),
                 file.path(RESDIR, paste0(GSE_ID, "_DESeq2_normalized_counts.csv")))

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
