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
