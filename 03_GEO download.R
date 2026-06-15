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
