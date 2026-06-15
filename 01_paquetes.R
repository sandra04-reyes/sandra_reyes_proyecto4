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
