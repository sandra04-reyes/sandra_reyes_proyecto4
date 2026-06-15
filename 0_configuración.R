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
