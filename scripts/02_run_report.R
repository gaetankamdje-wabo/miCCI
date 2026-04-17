# =============================================================================
# miCCI v1.0.0 — Stage 2 runner: BUILD REPORT (fast, ~minutes)
# -----------------------------------------------------------------------------
# Reads <INPUT_DIR>/predictions.parquet (written by 01_run_compute.R) and
# regenerates ALL tables and figures. Run as often as needed for cosmetic
# tweaks; the long compute stage is never re-executed.
# =============================================================================

INPUT_DIR  <- "A:/HLZ/Promotionen/Gaetan Kamdje Wabo/Results CCI Model v050"
OUTPUT_DIR <- INPUT_DIR    # write tables/figures next to predictions.parquet

# Bootstrap configuration
B_GLOBAL  <- 1000L         # global per-strategy CIs
B_STRATA  <- 500L          # per-stratum CIs (chapter / group / pattern)
TOP_PATS  <- 30L           # number of most-frequent multimorbid patterns
PARALLEL  <- FALSE         # set TRUE if future.apply is installed

if (requireNamespace("miCCI", quietly = TRUE)) {
  library(miCCI)
} else if (file.exists("DESCRIPTION")) {
  devtools::load_all(".")
} else {
  stop("miCCI is not installed and DESCRIPTION not found in the working ",
       "directory. Either run `devtools::install('.')` from the miCCI ",
       "folder, or setwd() into it and try again.")
}
library(data.table)

cat(strrep("=", 64), "\n")
cat("miCCI v1.0.0 — Stage 2: build_report\n")
cat("Input:  ", INPUT_DIR, "\n")
cat("Output: ", OUTPUT_DIR, "\n")
cat("Bootstrap: global=", B_GLOBAL, "  strata=", B_STRATA,
    "  parallel=", PARALLEL, "\n", sep = "")
cat(strrep("=", 64), "\n\n")

res <- build_report(
  input_dir    = INPUT_DIR,
  output_dir   = OUTPUT_DIR,
  B            = B_GLOBAL,
  B_strata     = B_STRATA,
  top_patterns = TOP_PATS,
  parallel     = PARALLEL
)

cat("\n=== Files in output directory ===\n")
print(list.files(OUTPUT_DIR, pattern = "^(tab|fig)_", full.names = FALSE))
cat("\nDone.\n")
