# =============================================================================
# miCCI — Main execution script
# =============================================================================

DATA_PATH  <- "A:/HLZ/Projekt/2025/comorbidity_gaetan/comorbidity_diagnoses_only_pseudonym.parquet"
OUTPUT_DIR <- "A:/HLZ/Promotionen/Gaetan Kamdje Wabo/Results CCI Model"
SAMPLE_SIZE <- NULL  # Set NULL for full run (~720k encounters)

if (file.exists("DESCRIPTION")) devtools::load_all(".") else library(miCCI)
library(data.table)

cat(strrep("=", 60), "\n")
cat("miCCI Validation Pipeline\n")
cat("Data:   ", DATA_PATH, "\n")
cat("Output: ", OUTPUT_DIR, "\n")
cat("Sample: ", ifelse(is.null(SAMPLE_SIZE), "FULL", as.character(SAMPLE_SIZE)), "\n")
cat(strrep("=", 60), "\n\n")

results <- run_pipeline(data_path = DATA_PATH, output_dir = OUTPUT_DIR,
                        sample_size = SAMPLE_SIZE)
cat("\n=== DONE ===\n")
