library(testthat)
library(miCCI)
library(data.table)

# Run all tests with verbose reporter so labels and content print
test_check("miCCI", reporter = testthat::ProgressReporter$new(show_praise = FALSE))
