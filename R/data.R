#' Dummy patient dataset for testing all miCCI strategies
#'
#' @description
#' A data frame with 39 synthetic inpatient records covering the full spectrum
#' of Charlson Comorbidity Index (CCI) values, including healthy patients
#' (CCI = 0), single-comorbidity cases (CCI = 1-2), multi-morbid patients
#' (CCI = 3-7), and high-burden cases (CCI >= 7). The dataset also contains
#' four rows specifically designed to test the \emph{depends_on} suppression
#' logic (e.g., dm_complicated suppresses dm_simple).
#'
#' All ICD codes are provided at the full 4-digit level in
#' \code{list_full_icd_code}. The \code{run_demo()} function automatically
#' truncates these to 3 characters, simulating the §21 KHEntgG anonymization
#' process, before running the six strategies.
#'
#' @format A data frame with 39 rows and 3 columns.
#' \describe{
#'   \item{p_id}{Character. Unique patient identifier (PAT_001 to PAT_039).}
#'   \item{list_full_icd_code}{Character. Pipe-separated full 4-digit
#'     ICD-10-GM codes (e.g., \code{"E11.40|N18.40|I21.09"}). These are the
#'     codes available BEFORE anonymization.}
#'   \item{los_in_days}{Integer. Length of hospital stay in days.}
#' }
#'
#' @examples
#' data(micci_patients)
#' head(micci_patients)
#' table(nchar(gsub("[^|]", "", micci_patients$list_full_icd_code)) + 1)
#'
#' @source
#' Synthetic data created for package testing. Not based on real patient records.
"micci_patients"
