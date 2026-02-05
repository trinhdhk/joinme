# Global variable bindings for dplyr NSE and other functions
# These are used in dplyr pipelines and prevent R CMD check warnings

#' @keywords internal
utils::globalVariables(c(
  # dplyr NSE variables used in predict and plot functions
  "time", "marker", "id", 
  "Estimate", "Median", "Est.Error", "L95", "U95",
  # Additional variables used in internal functions
  ".data"
))
