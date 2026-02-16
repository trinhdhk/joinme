# File overview:
# - Registers NSE symbols used in dplyr pipelines.
# - Prevents R CMD check notes for predict/plot helpers.

# Global variable bindings for dplyr NSE and other functions
# These are used in dplyr pipelines and prevent R CMD check warnings

#' @keywords internal
utils::globalVariables(c(
  # Avoid NOTES for NSE variables used in dplyr/pipeline helpers
  # dplyr NSE variables used in predict and plot functions
  "time", "marker", "id", 
  "Estimate", "Median", "Est.Error", "L95", "U95",
  # Additional variables used in internal functions
  ".data"
))


# .onLoad <- function(libname, pkgname){
  
# }