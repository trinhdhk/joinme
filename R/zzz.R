# File overview:
# - Registers NSE symbols used in tidy evaluation helpers.
# - Prevents R CMD check notes for predict/plot helpers.

# Global variable bindings for tidy evaluation helpers
# These are used in tidytable/dplyr-style code and prevent R CMD check warnings

#' @keywords internal
utils::globalVariables(c(
  # Avoid NOTES for NSE variables used in dplyr/pipeline helpers
  # dplyr NSE variables used in predict and plot functions
  "time", "marker", "id", 
  "variable", "metric", "iteration", "value", "chain", "stat",
  "Estimate", "Median", "Est.Error", "L95", "U95", "marker_idx",
  # Additional variables used in internal functions
  ".data", ".by", ".SD", ".", ":="
))
