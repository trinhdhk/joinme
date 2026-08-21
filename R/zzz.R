# File overview:
# - Registers NSE symbols used in tidy evaluation helpers.
# - Prevents R CMD check notes for predict/plot helpers.

# Global variable bindings for tidy evaluation helpers
# These are used in tidytable/dplyr-style code and prevent R CMD check warnings

#' @keywords internal
utils::globalVariables(c(
  # Avoid NOTES for NSE variables used in dplyr/pipeline helpers
  # dplyr NSE variables used in predict and plot functions
  "time",
  "marker",
  "id",
  "variable",
  "metric",
  "iteration",
  "value",
  "chain",
  "stat",
  "Estimate",
  "Median",
  "Est.Error",
  "L95",
  "U95",
  "marker_idx",
  # Additional variables used in internal functions
  ".data",
  ".by",
  ".SD",
  ".",
  ":="
))

.onAttach <- function(libname, pkgname) {
  cached_dir_files <- list.files(.stan_cache_dir(), full.names = TRUE)
  has_model_exe <- any(grepl("joinme_fit", cached_dir_files)) &&
    any(grepl("joinme_dynpred", cached_dir_files))
  cli::simple_theme()
  cli::cli_h2('Join Mixed-Effects Model')
  cli::cli_text('Version: {utils::packageVersion("joinme")}')
  # cli::cli_alert_success('RStan version: { stanmodels$.META$stan_version }')
  if (has_model_exe) {
    cli::cli_alert_success('CmdStanR model compiled')
  } else {
    cli::cli_alert_danger('CmdStanR model not compiled.')
    if (sample(c(TRUE, FALSE), size = 1, prob = c(0.2, 0.8))) {
      cli::cli_alert_info(
        'Call `joinme::precompile_cmdstanr_models()` to compile and cache the models for the current environment.'
      )
    }
  }
  if (sample(c(TRUE, FALSE), size = 1, prob = c(0.001, 0.999)) | .devmode()) {
      .oucru_logo_ascii()
  }
}

.devmode <- function() {
  Sys.getenv("JOINME_DEVMODE", unset = "0") == "1"
}

.oucru_logo_ascii <- function(){
  readLines(system.file('etc/oucru', package = 'joinme')) |>
    paste(collapse = "\n") |>
    cat()
}
