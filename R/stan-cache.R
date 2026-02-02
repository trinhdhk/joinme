#' @name joinme_stan_cache
#' @title Stan Model Caching Helpers
#'
#' @description
#' These helpers centralize compiled Stan model caching for both CmdStanR
#' and RStan. The goal is to keep compilation deterministic, reusable, and
#' aggressively documented for maintenance and debugging.
#'
#' @keywords internal
NULL

#' Resolve cache directory for compiled Stan models
#'
#' @description
#' Creates a dedicated cache directory under the user's R cache root, then
#' returns its absolute path. This is the single source of truth for cached
#' CmdStanR and RStan compiled objects.
#'
#' @return A character scalar with the cache directory path.
#'
#' @keywords internal
.stan_cache_dir <- function() {
  # Use the standardized per-package cache location for portability.
  dir <- tools::R_user_dir("joinme", "cache")
  # Ensure the directory exists before we return it.
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

#' Resolve cached CmdStanR model
#'
#' @description
#' Builds (or retrieves) a cached CmdStanR compiled model for a given Stan file.
#' This function also validates that CmdStan is installed and discoverable.
#'
#' @param stan_file Path to the Stan file.
#' @param cpp_options Optional C++ options for compilation.
#' @param force_recompile Logical; force recompilation even if cached.
#'
#' @return A CmdStanR model object.
#'
#' @keywords internal
.get_cmdstan_model <- function(stan_file, cpp_options = NULL, force_recompile = FALSE) {
  # Validate inputs aggressively.
  assertthat::assert_that(is.character(stan_file), length(stan_file) == 1)
  assertthat::assert_that(is.logical(force_recompile), length(force_recompile) == 1)

  if (!file.exists(stan_file)) {
    cli::cli_abort(c(
      x = "Stan file not found: {stan_file}.",
      i = "Check installation or restore inst/stan files."
    ))
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli::cli_abort(c(
      x = "Package {.pkg cmdstanr} is required for engine = 'cmdstanr'.",
      i = "Install with install.packages('cmdstanr') and set up CmdStan."
    ))
  }

  # Ensure CmdStan is installed and discoverable.
  ver <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NA)
  if (is.na(ver)) {
    cli::cli_abort(c(
      x = "CmdStan is not installed or not detected by cmdstanr.",
      i = "Run cmdstanr::install_cmdstan() and retry."
    ))
  }

  # Use a stable cache directory so the compiled model can be reused.
  cache_dir <- file.path(.stan_cache_dir(), "cmdstanr")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Compile or load the model using CmdStanR's model cache.
  cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = cpp_options,
    include_paths = dirname(stan_file),
    dir = cache_dir,
    force_recompile = force_recompile
  )
}

#' Resolve cached rstan model
#'
#' @description
#' Builds (or retrieves) a cached rstan compiled model for a given Stan file.
#' The compiled model is stored as an RDS file in the package cache.
#'
#' @param stan_file Path to the Stan file.
#' @param force_recompile Logical; force recompilation even if cached.
#'
#' @return A stanfit model object.
#'
#' @keywords internal
.get_rstan_model <- function(stan_file, force_recompile = FALSE) {
  # Validate inputs aggressively.
  assertthat::assert_that(is.character(stan_file), length(stan_file) == 1)
  assertthat::assert_that(is.logical(force_recompile), length(force_recompile) == 1)

  if (!file.exists(stan_file)) {
    cli::cli_abort(c(
      x = "Stan file not found: {stan_file}.",
      i = "Check installation or restore inst/stan files."
    ))
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    cli::cli_abort(c(
      x = "Package {.pkg rstan} is required for engine = 'rstan'.",
      i = "Install with install.packages('rstan') or use engine = 'cmdstanr'."
    ))
  }

  # Use a stable cache directory for compiled rstan models.
  cache_dir <- file.path(.stan_cache_dir(), "rstan")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  rds_path <- file.path(cache_dir, paste0(basename(stan_file), ".rds"))

  # Return cached model unless recompilation is explicitly requested.
  if (file.exists(rds_path) && !isTRUE(force_recompile)) {
    return(readRDS(rds_path))
  }

  # Compile and cache the model.
  rstan::rstan_options(auto_write = TRUE)
  old_wd <- getwd()
  setwd(dirname(stan_file))
  on.exit(setwd(old_wd), add = TRUE)
  mod <- rstan::stan_model(
    file = basename(stan_file),
    save_dso = TRUE,
    verbose = FALSE
  )
  saveRDS(mod, rds_path)
  mod
}

#' Precompile rstan models on load
#'
#' @description
#' Precompiles all packaged Stan models for rstan into the cache, allowing
#' faster model usage later. This runs on package load when enabled.
#'
#' @return Invisibly returns TRUE if any models were attempted.
#'
#' @keywords internal
.warm_rstan_cache <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  # Enumerate all packaged Stan model files we want to cache.
  stan_files <- c(
    system.file("stan/joinme_fit.stan", package = "joinme"),
    system.file("stan/joinme_fit_threading.stan", package = "joinme"),
    system.file("stan/joinme_dynpred.stan", package = "joinme"),
    system.file("stan/joinme_dynpred_threading.stan", package = "joinme")
  )
  stan_files <- stan_files[nzchar(stan_files) & file.exists(stan_files)]
  if (length(stan_files) == 0) {
    return(invisible(FALSE))
  }

  # Attempt compilation one-by-one to keep errors isolated.
  for (sf in stan_files) {
    try(.get_rstan_model(sf, force_recompile = FALSE), silent = TRUE)
  }
  invisible(TRUE)
}

#' Precompile CmdStanR models
#'
#' @description
#' Manually precompiles the packaged Stan models using CmdStanR. This is
#' useful when you want to force compilation ahead of time, rather than
#' waiting for the first call with engine = "cmdstanr".
#'
#' @param force_recompile Logical; recompile even if cached.
#' @param threads_per_chain Integer; number of threads for Stan (>= 1).
#'
#' @return Invisibly returns a named list of compiled CmdStanR models.
#'
#' @export
precompile_cmdstanr_models <- function(force_recompile = FALSE, threads_per_chain = 1L) {
  # Validate inputs with assertthat, then provide structured cli errors.
  assertthat::assert_that(is.logical(force_recompile), length(force_recompile) == 1)
  assertthat::assert_that(is.numeric(threads_per_chain), length(threads_per_chain) == 1)

  threads_per_chain <- as.integer(threads_per_chain)
  if (threads_per_chain < 1L) {
    cli::cli_abort(c(
      x = "{.arg threads_per_chain} must be >= 1.",
      i = "Use threads_per_chain = 1 for non-threaded compilation."
    ))
  }

  # Prepare the Stan file list in a deterministic order.
  stan_files <- c(
    joinme_fit = system.file("stan/joinme_fit.stan", package = "joinme"),
    joinme_fit_threading = system.file("stan/joinme_fit_threading.stan", package = "joinme"),
    joinme_dynpred = system.file("stan/joinme_dynpred.stan", package = "joinme"),
    joinme_dynpred_threading = system.file("stan/joinme_dynpred_threading.stan", package = "joinme")
  )
  stan_files <- stan_files[nzchar(stan_files) & file.exists(stan_files)]
  if (length(stan_files) == 0) {
    cli::cli_abort(c(
      x = "No packaged Stan files were found to compile.",
      i = "Reinstall the package or restore inst/stan files."
    ))
  }

  # Compile all models, enabling threading if requested.
  cpp_opts <- if (threads_per_chain > 1L) list(stan_threads = TRUE) else NULL
  models <- lapply(stan_files, function(sf) {
    .get_cmdstan_model(
      stan_file = sf,
      cpp_options = cpp_opts,
      force_recompile = force_recompile
    )
  })

  # Return an invisibly named list for programmatic use.
  names(models) <- names(stan_files)
  invisible(models)
}
