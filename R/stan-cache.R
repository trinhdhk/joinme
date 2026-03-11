#' @name joinme_stan_cache
#' @title Stan Model Caching Helpers
#'
#' @description
#' These helpers centralise compiled Stan model caching for both CmdStanR
#' and RStan. The goal is to keep compilation deterministic, reusable, and
#' aggressively documented for maintenance and debugging.
#'
#' @keywords internal
NULL

# File overview:
# - Resolve per-user cache location for compiled Stan models.
# - Provide CmdStanR/RStan model lookup with deterministic caching.

#' Find cmdstan exe file
#' @description #' Helper function to find the compiled CmdStanR model executable in the cache. #' This is used internally to ensure we are using the cached model for sampling. #' #' @param stan_file Path to the Stan file. #' #' @return The path to the compiled executable, or NULL if not found.
#'
#' @return A character scalar with the path to the compiled CmdStanR model executable, or NULL if not found. #' @keywords internal .find_cmdstan_exe <- function(stan_file) { cache_dir <- .stan_cache_dir() is_windows <- isTRUE(.Platform$OS.type == "windows") ext <- if (is_windows) ".exe" else "" exe_path <- file.path( cache_dir, paste0(tools::file_path_sans_ext(basename(stan_file)), ext) ) if (file.exists(exe_path)) { return(exe_path) } NULL }
#'
#' @keywords internal
.get_cmdstan_exe <- function(stan_file, cache_dir, warn_missing = TRUE, model_base_name = NULL) {
  # Rip-off from cmdstanr
  is_windows <- isTRUE(.Platform$OS.type == "windows")
  is_wsl <- is_windows &&
    (grepl('//wsl$/', tolower(cmdstanr::cmdstan_path()), fixed = TRUE) ||
      Sys.getenv("CMDSTANR_USE_WSL") == 1)
  # If on Windows and not using WSL, look for exe format
  # The problem is sometimes cmdstanr compiled on WSL but then
  # WSL is not on when running R the next session or vice versa. So we should look for both formats in the cache directory and use the one that exists.
  # We should safeguard this case by looking for both exe and non-exe formats in the cache directory
  # And warn out if cmdstanr run incompatible mode
  base_name <- model_base_name %||% tools::file_path_sans_ext(basename(stan_file))

  if (is_windows && !is_wsl) {
    # Search for native Windows
    exe_path <- file.path(
      cache_dir,
      paste0(base_name, '.exe')
    )
    if (!file.exists(exe_path)) {
      exe_path_wsl <- file.path(
        cache_dir,
        base_name
      )
      if (file.exists(exe_path_wsl) && warn_missing) {
        cli::cli_warn(c(
          x = "Cmdstan model is not compiled in native Windows mode or cached executable not found: {exe_path}.",
          i = "`cmdstanr` may have switched mode, WSL has been shut-down, or you did not compile the model when installing `joinme`. Call `joinme::precompile_cmdstanr_models()` to compile and cache the models for the current environment."
        ))
      }
    }
  } else {
    exe_path <- file.path(
      cache_dir,
      base_name
    )
  }
  exe_path
}

#' Resolve packaged/local Stan source file
#'
#' @param program One of "joinme_fit" or "joinme_dynpred".
#' @param threaded Logical; whether to use the threading Stan twin.
#'
#' @return Absolute or relative path to an existing Stan source file.
#' @keywords internal
.get_stan_file <- function(program = c("joinme_fit", "joinme_dynpred"), threaded = FALSE) {
  program <- match.arg(program)
  suffix <- if (isTRUE(threaded)) "_threading" else ""
  file_name <- paste0(program, suffix, ".stan")

  stan_candidates <- c(
    system.file(file.path("stan", file_name), package = "joinme"),
    file.path("inst", "stan", file_name),
    file.path("..", "inst", "stan", file_name),
    file.path("..", "..", "inst", "stan", file_name)
  )

  stan_file <- stan_candidates[file.exists(stan_candidates)][1]
  if (is.na(stan_file) || !nzchar(stan_file)) {
    cli::cli_abort(c(
      x = "Stan model file not found: {file_name}.",
      i = "Reinstall the package or restore inst/stan files."
    ))
  }
  stan_file
}

#' Resolve cached CmdStan executable path for a Stan source file
#'
#' @param stan_file Path to Stan source file.
#' @param cpp_options Optional C++ options that affect cache key.
#' @param warn_missing Logical; whether to warn if executable is missing.
#'
#' @return Character scalar executable path.
#' @keywords internal
.get_stan_exe <- function(stan_file, cpp_options = NULL, warn_missing = TRUE) {
  cache_dir <- .stan_cache_dir()
  cache_key <- .stan_cache_key(stan_file, cpp_options = cpp_options)
  model_base_name <- paste0(tools::file_path_sans_ext(basename(stan_file)), "-", cache_key)
  .get_cmdstan_exe(
    stan_file = stan_file,
    cache_dir = cache_dir,
    warn_missing = warn_missing,
    model_base_name = model_base_name
  )
}

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
  # Canonical per-user cache directory for compiled Stan models
  # Use the standardised per-package cache location for portability.
  dir <- tools::R_user_dir("joinme", "cache")
  # Ensure the directory exists before we return it.
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

#' Resolve recursive Stan include dependencies
#'
#' @param stan_file Path to root Stan file.
#' @param visited Internal recursion guard.
#'
#' @return Character vector of dependency file paths (including the root file).
#' @keywords internal
.stan_dependency_files <- function(stan_file, visited = character()) {
  if (!file.exists(stan_file)) return(character())

  stan_file <- normalizePath(stan_file, winslash = "/", mustWork = TRUE)
  if (stan_file %in% visited) return(character())

  visited <- c(visited, stan_file)
  out <- stan_file

  lines <- tryCatch(readLines(stan_file, warn = FALSE), error = function(e) character())
  if (length(lines) == 0) return(out)

  include_lines <- grep("^\\s*#include\\s+", lines, value = TRUE)
  if (length(include_lines) == 0) return(out)

  inc_tokens <- sub("^\\s*#include\\s+", "", include_lines)
  inc_tokens <- trimws(gsub("[\"<>]", "", inc_tokens))
  inc_tokens <- inc_tokens[nzchar(inc_tokens)]
  if (length(inc_tokens) == 0) return(out)

  for (inc in inc_tokens) {
    inc_file <- file.path(dirname(stan_file), inc)
    if (!file.exists(inc_file)) next
    out <- c(out, .stan_dependency_files(inc_file, visited = visited))
  }

  unique(out)
}

#' Build deterministic Stan cache key from full source graph
#'
#' @param stan_file Path to root Stan file.
#' @param cpp_options Optional CmdStan C++ options.
#'
#' @return Short hexadecimal cache key.
#' @keywords internal
.stan_cache_key <- function(stan_file, cpp_options = NULL) {
  deps <- unique(.stan_dependency_files(stan_file))
  deps <- deps[file.exists(deps)]
  deps <- sort(normalizePath(deps, winslash = "/", mustWork = TRUE))
  dep_md5 <- if (length(deps) > 0) tools::md5sum(deps) else character()

  root_file <- normalizePath(stan_file, winslash = "/", mustWork = TRUE)
  root_dir <- dirname(root_file)
  dep_ids <- vapply(deps, function(p) {
    prefix <- paste0(root_dir, "/")
    rel <- if (startsWith(p, prefix)) substr(p, nchar(prefix) + 1L, nchar(p)) else p
    if (!nzchar(rel) || identical(rel, p)) {
      rel <- basename(p)
    }
    rel
  }, character(1))

  cpp_sig <- ""
  if (!is.null(cpp_options) && length(cpp_options) > 0) {
    if (is.null(names(cpp_options))) names(cpp_options) <- rep("", length(cpp_options))
    ord <- order(names(cpp_options))
    cpp_parts <- vapply(ord, function(i) {
      paste0(names(cpp_options)[i], "=", paste(cpp_options[[i]], collapse = ","))
    }, character(1))
    cpp_sig <- paste(cpp_parts, collapse = ";")
  }

  payload <- c(
    "joinme-stan-cache-v2",
    paste(dep_ids, unname(dep_md5), sep = "="),
    paste0("cpp:", cpp_sig)
  )

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(payload, con = tmp, useBytes = TRUE)
  substr(unname(tools::md5sum(tmp)), 1, 16)
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
#' @param warn_missing Logical; whether to warn out if no exe file exists
#' @return A CmdStanR model object.
#'
#' @keywords internal
.get_cmdstan_model <- function(
  stan_file,
  cpp_options = NULL,
  force_recompile = FALSE,
  warn_missing = TRUE
) {
  assertthat::assert_that(is.character(stan_file), length(stan_file) == 1)
  assertthat::assert_that(
    is.logical(force_recompile),
    length(force_recompile) == 1
  )

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
  ver <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NA
  )
  if (is.na(ver)) {
    cli::cli_abort(c(
      x = "CmdStan is not installed or not detected by cmdstanr.",
      i = "Run cmdstanr::install_cmdstan() and retry."
    ))
  }

  # Use a stable cache directory so the compiled model can be reused.
  cache_dir <- .stan_cache_dir()


  exe_path <- .get_stan_exe(
    stan_file = stan_file,
    cpp_options = cpp_options,
    warn_missing = warn_missing
  )
  if (!file.exists(exe_path) && warn_missing) {
    cli::cli_inform(
      c(
        'Executable file not found {exe_path}. Compiling model. 
        To save time in the future, call `joinme::precompile_cmdstanr_models()` to compile and cache the models ahead of time.'
      )
    )
  }
  need_compile <- isTRUE(force_recompile) || !file.exists(exe_path)
  if (need_compile) {
    return(cmdstanr::cmdstan_model(
      stan_file = stan_file,
      exe_file = exe_path,
      compile = TRUE,
      cpp_options = cpp_options,
      include_paths = dirname(stan_file),
      dir = cache_dir,
      force_recompile = force_recompile
    ))
  }
  # writeLines(exe_path)
  # browser()
  mod <- cmdstanr::cmdstan_model(
    exe_file = exe_path,
    compile = FALSE,
    cpp_options = cpp_options,
    include_paths = dirname(stan_file),
    dir = cache_dir
  )
  # HACK: current parser parses cpp_options to UPPERCASE but sampler wants lowercase.
  # SEE: https://github.com/stan-dev/cmdstanr/blob/453084bb67996c64b0c49b363312f9a60181b391/R/cpp_opts.R#L53
  # TODO: FUTURE version of cmdstanr should fix it.
  names(mod$.__enclos_env__$private$cpp_options_) <- tolower(names(mod$.__enclos_env__$private$cpp_options_))
  mod
}

#' Resolve cached rstan model
#'
#' @description
#' Builds (or retrieves) a cached rstan compiled model for a given Stan file.
#' The compiled model is stored as an RDS file in the package cache.
#'
#' @param stan_file Path to the Stan file.
#'
#' @return A stanfit model object.
#'
#' @keywords internal
.get_rstan_model <- function(stan_file) {
  # Compile from current Stan source to avoid stale precompiled module/data-schema mismatches.
  if (!file.exists(stan_file)) {
    cli::cli_abort(c(
      x = "Stan file not found: {stan_file}",
      i = "Check the model path passed to {.fn .get_rstan_model}."
    ))
  }

  if (!exists(".joinme_rstan_model_cache", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".joinme_rstan_model_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  cache_env <- get(".joinme_rstan_model_cache", envir = .GlobalEnv, inherits = FALSE)

  file_info <- file.info(stan_file)
  cache_key <- paste0(normalizePath(stan_file, winslash = "/", mustWork = TRUE), "::", file_info$mtime)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  model_name <- tools::file_path_sans_ext(basename(stan_file))
  stanc_ret <- rstan::stanc(
    file = stan_file,
    model_name = model_name,
    allow_undefined = TRUE,
    verbose = FALSE
  )
  mod <- rstan::stan_model(
    stanc_ret = stanc_ret,
    auto_write = FALSE,
    verbose = FALSE
  )

  assign(cache_key, mod, envir = cache_env)
  mod
}

#' Precompile CmdStanR models
#'
#' @description
#' Manually precompiles the packaged Stan models using CmdStanR. This is
#' useful when you want to force compilation ahead of time, rather than
#' waiting for the first call with engine = "cmdstanr".
#'
#' @param force_recompile Logical; recompile even if cached.
#' @param cleanup Logical; whether to clean up old cached executables before compiling.
#' @return Invisibly returns a named list of compiled CmdStanR models.
#'
#' @export
precompile_cmdstanr_models <- function(force_recompile = TRUE, cleanup = TRUE) {
  # Precompile packaged Stan models into the cache
  # Validate inputs with assertthat, then provide structured cli errors.
  assertthat::assert_that(
    is.logical(force_recompile),
    length(force_recompile) == 1
  )

  # Prepare the Stan file list in a deterministic order.
  stan_files <- c(
    joinme_fit = system.file("stan/joinme_fit.stan", package = "joinme"),
    joinme_fit_threading = system.file(
      "stan/joinme_fit_threading.stan",
      package = "joinme"
    ),
    joinme_dynpred = system.file(
      "stan/joinme_dynpred.stan",
      package = "joinme"
    ),
    joinme_dynpred_threading = system.file(
      "stan/joinme_dynpred_threading.stan",
      package = "joinme"
    )
  )
  stan_files <- stan_files[nzchar(stan_files) & file.exists(stan_files)]
  if (length(stan_files) == 0) {
    cli::cli_abort(c(
      x = "No packaged Stan files were found to compile.",
      i = "Something has removed the stan files. Reinstall the package"
    ))
  }

  # Clean-up old cached executables that match the naming pattern to prevent stale models. We can be aggressive here since the cache key includes file hashes, so old executables won't be reused anyway. This ensures that if the source files change, we won't accidentally use an old cached executable that doesn't match the new source.
  cache_dir <- .stan_cache_dir()
  if (cleanup) {
    old_exes <- list.files(cache_dir, pattern = "^joinme_.*\\.(exe)?$", full.names = TRUE)
    if (length(old_exes) > 0) {
      unlink(old_exes, force = TRUE)
    }
  }
  
  # Compile all models, enabling threading if requested.
  models <- lapply(stan_files, function(sf) {
    cat(sf, '')
    cpp_opts <- list(stan_threads = grepl('threading', sf, fixed = TRUE))
    .get_cmdstan_model(
      stan_file = sf,
      cpp_options = cpp_opts,
      force_recompile = force_recompile,
      warn_missing = FALSE
    )
  })

  # Return an invisibly named list for programmatic use.
  names(models) <- names(stan_files)
  invisible(models)
}
