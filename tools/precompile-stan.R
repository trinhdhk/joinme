# Precompile Stan models at installation time (best-effort).
# This script is intentionally dependency-light and should never fail install.

quiet_require <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

# This will quack somewhat like message() or cli::cli_inform()
# Why the names? Why not.
quack <- function(msg) {
    if (quiet_require("cli")) {
        cli::cli_inform(msg)
    } else {
        message(msg)
    }
    invisible(NULL)
}

try_safely <- function(expr) {
    try(eval(expr), silent = FALSE)
}

`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

stan_dependency_files <- function(stan_file, visited = character()) {
    if (!file.exists(stan_file)) {
        return(character())
    }

    stan_file <- normalizePath(stan_file, winslash = "/", mustWork = TRUE)
    if (stan_file %in% visited) {
        return(character())
    }

    visited <- c(visited, stan_file)
    out <- stan_file

    lines <- tryCatch(readLines(stan_file, warn = FALSE), error = function(e) {
        character()
    })
    if (length(lines) == 0) {
        return(out)
    }

    include_lines <- grep("^\\s*#include\\s+", lines, value = TRUE)
    if (length(include_lines) == 0) {
        return(out)
    }

    inc_tokens <- sub("^\\s*#include\\s+", "", include_lines)
    inc_tokens <- trimws(gsub("[\"<>]", "", inc_tokens))
    inc_tokens <- inc_tokens[nzchar(inc_tokens)]
    if (length(inc_tokens) == 0) {
        return(out)
    }

    for (inc in inc_tokens) {
        inc_file <- file.path(dirname(stan_file), inc)
        if (!file.exists(inc_file)) {
            next
        }
        out <- c(out, stan_dependency_files(inc_file, visited = visited))
    }

    unique(out)
}

stan_cache_key <- function(stan_file, cpp_options = NULL) {
    deps <- unique(stan_dependency_files(stan_file))
    deps <- deps[file.exists(deps)]
    deps <- sort(normalizePath(deps, winslash = "/", mustWork = TRUE))
    dep_md5 <- if (length(deps) > 0) tools::md5sum(deps) else character()

    root_file <- normalizePath(stan_file, winslash = "/", mustWork = TRUE)
    root_dir <- dirname(root_file)
    dep_ids <- vapply(
        deps,
        function(p) {
            prefix <- paste0(root_dir, "/")
            rel <- if (startsWith(p, prefix)) {
                substr(p, nchar(prefix) + 1L, nchar(p))
            } else {
                p
            }
            if (!nzchar(rel) || identical(rel, p)) {
                rel <- basename(p)
            }
            rel
        },
        character(1)
    )

    cpp_sig <- ""
    if (!is.null(cpp_options) && length(cpp_options) > 0) {
        if (is.null(names(cpp_options))) {
            names(cpp_options) <- rep("", length(cpp_options))
        }
        ord <- order(names(cpp_options))
        cpp_parts <- vapply(
            ord,
            function(i) {
                paste0(
                    names(cpp_options)[i],
                    "=",
                    paste(cpp_options[[i]], collapse = ",")
                )
            },
            character(1)
        )
        cpp_sig <- paste(cpp_parts, collapse = ";")
    }

    payload <- c(
        "JoiNMe-stan-cache-v2",
        paste(dep_ids, unname(dep_md5), sep = "="),
        paste0("cpp:", cpp_sig)
    )

    tmp <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(payload, con = tmp, useBytes = TRUE)
    substr(unname(tools::md5sum(tmp)), 1, 16)
}

is_wsl_mode <- function() {
    if (!isTRUE(.Platform$OS.type == "windows")) {
        return(FALSE)
    }
    cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
    grepl("//wsl$/", tolower(cmdstan_path), fixed = TRUE) ||
        Sys.getenv("CMDSTANR_USE_WSL") == "1"
}

cached_exe_path <- function(cache_dir, model_base_name) {
    if (isTRUE(.Platform$OS.type == "windows") && !is_wsl_mode()) {
        file.path(cache_dir, paste0(model_base_name, ".exe"))
    } else {
        file.path(cache_dir, model_base_name)
    }
}

clean_cache_dir <- function(cache_dir) {
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        return(invisible(NULL))
    }

    cached_files <- list.files(
        cache_dir,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    if (length(cached_files) > 0) {
        unlink(cached_files, recursive = TRUE, force = TRUE)
    }
    invisible(NULL)
}

stan_dir <- if (dir.exists("stan")) "stan" else file.path("inst", "stan")
stan_files <- list.files(
    path = stan_dir,
    pattern = "_threading\\.stan$",
    full.names = TRUE,
    recursive = FALSE
)

# ---- RStan precompile --------------------------------------------
orig_zzz <- readLines("R/zzz.R", warn = FALSE)
on.exit(writeLines(orig_zzz, con = "R/zzz.R"), add = TRUE)

quack("- Precompiling Stan models with Rstan...")
try_safely(quote({
    old_stan_thread <- getOption("stan.thread")
    options(stan.thread = 1L)
    on.exit(options(stan.thread = old_stan_thread), add = TRUE)
    # on.exit(unlink('R/stanmodels.R', force = TRUE), add = TRUE)
    rstantools::rstan_config()
    stan_version <- rstan::stan_version()
    # cat(
    #     paste0('.METADATA$stan_version <- "', stan_version, '"\n'),
    #     file = 'R/zzz.R',
    #     append = TRUE
    # )

    # orig <- readLines("NAMESPACE", warn = FALSE)
    # on.exit(writeLines(orig, con = "NAMESPACE"), add = TRUE)
    # cat('useDynLib(joinme, .registration = TRUE)\n', file = 'NAMESPACE', append = TRUE)
    # invisible(lapply(stan_files, function(sf) {
    #     rstan::stan_model(file = sf, save_dso = FALSE, verbose = FALSE)
    # }))
}))


# ---- CmdStanR precompile (if available) -----------------------------------
compile_cmdstan <- identical(
    Sys.getenv("JOINME_COMPILE_CMDSTANR", unset = "0"),
    "1"
)
if (compile_cmdstan && quiet_require("cmdstanr")) {
    cache_dir <- tools::R_user_dir("joinme", "cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (quiet_require("cmdstanr")) {
        ver <- tryCatch(
            cmdstanr::cmdstan_version(error_on_NA = FALSE),
            error = function(e) NA
        )
        if (!is.na(ver)) {
            quack("- Cleaning Stan cache directory before recompilation...")
            clean_cache_dir(cache_dir)

            quack("- Precompiling Stan models with CmdStanR...")
            for (sf in stan_files) {
                cpp_options <- list(stan_threads = TRUE)

                key <- stan_cache_key(sf, cpp_options = cpp_options)
                model_base_name <- paste0(
                    tools::file_path_sans_ext(basename(sf)),
                    "-",
                    key
                )
                exe_file <- cached_exe_path(
                    cache_dir = cache_dir,
                    model_base_name = model_base_name
                )

                cmdstanr::cmdstan_model(
                    stan_file = sf,
                    exe_file = exe_file,
                    compile = TRUE,
                    cpp_options = cpp_options,
                    include_paths = dirname(sf),
                    dir = cache_dir,
                    force_recompile = FALSE
                )
                ver <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
        # if (!is.na(ver)) {
        #     cat(
        #         paste0('.METADATA$cmdstan_version <- "', ver, '"\n'),
        #         file = 'R/zzz.R',
        #         append = TRUE
        #     )
        # }
            }
        } else {
            quack(
                " -- CmdStanR detected but CmdStan not installed; skipping CmdStanR precompile."
            )
        }
    } else {
        quack(" -- CmdStanR not available; skipping CmdStanR precompile.")
    }
} else if (!compile_cmdstan && quiet_require("cmdstanr")) {
    quack(
        " -- Skipping CmdStanR precompile. If you want to force compiling in CmdStanR, set JOINME_COMPILE_CMDSTANR = 1; "
    )
}


quit(save = "no", status = 0)
