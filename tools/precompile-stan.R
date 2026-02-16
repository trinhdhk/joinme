# Precompile Stan models at installation time (best-effort).
# This script is intentionally dependency-light and should never fail install.

quiet_require <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}

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
    invisible(NULL)
}

stan_dir <- if (dir.exists("stan")) "stan" else file.path("inst", "stan")
stan_files <- list.files(
    path = stan_dir,
    pattern = "\\.stan$",
    full.names = TRUE,
    recursive = FALSE
)

# ---- CmdStanR precompile (if available) -----------------------------------   
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
        quack("- Precompiling Stan models with CmdStanR...")
        for (sf in stan_files) {
            cmdstanr::cmdstan_model(
                sf,
                cpp_options = list(STAN_THREADS = TRUE), #if (grepl('thread', sf, fixed=TRUE)) list(stan_threads = TRUE) else NULL,
                include_paths = dirname(sf),
                dir = cache_dir,
                force_recompile = TRUE
            )
        }
    } else {
        quack(
            " -- CmdStanR detected but CmdStan not installed; skipping CmdStanR precompile."
        )
    }
} else {
    quack(" -- CmdStanR not available; skipping CmdStanR precompile.")
}

# ---- RStan precompile (default) -------------------------------------------
if (quiet_require("rstan")) {
    quack(
        " - Precompiling Stan models with RStan... You shall expect chaotic output here."
    )
    options(stan.threads = 2)
    rstantools::rstan_config()
    # system("echo \"PKG_CXXFLAGS += -Wa,-mbig-obj\" >> ./src/Makevars.win")
    # system("echo \"PKG_CXXFLAGS += -Wa,-mbig-obj\" >> ./src/Makevars")
    options(stan.threads = NULL)
    #   rstan_cache <- file.path(system.file("stan", package = "joinme"), "pre-compiled")
    #   if (!dir.exists(rstan_cache)) {
    #     dir.create(rstan_cache, recursive = TRUE, showWarnings = FALSE)
    #   }

    #   for (sf in stan_files) {
    #     rds_path <- file.path(rstan_cache, paste0(basename(sf), ".rds"))
    #     if (file.exists(rds_path)) next

    #     try_safely(quote({
    #       old_wd <- getwd()
    #       setwd(dirname(sf))
    #       on.exit(setwd(old_wd), add = TRUE)
    #       rstan::rstan_options(auto_write = TRUE)
    #       mod <- rstan::stan_model(
    #         file = basename(sf),
    #         save_dso = TRUE,
    #         verbose = FALSE
    #       )
    #       saveRDS(mod, rds_path)
    #     }))
    #   }
}


quit(save = "no", status = 0)
