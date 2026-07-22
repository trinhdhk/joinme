#' Plot longitudinal trajectories from JoiNMe objects
#'
#' @param object A `JoiNMeFit` or `JoiNMeDynPred` object.
#' @param longitudinal_times Optional numeric vector of fitted trajectory times
#'   on the original study-time scale. This argument applies to `JoiNMeFit`
#'   objects; existing `JoiNMeDynPred` objects retain their stored time grid.
#' @param longitudinal_points Number of evenly spaced fitted trajectory times
#'   used for a `JoiNMeFit` object when `longitudinal_times = NULL`.
#' @param ... Additional arguments forwarded to [plot()].
#'
#' @return A `ggplot` object, a combined plot, or a named list of plots.
#' @export
longitudinal_plot <- function(object,
                              longitudinal_times = NULL,
                              longitudinal_points = 80L,
                              ...) {
  # Step 1: For a fitted model, evaluate posterior trajectories on the explicit
  # or default smooth time design before constructing the longitudinal display.
  if (inherits(object, "JoiNMeFit")) {
    return(plot.JoiNMeFit(
      object,
      type = "longitudinal",
      longitudinal_times = longitudinal_times,
      longitudinal_points = longitudinal_points,
      ...,
      .use_wrapper_dispatch = FALSE
    ))
  }

  # Step 2: For an existing dynamic-prediction object, retain the time design
  # on which its posterior trajectories were originally estimated.
  if (inherits(object, "JoiNMeDynPred")) {
    return(plot.JoiNMeDynPred(object, type = "longitudinal", ..., .use_wrapper_dispatch = FALSE))
  }

  # Step 3: Reject objects that do not contain fitted or predicted trajectories.
  cli::cli_abort("{.arg object} must inherit from JoiNMeFit or JoiNMeDynPred.")
}

#' Plot survival trajectories from JoiNMe objects
#'
#' @param object A `JoiNMeFit` or `JoiNMeDynPred` object.
#' @param ... Additional arguments forwarded to [plot()].
#'
#' @return A `ggplot` object, a combined plot, or a named list of plots.
#' @export
survival_plot <- function(object, ...) {
  if (inherits(object, "JoiNMeFit")) {
    return(plot.JoiNMeFit(object, type = "survival", ..., .use_wrapper_dispatch = FALSE))
  }
  if (inherits(object, "JoiNMeDynPred")) {
    return(plot.JoiNMeDynPred(object, type = "survival", ..., .use_wrapper_dispatch = FALSE))
  }
  cli::cli_abort("{.arg object} must inherit from JoiNMeFit or JoiNMeDynPred.")
}

#' Plot cumulative hazard trajectories from JoiNMe objects
#'
#' @param object A `JoiNMeFit` or `JoiNMeDynPred` object.
#' @param ... Additional arguments forwarded to [plot()].
#'
#' @return A `ggplot` object, a combined plot, or a named list of plots.
#' @export
cumhaz_plot <- function(object, ...) {
  if (inherits(object, "JoiNMeFit")) {
    return(plot.JoiNMeFit(object, type = "cumhaz", ..., .use_wrapper_dispatch = FALSE))
  }
  if (inherits(object, "JoiNMeDynPred")) {
    return(plot.JoiNMeDynPred(object, type = "cumhaz", ..., .use_wrapper_dispatch = FALSE))
  }
  cli::cli_abort("{.arg object} must inherit from JoiNMeFit or JoiNMeDynPred.")
}

#' Plot fitted association curves from a JoiNMe fit
#'
#' @param object A `JoiNMeFit` object.
#' @param ... Additional arguments forwarded to [plot()].
#'
#' @return A `ggplot` object, a combined plot, or a named list of plots.
#' @export
association_plot <- function(object, ...) {
  if (!inherits(object, "JoiNMeFit")) {
    cli::cli_abort("{.arg object} must be a JoiNMeFit object for association plotting.")
  }
  plot.JoiNMeFit(object, type = "association", ..., .use_wrapper_dispatch = FALSE)
}

#' Plot posterior diagnostics from a JoiNMe fit
#'
#' @param object A `JoiNMeFit` object.
#' @param type Diagnostic plot type.
#' @param ... Additional arguments forwarded to [plot()].
#'
#' @return A `ggplot` object.
#' @export
diagnostic_plot <- function(object,
                            type = c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_sd", "running_mean", "running_quantile"),
                            ...) {
  if (!inherits(object, "JoiNMeFit")) {
    cli::cli_abort("{.arg object} must be a JoiNMeFit object for diagnostic plotting.")
  }
  type <- match.arg(type)
  plot.JoiNMeFit(object, type = type, ..., .use_wrapper_dispatch = FALSE)
}

#' Plot posterior MCMC summaries using bayesplot
#'
#' @description
#' Provides a JoiNMe-friendly wrapper around `bayesplot::mcmc_*` functions. The
#' posterior draws are first relabelled with user-facing parameter names via
#' [draws()], after which the selected bayesplot geometry is applied.
#'
#' @param object A `JoiNMeFit` or `JoiNMeDynPred` object.
#' @param pars Deprecated alias of `variable`.
#' @param type Plot type. Supported values are `"intervals"`, `"areas"`,
#'   `"dens"`, `"dens_overlay"`, `"hist"`, `"trace"`, `"violin"`,
#'   `"acf"`, `"rhat"`, and `"neff"`.
#' @param variable Optional variable names after relabelling.
#' @param regex Logical; treat `variable` as a regular expression.
#' @param fixed Deprecated logical alias controlling exact versus regex matching
#'   when `pars` is supplied.
#' @param draws Optional number of posterior draws to keep.
#' @param seed Integer seed used when subsetting draws.
#' @param ... Additional arguments passed to the selected `bayesplot` function.
#'
#' @return A `ggplot` object.
#' @export
mcmc_plot <- function(object,
                      pars = NA,
                      type = c("intervals", "areas", "dens", "dens_overlay", "hist", "trace", "violin", "acf", "rhat", "neff"),
                      variable = NULL,
                      regex = FALSE,
                      fixed = FALSE,
                      draws = NULL,
                      seed = 1,
                      ...) {
  if (!inherits(object, c("JoiNMeFit", "JoiNMeDynPred"))) {
    cli::cli_abort("{.arg object} must inherit from JoiNMeFit or JoiNMeDynPred.")
  }
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg bayesplot} is required for {.fn mcmc_plot}.")
  }

  type <- match.arg(type)
  if (is.null(variable) && !all(is.na(pars))) {
    variable <- pars
    if (isFALSE(fixed)) {
      regex <- TRUE
    }
  }

  draws_obj <- draws(object, format = "draws_array")
  all_vars <- posterior::variables(draws_obj)
  if (is.null(variable)) {
    variable <- setdiff(all_vars, "lp__")
    if (!length(variable)) {
      variable <- all_vars
    }
    variable <- utils::head(variable, 8L)
  }

  if (type %in% c("rhat", "neff")) {
    if (!inherits(object, "JoiNMeFit")) {
      cli::cli_abort("Sampler diagnostic MCMC plots are currently available only for JoiNMeFit objects.")
    }
    return(diagnostic_plot(
      object,
      type = if (identical(type, "rhat")) "rhat" else "ess_bulk",
      pars = if (!isTRUE(regex)) variable else NULL,
      regex_pars = if (isTRUE(regex)) variable else NULL,
      draws = draws,
      seed = seed,
      ...
    ))
  }

  plot_draws <- draws(object, variables = variable, regex = regex, draws = draws, seed = seed, format = "draws_array")

  plot_fun <- switch(
    type,
    intervals = bayesplot::mcmc_intervals,
    areas = bayesplot::mcmc_areas,
    dens = bayesplot::mcmc_dens,
    dens_overlay = bayesplot::mcmc_dens_overlay,
    hist = bayesplot::mcmc_hist,
    trace = bayesplot::mcmc_trace,
    violin = bayesplot::mcmc_violin,
    acf = bayesplot::mcmc_acf
  )

  plot_fun(plot_draws, ...)
}
