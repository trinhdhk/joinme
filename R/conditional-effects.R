#' Conditional effects for a fitted joint model
#'
#' @description
#' `conditional_effects()` evaluates model-implied changes in the longitudinal
#' and event processes while holding all predictors not named in `effects` at
#' explicit conditioning values. Its interface follows
#' [brms::conditional_effects()] where the joint-model structure permits a
#' direct correspondence.
#'
#' Unlike a univariate regression, a `JoiNMeFit` contains two statistical
#' processes. The `process` argument therefore selects the longitudinal
#' expected response, the event relative hazard, or both. Longitudinal effects
#' may contain only the fixed-effect contribution, add the fitted marker-level
#' deviation, or average marker-specific predictions over the selected markers.
#' Subject- and marker-by-subject deviations are excluded from all three
#' longitudinal estimands. Event effects describe the proportional-hazards
#' multiplier `exp(W gamma)` (or its logarithm), with the baseline hazard and
#' longitudinal association contribution held fixed. Consequently, event
#' contrasts isolate the part of the event process attributable to covariates
#' in `formulaEvent`.
#'
#' The returned data follow the `brms` conditional-effects convention. Each
#' effect table contains `estimate__`, `se__`, `lower__`, `upper__`, and
#' `cond__`, plus `effect1__` and, for two-predictor effects, `effect2__`.
#' JoiNMe adds `process__`, `marker__`, `longitudinal_estimand__`, and
#' `event_type__` where applicable. This permits JoiNMe to delegate the
#' graphical construction to the tested `brms` plotting method rather than
#' maintaining a second plotting grammar.
#'
#' @param x A fitted object of class `JoiNMeFit`.
#' @param effects Effects to evaluate. Supply a character vector such as
#'   `c("time", "time:treatment")` to use the same requests wherever they are
#'   valid, or a named list with `longitudinal` and `event` components to choose
#'   process-specific effects. Interactions may contain at most two predictors.
#'   If `NULL`, all model predictors and fitted two-way interactions are used.
#' @param conditions Optional data frame containing values of predictors on
#'   which to condition. One set of effects is evaluated for every row.
#'   `conditions$cond__`, when supplied, provides the facet label; otherwise row
#'   names are used. Tables returned by [make_conditions()] can be supplied
#'   directly.
#' @param int_conditions Optional named list controlling the evaluation values
#'   of predictors named in `effects`. Each element may be a vector or a
#'   function applied to the observed predictor. By default, the first numeric
#'   predictor spans its observed range, factors use all levels, and a second
#'   numeric predictor is evaluated at its mean and mean plus or minus one
#'   standard deviation. This mirrors the principal `brms` convention.
#' @param process Character vector selecting `"longitudinal"`, `"event"`, or
#'   both. The shorthand `"both"` and the default
#'   `c("longitudinal", "event")` both return the two processes.
#' @param prob Probability covered by the equal-tailed posterior uncertainty
#'   interval. The default is `0.95`.
#' @param robust Logical. If `TRUE`, posterior medians define `estimate__`; if
#'   `FALSE`, posterior means are used.
#' @param method Longitudinal posterior scale. Supported values are
#'   `"posterior_epred"` and `"posterior_linpred"`. The event process is
#'   controlled separately by `event_scale`.
#' @param longitudinal_estimand Longitudinal quantity to evaluate.
#'   `"population"` uses only the fixed-effect contribution, corresponding to
#'   exclusion of group-level effects in `brms`. `"marker"` adds the fitted
#'   marker-level deviation but excludes subject and marker-by-subject
#'   deviations. `"marginal_marker"` first calculates each selected marker's
#'   posterior prediction, including its marker-level deviation and fitted
#'   inverse link, and then takes an equally weighted mean across markers within
#'   every posterior draw. Thus, on the expected-response scale it estimates an
#'   average of marker responses rather than the response obtained from an
#'   averaged linear predictor.
#' @param event_scale Event-process estimand. `"hazard_ratio"` returns the
#'   covariate-specific proportional-hazards multiplier `exp(W gamma)`;
#'   `"log_hazard_ratio"` returns `W gamma`.
#' @param resolution Number of support points for the first continuous
#'   predictor. For a two-dimensional surface it is used on both axes.
#' @param surface Logical. If `TRUE`, two continuous predictors are evaluated
#'   over a rectangular surface. If `FALSE`, the second continuous predictor is
#'   represented by a small set of conditioning curves.
#' @param markers Optional character vector restricting longitudinal results to
#'   selected marker levels. By default all fitted markers are used. For
#'   `longitudinal_estimand = "marginal_marker"`, this argument defines the
#'   markers entering the equally weighted posterior average.
#' @param draws Optional positive integer limiting the posterior draws used in
#'   the calculation. `NULL` uses all available draws.
#' @param seed Integer seed used when posterior draws are subsampled.
#' @param plot Logical. If `TRUE`, construct and display conditional-effects
#'   plots. If `FALSE`, return the underlying `JoiNMeConditionalEffects` data
#'   object. The data object can subsequently be plotted with `plot()`.
#' @param ... Additional arguments passed to the plotting method when
#'   `plot = TRUE`, for example `points`, `rug`, `stype`, or `theme`; see
#'   [brms::conditional_effects()] for the corresponding graphical arguments.
#'
#' @return If `plot = FALSE`, a `JoiNMeConditionalEffects` object: a named list
#'   with one component per requested process. Each process component is a
#'   `brms_conditional_effects`-compatible named list containing one data frame
#'   per effect. Longitudinal tables identify their estimand in
#'   `longitudinal_estimand__`. If `plot = TRUE`, a similarly nested named list
#'   of `ggplot` objects is returned invisibly after plotting.
#'
#' @details
#' Numeric predictors not named in an effect are held at their observed mean;
#' factors are held at their first level. Values supplied in `conditions`
#' override those defaults. Group-level subject deviations are excluded, as in
#' the default `re_formula = NA` behavior of `brms::conditional_effects()`.
#' Marker-level deviations are also excluded for the `"population"` estimand,
#' retained for the `"marker"` estimand, and integrated by finite averaging for
#' the `"marginal_marker"` estimand.
#'
#' The event result is a relative-hazard effect rather than a dynamic survival
#' prediction. Dynamic survival probabilities depend on a subject's observed
#' marker history and conditioning time and therefore remain the estimand of
#' [predict.JoiNMeFit()]. Keeping these two estimands separate prevents a table
#' of baseline covariate effects from being mislabeled as subject-specific
#' dynamic prediction.
#'
#' @examples
#' \dontrun{
#' # Default longitudinal and event effects, returned as plot-ready data.
#' cond_eff <- conditional_effects(fit, plot = FALSE)
#'
#' # Compare treatment profiles for a longitudinal time effect and an event
#' # age effect, using labels generated by brms::make_conditions().
#' profiles <- make_conditions(fit$dataEvent, vars = "treatment")
#' cond_eff <- conditional_effects(
#'   fit,
#'   effects = list(longitudinal = "time", event = "age"),
#'   conditions = profiles,
#'   process = c("longitudinal", "event"),
#'   plot = FALSE
#' )
#' plot(cond_eff, ask = FALSE)
#'
#' # Request marker-specific longitudinal expected responses.
#' conditional_effects(
#'   fit,
#'   effects = "time:treatment",
#'   process = "longitudinal",
#'   longitudinal_estimand = "marker",
#'   int_conditions = list(treatment = c("control", "active"))
#' )
#'
#' # Average marker-specific expected responses draw by draw.
#' conditional_effects(
#'   fit,
#'   effects = "time",
#'   process = "longitudinal",
#'   longitudinal_estimand = "marginal_marker",
#'   markers = c("marker_1", "marker_2"),
#'   plot = FALSE
#' )
#' }
#'
#' @seealso [brms::conditional_effects()], [make_conditions()],
#'   [predict.JoiNMeFit()]
#' @export
conditional_effects.JoiNMeFit <- function(
    x,
    effects = NULL,
    conditions = NULL,
    int_conditions = NULL,
    process = c("longitudinal", "event"),
    prob = 0.95,
    robust = TRUE,
    method = c("posterior_epred", "posterior_linpred"),
    longitudinal_estimand = c("population", "marker", "marginal_marker"),
    event_scale = c("hazard_ratio", "log_hazard_ratio"),
    resolution = 100L,
    surface = FALSE,
    markers = NULL,
    draws = NULL,
    seed = 1,
    plot = TRUE,
    ...) {
  # Step 1: Validate the joint-model object and all scalar controls before any
  # posterior calculation. Early validation keeps errors attached to the
  # scientific request rather than to a later model-matrix operation.
  if (!inherits(x, "JoiNMeFit")) {
    cli::cli_abort("{.arg x} must be a fitted {.cls JoiNMeFit} object.")
  }
  process <- match.arg(process, c("both", "longitudinal", "event"), several.ok = TRUE)
  if ("both" %in% process) process <- c("longitudinal", "event")
  process <- unique(process)
  method <- match.arg(method)
  longitudinal_estimand <- match.arg(longitudinal_estimand)
  event_scale <- match.arg(event_scale)

  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1) {
    cli::cli_abort("{.arg prob} must be a single number strictly between 0 and 1.")
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    cli::cli_abort("{.arg robust} must be either TRUE or FALSE.")
  }
  if (!is.logical(surface) || length(surface) != 1L || is.na(surface)) {
    cli::cli_abort("{.arg surface} must be either TRUE or FALSE.")
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    cli::cli_abort("{.arg plot} must be either TRUE or FALSE.")
  }
  resolution <- as.integer(resolution)
  if (length(resolution) != 1L || !is.finite(resolution) || resolution < 2L) {
    cli::cli_abort("{.arg resolution} must be an integer greater than or equal to 2.")
  }
  if (!is.null(draws)) {
    draws <- as.integer(draws)
    if (length(draws) != 1L || !is.finite(draws) || draws < 1L) {
      cli::cli_abort("{.arg draws} must be NULL or a positive integer.")
    }
  }
  if (!is.list(int_conditions) && !is.null(int_conditions)) {
    cli::cli_abort("{.arg int_conditions} must be NULL or a named list.")
  }

  # Step 2: Recover process-specific fixed-effect predictors and default effect
  # sets. The event response variables on the left of Surv() are intentionally
  # excluded because the event estimand here is the covariate hazard multiplier.
  process_specs <- list(
    longitudinal = .conditional_effect_process_spec(
      formula = reformulas::nobars(x$formulaLong),
      data = x$dataLong,
      process = "longitudinal"
    ),
    event = .conditional_effect_process_spec(
      formula = x$formulaEvent,
      data = x$dataEvent,
      process = "event"
    )
  )
  effects_by_process <- .conditional_effect_requests(
    effects = effects,
    selected_processes = process,
    process_specs = process_specs
  )

  # Step 3: Validate and label the shared conditioning rows. A single condition
  # table is deliberately shared across processes so longitudinal and event
  # panels describe the same named scientific profiles.
  all_model_variables <- unique(unlist(lapply(process_specs[process], `[[`, "variables"), use.names = FALSE))
  condition_rows <- .conditional_effect_conditions(
    conditions = conditions,
    model_variables = all_model_variables
  )

  # Step 4: Build each requested process independently. Both builders return
  # the exact table contract expected by plot.brms_conditional_effects().
  out <- setNames(vector("list", length(process)), process)
  if ("longitudinal" %in% process) {
    out$longitudinal <- .conditional_effects_longitudinal(
      object = x,
      effect_requests = effects_by_process$longitudinal,
      process_spec = process_specs$longitudinal,
      condition_rows = condition_rows,
      int_conditions = int_conditions,
      probability = prob,
      robust_center = robust,
      prediction_method = method,
      longitudinal_estimand = longitudinal_estimand,
      grid_resolution = resolution,
      draw_surface = surface,
      marker_selection = markers,
      posterior_draws = draws,
      random_seed = seed
    )
  }
  if ("event" %in% process) {
    out$event <- .conditional_effects_event(
      object = x,
      effect_requests = effects_by_process$event,
      process_spec = process_specs$event,
      condition_rows = condition_rows,
      int_conditions = int_conditions,
      probability = prob,
      robust_center = robust,
      event_estimand = event_scale,
      grid_resolution = resolution,
      draw_surface = surface,
      posterior_draws = draws,
      random_seed = seed
    )
  }

  # Step 5: Attach a compact description of the estimands and calculation so a
  # stored result remains self-describing when it is plotted or inspected later.
  out <- structure(
    out,
    class = c("JoiNMeConditionalEffects", "list"),
    process = process,
    probability = prob,
    robust = robust,
    method = method,
    longitudinal_estimand = longitudinal_estimand,
    event_scale = event_scale,
    call = match.call()
  )

  # Step 6: Match the requested convenience behavior. plot = FALSE exposes the
  # brms-compatible effect tables; plot = TRUE delegates their graphical
  # representation to the class method below.
  if (!isTRUE(plot)) {
    return(out)
  }
  invisible(plot(out, plot = TRUE, ...))
}

#' Plot JoiNMe conditional effects
#'
#' @description
#' Plots the process-specific effect tables returned by
#' `conditional_effects(..., plot = FALSE)`. Each process is delegated to
#' [brms::conditional_effects()]'s plotting method, retaining its line,
#' interval, surface, point, rug, and theme conventions.
#'
#' @param x A `JoiNMeConditionalEffects` object.
#' @param plot Logical. If `TRUE`, draw each plot in the active graphics device.
#'   If `FALSE`, return the plots without drawing them.
#' @param ask Logical. Whether to prompt before drawing a subsequent plot.
#' @param ... Arguments passed to the `brms_conditional_effects` plot method,
#'   including `points`, `rug`, `stype`, `line_args`, `surface_args`, and
#'   `theme`.
#'
#' @return Invisibly, a named list by process and effect containing `ggplot`
#'   objects.
#' @export
plot.JoiNMeConditionalEffects <- function(x, plot = TRUE, ask = FALSE, ...) {
  # Validate the stored conditional-effect contract before delegating.
  if (!inherits(x, "JoiNMeConditionalEffects")) {
    cli::cli_abort("{.arg x} must inherit from {.cls JoiNMeConditionalEffects}.")
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    cli::cli_abort("{.arg plot} must be either TRUE or FALSE.")
  }

  # Delegate construction and optional drawing to brms for each process.
  # Passing plot and ask through preserves the established brms behavior. Empty
  # process components are retained so the requested process selection remains
  # visible in the returned structure.
  plots <- lapply(x, function(process_tables) {
    if (!length(process_tables)) {
      return(list())
    }
    plot(process_tables, plot = plot, ask = ask, ...)
  })
  invisible(plots)
}

#' Print a JoiNMe conditional-effects data object
#'
#' @param x A `JoiNMeConditionalEffects` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.JoiNMeConditionalEffects <- function(x, ...) {
  process_names <- names(x)
  effect_counts <- vapply(x, length, integer(1))
  .cli_summary_heading("Conditional effects for JoiNMe model\n")
  longitudinal_estimand <- attr(x, "longitudinal_estimand", exact = TRUE)
  if ("longitudinal" %in% process_names && !is.null(longitudinal_estimand)) {
    cli::cli_bullets('*' = sprintf(" Longitudinal estimand: %s\n", longitudinal_estimand))
  }
  for (i in seq_along(process_names)) {
    cli::cli_bullets('*' = sprintf("%s process: %d effect%s\n", toTitleCase(process_names[[i]]), effect_counts[[i]], if (effect_counts[[i]] == 1L) "" else "s"))
  }
  invisible(x)
}

#' Recover predictors and default conditional-effect terms for one process
#'
#' @param formula Model formula for one process.
#' @param data Observed process data.
#' @param process Process name.
#'
#' @return A list containing predictor variables and default effects.
#' @keywords internal
#' @noRd
.conditional_effect_process_spec <- function(formula, data, process) {
  # Step 1: Remove the response and group-level terms so the effect vocabulary
  # contains only population-level predictors that can be varied in new rows.
  formula_use <- if (identical(process, "event")) {
    stats::delete.response(stats::terms(formula))
  } else {
    stats::delete.response(stats::terms(reformulas::nobars(formula)))
  }
  variables <- intersect(all.vars(formula_use), names(data))

  # Step 2: Translate fitted term labels to one- or two-predictor requests.
  # Transformations such as splines remain represented by their underlying
  # predictor, while interactions retain the familiar colon notation.
  term_labels <- attr(formula_use, "term.labels") %||% character(0)
  interaction_effects <- lapply(term_labels, function(term_label) {
    term_variables <- tryCatch(
      all.vars(stats::as.formula(paste("~", term_label))),
      error = function(e) character(0)
    )
    term_variables <- unique(intersect(term_variables, variables))
    if (length(term_variables) == 2L) paste(term_variables, collapse = ":") else character(0)
  })
  defaults <- unique(c(variables, unlist(interaction_effects, use.names = FALSE)))

  list(
    process = process,
    formula = formula,
    variables = variables,
    effects = defaults,
    data = data
  )
}

#' Resolve process-specific effect requests
#'
#' @param effects User effect specification.
#' @param selected_processes Requested processes.
#' @param process_specs Process metadata.
#'
#' @return Named list of validated effects by process.
#' @keywords internal
#' @noRd
.conditional_effect_requests <- function(effects, selected_processes, process_specs) {
  # Step 1: Interpret a named list as an explicit process-specific request and a
  # character vector as a shared request. Missing list entries use defaults.
  if (is.null(effects)) {
    requested <- lapply(process_specs, `[[`, "effects")
  } else if (is.list(effects)) {
    invalid_names <- setdiff(names(effects), c("longitudinal", "event"))
    if (length(invalid_names)) {
      cli::cli_abort("Unknown process names in {.arg effects}: {paste(invalid_names, collapse = ', ')}.")
    }
    requested <- lapply(names(process_specs), function(process_name) {
      if (process_name %in% names(effects)) effects[[process_name]] else process_specs[[process_name]]$effects
    })
    names(requested) <- names(process_specs)
  } else {
    requested <- setNames(rep(list(as.character(effects)), length(process_specs)), names(process_specs))
  }

  # Step 2: Validate each effect through its underlying predictor names. Manual
  # two-way combinations are allowed when both predictors occur in the process,
  # matching brms behavior even when the interaction was not fitted explicitly.
  any_valid <- FALSE
  out <- setNames(vector("list", length(process_specs)), names(process_specs))
  for (process_name in names(process_specs)) {
    process_requests <- unique(as.character(requested[[process_name]] %||% character(0)))
    process_requests <- process_requests[nzchar(process_requests)]
    valid <- vapply(process_requests, function(effect_name) {
      effect_variables <- unique(strsplit(effect_name, ":", fixed = TRUE)[[1L]])
      length(effect_variables) %in% c(1L, 2L) && all(effect_variables %in% process_specs[[process_name]]$variables)
    }, logical(1))
    out[[process_name]] <- process_requests[valid]
    if (process_name %in% selected_processes && any(valid)) any_valid <- TRUE
  }

  # Step 3: Fail only when the explicit request is invalid for every selected
  # process. A predictor belonging to one process need not be present in the
  # other process of the joint model.
  if (!is.null(effects) && !any_valid) {
    valid_summary <- vapply(selected_processes, function(process_name) {
      paste(process_specs[[process_name]]$variables, collapse = ", ")
    }, character(1))
    cli::cli_abort(c(
      x = "None of the requested {.arg effects} are valid for the selected processes.",
      i = "Available predictors: {paste(paste0(selected_processes, ': ', valid_summary), collapse = '; ')}."
    ))
  }
  out
}

#' Validate and label conditioning rows
#'
#' @param conditions Optional condition data frame.
#' @param model_variables Variables used by either selected process.
#'
#' @return A unique data frame with a `cond__` label column.
#' @keywords internal
#' @noRd
.conditional_effect_conditions <- function(conditions, model_variables) {
  # Step 1: Represent the default profile by one row. Predictor defaults are
  # filled later from the process-specific observed data.
  if (is.null(conditions)) {
    return(data.frame(cond__ = factor("1", levels = "1"), stringsAsFactors = FALSE))
  }
  conditions <- as.data.frame(conditions)
  if (!nrow(conditions)) {
    cli::cli_abort("{.arg conditions} must contain at least one row.")
  }
  conditions <- unique(conditions)

  # Step 2: Follow brms label precedence: an explicit cond__ column takes
  # priority, followed by informative row names, followed by numbered labels.
  if ("cond__" %in% names(conditions)) {
    condition_labels <- as.character(conditions$cond__)
  } else {
    condition_labels <- rownames(conditions)
    automatic_row_names <- is.null(condition_labels) || identical(condition_labels, as.character(seq_len(nrow(conditions))))
    if (automatic_row_names) condition_labels <- paste("condition", seq_len(nrow(conditions)))
  }
  if (anyNA(condition_labels) || any(!nzchar(condition_labels)) || anyDuplicated(condition_labels)) {
    cli::cli_abort("Condition labels in {.code cond__} or row names must be non-missing and unique.")
  }
  conditions$cond__ <- factor(condition_labels, levels = condition_labels)

  # Step 3: Retain unused columns for informative output, but alert the analyst
  # because those values cannot alter either process.
  unused <- setdiff(names(conditions), c(model_variables, "cond__"))
  if (length(unused)) {
    cli::cli_warn("Condition variables not used by either selected process: {paste(unused, collapse = ', ')}.")
  }
  conditions
}

#' Construct a brms-style conditional evaluation grid
#'
#' @param process_spec Process metadata.
#' @param effect_name One requested effect.
#' @param condition_rows Labelled conditions.
#' @param int_conditions Optional custom focal values.
#' @param grid_resolution Continuous grid resolution.
#' @param draw_surface Whether two numeric predictors form a surface.
#'
#' @return Evaluation data with effect metadata attributes.
#' @keywords internal
#' @noRd
.conditional_effect_grid <- function(process_spec,
                                     effect_name,
                                     condition_rows,
                                     int_conditions,
                                     grid_resolution,
                                     draw_surface) {
  effect_variables <- unique(strsplit(effect_name, ":", fixed = TRUE)[[1L]])
  source_data <- process_spec$data

  # Step 1: Place a numeric predictor on the first axis when an interaction
  # contains one numeric and one factor predictor. This reproduces the plotting
  # orientation used by brms and gives line plots their conventional x axis.
  effect_is_factor <- vapply(effect_variables, function(variable) {
    value <- source_data[[variable]]
    is.factor(value) || is.character(value) || is.logical(value)
  }, logical(1))
  if (length(effect_variables) == 2L) {
    ordering <- order(effect_is_factor, decreasing = FALSE)
    effect_variables <- effect_variables[ordering]
    effect_is_factor <- effect_is_factor[ordering]
  }

  # Step 2: Determine focal support. Custom int_conditions take precedence;
  # otherwise factors use all fitted levels and numeric variables follow the
  # range/representative-value rules documented above.
  focal_values <- vector("list", length(effect_variables))
  names(focal_values) <- effect_variables
  for (j in seq_along(effect_variables)) {
    variable <- effect_variables[[j]]
    observed <- source_data[[variable]]
    custom_values <- int_conditions[[variable]] %||% NULL
    if (is.function(custom_values)) custom_values <- custom_values(observed)

    if (!is.null(custom_values)) {
      values <- custom_values
    } else if (effect_is_factor[[j]]) {
      values <- if (is.factor(observed)) levels(observed) else unique(stats::na.omit(observed))
    } else if (j == 1L || isTRUE(draw_surface)) {
      observed_numeric <- as.numeric(observed)
      observed_range <- range(observed_numeric[is.finite(observed_numeric)])
      values <- if (diff(observed_range) > 0) {
        seq(observed_range[[1L]], observed_range[[2L]], length.out = grid_resolution)
      } else {
        observed_range[[1L]]
      }
    } else {
      observed_numeric <- as.numeric(observed)
      center <- mean(observed_numeric, na.rm = TRUE)
      spread <- stats::sd(observed_numeric, na.rm = TRUE)
      if (!is.finite(spread)) spread <- 0
      values <- unique(center + c(-1, 0, 1) * spread)
    }
    if (!length(values)) {
      cli::cli_abort("No finite evaluation values are available for effect predictor {.val {variable}}.")
    }
    focal_values[[variable]] <- .conditional_effect_cast(values, observed)
  }

  # Step 3: Form the focal Cartesian product and replicate it for each named
  # condition. The process model variables not in the effect are then filled by
  # the supplied condition or their observed reference value.
  focal_grid <- do.call(expand.grid, c(focal_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  grids <- vector("list", nrow(condition_rows))
  for (condition_index in seq_len(nrow(condition_rows))) {
    grid_i <- focal_grid
    for (variable in process_spec$variables) {
      if (variable %in% effect_variables) next
      if (variable %in% names(condition_rows) && !is.na(condition_rows[[variable]][condition_index])) {
        value <- condition_rows[[variable]][condition_index]
      } else {
        value <- .conditional_effect_reference(process_spec$data[[variable]])
      }
      grid_i[[variable]] <- .conditional_effect_cast(rep(value, nrow(grid_i)), process_spec$data[[variable]])
    }
    extra_condition_variables <- setdiff(names(condition_rows), c(names(grid_i), "cond__"))
    for (variable in extra_condition_variables) {
      grid_i[[variable]] <- rep(condition_rows[[variable]][condition_index], nrow(grid_i))
    }
    grid_i$cond__ <- as.character(condition_rows$cond__[[condition_index]])
    grids[[condition_index]] <- grid_i
  }
  evaluation_data <- do.call(rbind, grids)
  rownames(evaluation_data) <- NULL
  evaluation_data$cond__ <- factor(evaluation_data$cond__, levels = as.character(condition_rows$cond__))

  # Step 4: Add the special plotting columns used by brms. A second continuous
  # predictor is converted to an ordered factor for separate curves unless a
  # surface was explicitly requested.
  evaluation_data$effect1__ <- evaluation_data[[effect_variables[[1L]]]]
  if (length(effect_variables) == 2L) {
    effect2 <- evaluation_data[[effect_variables[[2L]]]]
    if (is.numeric(effect2) && !isTRUE(draw_surface)) {
      rounded <- signif(effect2, 6L)
      custom_labels <- names(focal_values[[effect_variables[[2L]]]])
      effect2 <- factor(rounded, levels = sort(unique(rounded), decreasing = TRUE))
      if (length(custom_labels) == nlevels(effect2)) levels(effect2) <- custom_labels
    }
    evaluation_data$effect2__ <- effect2
  }

  attr(evaluation_data, "effects") <- effect_variables
  attr(evaluation_data, "types") <- ifelse(effect_is_factor, "factor", "numeric")
  evaluation_data
}

#' Reference value for a conditioned predictor
#'
#' @param x Observed predictor.
#'
#' @return A scalar reference value.
#' @keywords internal
#' @noRd
.conditional_effect_reference <- function(x) {
  if (is.factor(x)) return(levels(x)[[1L]])
  if (is.character(x)) return(unique(stats::na.omit(x))[[1L]])
  if (is.logical(x)) return(FALSE)
  mean(as.numeric(x), na.rm = TRUE)
}

#' Cast conditional values to the observed predictor type
#'
#' @param values Values to cast.
#' @param observed Observed predictor defining the target type and levels.
#'
#' @return Values compatible with the fitted model matrix.
#' @keywords internal
#' @noRd
.conditional_effect_cast <- function(values, observed) {
  if (is.factor(observed)) {
    cast_values <- factor(as.character(values), levels = levels(observed), ordered = is.ordered(observed))
    observed_contrasts <- attr(observed, "contrasts", exact = TRUE)
    if (!is.null(observed_contrasts)) contrasts(cast_values) <- observed_contrasts
    return(cast_values)
  }
  if (is.character(observed)) {
    # Character predictors are converted to factors by model.matrix(). Retain
    # every observed level here so a profile containing only one category still
    # produces the full fitted dummy-variable design.
    observed_levels <- sort(unique(as.character(stats::na.omit(observed))))
    return(factor(as.character(values), levels = observed_levels))
  }
  if (is.logical(observed)) return(as.logical(values))
  if (inherits(observed, "Date")) return(as.Date(values, origin = "1970-01-01"))
  as.numeric(values)
}

#' Summarise posterior conditional predictions
#'
#' @param posterior_matrix Draw-by-row posterior matrix.
#' @param probability Interval probability.
#' @param robust_center Whether to use the median.
#'
#' @return Four-column posterior summary data frame.
#' @keywords internal
#' @noRd
.conditional_effect_summary <- function(posterior_matrix, probability, robust_center) {
  lower_probability <- (1 - probability) / 2
  upper_probability <- 1 - lower_probability
  estimate <- if (isTRUE(robust_center)) {
    apply(posterior_matrix, 2L, stats::median, na.rm = TRUE)
  } else {
    colMeans(posterior_matrix, na.rm = TRUE)
  }
  data.frame(
    estimate__ = as.numeric(estimate),
    se__ = apply(posterior_matrix, 2L, stats::sd, na.rm = TRUE),
    lower__ = apply(posterior_matrix, 2L, stats::quantile, probs = lower_probability, names = FALSE, na.rm = TRUE),
    upper__ = apply(posterior_matrix, 2L, stats::quantile, probs = upper_probability, names = FALSE, na.rm = TRUE),
    check.names = FALSE
  )
}

#' Construct observed points for a longitudinal conditional-effect table
#'
#' @param object Fitted model.
#' @param effect_table Completed conditional-effect table.
#' @param condition_rows Labelled conditioning profiles.
#' @param marker_levels Marker strata retained in the table.
#' @param draw_surface Whether a two-numeric-predictor surface was requested.
#'
#' @return A data frame following the point-data contract used by the brms
#'   conditional-effects plot method.
#' @keywords internal
#' @noRd
.conditional_effect_longitudinal_points <- function(object,
                                                    effect_table,
                                                    condition_rows,
                                                    marker_levels,
                                                    draw_surface) {
  effect_variables <- attr(effect_table, "effects")
  marker_variable <- .JoiNMefit_call_arg_chr(object$call, "marker_var", "marker")
  response_variable <- .resolve_response_var(object$formulaLong, object$dataLong, context = "conditional_effects()")
  observed_data <- object$dataLong

  # Step 1: Reproduce every condition-marker facet represented in the posterior
  # table. Categorical conditioning variables are matched exactly; continuous
  # conditions retain all observations, consistent with brms until an explicit
  # point-distance rule is requested.
  point_blocks <- list()
  for (condition_index in seq_len(nrow(condition_rows))) {
    for (marker_level in marker_levels) {
      points <- observed_data[as.character(observed_data[[marker_variable]]) == marker_level, , drop = FALSE]
      categorical_conditions <- intersect(
        setdiff(names(condition_rows), c(effect_variables, "cond__")),
        names(points)
      )
      for (variable in categorical_conditions) {
        observed <- points[[variable]]
        if (is.factor(observed) || is.character(observed) || is.logical(observed)) {
          requested <- condition_rows[[variable]][condition_index]
          if (!is.na(requested)) points <- points[as.character(observed) == as.character(requested), , drop = FALSE]
        }
      }
      if (!nrow(points)) next

      # Step 2: Add the special effect and response columns consumed by brms.
      # Numeric second predictors are assigned to their nearest displayed curve;
      # surfaces retain their original numeric coordinate.
      points$effect1__ <- points[[effect_variables[[1L]]]]
      if (length(effect_variables) == 2L) {
        second_variable <- effect_variables[[2L]]
        second_observed <- points[[second_variable]]
        if (is.numeric(second_observed) && !isTRUE(draw_surface)) {
          focal_rows <- !duplicated(effect_table[[second_variable]])
          focal_values <- as.numeric(effect_table[[second_variable]][focal_rows])
          focal_labels <- as.character(effect_table$effect2__[focal_rows])
          nearest <- vapply(as.numeric(second_observed), function(value) {
            which.min(abs(focal_values - value))
          }, integer(1))
          points$effect2__ <- factor(focal_labels[nearest], levels = levels(effect_table$effect2__))
        } else if (is.factor(effect_table$effect2__)) {
          points$effect2__ <- factor(as.character(second_observed), levels = levels(effect_table$effect2__))
        } else {
          points$effect2__ <- second_observed
        }
      }
      points$resp__ <- points[[response_variable]]
      points$marker__ <- marker_level
      condition_label <- as.character(condition_rows$cond__[[condition_index]])
      points$cond__ <- if (length(marker_levels) > 1L) paste0(condition_label, " | marker: ", marker_level) else condition_label
      point_blocks[[length(point_blocks) + 1L]] <- points
    }
  }

  # Step 3: Return a typed zero-row frame when no observation matches a profile.
  # This keeps points = TRUE safe without implying that unmatched observations
  # belong to a conditional panel.
  if (!length(point_blocks)) {
    columns <- c("effect1__", if (length(effect_variables) == 2L) "effect2__", "resp__", "cond__")
    empty <- effect_table[0, intersect(columns, names(effect_table)), drop = FALSE]
    if (!("resp__" %in% names(empty))) empty$resp__ <- numeric(0)
    return(empty)
  }
  points <- do.call(rbind, point_blocks)
  rownames(points) <- NULL
  points$cond__ <- factor(points$cond__, levels = levels(effect_table$cond__))
  points
}

#' Construct a typed empty point frame
#'
#' @param effect_table Conditional-effect table.
#'
#' @return A zero-row frame safe for the brms plot method.
#' @keywords internal
#' @noRd
.conditional_effect_empty_points <- function(effect_table) {
  special_columns <- intersect(c("effect1__", "effect2__", "cond__"), names(effect_table))
  points <- effect_table[0, special_columns, drop = FALSE]
  points$resp__ <- numeric(0)
  points
}

#' Extract a requested posterior matrix without exposing backend differences
#'
#' @param object Fitted model.
#' @param variables Posterior variable names.
#' @param posterior_draws Optional draw limit.
#' @param random_seed Subsampling seed.
#'
#' @return A posterior draw matrix, or a zero-column matrix when unavailable.
#' @keywords internal
#' @noRd
.conditional_effect_draw_matrix <- function(object, variables, posterior_draws, random_seed) {
  tryCatch(
    .get_draws_matrix(object$fit, variables = variables, draws = posterior_draws, seed = random_seed),
    error = function(e) matrix(numeric(0), nrow = 0L, ncol = 0L)
  )
}

#' Evaluate longitudinal conditional draws for a declared estimand
#'
#' @param object Fitted model.
#' @param evaluation_data New longitudinal rows.
#' @param prediction_method Posterior scale.
#' @param longitudinal_estimand Whether to retain fixed effects only or add the
#'   fitted marker-level deviation.
#' @param posterior_draws Optional draw count.
#' @param random_seed Subsampling seed.
#'
#' @return Draw-by-row posterior prediction matrix.
#' @keywords internal
#' @noRd
.conditional_effect_longitudinal_draws <- function(object,
                                                   evaluation_data,
                                                   prediction_method,
                                                   longitudinal_estimand,
                                                   posterior_draws,
                                                   random_seed) {
  stan_data <- object$stan_data
  time_variable <- .JoiNMefit_call_arg_chr(object$call, "time_var", stan_data$time_var %||% "time")
  marker_variable <- .JoiNMefit_call_arg_chr(object$call, "marker_var", "marker")

  # Step 1: Reuse the fitted-trajectory design builder so conditional effects
  # inherit the exact fixed- and marker-effect bases used by plot.JoiNMeFit().
  # The builder applies the stored time transformation and validates every
  # model-matrix dimension against the fitted Stan data.
  longitudinal_design <- .JoiNMefit_longitudinal_design_matrices(
    fitted_model = object,
    longitudinal_evaluation_data = evaluation_data,
    time_variable = time_variable,
    marker_variable = marker_variable
  )
  fixed_design <- longitudinal_design$fixed
  marker_design <- longitudinal_design$marker
  expected_coefficients <- as.integer(stan_data$P %||% ncol(fixed_design))
  marker_coefficient_count <- as.integer(stan_data$R_mk %||% ncol(marker_design))
  fitted_marker_levels <- as.character(
    stan_data$marker_levels %||%
      levels(object$dataLong[[marker_variable]]) %||%
      unique(as.character(object$dataLong[[marker_variable]]))
  )
  include_marker_deviation <- longitudinal_estimand %in% c("marker", "marginal_marker")

  # Step 2: Request fixed- and, when required, marker-level draws together. A
  # joint request guarantees that every coefficient block refers to the same
  # sampled MCMC iterations when a draw limit is used.
  scaled_names <- paste0("beta_scaled[", seq_len(expected_coefficients), "]")
  marker_names <- if (include_marker_deviation && marker_coefficient_count > 0L) {
    as.vector(outer(
      seq_along(fitted_marker_levels),
      seq_len(marker_coefficient_count),
      function(marker_index, coefficient_index) {
        paste0("v_marker[", marker_index, ",", coefficient_index, "]")
      }
    ))
  } else {
    character(0)
  }
  requested_draws <- .conditional_effect_draw_matrix(
    object,
    c(scaled_names, marker_names),
    posterior_draws,
    random_seed
  )
  scaled_positions <- match(scaled_names, colnames(requested_draws))
  if (!anyNA(scaled_positions)) {
    coefficient_draws <- requested_draws[, scaled_positions, drop = FALSE]
  } else {
    # Step 3: Retain the established raw-coefficient reconstruction for fitted
    # objects that do not store beta_scaled. Marker draws are requested in the
    # same call so posterior iteration alignment is preserved in this branch.
    raw_names <- paste0("beta[", seq_len(expected_coefficients), "]")
    requested_draws <- .conditional_effect_draw_matrix(
      object,
      c(raw_names, marker_names),
      posterior_draws,
      random_seed
    )
    raw_positions <- match(raw_names, colnames(requested_draws))
    coefficient_draws <- if (!anyNA(raw_positions)) {
      requested_draws[, raw_positions, drop = FALSE]
    } else {
      matrix(numeric(0), nrow = nrow(requested_draws), ncol = 0L)
    }
    time_indices <- as.integer(stan_data$idx_time_beta %||% integer(0))
    time_indices <- time_indices[time_indices >= 1L & time_indices <= ncol(coefficient_draws)]
    if (length(time_indices)) {
      time_scale <- as.numeric(object$tmax %||% object$config$tmax %||% stan_data$tmax %||% 1)
      coefficient_draws[, time_indices] <- coefficient_draws[, time_indices, drop = FALSE] * time_scale
    }
  }
  if (!nrow(coefficient_draws) || ncol(coefficient_draws) != expected_coefficients) {
    cli::cli_abort("Could not extract fitted fixed-effect draws for longitudinal conditional effects.")
  }

  # Step 4: Form the fixed-effect linear predictor. Subject and
  # marker-by-subject random effects are intentionally absent for every
  # conditional-effects estimand because no observed individual is conditioned
  # upon.
  linear_predictor <- coefficient_draws %*% t(fixed_design)

  # Step 5: For marker-specific and marker-marginal estimands, add the fitted
  # marker-level deviation at every evaluation row. This is the same v_marker
  # contribution used by fitted trajectory plots, without either subject-level
  # contribution.
  if (include_marker_deviation && marker_coefficient_count > 0L) {
    marker_positions <- match(marker_names, colnames(requested_draws))
    if (anyNA(marker_positions)) {
      cli::cli_abort("Could not extract fitted marker-level draws for longitudinal conditional effects.")
    }
    marker_draws <- requested_draws[, marker_positions, drop = FALSE]
    marker_index <- match(as.character(evaluation_data[[marker_variable]]), fitted_marker_levels)
    if (anyNA(marker_index)) {
      cli::cli_abort("The longitudinal conditional grid contains an unknown marker level.")
    }
    draw_count <- nrow(coefficient_draws)
    marker_count <- length(fitted_marker_levels)
    for (coefficient_index in seq_len(marker_coefficient_count)) {
      block_columns <- (coefficient_index - 1L) * marker_count + seq_len(marker_count)
      coefficient_by_row <- marker_draws[, block_columns, drop = FALSE][, marker_index, drop = FALSE]
      design_by_draw <- matrix(
        rep(marker_design[, coefficient_index], each = draw_count),
        nrow = draw_count
      )
      linear_predictor <- linear_predictor + coefficient_by_row * design_by_draw
    }
  }

  if (identical(prediction_method, "posterior_linpred")) return(linear_predictor)

  # Step 6: Apply each marker's fitted inverse link before any requested marker
  # average. Performing this transformation at the draw-marker level is what
  # distinguishes an average response from a response evaluated at an average
  # linear predictor.
  marker_index <- match(as.character(evaluation_data[[marker_variable]]), fitted_marker_levels)
  expected_response <- matrix(NA_real_, nrow = nrow(linear_predictor), ncol = ncol(linear_predictor))
  for (marker_number in seq_along(fitted_marker_levels)) {
    columns <- which(marker_index == marker_number)
    if (!length(columns)) next
    eta <- linear_predictor[, columns, drop = FALSE]
    link_code <- as.integer(stan_data$link_long[[marker_number]] %||% 1L)
    if (link_code == 1L) {
      expected_response[, columns] <- eta
    } else if (link_code %in% c(2L, 5L)) {
      expected_response[, columns] <- exp(eta)
    } else if (link_code == 3L) {
      expected_response[, columns] <- stats::plogis(eta)
    } else if (link_code == 4L) {
      expected_response[, columns] <- stats::pnorm(eta)
    } else {
      operation_count <- as.integer(stan_data$inv_link_n_ops[[marker_number]] %||% 0L)
      constant_count <- as.integer(stan_data$inv_link_n_const[[marker_number]] %||% 0L)
      bytecode <- if (operation_count > 0L) as.integer(stan_data$inv_link_ops[marker_number, seq_len(operation_count)]) else integer(0)
      constants <- if (constant_count > 0L) as.numeric(stan_data$inv_link_const[marker_number, seq_len(constant_count)]) else numeric(0)
      for (column in columns) {
        expected_response[, column] <- eval_bytecode_vector(linear_predictor[, column], bytecode = bytecode, const_data = constants)
      }
    }
  }
  expected_response
}

#' Average marker-specific posterior predictions within each focal grid row
#'
#' @param posterior_matrix Draw-by-row marker-specific posterior predictions.
#' @param evaluation_data Marker-expanded conditional grid.
#' @param marker_levels Marker levels included in the marginal estimand.
#'
#' @return A list containing the draw-by-grid-row marginal posterior matrix and
#'   one evaluation row per original focal grid row.
#' @keywords internal
#' @noRd
.conditional_effect_marginalize_markers <- function(posterior_matrix,
                                                     evaluation_data,
                                                     marker_levels) {
  grid_row_variable <- ".conditional_grid_row__"
  base_condition_variable <- ".conditional_base_condition__"
  if (!(grid_row_variable %in% names(evaluation_data)) ||
      !(base_condition_variable %in% names(evaluation_data))) {
    cli::cli_abort("The marker-expanded conditional grid is missing its averaging index.")
  }

  # Step 1: Verify that every focal covariate row contains exactly one posterior
  # prediction for each selected marker. This establishes a finite, equally
  # weighted marker distribution at every point of the effect curve or surface.
  grid_rows <- unique(evaluation_data[[grid_row_variable]])
  marker_counts <- vapply(grid_rows, function(grid_row) {
    sum(evaluation_data[[grid_row_variable]] == grid_row)
  }, integer(1))
  if (any(marker_counts != length(marker_levels))) {
    cli::cli_abort("Each longitudinal focal grid row must contain every selected marker before marginalisation.")
  }

  # Step 2: Average expected responses (or linear predictors) inside every MCMC
  # draw. Posterior summaries are deliberately deferred until after averaging,
  # thereby retaining dependence among marker-specific fitted quantities.
  marginal_draws <- do.call(cbind, lapply(grid_rows, function(grid_row) {
    columns <- which(evaluation_data[[grid_row_variable]] == grid_row)
    rowMeans(posterior_matrix[, columns, drop = FALSE])
  }))

  # Step 3: Collapse the expanded data to one row per focal covariate value and
  # restore the original profile label, which no longer has a marker facet.
  representative_rows <- match(grid_rows, evaluation_data[[grid_row_variable]])
  marginal_data <- evaluation_data[representative_rows, , drop = FALSE]
  marginal_data$cond__ <- marginal_data[[base_condition_variable]]
  marginal_data$marker__ <- "Marginal over markers"
  rownames(marginal_data) <- NULL

  list(draws = marginal_draws, data = marginal_data)
}

#' Build longitudinal conditional-effect tables
#'
#' @param object Fitted model.
#' @param effect_requests Longitudinal focal effects to evaluate.
#' @param process_spec Longitudinal fixed-effect specification.
#' @param condition_rows Labelled conditioning profiles.
#' @param int_conditions Optional focal predictor values.
#' @param probability Posterior interval probability.
#' @param robust_center Whether posterior medians define point estimates.
#' @param prediction_method Longitudinal posterior scale.
#' @param longitudinal_estimand Population, marker-specific, or marker-marginal
#'   longitudinal quantity.
#' @param grid_resolution Number of continuous focal values.
#' @param draw_surface Whether two continuous predictors form a surface.
#' @param marker_selection Optional fitted marker levels to retain or average.
#' @param posterior_draws Optional number of posterior draws.
#' @param random_seed Seed used when posterior draws are subsampled.
#'
#' @return A `brms_conditional_effects`-compatible list of longitudinal effect
#'   tables.
#' @keywords internal
#' @noRd
.conditional_effects_longitudinal <- function(object,
                                              effect_requests,
                                              process_spec,
                                              condition_rows,
                                              int_conditions,
                                              probability,
                                              robust_center,
                                              prediction_method,
                                              longitudinal_estimand,
                                              grid_resolution,
                                              draw_surface,
                                              marker_selection,
                                              posterior_draws,
                                              random_seed) {
  marker_variable <- .JoiNMefit_call_arg_chr(object$call, "marker_var", "marker")
  marker_levels <- as.character(object$stan_data$marker_levels %||% levels(object$dataLong[[marker_variable]]) %||% unique(as.character(object$dataLong[[marker_variable]])))

  # Step 1: Restrict marker strata before expanding any evaluation grid.
  if (!is.null(marker_selection)) {
    unknown_markers <- setdiff(as.character(marker_selection), marker_levels)
    if (length(unknown_markers)) {
      cli::cli_abort("Unknown marker levels: {paste(unknown_markers, collapse = ', ')}.")
    }
    marker_levels <- marker_levels[marker_levels %in% as.character(marker_selection)]
  }
  if (!length(marker_levels)) {
    cli::cli_abort("At least one fitted marker must be retained for longitudinal conditional effects.")
  }

  # Step 2: Evaluate every requested effect and retain one brms-compatible table
  # per effect. Population and marker-specific estimands retain marker strata
  # because fitted inverse links may differ by marker. The marker-marginal
  # estimand removes those strata only after draw-level response calculation.
  tables <- setNames(vector("list", length(effect_requests)), effect_requests)
  for (effect_name in effect_requests) {
    grid <- .conditional_effect_grid(
      process_spec = process_spec,
      effect_name = effect_name,
      condition_rows = condition_rows,
      int_conditions = int_conditions,
      grid_resolution = grid_resolution,
      draw_surface = draw_surface
    )
    grid_effects <- attr(grid, "effects")
    grid_types <- attr(grid, "types")
    marker_is_focal <- marker_variable %in% grid_effects
    if (identical(longitudinal_estimand, "marginal_marker") && marker_is_focal) {
      cli::cli_abort(c(
        x = "A marker-marginal conditional effect cannot use the marker variable as a focal predictor.",
        i = "Choose another focal effect or use {.arg longitudinal_estimand} = {.val marker}."
      ))
    }

    # Step 3: Give each unexpanded focal row a stable averaging index and retain
    # its original condition label. These internal columns ensure the marker
    # average pairs identical predictor profiles even when the grid contains
    # interactions or multiple declared conditions.
    grid$.conditional_grid_row__ <- seq_len(nrow(grid))
    grid$.conditional_base_condition__ <- as.character(grid$cond__)

    if (marker_is_focal) {
      # When marker itself is a fitted fixed effect, retain its focal grid rather
      # than duplicating and overwriting it as an external outcome stratum.
      evaluation_data <- grid[as.character(grid[[marker_variable]]) %in% marker_levels, , drop = FALSE]
      evaluation_data$marker__ <- as.character(evaluation_data[[marker_variable]])
      base_condition <- evaluation_data$.conditional_base_condition__
      evaluation_data$cond__ <- if (length(marker_levels) > 1L) {
        paste0(base_condition, " | marker: ", evaluation_data$marker__)
      } else {
        base_condition
      }
    } else {
      marker_grids <- lapply(marker_levels, function(marker_level) {
        marker_grid <- grid
        marker_grid[[marker_variable]] <- factor(marker_level, levels = as.character(object$stan_data$marker_levels %||% marker_levels))
        marker_grid$marker__ <- marker_level
        base_condition <- marker_grid$.conditional_base_condition__
        marker_grid$cond__ <- if (length(marker_levels) > 1L) paste0(base_condition, " | marker: ", marker_level) else base_condition
        marker_grid
      })
      evaluation_data <- do.call(rbind, marker_grids)
    }
    rownames(evaluation_data) <- NULL

    # Step 4: Calculate marker-indexed posterior quantities before any summary.
    # The draw evaluator includes marker effects only for the marker-specific
    # and marker-marginal estimands and applies inverse links marker by marker.
    posterior_values <- .conditional_effect_longitudinal_draws(
      object = object,
      evaluation_data = evaluation_data,
      prediction_method = prediction_method,
      longitudinal_estimand = longitudinal_estimand,
      posterior_draws = posterior_draws,
      random_seed = random_seed
    )

    # Step 5: For the marker-marginal estimand, average marker-specific values
    # within each draw and predictor profile. This order yields an average of
    # expected marker responses on posterior_epred rather than an inverse-link
    # transformation of an averaged marker effect.
    if (identical(longitudinal_estimand, "marginal_marker")) {
      marginal_result <- .conditional_effect_marginalize_markers(
        posterior_matrix = posterior_values,
        evaluation_data = evaluation_data,
        marker_levels = marker_levels
      )
      posterior_values <- marginal_result$draws
      evaluation_data <- marginal_result$data
    }

    # Step 6: Remove calculation-only matching columns, restore stable facet
    # levels, and then summarise the declared posterior estimand.
    internal_columns <- c(".conditional_grid_row__", ".conditional_base_condition__")
    evaluation_data <- evaluation_data[, setdiff(names(evaluation_data), internal_columns), drop = FALSE]
    if (identical(longitudinal_estimand, "marginal_marker") && !marker_is_focal) {
      evaluation_data[[marker_variable]] <- NULL
    }
    evaluation_data$cond__ <- factor(
      evaluation_data$cond__,
      levels = unique(as.character(evaluation_data$cond__))
    )
    summary_data <- .conditional_effect_summary(posterior_values, probability, robust_center)
    table <- cbind(evaluation_data, summary_data)
    table$process__ <- "longitudinal"
    table$longitudinal_estimand__ <- longitudinal_estimand

    # Step 7: Attach the attributes consumed by the brms plotting method. Raw
    # marker observations are available for population and marker-specific
    # panels. No individual observation represents an across-marker mean, so
    # the marker-marginal table receives a typed empty point data frame.
    attr(table, "effects") <- grid_effects
    response_prefix <- switch(
      longitudinal_estimand,
      population = "Population",
      marker = "Marker-specific",
      marginal_marker = "Marker-marginal"
    )
    attr(table, "response") <- if (identical(prediction_method, "posterior_linpred")) {
      paste(response_prefix, "longitudinal linear predictor")
    } else {
      paste(response_prefix, "expected longitudinal response")
    }
    attr(table, "surface") <- isTRUE(draw_surface) && length(grid_effects) == 2L && all(grid_types == "numeric")
    attr(table, "categorical") <- FALSE
    attr(table, "catscale") <- NULL
    attr(table, "ordinal") <- FALSE
    attr(table, "spaghetti") <- NULL
    attr(table, "longitudinal_estimand") <- longitudinal_estimand
    attr(table, "points") <- if (identical(longitudinal_estimand, "marginal_marker")) {
      .conditional_effect_empty_points(table)
    } else {
      .conditional_effect_longitudinal_points(
        object = object,
        effect_table = table,
        condition_rows = condition_rows,
        marker_levels = marker_levels,
        draw_surface = draw_surface
      )
    }
    tables[[effect_name]] <- table
  }
  structure(tables, class = c("brms_conditional_effects", "list"))
}

#' Evaluate event-regression conditional draws
#'
#' @param object Fitted model.
#' @param evaluation_data New event rows.
#' @param event_estimand Event scale.
#' @param posterior_draws Optional draw count.
#' @param random_seed Subsampling seed.
#'
#' @return List of draw matrices, one per event type.
#' @keywords internal
#' @noRd
.conditional_effect_event_draws <- function(object,
                                            evaluation_data,
                                            event_estimand,
                                            posterior_draws,
                                            random_seed) {
  stan_data <- object$stan_data
  event_design <- .mm_event(object$formulaEvent, evaluation_data)
  expected_columns <- as.character(stan_data$w_cols %||% colnames(event_design))

  # Step 1: Align the event design exactly to the fitted coefficient order.
  missing_columns <- setdiff(expected_columns, colnames(event_design))
  for (column in missing_columns) event_design <- cbind(event_design, setNames(data.frame(rep(0, nrow(event_design))), column))
  event_design <- as.matrix(event_design[, expected_columns, drop = FALSE])
  coefficient_count <- as.integer(stan_data$p_w %||% ncol(event_design))
  event_type_count <- as.integer(stan_data$K_event %||% 1L)

  # Step 2: Extract the event-type-specific gamma vectors. Stan arrays are
  # represented by two indices; the one-event fallback accommodates posterior
  # backends that simplify a length-one outer array.
  output <- vector("list", event_type_count)
  for (event_type in seq_len(event_type_count)) {
    matrix_names <- paste0("gamma_w[", event_type, ",", seq_len(coefficient_count), "]")
    gamma_draws <- .conditional_effect_draw_matrix(object, matrix_names, posterior_draws, random_seed)
    if (event_type_count == 1L && ncol(gamma_draws) != coefficient_count) {
      vector_names <- paste0("gamma_w[", seq_len(coefficient_count), "]")
      gamma_draws <- .conditional_effect_draw_matrix(object, vector_names, posterior_draws, random_seed)
    }
    if (!nrow(gamma_draws) || ncol(gamma_draws) != coefficient_count) {
      cli::cli_abort("Could not extract fitted event-regression draws for event type {event_type}.")
    }

    # Step 3: Calculate the direct event covariate contribution. Exponentiation
    # converts the log-relative hazard to its proportional-hazards multiplier.
    event_linear_predictor <- gamma_draws %*% t(event_design)
    output[[event_type]] <- if (identical(event_estimand, "hazard_ratio")) exp(event_linear_predictor) else event_linear_predictor
  }
  output
}

#' Build event conditional-effect tables
#'
#' @keywords internal
#' @noRd
.conditional_effects_event <- function(object,
                                       effect_requests,
                                       process_spec,
                                       condition_rows,
                                       int_conditions,
                                       probability,
                                       robust_center,
                                       event_estimand,
                                       grid_resolution,
                                       draw_surface,
                                       posterior_draws,
                                       random_seed) {
  # Step 1: Retain an explicit empty event component for intercept-only event
  # models. This records that the process was requested even though no event
  # covariate effect exists to vary.
  if (!length(effect_requests)) {
    return(structure(list(), class = c("brms_conditional_effects", "list")))
  }

  event_type_count <- as.integer(object$stan_data$K_event %||% 1L)
  tables <- setNames(vector("list", length(effect_requests)), effect_requests)
  for (effect_name in effect_requests) {
    grid <- .conditional_effect_grid(
      process_spec = process_spec,
      effect_name = effect_name,
      condition_rows = condition_rows,
      int_conditions = int_conditions,
      grid_resolution = grid_resolution,
      draw_surface = draw_surface
    )
    grid_effects <- attr(grid, "effects")
    grid_types <- attr(grid, "types")

    # Step 2: Calculate posterior effects once per event type and bind them with
    # event-type labels. Competing events therefore remain statistically and
    # graphically distinct.
    draws_by_event_type <- .conditional_effect_event_draws(
      object = object,
      evaluation_data = grid,
      event_estimand = event_estimand,
      posterior_draws = posterior_draws,
      random_seed = random_seed
    )
    event_tables <- lapply(seq_len(event_type_count), function(event_type) {
      event_grid <- grid
      event_grid$event_type__ <- event_type
      base_condition <- as.character(event_grid$cond__)
      event_grid$cond__ <- if (event_type_count > 1L) paste0(base_condition, " | event type: ", event_type) else base_condition
      cbind(event_grid, .conditional_effect_summary(draws_by_event_type[[event_type]], probability, robust_center))
    })
    table <- do.call(rbind, event_tables)
    rownames(table) <- NULL
    table$cond__ <- factor(table$cond__, levels = unique(as.character(table$cond__)))
    table$process__ <- "event"

    # Step 3: Attach the standard brms conditional-effects metadata so its plot
    # method can render numeric, categorical, interaction, and surface effects.
    attr(table, "effects") <- grid_effects
    attr(table, "response") <- if (identical(event_estimand, "hazard_ratio")) "Relative event hazard" else "Log relative event hazard"
    attr(table, "surface") <- isTRUE(draw_surface) && length(grid_effects) == 2L && all(grid_types == "numeric")
    attr(table, "categorical") <- FALSE
    attr(table, "catscale") <- NULL
    attr(table, "ordinal") <- FALSE
    attr(table, "spaghetti") <- NULL
    # Event observations are binary/censored outcomes and are not observations
    # of the relative-hazard estimand; retain a typed empty point frame rather
    # than plotting them on an incompatible vertical scale.
    attr(table, "points") <- .conditional_effect_empty_points(table)
    tables[[effect_name]] <- table
  }
  structure(tables, class = c("brms_conditional_effects", "list"))
}
