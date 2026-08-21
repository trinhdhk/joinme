#' Conditional contrasts for fitted JoiNMe models
#'
#' @description
#' `conditional_contrast()` compares two named covariate profiles within every
#' posterior draw and only then summarises the resulting contrast. This ordering
#' retains posterior dependence between the two counterfactual predictions and
#' gives the appropriate uncertainty for their paired difference or hazard
#' ratio.
#'
#' `groupA` and `groupB` describe the covariate values that distinguish the two
#' profiles. Other model predictors are taken from each row of `conditions`,
#' when supplied, or from reference values of the fitting data. Group values
#' take precedence when a variable also appears in `conditions`. A named list is
#' useful when values have different types; a named atomic vector is accepted
#' for concise specifications such as
#' `c(treatment = "active", sex = "female")`.
#'
#' For the longitudinal process, the reported contrast is group A minus group B
#' on the scale selected by `method`. For the event process,
#' `event_scale = "hazard_ratio"` reports the direct covariate hazard ratio of
#' group A relative to group B, whereas `"log_hazard_ratio"` reports its
#' logarithm. As in [conditional_effects.JoiNMeFit()], the event calculation
#' concerns the covariate component of `formulaEvent`; the baseline hazard and
#' longitudinal association contribution are held common between profiles.
#'
#' @param x A fitted object inheriting from `JoiNMeFit`.
#' @param groupA Named atomic vector or named list containing one scalar value
#'   for every covariate fixed to define group A.
#' @param groupB Named atomic vector or named list containing one scalar value
#'   for every covariate fixed to define group B.
#' @param conditions Optional data frame containing the common predictor values
#'   at which the contrast is evaluated. There is one posterior contrast for
#'   every row. A `cond__` column supplies display labels; otherwise informative
#'   row names or numbered labels are used. Values not supplied by either a
#'   group or a condition row are derived from the fitting data in the same way
#'   as [conditional_effects.JoiNMeFit()].
#' @param process Character vector selecting `"longitudinal"`, `"event"`, or
#'   both. `"survival"` is accepted as an alias for `"event"`, and `"both"`
#'   selects both fitted processes.
#' @param prob Probability covered by the equal-tailed posterior uncertainty
#'   interval. The default is `0.95`.
#' @param robust Logical. If `TRUE`, posterior medians define `estimate__`; if
#'   `FALSE`, posterior means are used.
#' @param method Longitudinal posterior scale. `"posterior_epred"` compares
#'   expected responses after applying the marker-specific inverse link;
#'   `"posterior_linpred"` compares longitudinal linear predictors.
#' @param longitudinal_estimand Longitudinal quantity entering the comparison.
#'   `"population"` uses population coefficients only and averages the selected
#'   marker trajectories within each draw. `"marker"` also includes the fitted
#'   marker-level deviation and retains a contrast for each selected marker.
#'   `"marginal_marker"` includes that deviation and averages the selected
#'   marker-specific predictions within each draw. Subject and marker-by-subject
#'   deviations are excluded unless `reuse_fitted_re = TRUE` selects a fitted
#'   subject.
#' @param event_scale Event-process contrast. `"hazard_ratio"` reports
#'   `exp(eta_A - eta_B)` and `"log_hazard_ratio"` reports `eta_A - eta_B`, where
#'   `eta` is the direct event-regression linear predictor.
#' @param markers Optional character vector restricting longitudinal marker
#'   levels. With `longitudinal_estimand = "population"` or
#'   `"marginal_marker"`, these are the marker trajectories averaged within
#'   each posterior draw.
#' @param draws Optional positive integer limiting the posterior draws used in
#'   the paired calculation. `NULL` uses all available draws.
#' @param summary Logical. If `TRUE`, return posterior centres, standard errors,
#'   and credible intervals. If `FALSE`, each process contains a tidy data frame
#'   with one row per posterior draw and estimand. The `.draw` and `.value`
#'   columns identify the retained draw and its value. Draw-level results are
#'   returned without plotting.
#' @param reuse_fitted_re Logical. If `TRUE`, `conditions` must contain a fitted
#'   subject identifier. The paired comparison then includes that subject's
#'   posterior subject and marker-by-subject effects without fitting new random
#'   effects. The identifier cannot be changed between `groupA` and `groupB`.
#' @param seed Integer seed used when posterior draws are subsampled.
#' @param plot Logical. If `TRUE` and `summary = TRUE`, construct and display
#'   contrast plots. If `FALSE`, return the `JoiNMeConditionalContrasts`
#'   object. `summary = FALSE` always returns draw-level results without
#'   plotting.
#' @param ... Additional arguments passed to
#'   [plot.JoiNMeConditionalContrasts()] when `plot = TRUE`.
#'
#' @return With `summary = TRUE` and `plot = FALSE`, a
#'   `JoiNMeConditionalContrasts` object containing
#'   one data frame per requested process. Every table includes `estimate__`,
#'   `se__`, `lower__`, and `upper__`, together with condition, group, marker,
#'   event-type, and estimand labels where applicable. With `summary = FALSE`,
#'   every process is instead a tidy draw-level data frame containing `.draw`,
#'   `.value`, `estimand__`, and the applicable scientific descriptors. With
#'   `summary = TRUE` and `plot = TRUE`, the plotting method invisibly returns
#'   either the named list of `ggplot` objects or an arranged `patchwork`
#'   display.
#'
#' @details
#' Factor levels, transformations, spline bases, time scaling, and coefficient
#' ordering are recovered from the fitted model. The two profile designs are
#' evaluated in one posterior calculation so their columns refer to precisely
#' the same MCMC draws. The contrast is then formed draw by draw. Credible
#' intervals therefore describe the posterior distribution of the contrast,
#' rather than a subtraction of separately summarised intervals.
#'
#' If `conditions` is omitted, numeric predictors use their observed mean and
#' categorical predictors use their first fitted level. Multiple condition rows
#' may describe distinct profiles or a trajectory over a continuous variable
#' such as follow-up time. The plotting method can select such a variable with
#' `condition_variable`.
#'
#' @examples
#' \dontrun{
#' treatment_contrast <- conditional_contrast(
#'   fit,
#'   groupA = c(treatment = "active"),
#'   groupB = c(treatment = "control"),
#'   process = c("longitudinal", "event"),
#'   plot = FALSE
#' )
#' print(treatment_contrast)
#' plot(treatment_contrast, arrange = "row", guides = "collect")
#'
#' time_profiles <- data.frame(
#'   time = seq(0, 24, length.out = 50),
#'   cond__ = paste0("month ", seq(0, 24, length.out = 50))
#' )
#' treatment_over_time <- conditional_contrast(
#'   fit,
#'   groupA = list(treatment = "active", sex = "female"),
#'   groupB = list(treatment = "control", sex = "female"),
#'   conditions = time_profiles,
#'   process = "longitudinal",
#'   method = "posterior_epred",
#'   longitudinal_estimand = "marginal_marker",
#'   plot = FALSE
#' )
#' plot(treatment_over_time, condition_variable = "time")
#' }
#'
#' @seealso [conditional_effects.JoiNMeFit()], [make_conditions()]
#' @export
conditional_contrast <- function(x, ...) {
  UseMethod("conditional_contrast")
}

#' @rdname conditional_contrast
#' @export
conditional_contrast.JoiNMeFit <- function(
    x,
    groupA,
    groupB,
    conditions = NULL,
    process = c("longitudinal", "event"),
    prob = 0.95,
    robust = TRUE,
    method = c("posterior_epred", "posterior_linpred"),
    longitudinal_estimand = c("population", "marker", "marginal_marker"),
    event_scale = c("hazard_ratio", "log_hazard_ratio"),
    markers = NULL,
    draws = NULL,
    summary = TRUE,
    reuse_fitted_re = FALSE,
    seed = 1,
    plot = TRUE,
    ...) {
  # Validate scalar controls before forming either counterfactual profile. This
  # keeps malformed scientific requests separate from later design-matrix
  # errors and ensures both profiles are governed by the same estimand.
  if (!inherits(x, "JoiNMeFit")) {
    cli::cli_abort("{.arg x} must be a fitted {.cls JoiNMeFit} object.")
  }
  process_was_missing <- missing(process) # whether the default may adapt to a longitudinal-only fit
  process <- .conditional_contrast_processes(process) # canonical fitted-process names requested by the analyst
  method <- match.arg(method) # longitudinal posterior scale used before subtraction
  longitudinal_estimand <- match.arg(longitudinal_estimand) # random-effect level entering both longitudinal predictions
  event_scale <- match.arg(event_scale) # survival comparison reported after paired predictor evaluation
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1) {
    cli::cli_abort("{.arg prob} must be a single number strictly between 0 and 1.")
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    cli::cli_abort("{.arg robust} must be either TRUE or FALSE.")
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    cli::cli_abort("{.arg plot} must be either TRUE or FALSE.")
  }
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) {
    cli::cli_abort("{.arg summary} must be either TRUE or FALSE.")
  }
  if (!is.logical(reuse_fitted_re) || length(reuse_fitted_re) != 1L || is.na(reuse_fitted_re)) {
    cli::cli_abort("{.arg reuse_fitted_re} must be either TRUE or FALSE.")
  }
  if (!is.null(draws)) {
    draws <- as.integer(draws) # optional common number of paired posterior draws
    if (length(draws) != 1L || !is.finite(draws) || draws < 1L) {
      cli::cli_abort("{.arg draws} must be NULL or a positive integer.")
    }
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    cli::cli_abort("{.arg seed} must be a single finite number.")
  }

  # A longitudinal-only fit has no event estimand. The untouched default drops
  # that unavailable process, whereas an explicit event request receives a
  # direct explanation rather than returning a prior-only Stan quantity.
  has_survival_process <- .fit_includes_survival(x) # whether observed event data contributed to this fit
  if (!has_survival_process && "event" %in% process) {
    if (isTRUE(process_was_missing)) {
      process <- setdiff(process, "event")
    } else {
      .require_fitted_survival_process(x, "Conditional event contrasts")
    }
  }

  # Recover the covariate vocabulary separately for each fitted process. Event
  # responses and longitudinal grouping expressions are not counterfactual
  # covariates and are consequently excluded.
  process_specs <- list(
    longitudinal = .conditional_effect_process_spec(
      formula = reformulas::nobars(x$formulaLong),
      data = x$dataLong,
      process = "longitudinal"
    )
  ) # formula, data, and reference values needed by each selected process
  if (has_survival_process) {
    process_specs$event <- .conditional_effect_process_spec(
      formula = x$formulaEvent,
      data = x$dataEvent,
      process = "event"
    )
  }
  model_variables <- unique(unlist(lapply(process_specs[process], `[[`, "variables"), use.names = FALSE)) # admissible profile covariates
  id_variable <- .get_call_args(x$call, "id_var", "id") # fitted subject grouping variable, held common across a contrast
  if (isTRUE(reuse_fitted_re)) {
    model_variables <- union(model_variables, id_variable)
    process_specs$longitudinal$variables <- union(process_specs$longitudinal$variables, id_variable)
  }
  marker_variable <- .get_call_args(x$call, "marker_var", "marker") # multivariate outcome index controlled by markers
  group_a <- .conditional_contrast_group(groupA, "groupA") # checked scalar assignments defining profile A
  group_b <- .conditional_contrast_group(groupB, "groupB") # checked scalar assignments defining profile B
  group_variables <- union(names(group_a), names(group_b)) # covariates intervened upon in either profile
  unknown_group_variables <- setdiff(group_variables, model_variables) # assignments unable to alter a selected process
  if (length(unknown_group_variables)) {
    cli::cli_abort(c(
      x = "Unknown or unavailable contrast variables: {paste(unknown_group_variables, collapse = ', ')}.",
      i = "Available predictors for the selected processes: {paste(model_variables, collapse = ', ')}."
    ))
  }
  if (marker_variable %in% group_variables) {
    cli::cli_abort(c(
      x = "The marker index cannot be assigned through {.arg groupA} or {.arg groupB}.",
      i = "Use {.arg markers} to select marker-specific or marker-marginal longitudinal contrasts."
    ))
  }
  if (id_variable %in% group_variables) {
    cli::cli_abort(c(
      x = "The fitted subject identifier cannot differ between {.arg groupA} and {.arg groupB}.",
      i = "Supply {.field {id_variable}} in {.arg conditions} so both profiles reuse the same fitted effects."
    ))
  }
  if (.conditional_contrast_groups_identical(group_a, group_b)) {
    cli::cli_abort("{.arg groupA} and {.arg groupB} must describe different covariate profiles.")
  }

  # Conditions provide the common covariate setting for both groups. The shared
  # conditional-effects validator supplies stable labels and preserves tables
  # created by make_conditions().
  if (isTRUE(reuse_fitted_re) && (is.null(conditions) || !(id_variable %in% names(conditions)))) {
    cli::cli_abort(c(
      x = "{.arg conditions} must contain the fitted subject identifier when {.arg reuse_fitted_re = TRUE}.",
      i = "Add a {.field {id_variable}} column containing one or more identifiers represented during fitting."
    ))
  }
  if (isTRUE(reuse_fitted_re)) .fitted_subject_indices(x, conditions[[id_variable]], id_variable)
  condition_rows <- .conditional_effect_conditions(
    conditions = conditions,
    model_variables = model_variables
  ) # one labelled common covariate profile per requested contrast
  group_a_label <- .conditional_contrast_group_label(group_a) # concise description retained in every output row
  group_b_label <- .conditional_contrast_group_label(group_b) # corresponding description of the reference profile

  # Both process builders concatenate profile A and profile B before extracting
  # draws, thereby guaranteeing paired posterior arithmetic when draws are
  # sampled from the fitted object.
  output <- setNames(vector("list", length(process)), process) # process-indexed contrast tables
  if ("longitudinal" %in% process) {
    output$longitudinal <- .conditional_contrast_longitudinal(
      object = x,
      process_spec = process_specs$longitudinal,
      condition_rows = condition_rows,
      group_a = group_a,
      group_b = group_b,
      group_a_label = group_a_label,
      group_b_label = group_b_label,
      probability = prob,
      robust_center = robust,
      prediction_method = method,
      longitudinal_estimand = longitudinal_estimand,
      marker_selection = markers,
      posterior_draws = draws,
      summarise = summary,
      reuse_fitted_re = reuse_fitted_re,
      random_seed = seed
    )
  }
  if ("event" %in% process) {
    output$event <- .conditional_contrast_event(
      object = x,
      process_spec = process_specs$event,
      condition_rows = condition_rows,
      group_a = group_a,
      group_b = group_b,
      group_a_label = group_a_label,
      group_b_label = group_b_label,
      probability = prob,
      robust_center = robust,
      event_estimand = event_scale,
      posterior_draws = draws,
      summarise = summary,
      random_seed = seed
    )
  }

  # Store the estimand and profile descriptions so a saved result remains
  # interpretable and plottable without consulting the fitted object.
  output <- structure(
    output,
    class = c("JoiNMeConditionalContrasts", "list"),
    process = process,
    probability = prob,
    robust = robust,
    summary = summary,
    method = method,
    longitudinal_estimand = longitudinal_estimand,
    reuse_fitted_re = reuse_fitted_re,
    event_scale = event_scale,
    groupA = group_a,
    groupB = group_b,
    groupA_label = group_a_label,
    groupB_label = group_b_label,
    condition_variables = setdiff(names(condition_rows), "cond__"),
    time_variable = .get_call_args(x$call, "time_var", x$stan_data$time_var %||% "time"),
    call = match.call()
  ) # self-describing posterior contrast result

  if (!isTRUE(summary) || !isTRUE(plot)) return(output)
  invisible(plot(output, plot = TRUE, ...))
}

#' Plot JoiNMe conditional contrasts
#'
#' @description
#' Displays the posterior centre and credible interval for each conditional
#' contrast. A numeric or date-valued condition may be used as a trajectory
#' axis; otherwise labelled conditions are shown as point intervals. Markers and
#' competing event types occupy separate panels.
#'
#' @param x A `JoiNMeConditionalContrasts` object.
#' @param plot Logical. If `TRUE`, draw each requested process.
#' @param ask Logical. Whether to prompt before drawing a subsequent process.
#' @param condition_variable Optional single column from `conditions` to place
#'   on the horizontal axis. If `NULL`, fitted time is preferred when it varies,
#'   followed by another varying numeric condition and then `cond__`.
#' @param arrange Arrangement of multiple process plots. `"separate"` retains
#'   a named list and draws one figure at a time; `"grid"` uses an automatically
#'   sized grid; `"row"` or `"column"` uses a single row or column; and
#'   `"design"` follows the layout supplied through `design`.
#' @param ncol,nrow Optional positive whole numbers giving the grid dimensions
#'   when `arrange = "grid"`.
#' @param design A patchwork design string or area description used when
#'   `arrange = "design"`. Panels follow the process order in `x`.
#' @param widths,heights Optional positive numeric vectors giving relative
#'   column widths and row heights in the arranged display.
#' @param guides How legends are treated across an arranged display. One of
#'   `"keep"`, `"collect"`, or `"auto"`.
#' @param ... Unused and reserved for graphical extensions.
#'
#' @return Invisibly, a named list containing one `ggplot` per process when
#'   `arrange = "separate"`, or one `patchwork` display for every other
#'   arrangement.
#' @export
plot.JoiNMeConditionalContrasts <- function(
    x,
    plot = TRUE,
    ask = FALSE,
    condition_variable = NULL,
    arrange = c("separate", "grid", "row", "column", "design"),
    ncol = NULL,
    nrow = NULL,
    design = NULL,
    widths = NULL,
    heights = NULL,
    guides = c("keep", "collect", "auto"),
    ...) {
  if (!inherits(x, "JoiNMeConditionalContrasts")) {
    cli::cli_abort("{.arg x} must inherit from {.cls JoiNMeConditionalContrasts}.")
  }
  if (identical(attr(x, "summary", exact = TRUE), FALSE)) {
    cli::cli_abort(c(
      x = "Draw-level conditional contrasts cannot be plotted directly.",
      i = "Recalculate with {.arg summary = TRUE} to obtain interval plots."
    ))
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    cli::cli_abort("{.arg plot} must be either TRUE or FALSE.")
  }
  if (!is.logical(ask) || length(ask) != 1L || is.na(ask)) {
    cli::cli_abort("{.arg ask} must be either TRUE or FALSE.")
  }
  arrange <- match.arg(arrange) # requested relationship among longitudinal and event figures
  guides <- match.arg(guides) # requested treatment of repeated legends in an arranged display
  unused_arguments <- names(list(...)) # unrecognised named graphical arguments
  if (length(unused_arguments)) {
    cli::cli_abort("unused argument{?s}: {paste(unused_arguments, collapse = ', ')}")
  }
  if (!is.null(condition_variable) &&
      (!is.character(condition_variable) || length(condition_variable) != 1L || !nzchar(condition_variable))) {
    cli::cli_abort("{.arg condition_variable} must be NULL or one column name from {.arg conditions}.")
  }

  # Build figures independently so longitudinal differences and event hazard
  # ratios retain their distinct null values and vertical scales.
  plots <- setNames(vector("list", length(x)), names(x)) # process-indexed ggplot objects returned invisibly
  for (process_index in seq_along(x)) {
    process_name <- names(x)[[process_index]] # submodel represented by this contrast table
    plots[[process_name]] <- .conditional_contrast_plot(
      contrast_table = x[[process_name]],
      process_name = process_name,
      condition_variable = condition_variable,
      preferred_time_variable = attr(x, "time_variable", exact = TRUE),
      condition_variables = attr(x, "condition_variables", exact = TRUE) %||% character(0),
      event_estimand = attr(x, "event_scale", exact = TRUE)
    )
  }
  if (!identical(arrange, "separate")) {
    arranged_plot <- .arrange_conditional_plots(
      plots = plots,
      arrange = arrange,
      ncol = ncol,
      nrow = nrow,
      design = design,
      widths = widths,
      heights = heights,
      guides = guides
    ) # one display retaining the distinct scales of its process-specific panels
    if (isTRUE(plot)) print(arranged_plot)
    return(invisible(arranged_plot))
  }

  # Separate drawing retains the established prompt behaviour. Layout
  # arguments are rejected here because they otherwise have no visible effect.
  .validate_separate_conditional_layout(ncol, nrow, design, widths, heights)
  if (isTRUE(plot)) {
    for (plot_index in seq_along(plots)) {
      if (plot_index > 1L && isTRUE(ask) && interactive()) {
        readline("Press <Return> to display the next conditional-contrast plot: ")
      }
      print(plots[[plot_index]])
    }
  }
  invisible(plots)
}

#' Print a JoiNMe conditional-contrast object
#'
#' @param x A `JoiNMeConditionalContrasts` object.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.JoiNMeConditionalContrasts <- function(x, ...) {
  .cli_summary_heading("Conditional contrasts for JoiNMe model\n")
  draw_level_output <- identical(attr(x, "summary", exact = TRUE), FALSE) # whether terminal process entries contain raw sample matrices
  if (draw_level_output) {
    cli::cli_bullets(c("*" = "Output: draw-level MCMC samples"))
  }
  group_a_label <- attr(x, "groupA_label", exact = TRUE) %||% "group A" # displayed numerator or minuend profile
  group_b_label <- attr(x, "groupB_label", exact = TRUE) %||% "group B" # displayed denominator or subtrahend profile
  cli::cli_bullets(c("*" = sprintf("Comparison: [%s] versus [%s]", group_a_label, group_b_label)))
  if ("longitudinal" %in% names(x)) {
    cli::cli_bullets(c("*" = sprintf(
      "Longitudinal contrast: A - B; estimand: %s; method: %s",
      attr(x, "longitudinal_estimand", exact = TRUE),
      attr(x, "method", exact = TRUE)
    )))
  }
  if ("event" %in% names(x)) {
    event_description <- if (identical(attr(x, "event_scale", exact = TRUE), "hazard_ratio")) {
      "hazard ratio A / B"
    } else {
      "log hazard ratio A - B"
    } # direction and scale of the direct event-regression comparison
    cli::cli_bullets(c("*" = sprintf("Event contrast: %s", event_description)))
  }
  for (process_name in names(x)) {
    process_data <- if (draw_level_output) x[[process_name]]$data else x[[process_name]] # estimand description table for the process
    condition_count <- length(unique(process_data$condition_index__)) # common profiles counted once across markers or event types
    cli::cli_bullets(c("*" = sprintf(
      "%s process: %d condition%s",
      tools::toTitleCase(process_name),
      condition_count,
      if (condition_count == 1L) "" else "s"
    )))
  }
  invisible(x)
}

#' Canonicalise requested conditional-contrast processes
#'
#' @param process User process selection.
#'
#' @return Unique character vector containing `longitudinal` and/or `event`.
#' @keywords internal
#' @noRd
.conditional_contrast_processes <- function(process) {
  process <- as.character(process) # user-facing names before aliases are resolved
  if (!length(process) || anyNA(process) || any(!nzchar(process))) {
    cli::cli_abort("{.arg process} must select {.val longitudinal}, {.val event}, {.val survival}, or {.val both}.")
  }
  invalid_processes <- setdiff(process, c("both", "longitudinal", "event", "survival")) # entries outside the documented vocabulary
  if (length(invalid_processes)) {
    cli::cli_abort("Unknown {.arg process}: {paste(invalid_processes, collapse = ', ')}.")
  }
  process[process == "survival"] <- "event"
  if ("both" %in% process) process <- c("longitudinal", "event")
  unique(process)
}

#' Validate one named conditional-contrast profile
#'
#' @param group Named vector or list of scalar covariate assignments.
#' @param argument_name Argument name used in diagnostic messages.
#'
#' @return A named list retaining supplied scalar value types.
#' @keywords internal
#' @noRd
.conditional_contrast_group <- function(group, argument_name) {
  if (!(is.atomic(group) || is.list(group)) || is.data.frame(group) || !length(group)) {
    cli::cli_abort("{.arg {argument_name}} must be a non-empty named vector or named list.")
  }
  group_names <- names(group) # predictors assigned by this counterfactual profile
  if (is.null(group_names) || anyNA(group_names) || any(!nzchar(group_names))) {
    cli::cli_abort("Every value in {.arg {argument_name}} must have a predictor name.")
  }
  if (anyDuplicated(group_names)) {
    duplicated_names <- unique(group_names[duplicated(group_names)]) # ambiguous repeated assignments in one profile
    cli::cli_abort(sprintf(
      "Duplicated variables in `%s`: %s.",
      argument_name,
      paste(duplicated_names, collapse = ", ")
    ))
  }
  group_list <- as.list(group) # type-preserving scalar assignments used to construct fitted designs
  invalid_values <- vapply(group_list, function(value) {
    length(value) != 1L || is.list(value) || is.na(value)
  }, logical(1))
  if (any(invalid_values)) {
    cli::cli_abort("Every value in {.arg {argument_name}} must be a single non-missing value.")
  }
  group_list
}

#' Determine whether two supplied contrast profiles are identical
#'
#' @param group_a Checked group-A assignments.
#' @param group_b Checked group-B assignments.
#'
#' @return A single logical value.
#' @keywords internal
#' @noRd
.conditional_contrast_groups_identical <- function(group_a, group_b) {
  if (!setequal(names(group_a), names(group_b))) return(FALSE)
  all(vapply(sort(names(group_a)), function(variable) {
    identical(as.character(group_a[[variable]]), as.character(group_b[[variable]]))
  }, logical(1)))
}

#' Format a counterfactual profile for tables and printed output
#'
#' @param group Checked named group assignments.
#'
#' @return A single compact character label.
#' @keywords internal
#' @noRd
.conditional_contrast_group_label <- function(group) {
  assignments <- vapply(names(group), function(variable) {
    sprintf("%s=%s", variable, paste(as.character(group[[variable]]), collapse = ","))
  }, character(1)) # ordered variable-value descriptions retaining the supplied profile order
  paste(assignments, collapse = ", ")
}

#' Construct one counterfactual model-data profile per condition
#'
#' @param process_spec Process formula, fitting data, and predictor names.
#' @param condition_rows Labelled common conditioning profiles.
#' @param group Checked group-specific assignments.
#'
#' @return Data frame ready for a fitted process design matrix.
#' @keywords internal
#' @noRd
.conditional_contrast_evaluation_data <- function(process_spec, condition_rows, group) {
  # Work with an ordinary data frame even when a fitted object contains a
  # matrix-like copy of its analysis data. Matrix `$<-` assignment coerces the
  # object to an unnamed list; the ensuing model matrix then sees missing or
  # zero-length predictors. This explicit boundary preserves column names and
  # permits logical contrasts such as `c(arm = TRUE)`.
  fitted_data <- as.data.frame(process_spec$data, stringsAsFactors = FALSE) # named fitted columns defining types, levels, and reference values
  evaluation_data <- fitted_data[rep(1L, nrow(condition_rows)), process_spec$variables, drop = FALSE] # correctly typed shell for one row per condition
  for (variable in process_spec$variables) {
    observed <- fitted_data[[variable]] # fitting-data column defining levels and storage type
    reference_value <- .conditional_effect_reference(observed) # empirical reference when no group or condition value is supplied
    if (variable %in% names(group)) {
      requested_values <- rep(group[[variable]], nrow(condition_rows)) # intervention shared across all common conditions
    } else if (variable %in% names(condition_rows)) {
      requested_values <- condition_rows[[variable]] # condition-specific common values for both groups
      missing_values <- is.na(requested_values) # entries that defer to the training-data reference
      if (any(missing_values)) requested_values[missing_values] <- reference_value
    } else {
      requested_values <- rep(reference_value, nrow(condition_rows)) # reference profile derived from fitting data
    }
    evaluation_data[[variable]] <- .conditional_effect_cast(requested_values, observed)
    if (anyNA(evaluation_data[[variable]])) {
      cli::cli_abort(c(
        x = "The requested value for {.val {variable}} is incompatible with its fitted levels or data type.",
        i = "Use a value represented in the model-fitting data."
      ))
    }
  }
  evaluation_data$cond__ <- condition_rows$cond__ # stable condition label retained in the result
  evaluation_data$condition_index__ <- seq_len(nrow(condition_rows)) # pairing key shared by A and B designs
  evaluation_data
}

#' Retain condition and group columns in a contrast result
#'
#' @param condition_rows Labelled common conditioning profiles.
#' @param group_variables Covariates overridden by either group.
#' @param group_a_label Display label for group A.
#' @param group_b_label Display label for group B.
#'
#' @return One descriptive row per supplied condition.
#' @keywords internal
#' @noRd
.conditional_contrast_output_rows <- function(
    condition_rows,
    group_variables,
    group_a_label,
    group_b_label) {
  retained_columns <- setdiff(names(condition_rows), group_variables) # common conditions not superseded by group assignments
  output_rows <- condition_rows[, retained_columns, drop = FALSE] # analyst-supplied common profile values
  output_rows$condition_index__ <- seq_len(nrow(output_rows)) # explicit pairing index for verbose condition labels
  output_rows$groupA__ <- group_a_label # profile forming the minuend or ratio numerator
  output_rows$groupB__ <- group_b_label # profile forming the subtrahend or ratio denominator
  output_rows
}

#' Calculate paired longitudinal posterior contrasts
#'
#' @param summarise Whether to replace the paired contrast draws by posterior
#'   summaries.
#'
#' @return A summary data frame when `summarise` is true, otherwise a tidy
#'   draw-level data frame.
#' @keywords internal
#' @noRd
.conditional_contrast_longitudinal <- function(
    object,
    process_spec,
    condition_rows,
    group_a,
    group_b,
    group_a_label,
    group_b_label,
    probability,
    robust_center,
    prediction_method,
    longitudinal_estimand,
    marker_selection,
    posterior_draws,
    summarise = TRUE,
    reuse_fitted_re = FALSE,
    random_seed = 1) {
  marker_variable <- .get_call_args(object$call, "marker_var", "marker") # multivariate longitudinal outcome index
  marker_levels <- as.character(
    object$stan_data$marker_levels %||%
      levels(object$dataLong[[marker_variable]]) %||%
      unique(as.character(object$dataLong[[marker_variable]]))
  ) # fitted marker strata available for evaluation and finite marginalisation
  if (!is.null(marker_selection)) {
    unknown_markers <- setdiff(as.character(marker_selection), marker_levels) # requested outcomes absent from the fitted process
    if (length(unknown_markers)) {
      cli::cli_abort("Unknown marker levels: {paste(unknown_markers, collapse = ', ')}.")
    }
    marker_levels <- marker_levels[marker_levels %in% as.character(marker_selection)]
  }
  if (!length(marker_levels)) {
    cli::cli_abort("At least one fitted marker must be retained for a longitudinal contrast.")
  }

  # Construct both covariate profiles and repeat each over exactly the same
  # ordered marker strata. Their common order gives an unambiguous column
  # pairing after posterior prediction.
  profile_a <- .conditional_contrast_evaluation_data(process_spec, condition_rows, group_a) # group-A fixed-effect rows before marker expansion
  profile_b <- .conditional_contrast_evaluation_data(process_spec, condition_rows, group_b) # group-B fixed-effect rows before marker expansion
  expand_over_markers <- function(profile) {
    expanded_profiles <- lapply(marker_levels, function(marker_level) {
      marker_profile <- profile # all common conditions for one fitted marker
      marker_profile[[marker_variable]] <- factor(marker_level, levels = as.character(object$stan_data$marker_levels %||% marker_levels))
      marker_profile$marker__ <- marker_level # explicit outcome label retained through calculations
      marker_profile
    })
    expanded <- do.call(rbind, expanded_profiles) # marker-major order shared by A and B
    rownames(expanded) <- NULL
    expanded
  }
  evaluation_a <- expand_over_markers(profile_a) # A rows for all condition-marker combinations
  evaluation_b <- expand_over_markers(profile_b) # paired B rows in precisely the same order
  paired_row_count <- nrow(evaluation_a) # posterior columns belonging to either profile

  # Both groups pass through a single draw extraction and matrix multiplication.
  # A limited draw request therefore cannot independently sample iterations for
  # A and B, and their posterior covariance is preserved exactly.
  paired_evaluation_data <- rbind(evaluation_a, evaluation_b) # concatenated A then B evaluation rows
  paired_values <- .conditional_effect_longitudinal_draws(
    object = object,
    evaluation_data = paired_evaluation_data,
    prediction_method = prediction_method,
    longitudinal_estimand = longitudinal_estimand,
    posterior_draws = posterior_draws,
    reuse_fitted_re = reuse_fitted_re,
    random_seed = random_seed
  ) # posterior draws by concatenated evaluation row
  group_a_values <- paired_values[, seq_len(paired_row_count), drop = FALSE] # draw-aligned predictions under A
  group_b_values <- paired_values[, paired_row_count + seq_len(paired_row_count), drop = FALSE] # draw-aligned predictions under B
  contrast_values <- group_a_values - group_b_values # A minus B before posterior summary

  # Population and marker-marginal contrasts average each draw over selected
  # outcomes, leaving one contrast per condition. The former excludes marker
  # deviations and the latter includes them. Only the marker estimand retains a
  # separate contrast for each selected marker.
  output_base <- .conditional_contrast_output_rows(
    condition_rows = condition_rows,
    group_variables = union(names(group_a), names(group_b)),
    group_a_label = group_a_label,
    group_b_label = group_b_label
  ) # condition and group descriptors before marker replication
  marginalise_marker_trajectories <- longitudinal_estimand %in% c("population", "marginal_marker") # whether selected marker contrasts collapse within each posterior draw
  if (isTRUE(marginalise_marker_trajectories)) {
    condition_indices <- evaluation_a$condition_index__ # repeated condition keys over selected markers
    contrast_values <- do.call(cbind, lapply(seq_len(nrow(condition_rows)), function(condition_index) {
      rowMeans(contrast_values[, condition_indices == condition_index, drop = FALSE])
    })) # equally weighted draw-level marker average of A-minus-B
    output_rows <- output_base # one marginal contrast per common condition
    output_rows$marker__ <- "Marginal over markers"
  } else {
    output_rows <- output_base[evaluation_a$condition_index__, , drop = FALSE] # descriptions repeated in posterior-column order
    output_rows$marker__ <- evaluation_a$marker__
    rownames(output_rows) <- NULL
  }

  output_rows$process__ <- "longitudinal"
  output_rows$longitudinal_estimand__ <- longitudinal_estimand
  output_rows$method__ <- prediction_method
  output_rows$contrast__ <- "groupA - groupB"
  response_label <- if (identical(prediction_method, "posterior_linpred")) {
    "Longitudinal linear-predictor contrast (A - B)"
  } else {
    "Expected longitudinal-response contrast (A - B)"
  } # human-readable quantity shared by summarised and draw-level results
  if (!isTRUE(summarise)) {
    return(.conditional_mcmc_samples(
      posterior_matrix = contrast_values,
      estimand_data = output_rows,
      response_label = response_label
    ))
  }
  summary_data <- .conditional_effect_summary(contrast_values, probability, robust_center) # paired posterior centre, spread, and interval
  table <- cbind(output_rows, summary_data) # self-describing longitudinal contrast table
  attr(table, "response") <- response_label
  table
}

#' Calculate paired event-regression posterior contrasts
#'
#' @param summarise Whether to replace the paired event contrast draws by
#'   posterior summaries.
#'
#' @return A summary data frame when `summarise` is true, otherwise a tidy
#'   draw-level data frame.
#' @keywords internal
#' @noRd
.conditional_contrast_event <- function(
    object,
    process_spec,
    condition_rows,
    group_a,
    group_b,
    group_a_label,
    group_b_label,
    probability,
    robust_center,
    event_estimand,
    posterior_draws,
    summarise = TRUE,
    random_seed) {
  profile_a <- .conditional_contrast_evaluation_data(process_spec, condition_rows, group_a) # event rows under A
  profile_b <- .conditional_contrast_evaluation_data(process_spec, condition_rows, group_b) # event rows under B
  condition_count <- nrow(condition_rows) # posterior columns belonging to either profile
  paired_profiles <- rbind(profile_a, profile_b) # A then B evaluated against one coefficient sample
  event_type_count <- as.integer(object$stan_data$K_event %||% 1L) # fitted event causes requiring separate contrasts
  coefficient_count <- as.integer(object$stan_data$p_w %||% 0L) # direct event coefficients in formulaEvent

  # Without direct event covariates the shared baseline cancels and the contrast
  # is deterministically null. Otherwise the two linear predictors are evaluated
  # together and compared within every posterior draw.
  if (coefficient_count > 0L) {
    linear_predictors <- .conditional_effect_event_draws(
      object = object,
      evaluation_data = paired_profiles,
      event_estimand = "log_hazard_ratio",
      posterior_draws = posterior_draws,
      random_seed = random_seed
    ) # draw-by-profile direct event predictors for each cause
  } else {
    retained_draw_count <- .conditional_effect_draw_count(
      object = object,
      posterior_draws = posterior_draws,
      random_seed = random_seed
    ) # posterior rows replicated by the deterministic shared-baseline contrast
    linear_predictors <- rep(
      list(matrix(0, nrow = retained_draw_count, ncol = 2L * condition_count)),
      event_type_count
    ) # deterministic shared-baseline result for an intercept-only event model
  }

  output_base <- .conditional_contrast_output_rows(
    condition_rows = condition_rows,
    group_variables = union(names(group_a), names(group_b)),
    group_a_label = group_a_label,
    group_b_label = group_b_label
  ) # one descriptive row per common condition
  event_values <- vector("list", event_type_count) # cause-specific draw matrices retained before optional summary
  event_rows <- vector("list", event_type_count) # cause-specific estimand descriptions in matching order
  for (event_type in seq_len(event_type_count)) {
    event_draws <- linear_predictors[[event_type]] # paired predictors for one event cause
    log_contrast <- event_draws[, seq_len(condition_count), drop = FALSE] -
      event_draws[, condition_count + seq_len(condition_count), drop = FALSE] # eta_A minus eta_B by draw
    contrast_values <- if (identical(event_estimand, "hazard_ratio")) exp(log_contrast) else log_contrast # requested comparison before summary
    event_values[[event_type]] <- contrast_values
    cause_rows <- output_base # common conditions for this event cause
    cause_rows$event_type__ <- event_type
    event_rows[[event_type]] <- cause_rows
  }
  estimand_data <- do.call(rbind, event_rows) # vertically combined cause-specific descriptions
  rownames(estimand_data) <- NULL
  estimand_data$process__ <- "event"
  estimand_data$event_scale__ <- event_estimand
  estimand_data$contrast__ <- if (identical(event_estimand, "hazard_ratio")) "groupA / groupB" else "groupA - groupB"
  response_label <- if (identical(event_estimand, "hazard_ratio")) {
    "Direct event hazard ratio (A / B)"
  } else {
    "Direct event log hazard ratio (A - B)"
  } # cause-specific event comparison represented by every sample column
  if (!isTRUE(summarise)) {
    return(.conditional_mcmc_samples(
      posterior_matrix = do.call(cbind, event_values),
      estimand_data = estimand_data,
      response_label = response_label
    ))
  }
  event_tables <- lapply(seq_len(event_type_count), function(event_type) {
    cbind(event_rows[[event_type]], .conditional_effect_summary(event_values[[event_type]], probability, robust_center))
  })
  table <- do.call(rbind, event_tables) # vertically combined cause-specific contrasts
  rownames(table) <- NULL
  table$process__ <- "event"
  table$event_scale__ <- event_estimand
  table$contrast__ <- if (identical(event_estimand, "hazard_ratio")) "groupA / groupB" else "groupA - groupB"
  attr(table, "response") <- response_label
  table
}

#' Construct one conditional-contrast figure
#'
#' @param contrast_table Summarised process-specific contrast table.
#' @param process_name Canonical process name.
#' @param condition_variable Optional requested horizontal-axis variable.
#' @param preferred_time_variable Fitted longitudinal time variable.
#' @param condition_variables Columns originating in `conditions`.
#' @param event_estimand Event contrast scale used to select its null value.
#'
#' @return A `ggplot` object.
#' @keywords internal
#' @noRd
.conditional_contrast_plot <- function(
    contrast_table,
    process_name,
    condition_variable,
    preferred_time_variable,
    condition_variables,
    event_estimand) {
  available_condition_variables <- intersect(condition_variables, names(contrast_table)) # retained common condition columns
  varying_condition_variables <- available_condition_variables[vapply(available_condition_variables, function(variable) {
    length(unique(contrast_table[[variable]])) > 1L
  }, logical(1))] # conditions capable of defining an axis or concurrent profiles

  # Prefer an explicit axis, then fitted time, then another varying numeric
  # condition. Labelled conditions are a safe discrete fallback.
  if (!is.null(condition_variable)) {
    if (!(condition_variable %in% names(contrast_table))) {
      cli::cli_abort("Unknown {.arg condition_variable} {.val {condition_variable}} for the {process_name} contrast table.")
    }
    axis_variable <- condition_variable
  } else if (!is.null(preferred_time_variable) && preferred_time_variable %in% varying_condition_variables) {
    axis_variable <- preferred_time_variable
  } else {
    numeric_candidates <- varying_condition_variables[vapply(varying_condition_variables, function(variable) {
      is.numeric(contrast_table[[variable]]) || inherits(contrast_table[[variable]], c("Date", "POSIXt"))
    }, logical(1))]
    axis_variable <- if (length(numeric_candidates)) numeric_candidates[[1L]] else "cond__"
  } # column represented on the horizontal axis

  plot_data <- contrast_table # local graphical data augmented without changing the statistical result
  continuous_axis <- is.numeric(plot_data[[axis_variable]]) || inherits(plot_data[[axis_variable]], c("Date", "POSIXt")) # whether trajectories are meaningful
  profile_variables <- if (continuous_axis) setdiff(varying_condition_variables, axis_variable) else character(0) # non-axis conditions distinguishing curves
  if (length(profile_variables)) {
    plot_data$profile__ <- apply(
      as.data.frame(lapply(plot_data[, profile_variables, drop = FALSE], as.character)),
      1L,
      function(values) paste(sprintf("%s=%s", profile_variables, values), collapse = ", ")
    ) # readable label assembled from the remaining varying conditions
  } else {
    plot_data$profile__ <- "Contrast" # one common profile requiring no colour legend
  }

  # Marker and event-type panels prevent quantities representing distinct
  # outcomes or causes from being overlaid. Remaining condition combinations
  # are distinguished by line colour within each panel.
  panel_variable <- if (identical(process_name, "longitudinal") &&
                        "marker__" %in% names(plot_data) &&
                        length(unique(plot_data$marker__)) > 1L) {
    "marker__"
  } else if (identical(process_name, "event") &&
             "event_type__" %in% names(plot_data) &&
             length(unique(plot_data$event_type__)) > 1L) {
    "event_type__"
  } else {
    NULL
  } # optional scientific stratum assigned to separate facets
  reference_value <- if (identical(process_name, "event") && identical(event_estimand, "hazard_ratio")) 1 else 0 # null on the displayed scale
  response_label <- attr(contrast_table, "response", exact = TRUE) %||% "Posterior contrast" # vertical-axis description

  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[axis_variable]],
      y = .data$estimate__,
      colour = .data$profile__,
      fill = .data$profile__,
      group = .data$profile__
    )
  ) +
    ggplot2::geom_hline(yintercept = reference_value, linetype = 2L, colour = "grey45")
  panel_membership <- if (is.null(panel_variable)) rep("all", nrow(plot_data)) else plot_data[[panel_variable]] # panel label paired with every plotted condition
  profile_counts <- table(plot_data$profile__, panel_membership) # support for each prospective curve
  may_draw_curves <- continuous_axis && length(unique(plot_data[[axis_variable]])) > 1L && all(profile_counts >= 2L) # enough ordered support for lines and ribbons
  if (isTRUE(may_draw_curves)) {
    figure <- figure +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower__, ymax = .data$upper__),
        alpha = 0.18,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.85) +
      ggplot2::geom_point(size = 1.4)
  } else {
    figure <- figure +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower__, ymax = .data$upper__),
        width = 0.15
      ) +
      ggplot2::geom_point(size = 2.2)
  }
  if (!is.null(panel_variable)) {
    figure <- figure + ggplot2::facet_wrap(stats::as.formula(paste("~", panel_variable)), scales = "free_y")
  }
  if (length(unique(plot_data$profile__)) == 1L) {
    figure <- figure + ggplot2::guides(colour = "none", fill = "none")
  }
  figure +
    ggplot2::labs(
      title = sprintf("%s conditional contrast", tools::toTitleCase(process_name)),
      subtitle = sprintf("[%s] versus [%s]", unique(plot_data$groupA__)[[1L]], unique(plot_data$groupB__)[[1L]]),
      x = if (identical(axis_variable, "cond__")) "Condition" else axis_variable,
      y = response_label,
      colour = "Condition profile",
      fill = "Condition profile"
    ) +
    ggplot2::theme_bw()
}
