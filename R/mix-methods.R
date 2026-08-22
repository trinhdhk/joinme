#' Posterior class membership for a latent-progress fit
#'
#' @description
#' Returns posterior class probabilities for the natural allocation domains of
#' a [joinme_mix()] model. Subject, correlation and covariance classes are
#' reported by subject; marker classes are reported by marker.
#'
#' @param object A `JoiNMeMixFit` object.
#' @param draws Optional number of posterior draws.
#' @param seed Seed used for reproducible posterior subsetting.
#' @param digits Number of decimal places used for posterior summaries.
#' @param summary If `TRUE`, return the posterior mean (`Estimate`), posterior
#'   standard deviation (`Est.Error`), interval bounds, available sampling
#'   diagnostics, and the maximum-probability class. If `FALSE`, return draw
#'   matrices.
#'
#' @return A named list with `subject` and/or `marker` entries.
#' @export
posterior_class <- function(object, ...) {
  UseMethod("posterior_class")
}

#' @rdname posterior_class
#' @export
posterior_class.JoiNMeMixFit <- function(
  object,
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE
) {
  .stop_experimental("posterior_class()")
  mixture <- object$mixture %||% object$config$mixture %||% list()
  number_classes <- as.integer(mixture$n_classes %||% 0L)
  if (number_classes < 2L) {
    cli::cli_abort("The fitted mixture metadata is unavailable.")
  }
  if (
    length(digits) != 1L ||
      !is.finite(digits) ||
      digits < 0
  ) {
    cli::cli_abort("{.arg digits} must be a non-negative number.")
  }

  allocation_domains <- mixture$allocation_domains %||% list()
  output <- list()
  if (length(allocation_domains$subject %||% character(0)) > 0L) {
    output$subject <- .mixture_membership_(
      object = object,
      domain = "subject",
      number_units = as.integer(object$stan_data$n_id),
      unit_labels = .id_labels(object, as.integer(object$stan_data$n_id)),
      number_classes = number_classes,
      draws = draws,
      seed = seed,
      digits = digits,
      summary = summary
    )
  }
  if (length(allocation_domains$marker %||% character(0)) > 0L) {
    output$marker <- .mixture_membership_(
      object = object,
      domain = "marker",
      number_units = as.integer(object$stan_data$D),
      unit_labels = as.character(
        object$stan_data$marker_levels %||%
          seq_len(as.integer(object$stan_data$D))
      ),
      number_classes = number_classes,
      draws = draws,
      seed = seed,
      digits = digits,
      summary = summary
    )
  }
  available_draws <- tryCatch(
    posterior::ndraws(.get_draws_obj(object$fit)),
    error = function(error) NA_integer_
  ) # total post-warm-up draws stored by the fitted backend
  retained_draws <- if (
    is.null(draws) || !is.finite(draws)
  ) {
    available_draws
  } else {
    min(as.integer(draws), available_draws, na.rm = TRUE)
  } # posterior draws represented in the returned probabilities
  structure(
    output,
    class = c("JoiNMePosteriorClass", "list"),
    mixture = mixture,
    summary = isTRUE(summary),
    metadata = list(
      summary = isTRUE(summary),
      draws = retained_draws,
      digits = as.integer(digits),
      n_classes = number_classes,
      class_type = mixture$class_type %||% character(0),
      source = "fitted model"
    )
  )
}

#' Conditional class membership after dynamic prediction
#'
#' @description
#' Summarises the class probabilities obtained after conditioning a new
#' subject's random effects on their observed longitudinal history.  For a
#' combined subject and covariance-regression mixture, the probabilities refer
#' to their one shared allocation.
#'
#' @inheritParams posterior_class
#' @param object A `JoiNMeMixDynPred` object returned by
#'   `predict.JoiNMeMixFit()`.
#'
#' @return A `JoiNMePosteriorClass` list.  With `summary = TRUE`, the
#'   `subject` and `marker` entries are posterior summary tables.  With
#'   `summary = FALSE`, the subject-indexed matrices and arrays are returned.
#' @export
posterior_class.JoiNMeMixDynPred <- function(
  object,
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE
) {
  .stop_experimental("posterior_class()")
  conditional_draws <- object$draws$posterior_class %||% list()
  if (length(conditional_draws) == 0L) {
    cli::cli_abort(c(
      x = "Conditional class probabilities are unavailable.",
      i = "The parent mixture may contain no dynamically sampled class-specific block."
    ))
  }

  number_classes <- as.integer(
    (object$mixture %||% object$metadata$mixture)$n_classes %||% 0L
  )
  if (number_classes < 2L) {
    cli::cli_abort("The prediction's mixture metadata is unavailable.")
  }
  if (
    length(digits) != 1L ||
      !is.finite(digits) ||
      digits < 0
  ) {
    cli::cli_abort("{.arg digits} must be a non-negative number.")
  }

  subset_first_dimension <- function(values, subject_offset) {
    number_available <- dim(values)[1L]
    if (is.null(draws) || as.integer(draws) >= number_available) {
      return(values)
    }
    requested <- as.integer(draws)
    if (
      length(requested) != 1L ||
      !is.finite(requested) ||
      requested < 1L
    ) {
      cli::cli_abort("{.arg draws} must be a positive integer.")
    }
    has_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    if (has_seed) {
      previous_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(
        assign(".Random.seed", previous_seed, envir = .GlobalEnv),
        add = TRUE
      )
    } else {
      on.exit(
        rm(".Random.seed", envir = .GlobalEnv),
        add = TRUE
      )
    }
    set.seed(.mixture_safe_seed(seed, subject_offset))
    selected <- sample.int(number_available, requested)
    dimensions <- dim(values)
    indices <- c(
      list(selected),
      rep(list(TRUE), length(dimensions) - 1L),
      list(drop = FALSE)
    )
    do.call(`[`, c(list(values), indices))
  }

  subject_names <- names(conditional_draws)
  if (is.null(subject_names)) {
    subject_names <- as.character(seq_along(conditional_draws))
  }
  selected_draws <- conditional_draws
  for (subject_index in seq_along(selected_draws)) {
    for (domain in names(selected_draws[[subject_index]])) {
      selected_draws[[subject_index]][[domain]] <-
        subset_first_dimension(
          selected_draws[[subject_index]][[domain]],
          subject_index
        )
    }
  }
  if (!isTRUE(summary)) {
    raw_output <- list()
    raw_subject <- lapply(selected_draws, `[[`, "subject")
    names(raw_subject) <- subject_names
    raw_subject <- Filter(Negate(is.null), raw_subject)
    if (length(raw_subject) > 0L) {
      raw_output$subject <- raw_subject
    }
    raw_marker <- lapply(selected_draws, `[[`, "marker")
    names(raw_marker) <- subject_names
    raw_marker <- Filter(Negate(is.null), raw_marker)
    if (length(raw_marker) > 0L) {
      raw_output$marker <- raw_marker
    }
    represented_draws <- NA_integer_ # retained conditional draws in the raw result
    if (length(selected_draws) > 0L) {
      first_subject_draws <- Filter(
        Negate(is.null),
        unname(selected_draws[[1L]])
      ) # available allocation-domain arrays for the first predicted subject
      if (length(first_subject_draws) > 0L) {
        represented_draws <- dim(first_subject_draws[[1L]])[1L]
      }
    }
    return(structure(
      raw_output,
      class = c("JoiNMePosteriorClass", "list"),
      mixture = object$mixture %||% object$metadata$mixture,
      summary = FALSE,
      metadata = list(
        summary = FALSE,
        draws = represented_draws,
        digits = as.integer(digits),
        n_classes = number_classes,
        class_type = (
          object$mixture %||% object$metadata$mixture
        )$class_type %||% character(0),
        source = "dynamic prediction"
      )
    ))
  }

  subject_rows <- list()
  marker_rows <- list()
  marker_labels <- as.character(
    object$metadata$marker_levels %||% character(0)
  )
  for (subject_index in seq_along(selected_draws)) {
    subject_name <- subject_names[subject_index]
    subject_probability <-
      selected_draws[[subject_index]]$subject
    if (!is.null(subject_probability)) {
      posterior_mean <- colMeans(subject_probability, na.rm = TRUE)
      assigned_class <- which.max(posterior_mean)
      for (group in seq_len(number_classes)) {
        values <- subject_probability[, group]
        bounds <- stats::quantile(
          values,
          c(0.025, 0.975),
          names = FALSE,
          na.rm = TRUE
        )
        subject_rows[[length(subject_rows) + 1L]] <- data.frame(
          unit = subject_name,
          class = paste0("class_", group),
          Estimate = mean(values, na.rm = TRUE),
          Est.Error = stats::sd(values, na.rm = TRUE),
          Q2.5 = bounds[1L],
          Q97.5 = bounds[2L],
          Rhat = NA_real_,
          ess_bulk = NA_real_,
          ess_tail = NA_real_,
          assigned_class = paste0("class_", assigned_class),
          stringsAsFactors = FALSE
        )
      }
    }

    marker_probability <- selected_draws[[subject_index]]$marker
    if (!is.null(marker_probability)) {
      number_markers <- dim(marker_probability)[2L]
      if (length(marker_labels) != number_markers) {
        marker_labels <- as.character(seq_len(number_markers))
      }
      for (marker_index in seq_len(number_markers)) {
        posterior_mean <- colMeans(
          marker_probability[, marker_index, , drop = FALSE],
          dims = 2L,
          na.rm = TRUE
        )
        assigned_class <- which.max(posterior_mean)
        for (group in seq_len(number_classes)) {
          values <- marker_probability[, marker_index, group]
          bounds <- stats::quantile(
            values,
            c(0.025, 0.975),
            names = FALSE,
            na.rm = TRUE
          )
          marker_rows[[length(marker_rows) + 1L]] <- data.frame(
            subject = subject_name,
            unit = marker_labels[marker_index],
            class = paste0("class_", group),
            Estimate = mean(values, na.rm = TRUE),
            Est.Error = stats::sd(values, na.rm = TRUE),
            Q2.5 = bounds[1L],
            Q97.5 = bounds[2L],
            Rhat = NA_real_,
            ess_bulk = NA_real_,
            ess_tail = NA_real_,
            assigned_class = paste0("class_", assigned_class),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  output <- list()
  if (length(subject_rows) > 0L) {
    output$subject <- .round_summary_table(
      do.call(rbind, subject_rows),
      digits = digits
    )
  }
  if (length(marker_rows) > 0L) {
    output$marker <- .round_summary_table(
      do.call(rbind, marker_rows),
      digits = digits
    )
  }
  represented_draws <- if (length(selected_draws) > 0L) {
    first_subject_draws <- Filter(
      Negate(is.null),
      unname(selected_draws[[1L]])
    ) # available allocation-domain arrays for the first predicted subject
    first_domain <- if (length(first_subject_draws) > 0L) {
      first_subject_draws[[1L]]
    } else {
      NULL
    } # first available probability array used only to report retained draws
    if (is.null(first_domain)) NA_integer_ else dim(first_domain)[1L]
  } else {
    NA_integer_
  } # retained prediction draws represented in the membership summaries
  structure(
    output,
    class = c("JoiNMePosteriorClass", "list"),
    mixture = object$mixture %||% object$metadata$mixture,
    summary = TRUE,
    metadata = list(
      summary = TRUE,
      draws = represented_draws,
      digits = as.integer(digits),
      n_classes = number_classes,
      class_type = (
        object$mixture %||% object$metadata$mixture
      )$class_type %||% character(0),
      source = "dynamic prediction"
    )
  )
}

#' Extract and summarise one allocation domain
#'
#' @keywords internal
#' @noRd
.mixture_membership_ <- function(
  object,
  domain,
  number_units,
  unit_labels,
  number_classes,
  draws,
  seed,
  digits,
  summary
) {
  variables <- as.vector(outer(
    seq_len(number_units),
    seq_len(number_classes),
    function(unit, group) {
      paste0(
        "posterior_class_probability_",
        domain,
        "[",
        unit,
        ",",
        group,
        "]"
      )
    }
  ))
  if (!isTRUE(summary)) {
    return(.get_draws_matrix(
      object$fit,
      variables = variables,
      draws = draws,
      seed = seed
    ))
  }

  probability_draws <- .get_draws_array(
    object$fit,
    variables = variables,
    draws = draws,
    seed = seed
  ) # iteration-by-chain draws retain the structure needed for R-hat and ESS
  probability_summary <- .assoc_summary_from_draw_array(
    probability_draws,
    term_labels = variables,
    digits = digits
  ) # one diagnostic-aware summary row for every unit-by-class probability
  probability_summary <- .mixture_summary_columns(
    probability_summary,
    identifier_columns = "term"
  ) # retain the compact reporting schema and drop alias columns (Mean, Median, SD)

  rows <- vector("list", number_units * number_classes)
  row_position <- 1L
  posterior_means <- matrix(
    NA_real_,
    nrow = number_units,
    ncol = number_classes
  )
  for (unit in seq_len(number_units)) {
    for (group in seq_len(number_classes)) {
      variable <- paste0(
        "posterior_class_probability_",
        domain,
        "[",
        unit,
        ",",
        group,
        "]"
      )
      summary_position <- match(
        variable,
        probability_summary$term
      ) # exact Stan generated-quantity row for this allocation probability
      probability_row <- probability_summary[
        summary_position,
        setdiff(names(probability_summary), "term"),
        drop = FALSE
      ]
      posterior_means[unit, group] <- probability_row$Estimate
      rows[[row_position]] <- cbind(data.frame(
        unit = as.character(unit_labels[unit]),
        class = paste0("class_", group),
        stringsAsFactors = FALSE
      ), probability_row)
      row_position <- row_position + 1L
    }
  }
  table <- do.call(rbind, rows)
  maximum_class <- max.col(posterior_means, ties.method = "first")
  table$assigned_class <- paste0(
    "class_",
    maximum_class[match(table$unit, unit_labels)]
  )
  .round_summary_table(
    table,
    digits = digits
  )
}

#' Print posterior class probabilities
#'
#' @param x A `JoiNMePosteriorClass` object returned by [posterior_class()].
#' @param max_rows Maximum number of unit-by-class probability rows printed for
#'   each allocation domain. The returned object always retains the full table.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.JoiNMePosteriorClass <- function(x, max_rows = 50L, ...) {
  metadata <- attr(x, "metadata") %||% list(
    summary = isTRUE(attr(x, "summary")),
    n_classes = (attr(x, "mixture") %||% list())$n_classes,
    class_type = (attr(x, "mixture") %||% list())$class_type
  ) # presentation metadata retained when the posterior probabilities were built
  .cli_summary_heading("Posterior class membership", level = 1L)
  .cli_print_bullets(c(
    paste0("Source: ", metadata$source %||% "fitted model"),
    paste0("Classes: ", metadata$n_classes %||% NA_integer_),
    if (length(metadata$class_type %||% character(0)) > 0L) {
      paste0(
        "Class types: ",
        paste(metadata$class_type, collapse = ", ")
      )
    } else {
      NULL
    },
    paste0("Posterior draws: ", metadata$draws %||% NA_integer_),
    if (isTRUE(metadata$summary)) {
      "Probabilities use Estimate, Est.Error, Q2.5 and Q97.5; chain diagnostics are reported when available."
    } else {
      "Raw posterior probability draws are retained."
    }
  ))
  for (domain in names(x)) {
    .cli_summary_heading(
      paste0(tools::toTitleCase(domain), " allocation"),
      level = 2L
    )
    if (is.data.frame(x[[domain]])) {
      domain_table <- x[[domain]] # complete unit-by-class posterior summary
      allocation_identifiers <- intersect(
        c("subject", "unit", "assigned_class"),
        names(domain_table)
      ) # columns identifying each natural allocation and its modal class
      if ("assigned_class" %in% allocation_identifiers) {
        allocation_table <- unique(domain_table[, allocation_identifiers, drop = FALSE])
        allocation_counts <- as.data.frame(table(
          allocation_table$assigned_class
        )) # compact count of maximum-posterior-probability allocations
        names(allocation_counts) <- c("assigned_class", "units")
        .cli_print_table_section(
          "Maximum-probability allocation counts",
          allocation_counts,
          level = 3L
        )
      }
      probability_table <- domain_table[, setdiff(
        names(domain_table),
        "assigned_class"
      ), drop = FALSE] # probability summaries without a repeated modal-class column
      printed_table <- utils::head(
        probability_table,
        as.integer(max_rows)
      ) # concise console view; the full returned object is unchanged
      .cli_print_table_section(
        "Unit-by-class probabilities",
        printed_table,
        level = 3L
      )
      if (nrow(probability_table) > nrow(printed_table)) {
        .cli_print_bullets(paste0(
          nrow(probability_table) - nrow(printed_table),
          " further rows are retained in the returned object."
        ))
      }
    } else if (is.list(x[[domain]])) {
      .cli_print_bullets(
        paste0(
          length(x[[domain]]),
          " subject-specific posterior class draw objects."
        )
      )
    } else {
      .cli_print_bullets(
        paste0(
          nrow(x[[domain]]),
          " posterior draws for ",
          ncol(x[[domain]]),
          " unit-class probabilities."
        )
      )
    }
  }
  invisible(x)
}

#' Plot posterior class probabilities
#'
#' @description
#' Displays posterior mean class probabilities and their central 95% intervals
#' with a common colour for each latent class. Subject and marker allocation
#' domains are returned separately because they have different natural
#' observational units.
#'
#' @param x A summary-form `JoiNMePosteriorClass` object.
#' @param domain Optional allocation domain, either `"subject"` or `"marker"`.
#' @param ... Unused.
#'
#' @return A ggplot when one domain is requested or available; otherwise a
#'   named list of ggplots.
#' @export
plot.JoiNMePosteriorClass <- function(x, domain = NULL, ...) {
  metadata <- attr(x, "metadata") %||% list(
    summary = isTRUE(attr(x, "summary"))
  ) # membership representation and fitted-class metadata
  if (!isTRUE(metadata$summary)) {
    cli::cli_abort(c(
      x = "Raw class-probability draws cannot be plotted directly.",
      i = "Call {.fn posterior_class} with {.arg summary = TRUE}."
    ))
  }
  available_domains <- names(x)[vapply(
    x,
    is.data.frame,
    logical(1)
  )] # allocation domains with plottable posterior summary tables
  if (!is.null(domain)) {
    domain <- match.arg(domain, available_domains)
    available_domains <- domain
  }
  if (length(available_domains) == 0L) {
    cli::cli_abort("No posterior class-probability summary is available.")
  }

  plots <- lapply(available_domains, function(domain_name) {
    plot_data <- x[[domain_name]] # full unit-by-class probability summary
    plot_data$entity <- if ("subject" %in% names(plot_data)) {
      paste(plot_data$subject, plot_data$unit, sep = ": ")
    } else {
      as.character(plot_data$unit)
    } # natural allocation label used for colour and grouping
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$entity,
        y = .data$Estimate,
        ymin = .data$Q2.5,
        ymax = .data$Q97.5,
        colour = .data$class,
        group = .data$class
      )
    ) +
      ggplot2::geom_pointrange(
        position = ggplot2::position_dodge(width = 0.45),
        linewidth = 0.4
      ) +
      ggplot2::labs(
        title = paste0(
          tools::toTitleCase(domain_name),
          " posterior class probabilities"
        ),
        x = tools::toTitleCase(domain_name),
        y = "Posterior probability",
        colour = "Latent class"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )
  })
  names(plots) <- available_domains
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  plots
}

#' Summarise a latent-progress mixture fit
#'
#' @description
#' Extends the ordinary JoiNMe summary with baseline class probabilities,
#' selected component locations and scales, class-membership regression, and a
#' compact posterior allocation-count distribution. The allocation table
#' contains one row per class and active allocation domain; use
#' [posterior_class()] when subject-by-class or marker-by-class probabilities
#' are required.
#'
#' @inheritParams summary.JoiNMeFit
#' @return A `summary_JoiNMeMixFit` object inheriting from
#'   `summary_JoiNMeFit`. 
#' @export
summary.JoiNMeMixFit <- function(
  object,
  draws = NULL,
  seed = .Random.seed[[1]],
  digits = 3,
  include_corr = TRUE,
  ...
) {
  cache_key <- paste0(
    "mixture_summary_draws=",
    draws,
    "_seed=",
    seed,
    "_digits=",
    digits,
    "_corr=",
    as.integer(include_corr)
  )
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  # Obtain the established model summary
  output <- summary.JoiNMeFit(
    object,
    draws = draws,
    seed = seed,
    digits = digits,
    include_corr = include_corr,
    ...
  )
  mixture <- object$mixture %||% object$config$mixture
  if (!isTRUE(mixture$include_survival)) {
    # The likelihood-neutral scaffold exists only to reuse the common design
    # and Stan programme.  Its baseline-hazard coefficients are prior-only
    # nuisance quantities and must not be presented as though an event process
    # had been fitted.
    output$tables$baseline_hazard <- NULL
    output$tables$survival_process <- NULL
    output$tables$assoc <- NULL
    output$tables$transform_parameters <- NULL
    output$tables$piecewise_ordinates <- NULL
    output$metadata$event_process <- "not fitted"
  }
  number_classes <- as.integer(mixture$n_classes)
  total_dimension <- as.integer(mixture$total_dimension)

  probability_variables <- paste0(
    "mix_probability[",
    seq_len(number_classes),
    "]"
  )
  probability_draws <- .get_draws_array(
    object$fit,
    variables = probability_variables,
    draws = draws,
    seed = seed
  )
  class_probability <- .assoc_summary_from_draw_array(
    probability_draws,
    term_labels = paste0("class_", seq_len(number_classes)),
    digits = digits
  )
  class_probability <- .mixture_summary_columns(
    class_probability,
    identifier_columns = "term"
  ) # explicit mean, median, spread, interval and convergence columns by class

  class_location <- NULL
  class_scale <- NULL
  if (total_dimension > 0L) {
    location_variables <- as.vector(outer(
      seq_len(number_classes),
      seq_len(total_dimension),
      function(group, coordinate) {
        paste0("mix_location[", group, ",", coordinate, "]")
      }
    ))
    scale_variables <- as.vector(outer(
      seq_len(number_classes),
      seq_len(total_dimension),
      function(group, coordinate) {
        paste0("mix_scale[", group, ",", coordinate, "]")
      }
    ))
    coordinate_labels <- .mixture_coordinate_labels(object)
    class_location <- .mixture_parameter_summary(
      object,
      variables = location_variables,
      number_classes = number_classes,
      coordinate_labels = coordinate_labels,
      parameter = "location",
      draws = draws,
      seed = seed,
      digits = digits
    )
    class_scale <- .mixture_parameter_summary(
      object,
      variables = scale_variables,
      number_classes = number_classes,
      coordinate_labels = coordinate_labels,
      parameter = "scale",
      draws = draws,
      seed = seed,
      digits = digits
    )
  }

  output$tables$class_probability <- class_probability
  output$tables$class_location <- class_location
  output$tables$class_scale <- class_scale
  output$tables$class_regression <- list(
    subject = .mixture_class_regression_summary(
      object = object,
      domain = "subject",
      number_classes = number_classes,
      draws = draws,
      seed = seed,
      digits = digits
    ),
    marker = .mixture_class_regression_summary(
      object = object,
      domain = "marker",
      number_classes = number_classes,
      draws = draws,
      seed = seed,
      digits = digits
    )
  )
  output$tables$class_regression <- Filter(
    Negate(is.null),
    output$tables$class_regression
  )
  full_class_membership <- posterior_class(
    object,
    draws = draws,
    seed = seed,
    summary = TRUE
  )
  full_class_membership_draws <- posterior_class(
    object,
    draws = draws,
    seed = seed,
    summary = FALSE
  ) # draw-level probabilities needed for posterior allocation-count intervals
  output$tables$class_membership <- .mixture_membership_overview(
    class_membership = full_class_membership,
    class_membership_draws = full_class_membership_draws,
    digits = digits
  )
  output$tables$diagnostics <- .summary_diagnostics_table(
    sampler_diagnostics = output$diagnostics %||% list(),
    reported_tables = output$tables[
      setdiff(names(output$tables), "diagnostics")
    ]
  ) # include mixture parameters in headline convergence counts and extrema
  output$metadata$mixture <- mixture
  class(output) <- unique(c(
    "summary_JoiNMeMixFit",
    "summary_JoiNMeFit",
    class(output)
  ))
  object$cache_set(cache_key, output)
  output
}

#' Summarise allocation probabilities without printing every unit-class pair
#'
#' @description
#' `posterior_class()` deliberately returns one row for every allocation unit
#' and latent class because those probabilities are useful for classification
#' and uncertainty assessment. Repeating that complete table inside
#' `summary()` makes an otherwise concise model summary grow in proportion to
#' the number of subjects or markers. This helper therefore reduces each
#' active allocation domain to one row per fitted class. For every posterior
#' draw, the expected class count is the sum of that draw's soft allocation
#' probabilities. The returned estimate, posterior standard deviation and
#' interval therefore describe posterior uncertainty in the expected count. The
#' assigned count separately uses each unit's maximum posterior-mean class.
#'
#' @param class_membership Named posterior-class tables returned by
#'   `posterior_class(summary = TRUE)`.
#' @param class_membership_draws Named draw matrices returned by
#'   `posterior_class(summary = FALSE)`.
#' @param digits Number of decimal places retained in posterior summaries.
#'
#' @return A named list containing one compact data frame per active allocation
#'   domain. Inactive domains are absent.
#' @keywords internal
#' @noRd
.mixture_membership_overview <- function(
  class_membership,
  class_membership_draws = NULL,
  digits = 3
) {
  output <- lapply(names(class_membership), function(domain) {
    membership <- class_membership[[domain]] # unit-by-class posterior summary for one active allocation domain
    membership_draws <- as.matrix(
      class_membership_draws[[domain]] %||% matrix(numeric(0), 0L, 0L)
    ) # draws by unit-class probabilities for the same allocation domain
    required_columns <- c(
      "unit",
      "class",
      "assigned_class"
    ) # fields required to form class-level assigned counts
    if (
      !is.data.frame(membership) ||
        !all(required_columns %in% names(membership)) ||
        nrow(membership) == 0L ||
        nrow(membership_draws) == 0L ||
        ncol(membership_draws) == 0L
    ) {
      return(NULL)
    }

    class_names <- unique(
      as.character(membership$class)
    ) # fitted class labels in their established reporting order
    draw_variable_names <- colnames(membership_draws) %||%
      character(0) # Stan probability names carrying unit and class indices
    draw_class_index <- suppressWarnings(as.integer(sub(
      "^.*\\[[0-9]+,([0-9]+)\\]$",
      "\\1",
      draw_variable_names
    ))) # fitted class index parsed from each unit-class probability variable
    if (
      length(draw_class_index) != ncol(membership_draws) ||
        any(!is.finite(draw_class_index))
    ) {
      cli::cli_abort(
        "Could not identify class indices in posterior allocation draws for {domain}."
      )
    }

    rows <- lapply(seq_along(class_names), function(class_index) {
      class_name <- class_names[[class_index]] # display label of this fitted class
      class_rows <- membership[
        as.character(membership$class) == class_name,
        ,
        drop = FALSE
      ] # exactly one posterior-probability row per unit for this class
      assigned_units <- unique(
        as.character(
          class_rows$unit[
            as.character(class_rows$assigned_class) == class_name
          ]
        )
      ) # units whose maximum posterior probability selects this class
      expected_count_by_draw <- rowSums(
        membership_draws[
          ,
          draw_class_index == class_index,
          drop = FALSE
        ]
      ) # soft number of allocation units belonging to this class in each draw
      count_summary <- .summarize_draw_col(
        expected_count_by_draw
      ) # posterior distribution of the draw-specific expected class count
      data.frame(
        class = class_name,
        Estimate = unname(count_summary[["Estimate"]]),
        Est.Error = unname(count_summary[["Est.Error"]]),
        Q2.5 = unname(count_summary[["Q2.5"]]),
        Q97.5 = unname(count_summary[["Q97.5"]]),
        assigned_units = length(assigned_units),
        stringsAsFactors = FALSE
      )
    })
    .round_summary_table(
      do.call(rbind, rows),
      digits = digits
    )
  })
  names(output) <- names(class_membership)
  Filter(Negate(is.null), output)
}

#' Summarise class-membership regression coefficients
#'
#' @keywords internal
#' @noRd
.mixture_class_regression_summary <- function(
  object,
  domain,
  number_classes,
  draws,
  seed,
  digits
) {
  mixture <- object$mixture %||% object$config$mixture
  design_record <- (
    mixture$class_design %||% list()
  )[[domain]] %||% list() # fitted class-design record for this allocation domain

  # The fitted Stan dimension is the authoritative coefficient count.  Older
  # serialized objects can carry placeholder column labels such as "class_:"
  # even when no class-regression coefficient was estimated.
  fitted_coefficient_count <- as.integer(
    object$stan_data[[paste0("P_class_", domain)]] %||% 0L
  )
  if (
    length(fitted_coefficient_count) != 1L ||
      !is.finite(fitted_coefficient_count) ||
      fitted_coefficient_count < 1L
  ) {
    return(NULL)
  }

  available_variables <- tryCatch(
    posterior::variables(.get_draws_obj(object$fit)),
    error = function(error) character(0)
  ) # saved posterior variable names for this fitted object
  coefficient_pattern <- paste0(
    "^mix_class_coefficient_",
    domain,
    "\\[(\\d+)\\]$"
  )
  coefficient_variables <- grep(
    coefficient_pattern,
    available_variables,
    value = TRUE
  )
  if (length(coefficient_variables) < 1L) {
    # Zero-dimensional vectors are omitted from Stan output. Returning NULL
    # keeps summaries aligned with the fitted posterior rather than metadata
    # placeholders.
    return(NULL)
  }
  coefficient_index <- suppressWarnings(as.integer(sub(
    coefficient_pattern,
    "\\1",
    coefficient_variables
  )))
  coefficient_index <- coefficient_index[is.finite(coefficient_index)]
  if (length(coefficient_index) < 1L) {
    return(NULL)
  }
  number_covariates <- min(
    fitted_coefficient_count,
    max(coefficient_index)
  ) # common coefficient count represented in both standata and draws
  if (number_covariates < 1L) {
    return(NULL)
  }

  variables <- paste0(
    "mix_class_coefficient_",
    domain,
    "[",
    seq_len(number_covariates),
    "]"
  ) # compact Stan coefficient names in concatenated formula order
  coefficient_draws <- .get_draws_array(
    object$fit,
    variables = variables,
    draws = draws,
    seed = seed
  ) # posterior draws of the requested domain coefficients

  coefficient_labels <- as.character(
    design_record$coefficient_covariate %||% character(0)
  )
  coefficient_class <- as.integer(
    design_record$coefficient_class %||% integer(0)
  )

  # Use fitted class-design labels when they are complete; otherwise retain a
  # neutral statistical label based on coefficient position.
  if (length(coefficient_labels) < number_covariates) {
    coefficient_labels <- paste0(
      "class_formula_covariate_",
      seq_len(number_covariates)
    )
  } else {
    coefficient_labels <- coefficient_labels[seq_len(number_covariates)]
  }
  if (length(coefficient_class) < number_covariates) {
    coefficient_class <- rep.int(1L, number_covariates)
  } else {
    coefficient_class <- coefficient_class[seq_len(number_covariates)]
  }

  output <- .assoc_summary_from_draw_array(
    coefficient_draws,
    term_labels = coefficient_labels,
    digits = digits
  ) # established posterior summary columns and diagnostics
  output$domain <- domain
  output$class <- paste0("class_", coefficient_class)
  output$covariate <- coefficient_labels
  .mixture_summary_columns(
    output,
    identifier_columns = c("domain", "class", "covariate")
  )
}

#' Select explicit posterior columns for latent-class reporting
#'
#' @description
#' This helper gives every mixture table the same compact schema as the
#' ordinary JoiNMe summary: posterior mean (`Estimate`), posterior standard
#' deviation (`Est.Error`), central interval and MCMC diagnostics. Duplicate
#' mean, median and standard-deviation aliases are intentionally excluded.
#'
#' @param table Posterior summary table returned by
#'   .assoc_summary_from_draw_array().
#' @param identifier_columns Columns identifying the scientific estimand.
#'
#' @return The mixture-specific reporting columns in a stable order.
#' @keywords internal
#' @noRd
.mixture_summary_columns <- function(table, identifier_columns) {
  if (is.null(table)) {
    return(NULL)
  }
  required_columns <- c(
    identifier_columns,
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat",
    "ess_bulk",
    "ess_tail"
  ) # compact inferential schema shown for every latent-class parameter
  missing_columns <- setdiff(
    required_columns,
    names(table)
  ) # absent columns indicate a malformed or obsolete summary source
  if (length(missing_columns) > 0L) {
    cli::cli_abort(
      "Mixture summary is missing column{?s}: {paste(missing_columns, collapse = ', ')}."
    )
  }
  table[, required_columns, drop = FALSE]
}

#' Label mixture coordinates in their original random-effect blocks
#'
#' @keywords internal
#' @noRd
.mixture_coordinate_labels <- function(object) {
  mixture <- object$mixture %||% object$config$mixture
  stan_data <- object$stan_data
  block_terms <- list(
    subject = as.character(
      stan_data$zid_cols %||%
        paste0("subject_", seq_len(stan_data$R_id))
    ),
    marker = as.character(
      stan_data$zmk_cols %||%
        paste0("marker_", seq_len(stan_data$R_mk))
    ),
    corr = .mixture_covariance_coordinate_labels(stan_data),
    vcov = .mixture_covariance_coordinate_labels(stan_data)
  )
  labels <- character(0)
  for (level in names(mixture$dimensions)) {
    indices <- mixture$dimensions[[level]]
    if (length(indices) > 0L) {
      available <- block_terms[[level]]
      labels <- c(
        labels,
        paste0(level, ": ", available[indices])
      )
    }
  }
  labels
}

#' @keywords internal
#' @noRd
.mixture_covariance_coordinate_labels <- function(stan_data) {
  q_dimension <- as.integer(stan_data$Q_idm %||% 0L)
  term_labels <- as.character(
    stan_data$zidm_cols %||%
      paste0("covariance_", seq_len(q_dimension))
  )
  if (q_dimension < 1L) {
    return(character(0))
  }
  if (as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L) {
    return(paste0("SD(", term_labels, ")"))
  }
  labels <- character(0)
  for (row in seq_len(q_dimension)) {
    for (column in seq_len(row)) {
      labels <- c(
        labels,
        if (row == column) {
          paste0("SD(", term_labels[row], ")")
        } else {
          paste0(
            "partial-correlation(",
            term_labels[row],
            ", ",
            term_labels[column],
            ")"
          )
        }
      )
    }
  }
  labels
}

#' Extract independent covariance-regression slopes for latent-class displays
#'
#' @param object Fitted latent-progress object.
#' @param draws Number of posterior draws.
#' @param seed Optional draw-selection seed.
#'
#' @return A list containing SD and correlation slope arrays, their reference
#'   covariate profiles, and a function which packs observed-covariate
#'   contributions into the established lower-triangular coordinate order.
#' @keywords internal
#' @noRd
.mixture_covariance_regression_draws <- function(object, draws, seed) {
  stan_data <- object$stan_data # fitted dimensions and the two subject-level design matrices
  q_dimension <- as.integer(stan_data$Q_idm %||% 0L) # number of marker-by-subject basis coordinates
  correlation_dimension <- if (
    q_dimension > 1L &&
      as.integer(stan_data$indep_idmarker_cov %||% 0L) == 0L
  ) as.integer(q_dimension * (q_dimension - 1L) / 2L) else 0L
  covariance_design <- .stored_vcov_design(stan_data) # split design or earlier common-design view
  k_sd <- covariance_design$k_sd # formulaVCov$sd slope count
  k_corr <- covariance_design$k_corr # formulaVCov$corr slope count
  if (isTRUE(covariance_design$shared_format) && q_dimension > 0L && k_sd > 0L) {
    packed_dimension <- if (
      as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
    ) q_dimension else as.integer(q_dimension * (q_dimension + 1L) / 2L) # former beta_L row count
    shared_beta <- .mixture_matrix_draws(
      object, "beta_L", packed_dimension, k_sd, draws, seed
    ) # posterior slopes stored by fits predating the split covariance formula
    diagonal_positions <- integer(q_dimension) # packed beta_L rows supplying SD predictors
    correlation_positions <- integer(correlation_dimension) # packed beta_L rows supplying partial-correlation predictors
    packed_position <- 1L
    correlation_position <- 1L
    for (row in seq_len(q_dimension)) {
      columns <- if (as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L) row else seq_len(row)
      for (column in columns) {
        if (row == column) {
          diagonal_positions[row] <- packed_position
        } else {
          correlation_positions[correlation_position] <- packed_position
          correlation_position <- correlation_position + 1L
        }
        packed_position <- packed_position + 1L
      }
    }
    beta_sd <- shared_beta[, diagonal_positions, , drop = FALSE]
    beta_corr <- if (correlation_dimension > 0L) {
      shared_beta[, correlation_positions, , drop = FALSE]
    } else array(0, dim = c(draws, 0L, k_corr))
  } else {
    beta_sd <- if (q_dimension > 0L && k_sd > 0L) {
      .mixture_matrix_draws(object, "beta_L_sd", q_dimension, k_sd, draws, seed)
    } else array(0, dim = c(draws, q_dimension, k_sd))
    beta_corr <- if (correlation_dimension > 0L && k_corr > 0L) {
      .mixture_matrix_draws(object, "beta_L_corr", correlation_dimension, k_corr, draws, seed)
    } else array(0, dim = c(draws, correlation_dimension, k_corr))
  }
  reference_sd <- if (k_sd > 0L) {
    colMeans(covariance_design$x_sd, na.rm = TRUE)
  } else numeric(0)
  reference_corr <- if (k_corr > 0L) {
    colMeans(covariance_design$x_corr, na.rm = TRUE)
  } else numeric(0)

  pack_contribution <- function(sd_covariates = reference_sd,
                                corr_covariates = reference_corr) {
    number_draws <- dim(beta_sd)[1L] %||% dim(beta_corr)[1L] # common posterior draw count
    covariance_dimension <- if (
      as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L
    ) q_dimension else as.integer(q_dimension * (q_dimension + 1L) / 2L)
    contribution <- matrix(0, nrow = number_draws, ncol = covariance_dimension)
    packed_coordinate <- 1L
    correlation_coordinate <- 1L
    for (row in seq_len(q_dimension)) {
      for (column in if (as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L) row else seq_len(row)) {
        if (row == column) {
          if (k_sd > 0L) contribution[, packed_coordinate] <- beta_sd[, row, ] %*% sd_covariates
        } else {
          if (k_corr > 0L) contribution[, packed_coordinate] <- beta_corr[, correlation_coordinate, ] %*% corr_covariates
          correlation_coordinate <- correlation_coordinate + 1L
        }
        packed_coordinate <- packed_coordinate + 1L
      }
    }
    contribution
  }

  list(
    beta_sd = beta_sd,
    beta_corr = beta_corr,
    reference_sd = reference_sd,
    reference_corr = reference_corr,
    x_sd = covariance_design$x_sd,
    x_corr = covariance_design$x_corr,
    k_sd = k_sd,
    k_corr = k_corr,
    pack_contribution = pack_contribution
  )
}

#' Identify the selected covariance-regression class representation
#'
#' @param mixture Fitted latent-progress metadata.
#'
#' @return `"corr"` or `"vcov"` when selected, otherwise an empty string.
#' @keywords internal
#' @noRd
.mixture_covariance_class_type <- function(mixture) {
  selected <- intersect(
    c("corr", "vcov"),
    mixture$class_type %||% character(0)
  ) # public covariance representation active in the shared Stan latent block
  if (length(selected) == 1L) selected[[1L]] else ""
}

#' @keywords internal
#' @noRd
.mixture_parameter_summary <- function(
  object,
  variables,
  number_classes,
  coordinate_labels,
  parameter,
  draws,
  seed,
  digits
) {
  draw_array <- .get_draws_array(
    object$fit,
    variables = variables,
    draws = draws,
    seed = seed
  )
  term_labels <- rep(
    coordinate_labels,
    each = number_classes
  )
  output <- .assoc_summary_from_draw_array(
    draw_array,
    term_labels = term_labels,
    digits = digits
  )
  index <- seq_len(nrow(output))
  output$class <- paste0(
    "class_",
    ((index - 1L) %% number_classes) + 1L
  )
  output$coordinate <- term_labels
  output$parameter <- parameter
  .mixture_summary_columns(
    output,
    identifier_columns = c("class", "coordinate", "parameter")
  )
}

#' @export
print.JoiNMeMixFit <- function(x, ...) {
  x$print(...)
}

#' @export
print.summary_JoiNMeMixFit <- function(x, ...) {
  print.summary_JoiNMeFit(x, ...)
  mixture <- x$metadata$mixture %||% list()
  .cli_summary_heading("Latent-progress mixture", level = 2L)
  .cli_print_bullets(c(
    paste0("Classes: ", mixture$n_classes %||% NA_integer_),
    paste0(
      "Class types: ",
      paste(mixture$class_type %||% character(0), collapse = ", ")
    ),
    paste0("Class ordering: ", mixture$ordering %||% "none"),
    paste0("Component family: ", mixture$distribution %||% "unknown")
  ))
  .cli_print_bullets(c(
    "Class locations are centres of the selected latent random-effect coordinates; they are not response-scale means.",
    "Within-class scales describe dispersion around those centres before the ordinary random-effect covariance transformation is applied."
  ))
  .cli_print_table_section(
    "Baseline class probabilities (class covariates equal zero)",
    x$tables$class_probability,
    level = 3L
  )
  .cli_print_table_section(
    "Class locations on the standardised latent random-effect scale",
    x$tables$class_location,
    level = 3L
  )
  .cli_print_table_section(
    "Within-class scales on the standardised latent random-effect scale",
    x$tables$class_scale,
    level = 3L
  )
  class_regression <- x$tables$class_regression %||% list()
  for (domain in names(class_regression)) {
    class_design <- mixture$class_design[[domain]] %||%
      list() # fitted formula structure for this allocation domain
    if (!isTRUE(class_design$class_specific)) {
      .cli_print_bullets(
        paste0(
          "For the shared ",
          domain,
          " class formula, class_",
          mixture$n_classes,
          " is the zero-coefficient reference and is intentionally absent from the regression table."
        )
      )
    }
    .cli_print_table_section(
      paste0("Class-membership regression: ", domain),
      class_regression[[domain]],
      level = 3L
    )
  }
  membership_overviews <- x$tables$class_membership %||% list()
  if (length(membership_overviews) > 0L) {
    .cli_print_bullets(
      "Allocation Estimate, Est.Error and interval columns describe the posterior expected class count; assigned_units is the maximum-probability hard count."
    )
  }
  for (domain in names(membership_overviews)) {
    .cli_print_table_section(
      paste0("Posterior allocation counts: ", domain),
      membership_overviews[[domain]],
      level = 3L
    )
  }
  invisible(x)
}

#' Predict from a latent-progress mixture
#'
#' @description
#' Runs the established dynamic prediction calculation while carrying the
#' fitted component probability, location and scale into the priors for newly
#' sampled latent effects.  The new subject's history therefore informs a
#' posterior probability vector over the fitted classes.  Combined compatible
#' blocks retain one shared allocation. With `reuse_fitted_re = TRUE`, the
#' fitted subject's realised effects and fitted allocation probabilities are
#' retained draw by draw instead of evaluating a new-subject allocation.
#'
#' @inheritParams predict.JoiNMeFit
#' @return A `JoiNMeMixDynPred` object inheriting from `JoiNMeDynPred`.
#'   Conditional allocation draws are retained in
#'   `draws$posterior_class` and may be summarised with
#'   [posterior_class()].
#' @export
predict.JoiNMeMixFit <- function(
  object,
  newdataLong,
  newdataEvent = NULL,
  process = c("longitudinal", "event"),
  reuse_fitted_re = FALSE,
  ...
) {
  mixture <- object$mixture %||% object$config$mixture
  if (!isTRUE(mixture$include_survival)) {
    if (missing(process)) {
      process <- "longitudinal"
    }
    if ("event" %in% process) {
      cli::cli_abort(c(
        x = "Event prediction is unavailable for a longitudinal-only mixture.",
        i = "Use {.arg process = 'longitudinal'}."
      ))
    }
    if (is.null(newdataEvent)) {
      newdataEvent <- .longitudinal_only_event_scaffold(
        data_long = newdataLong,
        id_variable = .get_call_args(object$call, "id_var", "id"),
        time_variable = .get_call_args(
          object$call,
          "time_var",
          "time"
        )
      )
    }
  } else if (is.null(newdataEvent)) {
    cli::cli_abort(
      "{.arg newdataEvent} is required because the fitted mixture includes a survival process."
    )
  }

  prediction <- predict.JoiNMeFit(
    object = object,
    newdataLong = newdataLong,
    newdataEvent = newdataEvent,
    process = process,
    reuse_fitted_re = reuse_fitted_re,
    ...
  )
  prediction$metadata$mixture <- mixture
  JoiNMeMixDynPred$new(
    predictions = prediction$predictions,
    quantiles = prediction$quantiles,
    draws = prediction$draws,
    data = prediction$data,
    metadata = prediction$metadata,
    call = prediction$call,
    tmax = prediction$tmax,
    n_samples = prediction$n_samples,
    mixture = mixture
  )
}

#' @export
print.JoiNMeMixDynPred <- function(x, ...) {
  print.JoiNMeDynPred(x, ...)
  mixture <- x$mixture %||% x$metadata$mixture
  cat(
    "Latent classes: ",
    mixture$n_classes %||% NA_integer_,
    " (",
    paste(mixture$class_type %||% character(0), collapse = ", "),
    ")\n",
    sep = ""
  )
  invisible(x)
}

#' @export
summary.JoiNMeMixDynPred <- function(object, ...) {
  output <- summary.JoiNMeDynPred(object, ...)
  output$metadata$mixture <- object$mixture %||% object$metadata$mixture
  if (length(object$draws$posterior_class %||% list()) > 0L) {
    output$tables$class_membership <- posterior_class(object)
  }
  class(output) <- unique(c(
    "summary_JoiNMeMixDynPred",
    "summary_JoiNMeDynPred",
    class(output)
  ))
  output
}

#' @export
print.summary_JoiNMeMixDynPred <- function(x, ...) {
  print.summary_JoiNMeDynPred(x, ...)
  mixture <- x$metadata$mixture %||% list()
  .cli_print_bullets(
    paste0(
      "Parent mixture: ",
      mixture$n_classes %||% NA_integer_,
      " classes over ",
      paste(mixture$class_type %||% character(0), collapse = ", ")
    )
  )
  memberships <- x$tables$class_membership %||% list()
  for (domain in names(memberships)) {
    .cli_print_table_section(
      paste0("Conditional class membership: ", domain),
      memberships[[domain]],
      level = 3L
    )
  }
  invisible(x)
}

#' Plot class-specific latent-progress summaries
#'
#' @description
#' Adds four mixture displays to the inherited JoiNMe plotting interface:
#'
#' - `type = "longitudinal"` with `estimand = "mean_per_class"` evaluates the
#'   fixed-effects trajectory plus the centre of every class-specific random-effect
#'   block;
#' - `estimand = "marginal_per_class"` averages the inverse-link trajectory over
#'   the estimated within-class random-effect distribution;
#' - `type = "class_membership"` shows posterior allocation probabilities;
#' - `type = "covariance_class"` shows class-specific covariance-regression
#'   curves, with all \(G\) classes in each panel.
#'
#' All other plot types are delegated to [plot.JoiNMeFit()].
#'
#' @param x A `JoiNMeMixFit` object.
#' @param type Plot type.
#' @param estimand Optional class-specific longitudinal estimand.
#' @param draws Number of posterior draws.
#' @param seed Posterior subsetting seed.
#' @param longitudinal_times Optional trajectory time grid.
#' @param longitudinal_points Number of default trajectory time points.
#' @param marginal_samples Number of within-component Monte Carlo values per
#'   posterior draw for `estimand = "marginal_per_class"`.
#' @param ci_levels Credible interval levels.
#' @param marker Optional marker subset.
#' @param theme_fn ggplot2 theme function.
#' @param ... Arguments forwarded to inherited plotting methods.
#'
#' @return A ggplot, a combined display, or a named list of ggplots.
#' @export
plot.JoiNMeMixFit <- function(
  x,
  type = "longitudinal",
  estimand = NULL,
  draws = 400,
  seed = 1,
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  marginal_samples = 32L,
  ci_levels = c(0.5, 0.95),
  marker = NULL,
  theme_fn = ggplot2::theme_bw,
  ...
) {
  # Plot wrappers carry this private flag solely to avoid a dispatch cycle.
  # Remove it here, then add it exactly once if the request is delegated to the
  # ordinary JoiNMe plotting method.
  dot_arguments <- list(...)
  dot_arguments$.use_wrapper_dispatch <- NULL
  type <- as.character(type)
  if (
    any(type %in% c("survival", "cumhaz")) &&
      !isTRUE((x$mixture %||% x$config$mixture)$include_survival)
  ) {
    cli::cli_abort(c(
      x = "Survival plots are unavailable for a longitudinal-only mixture.",
      i = "Fit with {.arg formulaEvent} and {.arg dataEvent} to model an event process."
    ))
  }
  if (
    any(type %in% "association") &&
      !isTRUE((x$mixture %||% x$config$mixture)$include_survival)
  ) {
    cli::cli_abort(c(
      x = "Association plots are unavailable for a longitudinal-only mixture.",
      i = "Fit with {.arg formulaEvent} and {.arg dataEvent} to model an event process."
    ))
  }
  if (
    length(type) == 1L &&
      identical(type, "class_membership")
  ) {
    return(.plot_mixture_membership(
      x,
      draws = draws,
      seed = seed,
      theme_fn = theme_fn
    ))
  }
  if (
    length(type) == 1L &&
      identical(type, "covariance_class")
  ) {
    return(.plot_mixture_covariance(
      x,
      draws = draws,
      seed = seed,
      marker = marker,
      ci_levels = ci_levels,
      theme_fn = theme_fn
    ))
  }
  if (
    length(type) == 1L &&
      identical(type, "association") &&
      !is.null(estimand)
  ) {
    return(do.call(
      .plot_mixture_association_trajectory,
      c(
        list(
          object = x,
          estimand = estimand,
          draws = draws,
          seed = seed,
          longitudinal_times = longitudinal_times,
          longitudinal_points = longitudinal_points,
          marginal_samples = marginal_samples,
          ci_levels = ci_levels,
          marker = marker,
          theme_fn = theme_fn
        ),
        dot_arguments
      )
    ))
  }
  if (
    length(type) == 1L &&
      identical(type, "longitudinal") &&
      !is.null(estimand)
  ) {
    estimand <- match.arg(
      estimand,
      c("mean_per_class", "marginal_per_class")
    )
    trajectory <- .mixture_class_trajectory(
      x,
      estimand = estimand,
      draws = draws,
      seed = seed,
      longitudinal_times = longitudinal_times,
      longitudinal_points = longitudinal_points,
      marginal_samples = marginal_samples,
      ci_levels = ci_levels,
      marker = marker
    )
    return(.plot_mixture_trajectory_data(
      trajectory$data,
      estimand = estimand,
      ci_levels = ci_levels,
      theme_fn = theme_fn
    ))
  }

  do.call(
    plot.JoiNMeFit,
    c(
      list(
        x = x,
        type = type,
        draws = draws,
        seed = seed,
        longitudinal_times = longitudinal_times,
        longitudinal_points = longitudinal_points,
        ci_levels = ci_levels,
        marker = if (is.null(marker)) NA else marker,
        theme_fn = theme_fn,
        .use_wrapper_dispatch = FALSE
      ),
      dot_arguments
    )
  )
}

#' @export
plot.JoiNMeMixDynPred <- function(x, ...) {
  dot_arguments <- list(...)
  requested_type <- as.character(dot_arguments$type %||% "combined")
  mixture <- x$mixture %||% x$metadata$mixture %||% list()
  if (
    length(requested_type) == 1L &&
    identical(requested_type, "class_membership")
  ) {
    dot_arguments$type <- NULL
    posterior_draws <- dot_arguments$draws %||% NULL
    random_seed <- dot_arguments$seed %||% 1L
    theme_function <- dot_arguments$theme_fn %||% ggplot2::theme_bw
    dot_arguments$draws <- NULL
    dot_arguments$seed <- NULL
    dot_arguments$theme_fn <- NULL
    membership <- posterior_class(
      x,
      draws = posterior_draws,
      seed = random_seed,
      summary = TRUE
    )
    plot_rows <- list()
    if (!is.null(membership$subject)) {
      table <- membership$subject
      table$subject <- table$unit
      table$entity <- "subject"
      table$domain <- "subject allocation"
      plot_rows[[length(plot_rows) + 1L]] <- table
    }
    if (!is.null(membership$marker)) {
      table <- membership$marker
      table$entity <- table$unit
      table$domain <- "marker allocation"
      plot_rows[[length(plot_rows) + 1L]] <- table
    }
    plot_data <- do.call(rbind, plot_rows)
    output <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$class,
        y = .data$Estimate,
        colour = .data$entity,
        group = .data$entity
      )
    ) +
      ggplot2::geom_pointrange(
        ggplot2::aes(
          ymin = .data$Q2.5,
          ymax = .data$Q97.5
        ),
        position = ggplot2::position_dodge(width = 0.35)
      ) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$domain),
        cols = ggplot2::vars(.data$subject),
        scales = "free_x"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        x = "Latent-progress class",
        y = "Conditional posterior probability",
        colour = "Allocation unit"
      ) +
      theme_function()
    return(output)
  }
  if (
    !isTRUE(mixture$include_survival) &&
    !("type" %in% names(dot_arguments))
  ) {
    # The inherited dynamic plot defaults to a combined longitudinal-survival
    # display.  A longitudinal-only mixture has no event panel, so its natural
    # no-argument display is the longitudinal prediction itself.
    return(plot.JoiNMeDynPred(x, type = "longitudinal", ...))
  }
  if (
    any(requested_type %in% c("survival", "cumhaz", "combined")) &&
      !isTRUE(mixture$include_survival)
  ) {
    cli::cli_abort(c(
      x = "Event-process prediction plots are unavailable for a longitudinal-only mixture.",
      i = "Use {.arg type = 'longitudinal'}."
    ))
  }
  plot.JoiNMeDynPred(x, ...)
}

#' Require a fitted survival process for event-based summaries
#'
#' @keywords internal
#' @noRd
.require_mixture_survival <- function(object, method) {
  mixture <- object$mixture %||% object$config$mixture %||% list()
  if (!isTRUE(mixture$include_survival)) {
    cli::cli_abort(c(
      x = "{.fn {method}} is unavailable for a longitudinal-only mixture.",
      i = "Fit with {.arg formulaEvent} and {.arg dataEvent} to use event-based validation."
    ))
  }
  invisible(TRUE)
}

#' Event discrimination for latent-progress mixtures
#'
#' @inheritParams concordance.JoiNMeFit
#' @export
concordance.JoiNMeMixFit <- function(object, ...) {
  .require_mixture_survival(object, "concordance")
  concordance.JoiNMeFit(object, ...)
}

#' @inheritParams tvROC.JoiNMeFit
#' @export
tvROC.JoiNMeMixFit <- function(object, ...) {
  .require_mixture_survival(object, "tvROC")
  tvROC.JoiNMeFit(object, ...)
}

#' @inheritParams tvAUC.JoiNMeFit
#' @export
tvAUC.JoiNMeMixFit <- function(object, ...) {
  .require_mixture_survival(object, "tvAUC")
  tvAUC.JoiNMeFit(object, ...)
}

#' Plot posterior membership probabilities
#'
#' @keywords internal
#' @noRd
.plot_mixture_membership <- function(
  object,
  draws,
  seed,
  theme_fn
) {
  membership <- posterior_class(
    object,
    draws = draws,
    seed = seed,
    summary = TRUE
  )
  plots <- lapply(names(membership), function(domain) {
    data <- membership[[domain]]
    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data$unit,
        y = .data$Estimate,
        fill = .data$class
      )
    ) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = paste("Posterior class probabilities:", domain),
        x = domain,
        y = "Posterior probability",
        fill = "Latent class"
      ) +
      theme_fn()
  })
  names(plots) <- names(membership)
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  .combine_plot_grid(plots, fallback = "input")
}

#' Extract a matrix-valued posterior variable into a draw array
#'
#' @keywords internal
#' @noRd
.mixture_matrix_draws <- function(
  object,
  variable,
  number_rows,
  number_columns,
  draws,
  seed
) {
  if (number_rows < 1L || number_columns < 1L) {
    return(array(
      numeric(0),
      dim = c(0L, number_rows, number_columns)
    ))
  }
  variables <- as.vector(outer(
    seq_len(number_rows),
    seq_len(number_columns),
    function(row, column) {
      paste0(variable, "[", row, ",", column, "]")
    }
  ))
  draw_matrix <- .get_draws_matrix(
    object$fit,
    variables = variables,
    draws = draws,
    seed = seed
  )
  output <- array(
    NA_real_,
    dim = c(nrow(draw_matrix), number_rows, number_columns)
  )
  for (row in seq_len(number_rows)) {
    for (column in seq_len(number_columns)) {
      output[, row, column] <- draw_matrix[
        ,
        paste0(variable, "[", row, ",", column, "]")
      ]
    }
  }
  output
}

#' Draw centred class trajectories from the fitted posterior
#'
#' @keywords internal
#' @noRd
.mixture_class_trajectory <- function(
  object,
  estimand,
  draws,
  seed,
  longitudinal_times,
  longitudinal_points,
  marginal_samples,
  ci_levels,
  marker = NULL,
  retain_association_samples = FALSE
) {
  stan_data <- object$stan_data
  mixture <- object$mixture %||% object$config$mixture
  number_classes <- as.integer(mixture$n_classes)
  marker_levels <- as.character(
    stan_data$marker_levels %||%
      seq_len(as.integer(stan_data$D))
  )
  if (!is.null(marker) && !all(is.na(marker))) {
    marker_levels <- intersect(marker_levels, as.character(marker))
  }
  if (length(marker_levels) == 0L) {
    cli::cli_abort("No fitted marker remains after applying {.arg marker}.")
  }

  id_variable <- .get_call_args(object$call, "id_var", "id")
  time_variable <- .get_call_args(
    object$call,
    "time_var",
    "time"
  )
  marker_variable <- .get_call_args(
    object$call,
    "marker_var",
    "marker"
  )

  # Select a covariate profile from the observed data for every marker.  The
  # profile is held fixed over the common time grid, so differences between
  # displayed curves arise from posterior class structure rather than changing
  # covariate composition.
  marker_rows <- lapply(marker_levels, function(marker_label) {
    candidates <- object$dataLong[
      as.character(object$dataLong[[marker_variable]]) == marker_label,
      ,
      drop = FALSE
    ]
    candidates[1L, , drop = FALSE]
  })
  reference_rows <- do.call(rbind, marker_rows)
  observed_time <- as.numeric(object$dataLong[[time_variable]])
  if (is.null(longitudinal_times)) {
    longitudinal_points <- as.integer(longitudinal_points)
    if (
      length(longitudinal_points) != 1L ||
        !is.finite(longitudinal_points) ||
        longitudinal_points < 2L
    ) {
      cli::cli_abort(
        "{.arg longitudinal_points} must be an integer of at least two."
      )
    }
    time_range <- range(observed_time[is.finite(observed_time)])
    longitudinal_times <- seq(
      time_range[1L],
      time_range[2L],
      length.out = longitudinal_points
    )
  } else {
    longitudinal_times <- sort(unique(as.numeric(longitudinal_times)))
  }

  evaluation_rows <- lapply(seq_len(nrow(reference_rows)), function(row) {
    output <- reference_rows[
      rep(row, length(longitudinal_times)),
      ,
      drop = FALSE
    ]
    output[[time_variable]] <- longitudinal_times
    output
  })
  evaluation_data <- do.call(rbind, evaluation_rows)
  evaluation_data[[marker_variable]] <- factor(
    evaluation_data[[marker_variable]],
    levels = as.character(stan_data$marker_levels)
  )
  evaluation_data$marker_int <- match(
    as.character(evaluation_data[[marker_variable]]),
    as.character(stan_data$marker_levels)
  )
  rownames(evaluation_data) <- NULL

  design <- ..longitudinal_design_matrices(
    fitted_model = object,
    longitudinal_evaluation_data = evaluation_data,
    time_variable = time_variable,
    marker_variable = marker_variable
  )
  beta_variables <- paste0("beta[", seq_len(as.integer(stan_data$P)), "]")
  beta_draws <- .get_draws_matrix(
    object$fit,
    variables = beta_variables,
    draws = draws,
    seed = seed
  )
  number_draws <- nrow(beta_draws)
  fixed_predictor <- beta_draws %*% t(design$fixed)

  subject_mean <- .mixture_matrix_draws(
    object,
    "class_mean_subject",
    number_classes,
    as.integer(stan_data$R_id),
    draws = number_draws,
    seed = seed
  )
  marker_mean <- .mixture_matrix_draws(
    object,
    "class_mean_marker",
    number_classes,
    as.integer(stan_data$R_mk),
    draws = number_draws,
    seed = seed
  )
  subject_cholesky <- .mixture_matrix_draws(
    object,
    "L_u",
    as.integer(stan_data$R_id),
    as.integer(stan_data$R_id),
    draws = number_draws,
    seed = seed
  )
  marker_cholesky <- .mixture_matrix_draws(
    object,
    "L_v",
    as.integer(stan_data$R_mk),
    as.integer(stan_data$R_mk),
    draws = number_draws,
    seed = seed
  )
  q_dimension <- as.integer(stan_data$Q_idm %||% 0L)
  covariance_dimension <- if (q_dimension > 0L) {
    length(.mixture_covariance_coordinate_labels(stan_data))
  } else {
    0L
  }
  covariance_mean <- .mixture_matrix_draws(
    object,
    "class_mean_covariance",
    number_classes,
    covariance_dimension,
    draws = number_draws,
    seed = seed
  )
  covariance_intercept <- if (covariance_dimension > 0L) {
    .get_draws_matrix(
      object$fit,
      variables = paste0(
        "alpha_L[",
        seq_len(covariance_dimension),
        "]"
      ),
      draws = number_draws,
      seed = seed
    )
  } else {
    matrix(0, nrow = number_draws, ncol = 0L)
  }
  covariance_loading <- if (covariance_dimension > 0L) {
    .get_draws_matrix(
      object$fit,
      variables = paste0(
        "lambda_L[",
        seq_len(covariance_dimension),
        "]"
      ),
      draws = number_draws,
      seed = seed
    )
  } else {
    matrix(0, nrow = number_draws, ncol = 0L)
  }
  covariance_regression <- .mixture_covariance_regression_draws(
    object,
    draws = number_draws,
    seed = seed
  ) # independent posterior SD/correlation slopes and their observed reference profiles
  covariance_reference_contribution <- covariance_regression$pack_contribution()
  marker_cross_loading <- if (
    q_dimension > 0L &&
      as.integer(stan_data$R_mk) > 0L &&
      as.integer(stan_data$allow_marker_crosscorr %||% 0L) == 1L
  ) {
    .mixture_matrix_draws(
      object,
      "B_cross",
      q_dimension,
      as.integer(stan_data$R_mk),
      draws = number_draws,
      seed = seed
    )
  } else {
    array(
      0,
      dim = c(
        number_draws,
        q_dimension,
        as.integer(stan_data$R_mk)
      )
    )
  }
  class_scale <- NULL
  class_location <- NULL
  marginal_subject_samples <- NULL
  marginal_marker_samples <- NULL
  marginal_covariance_samples <- NULL
  marginal_nested_noise <- NULL
  if (identical(estimand, "marginal_per_class")) {
    marginal_samples <- as.integer(marginal_samples)
    if (
      length(marginal_samples) != 1L ||
        is.na(marginal_samples) ||
        marginal_samples < 1L
    ) {
      cli::cli_abort(
        "{.arg marginal_samples} must be a positive integer."
      )
    }
    class_scale <- .mixture_matrix_draws(
      object,
      "mix_scale",
      number_classes,
      as.integer(mixture$total_dimension),
      draws = number_draws,
      seed = seed
    )
    class_location <- .mixture_matrix_draws(
      object,
      "mix_location",
      number_classes,
      as.integer(mixture$total_dimension),
      draws = number_draws,
      seed = seed
    )
    # Marginalisation is performed within every posterior draw.  Several
    # component realisations are transformed to the response scale and
    # averaged, leaving a posterior distribution of the response-scale
    # within-class expectation rather than a predictive distribution from one
    # arbitrary component draw.
    sampled_latent <- lapply(
      seq_len(marginal_samples),
      function(sample_index) {
        .mixture_within_class_samples(
          mixture,
          class_location,
          class_scale,
          seed = .mixture_safe_seed(
            seed,
            104729 * sample_index
          )
        )
      }
    )
    marginal_subject_samples <- lapply(
      seq_along(sampled_latent),
      function(sample_index) {
        .mixture_full_block_sample(
          block_centre = subject_mean,
          sampled_latent = sampled_latent[[sample_index]],
          mixture = mixture,
          level = "subject",
          seed = .mixture_safe_seed(
            seed,
            130363 * sample_index
          )
        )
      }
    )
    marginal_marker_samples <- lapply(
      seq_along(sampled_latent),
      function(sample_index) {
        .mixture_full_block_sample(
          block_centre = marker_mean,
          sampled_latent = sampled_latent[[sample_index]],
          mixture = mixture,
          level = "marker",
          seed = .mixture_safe_seed(
            seed,
            155921 * sample_index
          )
        )
      }
    )
    marginal_covariance_samples <- lapply(
      seq_along(sampled_latent),
      function(sample_index) {
        .mixture_full_block_sample(
          block_centre = covariance_mean,
          sampled_latent = sampled_latent[[sample_index]],
          mixture = mixture,
          level = {
            selected_covariance_type <-
              .mixture_covariance_class_type(mixture)
            if (nzchar(selected_covariance_type)) {
              selected_covariance_type
            } else {
              "vcov"
            }
          },
          seed = .mixture_safe_seed(
            seed,
            180289 * sample_index
          )
        )
      }
    )
    marginal_nested_noise <- lapply(
      seq_len(marginal_samples),
      function(sample_index) {
        .mixture_normal_array(
          dimensions = c(
            number_draws,
            number_classes,
            length(marker_levels),
            q_dimension
          ),
          seed = .mixture_safe_seed(
            seed,
            205759 * sample_index
          )
        )
      }
    )
  }

  predictor_by_class <- vector("list", number_classes)
  mean_predictor_by_class <- vector("list", number_classes)
  marker_predictor_by_class <- vector("list", number_classes)
  response_by_class <- vector("list", number_classes)
  association_parts_by_class <- vector("list", number_classes)
  marker_index <- as.integer(evaluation_data$marker_int)

  # Construct the three structural longitudinal parts for one class and one
  # set of standardised component values.
  class_predictor_parts <- function(
    subject_latent,
    marker_latent,
    covariance_latent,
    nested_noise,
    group
  ) {
    mean_predictor <- fixed_predictor
    if (as.integer(stan_data$R_id) > 0L) {
      realised_subject <- matrix(
        0,
        nrow = number_draws,
        ncol = as.integer(stan_data$R_id)
      )
      for (draw in seq_len(number_draws)) {
        realised_subject[draw, ] <-
          subject_cholesky[draw, , ] %*%
            subject_latent[draw, group, ]
      }
      mean_predictor <-
        mean_predictor + realised_subject %*% t(design$subject)
    }
    marker_predictor <- matrix(
      0,
      nrow = number_draws,
      ncol = nrow(evaluation_data)
    )
    if (as.integer(stan_data$R_mk) > 0L) {
      realised_marker <- matrix(
        0,
        nrow = number_draws,
        ncol = as.integer(stan_data$R_mk)
      )
      for (draw in seq_len(number_draws)) {
        realised_marker[draw, ] <-
          marker_cholesky[draw, , ] %*%
            marker_latent[draw, group, ]
      }
      marker_predictor <-
        marker_predictor + realised_marker %*% t(design$marker)
    } else {
      realised_marker <- matrix(
        0,
        nrow = number_draws,
        ncol = as.integer(stan_data$R_mk)
      )
    }

    # The nested marker effect is centred through `z_w_lat`, but marker
    # cross-correlation can give it a non-zero class centre whenever marker
    # effects have latent classes.  For a marginal trajectory, `nested_noise`
    # contributes one standard-Normal innovation per marker and posterior draw.
    if (q_dimension > 0L) {
      if (is.null(nested_noise)) {
        nested_noise <- array(
          0,
          dim = c(
            number_draws,
            number_classes,
            length(marker_levels),
            q_dimension
          )
        )
      }
      for (draw in seq_len(number_draws)) {
        covariance_predictor <- as.numeric(
          covariance_intercept[draw, ]
        )
        covariance_predictor <- covariance_predictor +
          covariance_reference_contribution[draw, ]
        covariance_predictor <- covariance_predictor +
          as.numeric(covariance_loading[draw, ]) *
            covariance_latent[draw, group, ]
        covariance_cholesky <-
            .mixture_cholesky_from_predictor(
              predictor = covariance_predictor,
              q_dimension = q_dimension,
              diagonal_only = as.integer(
                stan_data$indep_idmarker_cov %||% 0L
              ) == 1L,
              diagonal_link = as.integer(
                stan_data$vcov_diag_link %||% 0L
              )
            )
        cross_centre <- if (
          as.integer(stan_data$R_mk) > 0L &&
            as.integer(
              stan_data$allow_marker_crosscorr %||% 0L
            ) == 1L
        ) {
          matrix(
            marker_cross_loading[
              draw,
              ,
              ,
              drop = FALSE
            ],
            nrow = q_dimension,
            ncol = as.integer(stan_data$R_mk)
          ) %*%
            realised_marker[draw, ]
        } else {
          rep(0, q_dimension)
        }
        for (marker_position in seq_along(marker_levels)) {
          marker_columns <- which(
            marker_index == match(
              marker_levels[marker_position],
              as.character(stan_data$marker_levels)
            )
          )
          nested_effect <- covariance_cholesky %*% (
            as.numeric(cross_centre) +
              nested_noise[
                draw,
                group,
                marker_position,
                ,
                drop = TRUE
              ]
          )
          marker_predictor[draw, marker_columns] <-
            marker_predictor[draw, marker_columns] +
              as.numeric(
                design$subject_marker[
                  marker_columns,
                  ,
                  drop = FALSE
                ] %*% nested_effect
              )
        }
      }
    }
    list(
      mean = mean_predictor,
      marker = marker_predictor,
      total = mean_predictor + marker_predictor
    )
  }

  # Apply the correct marker-specific inverse link to an entire draws-by-row
  # predictor matrix.
  response_from_predictor <- function(predictor) {
    response <- predictor
    for (marker_number in unique(marker_index)) {
      columns <- which(marker_index == marker_number)
      n_operations <- as.integer(stan_data$inv_link_n_ops[marker_number])
      n_constants <- as.integer(stan_data$inv_link_n_const[marker_number])
      response[, columns] <- .apply_inverse_link_matrix(
        eta = predictor[, columns, drop = FALSE],
        link_code = as.integer(stan_data$link_long[marker_number]),
        bytecode = if (n_operations > 0L) {
          as.integer(
            stan_data$inv_link_ops[
              marker_number,
              seq_len(n_operations)
            ]
          )
        } else {
          integer(0)
        },
        const_data = if (n_constants > 0L) {
          as.numeric(
            stan_data$inv_link_const[
              marker_number,
              seq_len(n_constants)
            ]
          )
        } else {
          numeric(0)
        }
      )
    }
    response
  }

  for (group in seq_len(number_classes)) {
    # The event model distinguishes the population-plus-subject trajectory
    # from the marker-specific deviation and their total.  Keep each centre
    # part so an association plot can use exactly its fitted channel.
    centre_parts <- class_predictor_parts(
      subject_latent = subject_mean,
      marker_latent = marker_mean,
      covariance_latent = covariance_mean,
      nested_noise = NULL,
      group = group
    )
    response <- if (identical(estimand, "marginal_per_class")) {
      marginal_response <- matrix(
        0,
        nrow = number_draws,
        ncol = nrow(evaluation_data)
      )
      sampled_association_parts <- if (
        isTRUE(retain_association_samples)
      ) {
        vector("list", marginal_samples)
      } else {
        NULL
      }
      for (sample_index in seq_len(marginal_samples)) {
        sampled_parts <- class_predictor_parts(
          subject_latent =
            marginal_subject_samples[[sample_index]],
          marker_latent =
            marginal_marker_samples[[sample_index]],
          covariance_latent =
            marginal_covariance_samples[[sample_index]],
          nested_noise =
            marginal_nested_noise[[sample_index]],
          group = group
        )
        marginal_response <- marginal_response +
          response_from_predictor(sampled_parts$total)
        if (isTRUE(retain_association_samples)) {
          sampled_association_parts[[sample_index]] <- sampled_parts
        }
      }
      association_parts_by_class[[group]] <-
        sampled_association_parts
      marginal_response / marginal_samples
    } else {
      association_parts_by_class[[group]] <- list(centre_parts)
      response_from_predictor(centre_parts$total)
    }

    mean_predictor_by_class[[group]] <- centre_parts$mean
    marker_predictor_by_class[[group]] <- centre_parts$marker
    predictor_by_class[[group]] <- centre_parts$total
    response_by_class[[group]] <- response
  }

  probability_points <- .quantile_probs_from_ci_plot(ci_levels)
  quantile_names <- .quantile_colnames(probability_points)
  summary_rows <- lapply(seq_len(number_classes), function(group) {
    response <- response_by_class[[group]]
    quantiles <- t(apply(
      response,
      2,
      stats::quantile,
      probs = probability_points,
      na.rm = TRUE,
      names = FALSE
    ))
    output <- data.frame(
      class = factor(
        paste0("class_", group),
        levels = paste0("class_", seq_len(number_classes))
      ),
      marker = as.character(evaluation_data[[marker_variable]]),
      time = as.numeric(evaluation_data[[time_variable]]),
      mean = colMeans(response),
      stringsAsFactors = FALSE
    )
    output[quantile_names] <- quantiles
    output
  })
  list(
    data = do.call(rbind, summary_rows),
    mean_linpred = mean_predictor_by_class,
    marker_linpred = marker_predictor_by_class,
    linpred = predictor_by_class,
    epred = response_by_class,
    evaluation_data = evaluation_data,
    marker_index = marker_index,
    marker_levels = as.character(stan_data$marker_levels),
    time_variable = time_variable,
    number_draws = number_draws,
    association_parts = association_parts_by_class
  )
}

#' Draw standardised within-class effects
#'
#' @keywords internal
#' @noRd
.mixture_within_class_samples <- function(
  mixture,
  location,
  scale,
  seed
) {
  dimensions <- dim(location)
  if (length(dimensions) != 3L || dimensions[3L] == 0L) {
    return(location)
  }
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(
      assign(".Random.seed", previous_seed, envir = .GlobalEnv),
      add = TRUE
    )
  } else {
    on.exit(
      rm(".Random.seed", envir = .GlobalEnv),
      add = TRUE
    )
  }
  set.seed(.mixture_safe_seed(seed, 319))
  component_family <- as.character(mixture$distribution %||% "student_t") # family retained from jm_prior(class$family)
  component_df <- as.numeric(mixture$distribution_df %||% 6) # fixed Student-t degrees of freedom, ignored by other component families
  noise <- switch(
    component_family,
    student_t = array(
      stats::rt(prod(dimensions), df = component_df),
      dim = dimensions
    ),
    laplace = {
      signs <- sample(c(-1, 1), prod(dimensions), replace = TRUE)
      array(
        signs * stats::rexp(prod(dimensions), rate = 1),
        dim = dimensions
      )
    },
    array(stats::rnorm(prod(dimensions)), dim = dimensions)
  )
  location + scale * noise
}

#' Convert a user seed and offset to a valid R integer seed
#'
#' @keywords internal
#' @noRd
.mixture_safe_seed <- function(seed, offset = 0) {
  seed <- as.numeric(seed %||% 1)[1L]
  offset <- as.numeric(offset %||% 0)[1L]
  if (!is.finite(seed)) {
    seed <- 1
  }
  if (!is.finite(offset)) {
    offset <- 0
  }
  as.integer(
    (abs(seed + offset) %% (.Machine$integer.max - 1)) + 1
  )
}

#' Draw a reproducible standard-Normal array without changing the R session
#'
#' @keywords internal
#' @noRd
.mixture_normal_array <- function(dimensions, seed) {
  dimensions <- as.integer(dimensions)
  if (length(dimensions) == 0L) {
    return(numeric(0))
  }
  has_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  if (has_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(
      assign(
        ".Random.seed",
        previous_seed,
        envir = .GlobalEnv
      ),
      add = TRUE
    )
  } else {
    on.exit(
      rm(".Random.seed", envir = .GlobalEnv),
      add = TRUE
    )
  }
  set.seed(.mixture_safe_seed(seed))
  array(
    stats::rnorm(prod(dimensions)),
    dim = dimensions
  )
}

#' Draw a complete standardised random-effect block within a class
#'
#' @description
#' Selected coordinates use their fitted component draw. Coordinates outside
#' the progress plane retain the ordinary standard-Normal distribution. This
#' distinction is required for a genuinely marginal response-scale trajectory;
#' setting unselected coordinates to zero would instead condition on their
#' centres.
#'
#' @keywords internal
#' @noRd
.mixture_full_block_sample <- function(
  block_centre,
  sampled_latent,
  mixture,
  level,
  seed
) {
  dimensions <- dim(block_centre)
  if (
    length(dimensions) != 3L ||
      any(dimensions == 0L)
  ) {
    return(block_centre)
  }
  full_sample <- .mixture_normal_array(
    dimensions = dimensions,
    seed = seed
  )
  .replace_mixture_block_samples(
    block = full_sample,
    sampled_latent = sampled_latent,
    mixture = mixture,
    level = level
  )
}

#' Insert sampled coordinates into a full random-effect block
#'
#' @keywords internal
#' @noRd
.replace_mixture_block_samples <- function(
  block,
  sampled_latent,
  mixture,
  level
) {
  indices <- as.integer(mixture$dimensions[[level]] %||% integer(0))
  if (length(indices) == 0L) {
    return(block)
  }
  start <- as.integer(mixture$starts[[level]])
  source_indices <- start + seq_along(indices) - 1L
  block[, , indices] <- sampled_latent[, , source_indices, drop = FALSE]
  block
}

#' Draw the class trajectory summary
#'
#' @keywords internal
#' @noRd
.plot_mixture_trajectory_data <- function(
  data,
  estimand,
  ci_levels,
  theme_fn
) {
  widest_interval <- max(ci_levels)
  interval_names <- .quantile_names_from_ci(widest_interval)
  median_name <- .quantile_name_from_prob(0.5)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$time,
      y = .data[[median_name]],
      colour = .data$class,
      fill = .data$class,
      group = .data$class
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data[[interval_names[1L]]],
        ymax = .data[[interval_names[2L]]]
      ),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::facet_wrap(ggplot2::vars(.data$marker), scales = "free_y") +
    ggplot2::labs(
      title = switch(
        estimand,
        mean_per_class = "Longitudinal trajectory at each class centre",
        marginal_per_class = "Marginal longitudinal trajectory within each class"
      ),
      x = "Time",
      y = "Expected marker value",
      colour = "Latent class",
      fill = "Latent class"
    ) +
    theme_fn()
}

#' Plot association contributions along class trajectories
#'
#' @keywords internal
#' @noRd
.plot_mixture_association_trajectory <- function(
  object,
  estimand,
  draws,
  seed,
  longitudinal_times,
  longitudinal_points,
  marginal_samples,
  ci_levels,
  marker,
  theme_fn,
  ...
) {
  dot_arguments <- list(...)
  association_term <- dot_arguments$association_term %||%
    dot_arguments$term %||%
    .available_association_terms(object)[1L]
  if (is.na(association_term) || !nzchar(association_term)) {
    cli::cli_abort("No fitted association term is available.")
  }
  term_key <- sub("\\[.*$", "", association_term)
  if (term_key %in% c("corr", "vcov")) {
    return(.plot_mixture_covariance(
      object,
      draws = draws,
      seed = seed,
      marker = marker,
      ci_levels = ci_levels,
      theme_fn = theme_fn,
      association_term = association_term,
      association_points = longitudinal_points
    ))
  }

  trajectory <- .mixture_class_trajectory(
    object,
    estimand = match.arg(
      estimand,
      c("mean_per_class", "marginal_per_class")
    ),
    draws = draws,
    seed = seed,
    longitudinal_times = longitudinal_times,
    longitudinal_points = longitudinal_points,
    marginal_samples = marginal_samples,
    ci_levels = ci_levels,
    marker = marker,
    retain_association_samples = TRUE
  )
  association_data <- .get_association_plot_data(object, seed = seed)
  association_data <- .mixture_subset_association_data(
    association_data,
    object = object,
    number_draws = trajectory$number_draws,
    seed = seed
  )
  coefficient_draws <- .assoc_coeff_draws(
    object,
    term = association_term,
    data = association_data,
    seed = seed
  )
  coefficient_draws <- coefficient_draws[
    seq_len(min(length(coefficient_draws), trajectory$number_draws))
  ]
  number_draws <- length(coefficient_draws)

  rows <- list()
  row_position <- 1L
  times <- as.numeric(
    trajectory$evaluation_data[[trajectory$time_variable]]
  )
  markers <- as.character(
    trajectory$evaluation_data[[
      .get_call_args(object$call, "marker_var", "marker")
    ]]
  )
  marker_levels <- unique(markers)
  common_times <- sort(unique(times))
  for (group in seq_along(trajectory$linpred)) {
    association_part_samples <-
      trajectory$association_parts[[group]]
    if (length(association_part_samples) == 0L) {
      association_part_samples <- list(list(
        mean = trajectory$mean_linpred[[group]],
        marker = trajectory$marker_linpred[[group]],
        total = trajectory$linpred[[group]]
      ))
    }
    transformed_curve <- matrix(
      0,
      nrow = number_draws,
      ncol = length(common_times)
    )
    weight_draws <- if (
      !term_key %in% c("cv_mean", "cs_mean")
    ) {
      .mixture_association_weight_draws(
        object = object,
        term_key = term_key,
        group = group,
        number_draws = number_draws,
        seed = seed,
        association_data = association_data,
        marker_levels = marker_levels
      )
    } else {
      NULL
    }

    # Evaluate every within-class realisation through the fitted transform
    # before averaging.  For a nonlinear transform, transforming the average
    # latent trajectory would be a different and generally biased estimand.
    for (sample_parts in association_part_samples) {
      source_draws <- switch(
        sub("^(cv|cs)_", "", term_key),
        mean = sample_parts$mean,
        marker = sample_parts$marker,
        total = sample_parts$total,
        sample_parts$total
      )
      source_draws <- source_draws[
        seq_len(number_draws),
        ,
        drop = FALSE
      ]

      # Stan transforms each marker contribution before it forms the weighted
      # marker average.  The mean channel is the only exception: it uses the
      # population-plus-subject trajectory once, without marker weighting.
      transformed_by_marker <- vector("list", length(marker_levels))
      for (marker_position in seq_along(marker_levels)) {
        marker_label <- marker_levels[marker_position]
        columns <- which(markers == marker_label)
        marker_times <- times[columns]
        ordered <- order(marker_times)
        raw_marker <- source_draws[, columns[ordered], drop = FALSE]
        if (grepl("^cs_", term_key)) {
          raw_marker <- t(apply(raw_marker, 1, function(value) {
            finite_difference <- diff(value) /
              pmax(
                diff(marker_times[ordered]),
                sqrt(.Machine$double.eps)
              )
            c(
              finite_difference,
              finite_difference[length(finite_difference)]
            )
          }))
        }
        transformed_by_marker[[marker_position]] <-
          .mixture_transform_draw_values(
            object = object,
            term_key = term_key,
            term = association_term,
            raw_values = raw_marker,
            seed = seed,
            association_data = association_data
          )
      }

      transformed_sample <- if (
        term_key %in% c("cv_mean", "cs_mean")
      ) {
        transformed_by_marker[[1L]]
      } else {
        weighted_average <- matrix(
          0,
          nrow = number_draws,
          ncol = length(common_times)
        )
        for (marker_position in seq_along(marker_levels)) {
          weighted_average <- weighted_average +
            .mixture_weight_draw_trajectory(
              transformed_by_marker[[marker_position]],
              weight_draws[, marker_position]
            ) # multiply each draw's complete time curve by its matching weight
        }
        weighted_average / length(marker_levels)
      }
      transformed_curve <- transformed_curve + transformed_sample
    }
    transformed_curve <-
      transformed_curve / length(association_part_samples)

    curve_draws <- transformed_curve * coefficient_draws
    for (time_position in seq_along(common_times)) {
      curve <- curve_draws[, time_position]
      rows[[row_position]] <- data.frame(
        class = paste0("class_", group),
        time = common_times[time_position],
        mean = mean(curve, na.rm = TRUE),
        sd = stats::sd(curve, na.rm = TRUE),
        lower = stats::quantile(
          curve,
          (1 - max(ci_levels)) / 2,
          na.rm = TRUE,
          names = FALSE
        ),
        median = stats::median(curve, na.rm = TRUE),
        upper = stats::quantile(
          curve,
          1 - (1 - max(ci_levels)) / 2,
          na.rm = TRUE,
          names = FALSE
        ),
        Q2.5 = stats::quantile(
          curve,
          0.025,
          na.rm = TRUE,
          names = FALSE
        ),
        Q97.5 = stats::quantile(
          curve,
          0.975,
          na.rm = TRUE,
          names = FALSE
        ),
        stringsAsFactors = FALSE
      )
      row_position <- row_position + 1L
    }
  }
  plot_data <- do.call(rbind, rows)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$time,
      y = .data$median,
      colour = .data$class,
      fill = .data$class,
      group = .data$class
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::labs(
      title = paste(
        "Association trajectory by latent class:",
        association_term
      ),
      x = "Time",
      y = "Contribution to log hazard",
      colour = "Latent class",
      fill = "Latent class"
    ) +
    theme_fn()
}

#' Evaluate one raw association value per posterior draw
#'
#' @description
#' Ordinary association plots evaluate every draw on a shared horizontal grid.
#' A class trajectory instead supplies a different raw value for every draw and
#' time.  The existing transform programme is evaluated once on a fine common
#' grid and each posterior curve is interpolated at its corresponding raw
#' value. 
#'
#' @keywords internal
#' @noRd
.mixture_transform_draw_values <- function(
  object,
  term_key,
  term,
  raw_values,
  seed,
  association_data
) {
  raw_values <- as.matrix(raw_values)
  number_draws <- nrow(raw_values)
  if (number_draws < 1L || ncol(raw_values) < 1L) {
    return(raw_values)
  }
  transform_specification <- .transform_spec_for_term(
    object,
    term_key,
    term = term
  )
  transform_type <- .canonicalise_transform_type(
    transform_specification$type %||% "identity"
  )
  if (identical(transform_type, "identity")) {
    return(raw_values)
  }

  finite_values <- raw_values[is.finite(raw_values)]
  if (length(finite_values) == 0L) {
    return(raw_values)
  }
  raw_range <- range(finite_values)
  evaluation_grid <- if (
    diff(raw_range) <= sqrt(.Machine$double.eps)
  ) {
    raw_range[1L]
  } else {
    seq(raw_range[1L], raw_range[2L], length.out = 321L)
  }
  transform_grid <- .association_transform_matrix(
    object,
    term_key = term_key,
    term = term,
    x_grid = evaluation_grid,
    n_draws = number_draws,
    seed = seed,
    data = association_data
  )
  output <- matrix(
    NA_real_,
    nrow = number_draws,
    ncol = ncol(raw_values)
  )
  for (draw in seq_len(number_draws)) {
    valid <- is.finite(raw_values[draw, ])
    if (any(valid)) {
      output[draw, valid] <- if (length(evaluation_grid) == 1L) {
        rep(transform_grid[draw, 1L], sum(valid))
      } else {
        stats::approx(
          x = evaluation_grid,
          y = transform_grid[draw, ],
          xout = raw_values[draw, valid],
          rule = 2,
          ties = "ordered"
        )$y
      }
    }
  }
  output
}

#' Apply one posterior marker weight to each posterior trajectory
#'
#' @description
#' `posterior::draws_matrix` deliberately preserves matrix dimensions when one
#' column is selected. Base array multiplication then rejects a draws-by-time
#' trajectory and a draws-by-one weight matrix as non-conformable. Coercing the
#' weights to a vector and using an explicit row-wise sweep preserves the
#' intended draw pairing for both ordinary matrices and posterior draw classes.
#'
#' @param trajectory A draws-by-time numeric matrix.
#' @param weight One marker weight per posterior draw.
#'
#' @return A draws-by-time matrix with draw-matched marker weighting.
#' @keywords internal
#' @noRd
.mixture_weight_draw_trajectory <- function(trajectory, weight) {
  trajectory <- as.matrix(
    trajectory
  ) # posterior association feature evaluated for every draw and time
  weight <- as.numeric(
    weight
  ) # draw-matched marker weights stripped of backend-specific matrix classes
  if (nrow(trajectory) != length(weight)) {
    cli::cli_abort(
      "The marker-weight draws do not align with the class trajectory draws."
    )
  }
  sweep(
    trajectory,
    MARGIN = 1L,
    STATS = weight,
    FUN = "*"
  )
}

#' Align cached association quantities with trajectory posterior draws
#'
#' @description
#' The ordinary association cache retains every posterior draw because its
#' usual plots evaluate a marginal curve.  Class trajectories jointly combine
#' fixed effects, component locations, transforms, marker weights and
#' association coefficients, so their rows must refer to the same posterior
#' draws.  This helper repeats the deterministic draw selection used by
#' `.get_draws_matrix()` and applies it recursively to cached draw matrices.
#'
#' @keywords internal
#' @noRd
.mixture_subset_association_data <- function(
  association_data,
  object,
  number_draws,
  seed
) {
  if (is.null(association_data)) {
    return(NULL)
  }
  total_draws <- posterior::ndraws(.get_draws_obj(object$fit))
  number_draws <- min(as.integer(number_draws), total_draws)
  draw_indices <- seq_len(total_draws)
  if (number_draws < total_draws) {
    has_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    if (has_seed) {
      previous_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(
        assign(
          ".Random.seed",
          previous_seed,
          envir = .GlobalEnv
        ),
        add = TRUE
      )
    } else {
      on.exit(
        rm(".Random.seed", envir = .GlobalEnv),
        add = TRUE
      )
    }
    set.seed(seed)
    draw_indices <- sample.int(
      total_draws,
      size = number_draws
    )
  }

  subset_draw_container <- function(value) {
    if (is.data.frame(value) || is.null(value)) {
      return(value)
    }
    if (is.matrix(value) && nrow(value) == total_draws) {
      return(value[draw_indices, , drop = FALSE])
    }
    if (
      is.array(value) &&
        length(dim(value)) > 0L &&
        dim(value)[1L] == total_draws
    ) {
      array_indices <- c(
        list(draw_indices),
        rep(list(TRUE), length(dim(value)) - 1L),
        list(drop = FALSE)
      )
      return(do.call(`[`, c(list(value), array_indices)))
    }
    if (is.list(value)) {
      return(lapply(value, subset_draw_container))
    }
    value
  }
  subset_draw_container(association_data)
}

#' Obtain effective marker weights for one latent class
#'
#' @description
#' Marker weights are not a latent-progress class type. The established
#' posterior effective-weight extractor is therefore used unchanged for every
#' class, preserving the ordinary JoiNMe association definition.
#'
#' @keywords internal
#' @noRd
.mixture_association_weight_draws <- function(
  object,
  term_key,
  group,
  number_draws,
  seed,
  association_data,
  marker_levels = NULL
) {
  stan_data <- object$stan_data
  fitted_marker_levels <- as.character(stan_data$marker_levels)
  marker_levels <- as.character(marker_levels %||% fitted_marker_levels)
  marker_positions <- match(marker_levels, fitted_marker_levels)
  if (anyNA(marker_positions)) {
    cli::cli_abort("A requested marker was not present in the fitted model.")
  }
  posterior_weights <- .marker_weight_draws(
    object,
    term_key = term_key,
    n_draws = number_draws,
    seed = seed,
    data = association_data
  )
  posterior_weights[, marker_positions, drop = FALSE]
}

#' Plot covariance-regression profiles at the class centres
#'
#' @keywords internal
#' @noRd
.plot_mixture_covariance <- function(
  object,
  draws,
  seed,
  marker,
  ci_levels,
  theme_fn,
  association_term = NULL,
  association_points = 80L
) {
  if (!is.null(association_term)) {
    return(.plot_mixture_covariance_association(
      object = object,
      association_term = association_term,
      draws = draws,
      seed = seed,
      ci_levels = ci_levels,
      theme_fn = theme_fn,
      points = association_points
    ))
  }
  mixture <- object$mixture %||% object$config$mixture
  if (!nzchar(.mixture_covariance_class_type(mixture))) {
    cli::cli_abort(c(
      x = "The covariance block was not assigned latent classes.",
      i = "Fit with {.arg class_type = 'corr'} or {.arg class_type = 'vcov'} to request class-specific covariance curves."
    ))
  }
  stan_data <- object$stan_data
  q_dimension <- as.integer(stan_data$Q_idm)
  covariance_dimension <- length(
    .mixture_covariance_coordinate_labels(stan_data)
  )
  number_classes <- as.integer(mixture$n_classes)
  class_mean <- .mixture_matrix_draws(
    object,
    "class_mean_covariance",
    number_classes,
    covariance_dimension,
    draws = draws,
    seed = seed
  )
  number_draws <- dim(class_mean)[1L]
  alpha <- .get_draws_matrix(
    object$fit,
    variables = paste0("alpha_L[", seq_len(covariance_dimension), "]"),
    draws = number_draws,
    seed = seed
  )
  lambda <- .get_draws_matrix(
    object$fit,
    variables = paste0("lambda_L[", seq_len(covariance_dimension), "]"),
    draws = number_draws,
    seed = seed
  )
  covariance_regression <- .mixture_covariance_regression_draws(
    object, number_draws, seed
  ) # separate formulaVCov$sd and formulaVCov$corr posterior slopes
  k_covariates <- max(covariance_regression$k_sd, covariance_regression$k_corr)
  x_grid <- if (k_covariates > 0L) {
    observed <- c(
      if (covariance_regression$k_sd > 0L) as.numeric(covariance_regression$x_sd[, 1L]) else numeric(0),
      if (covariance_regression$k_corr > 0L) as.numeric(covariance_regression$x_corr[, 1L]) else numeric(0)
    )
    seq(min(observed), max(observed), length.out = 60L)
  } else {
    seq(0, 1, length.out = 60L)
  }
  matrix_labels <- as.character(
    stan_data$zidm_cols %||%
      paste0("effect_", seq_len(q_dimension))
  )
  component_labels <- c(
    paste0("variance: ", matrix_labels),
    if (q_dimension > 1L) {
      unlist(lapply(2:q_dimension, function(row) {
        paste0(
          "correlation: ",
          matrix_labels[row],
          " with ",
          matrix_labels[seq_len(row - 1L)]
        )
      }))
    } else {
      character(0)
    }
  )

  result_rows <- list()
  result_position <- 1L
  interval <- max(ci_levels)
  for (group in seq_len(number_classes)) {
    for (x_position in seq_along(x_grid)) {
      component_draws <- matrix(
        NA_real_,
        nrow = number_draws,
        ncol = length(component_labels)
      )
      for (draw in seq_len(number_draws)) {
        sd_profile <- covariance_regression$reference_sd
        corr_profile <- covariance_regression$reference_corr
        if (length(sd_profile) > 0L) sd_profile[1L] <- x_grid[x_position]
        if (length(corr_profile) > 0L) corr_profile[1L] <- x_grid[x_position]
        observed_contribution <- covariance_regression$pack_contribution(
          sd_covariates = sd_profile,
          corr_covariates = corr_profile
        )
        linear_predictor <- as.numeric(alpha[draw, ]) + observed_contribution[draw, ]
        linear_predictor <- linear_predictor +
          as.numeric(lambda[draw, ]) *
            class_mean[draw, group, ]
        cholesky <- .mixture_cholesky_from_predictor(
          linear_predictor,
          q_dimension = q_dimension,
          diagonal_only = as.integer(
            stan_data$indep_idmarker_cov %||% 0L
          ) == 1L,
          diagonal_link = as.integer(stan_data$vcov_diag_link %||% 0L)
        )
        covariance <- tcrossprod(cholesky)
        standard_deviation <- sqrt(diag(covariance))
        correlation <- covariance /
          outer(standard_deviation, standard_deviation)
        component_draws[draw, ] <- c(
          diag(covariance),
          if (q_dimension > 1L) {
            unlist(lapply(2:q_dimension, function(row) {
              correlation[row, seq_len(row - 1L)]
            }))
          } else {
            numeric(0)
          }
        )
      }
      for (component in seq_along(component_labels)) {
        values <- component_draws[, component]
        result_rows[[result_position]] <- data.frame(
          class = paste0("class_", group),
          x = x_grid[x_position],
          component = component_labels[component],
          median = stats::median(values),
          lower = stats::quantile(
            values,
            (1 - interval) / 2,
            names = FALSE
          ),
          upper = stats::quantile(
            values,
            1 - (1 - interval) / 2,
            names = FALSE
          ),
          stringsAsFactors = FALSE
        )
        result_position <- result_position + 1L
      }
    }
  }
  curve_data <- do.call(rbind, result_rows)
  marker_levels <- as.character(stan_data$marker_levels)
  if (!is.null(marker) && !all(is.na(marker))) {
    marker_levels <- intersect(marker_levels, as.character(marker))
  }

  # The fitted covariance regression is common to marker-specific deviations.
  # Separate marker plots make that sharing explicit while satisfying the
  # multivariate reporting convention: within each plot every one of the G
  # curves is shown together.
  plots <- lapply(marker_levels, function(marker_label) {
    ggplot2::ggplot(
      curve_data,
      ggplot2::aes(
        x = .data$x,
        y = .data$median,
        colour = .data$class,
        fill = .data$class,
        group = .data$class
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
        alpha = 0.1,
        colour = NA
      ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$component),
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste("Marker-by-subject covariance:", marker_label),
        subtitle = "All latent classes are displayed in every component panel",
        x = if (k_covariates > 0L) {
          paste(
            c(
              if (covariance_regression$k_sd > 0L) paste0("SD: ", colnames(covariance_regression$x_sd)[1L]) else NULL,
              if (covariance_regression$k_corr > 0L) paste0("correlation: ", colnames(covariance_regression$x_corr)[1L]) else NULL
            ),
            collapse = "; "
          )
        } else {
          "Reference profile"
        },
        y = "Covariance summary",
        colour = "Latent class",
        fill = "Latent class"
      ) +
      theme_fn()
  })
  names(plots) <- marker_levels
  if (length(plots) == 1L) {
    return(plots[[1L]])
  }
  plots
}

#' Reconstruct a covariance matrix from one covariance-regression predictor
#'
#' @keywords internal
#' @noRd
.mixture_covariance_from_predictor <- function(
  predictor,
  q_dimension,
  diagonal_only,
  diagonal_link
) {
  cholesky <- .mixture_cholesky_from_predictor(
    predictor = predictor,
    q_dimension = q_dimension,
    diagonal_only = diagonal_only,
    diagonal_link = diagonal_link
  )
  tcrossprod(cholesky)
}

#' Reconstruct a covariance Cholesky factor from one regression predictor
#'
#' @description
#' This mirrors `include/submodels/longitudinal/transformed_parameters/fit.stan`:
#' diagonal predictors determine
#' row standard deviations, whilst off-diagonal predictors are sequential
#' partial correlations.  Keeping the factor available is necessary because
#' `corr` and `vcov` associations use Cholesky-correlation and effective-scale
#' features, not entries of the covariance matrix itself.
#'
#' @keywords internal
#' @noRd
.mixture_cholesky_from_predictor <- function(
  predictor,
  q_dimension,
  diagonal_only,
  diagonal_link
) {
  cholesky <- matrix(0, q_dimension, q_dimension)
  standard_deviation <- rep(1, q_dimension)
  position <- 1L
  if (isTRUE(diagonal_only)) {
    for (row in seq_len(q_dimension)) {
      standard_deviation[row] <- if (diagonal_link == 1L) {
        exp(predictor[position])
      } else {
        softplus(predictor[position])
      }
      cholesky[row, row] <- standard_deviation[row]
      position <- position + 1L
    }
  } else {
    for (row in seq_len(q_dimension)) {
      diagonal_position <- position + row - 1L
      standard_deviation[row] <- if (diagonal_link == 1L) {
        exp(predictor[diagonal_position])
      } else {
        softplus(predictor[diagonal_position])
      }
      cholesky[row, row] <- standard_deviation[row]
      for (column in seq_len(row)) {
        if (column < row) {
          partial_correlation <- tanh(predictor[position])
          cholesky[row, column] <-
            standard_deviation[row] * partial_correlation
        }
        position <- position + 1L
      }
    }
  }
  cholesky
}

#' Plot a covariance-style hazard contribution by latent class
#'
#' @description
#' Reconstructs the exact raw `corr` or `vcov` feature used by the Stan
#' likelihood, applies its fitted draw-specific transform, subtracts the
#' transform at zero, and multiplies by the matching association coefficient.
#' Covariance regression is usually time-constant, so the horizontal axis is
#' its first fitted covariate (or a reference profile for an intercept-only
#' model).
#'
#' @keywords internal
#' @noRd
.plot_mixture_covariance_association <- function(
  object,
  association_term,
  draws,
  seed,
  ci_levels,
  theme_fn,
  points = 80L
) {
  mixture <- object$mixture %||% object$config$mixture
  if (!nzchar(.mixture_covariance_class_type(mixture))) {
    cli::cli_abort(c(
      x = "The covariance block was not assigned latent classes.",
      i = "Fit with {.arg class_type = 'corr'} or {.arg class_type = 'vcov'} for class-specific covariance associations."
    ))
  }
  term_key <- sub("\\[.*$", "", association_term)
  if (!term_key %in% c("corr", "vcov")) {
    cli::cli_abort(
      "{.arg association_term} must identify a corr or vcov component."
    )
  }
  stan_data <- object$stan_data
  active_flag <- as.integer(
    stan_data[[paste0("assoc_", term_key)]] %||% 0L
  )
  if (active_flag != 1L) {
    cli::cli_abort(
      "Association term {.val {term_key}} was not fitted."
    )
  }

  q_dimension <- as.integer(stan_data$Q_idm)
  diagonal_only <- as.integer(
    stan_data$indep_idmarker_cov %||% 0L
  ) == 1L
  number_components <- .assoc_transform_component_count(
    term_key,
    q_idm = q_dimension,
    diagonal_only = diagonal_only
  )
  component_index <- .assoc_component_index(association_term)
  if (
    component_index < 1L ||
      component_index > number_components
  ) {
    cli::cli_abort(c(
      x = "Association component {.val {association_term}} is unavailable.",
      i = "Choose a component between 1 and {number_components}."
    ))
  }

  covariance_dimension <- length(
    .mixture_covariance_coordinate_labels(stan_data)
  )
  number_classes <- as.integer(mixture$n_classes)
  class_mean <- .mixture_matrix_draws(
    object,
    "class_mean_covariance",
    number_classes,
    covariance_dimension,
    draws = draws,
    seed = seed
  )
  number_draws <- dim(class_mean)[1L]
  alpha <- .get_draws_matrix(
    object$fit,
    variables = paste0(
      "alpha_L[",
      seq_len(covariance_dimension),
      "]"
    ),
    draws = number_draws,
    seed = seed
  )
  lambda <- .get_draws_matrix(
    object$fit,
    variables = paste0(
      "lambda_L[",
      seq_len(covariance_dimension),
      "]"
    ),
    draws = number_draws,
    seed = seed
  )
  covariance_regression <- .mixture_covariance_regression_draws(
    object, number_draws, seed
  ) # independent SD and correlation regression slopes used for association curves
  k_covariates <- max(covariance_regression$k_sd, covariance_regression$k_corr)
  points <- max(
    2L,
    as.integer(points)
  ) # number of covariance-covariate reference values used in this display
  x_grid <- if (k_covariates > 0L) {
    observed <- c(
      if (covariance_regression$k_sd > 0L) as.numeric(covariance_regression$x_sd[, 1L]) else numeric(0),
      if (covariance_regression$k_corr > 0L) as.numeric(covariance_regression$x_corr[, 1L]) else numeric(0)
    )
    seq(min(observed), max(observed), length.out = points)
  } else {
    seq(0, 1, length.out = points)
  }

  association_data <- .get_association_plot_data(
    object,
    seed = seed
  )
  association_data <- .mixture_subset_association_data(
    association_data,
    object = object,
    number_draws = number_draws,
    seed = seed
  )
  coefficient_draws <- .assoc_coeff_draws(
    object,
    term = association_term,
    data = association_data,
    seed = seed
  )
  number_draws <- min(number_draws, length(coefficient_draws))
  coefficient_draws <- coefficient_draws[seq_len(number_draws)]
  class_mean <- class_mean[seq_len(number_draws), , , drop = FALSE]
  alpha <- alpha[seq_len(number_draws), , drop = FALSE]
  lambda <- lambda[seq_len(number_draws), , drop = FALSE]

  interval <- max(ci_levels)
  result_rows <- list()
  result_position <- 1L
  for (group in seq_len(number_classes)) {
    raw_feature <- matrix(
      NA_real_,
      nrow = number_draws,
      ncol = length(x_grid)
    )
    for (x_position in seq_along(x_grid)) {
      for (draw in seq_len(number_draws)) {
        sd_profile <- covariance_regression$reference_sd
        corr_profile <- covariance_regression$reference_corr
        if (length(sd_profile) > 0L) sd_profile[1L] <- x_grid[x_position]
        if (length(corr_profile) > 0L) corr_profile[1L] <- x_grid[x_position]
        observed_contribution <- covariance_regression$pack_contribution(
          sd_covariates = sd_profile,
          corr_covariates = corr_profile
        )
        linear_predictor <- as.numeric(alpha[draw, ]) + observed_contribution[draw, ]
        linear_predictor <- linear_predictor +
          as.numeric(lambda[draw, ]) *
            class_mean[draw, group, ]
        cholesky <- .mixture_cholesky_from_predictor(
          predictor = linear_predictor,
          q_dimension = q_dimension,
          diagonal_only = diagonal_only,
          diagonal_link = as.integer(
            stan_data$vcov_diag_link %||% 0L
          )
        )
        feature_vector <- if (identical(term_key, "corr")) {
          .assoc_corr_features_from_chol(cholesky)
        } else {
          .assoc_vcov_features_from_chol(
            cholesky,
            diagonal_only = diagonal_only
          )
        }
        raw_feature[draw, x_position] <-
          feature_vector[component_index]
      }
    }

    transformed <- .mixture_transform_draw_values(
      object = object,
      term_key = term_key,
      term = association_term,
      raw_values = raw_feature,
      seed = seed,
      association_data = association_data
    )
    reference <- .association_transform_matrix(
      object,
      term_key = term_key,
      term = association_term,
      x_grid = 0,
      n_draws = number_draws,
      seed = seed,
      data = association_data
    )[, 1L]
    contribution <- (
      transformed - reference
    ) * coefficient_draws
    for (x_position in seq_along(x_grid)) {
      values <- contribution[, x_position]
      result_rows[[result_position]] <- data.frame(
        class = factor(
          paste0("class_", group),
          levels = paste0(
            "class_",
            seq_len(number_classes)
          )
        ),
        x = x_grid[x_position],
        mean = mean(values),
        sd = stats::sd(values),
        median = stats::median(values),
        lower = stats::quantile(
          values,
          (1 - interval) / 2,
          names = FALSE
        ),
        upper = stats::quantile(
          values,
          1 - (1 - interval) / 2,
          names = FALSE
        ),
        Q2.5 = stats::quantile(values, 0.025, names = FALSE),
        Q97.5 = stats::quantile(values, 0.975, names = FALSE),
        stringsAsFactors = FALSE
      )
      result_position <- result_position + 1L
    }
  }
  plot_data <- do.call(rbind, result_rows)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$x,
      y = .data$median,
      colour = .data$class,
      fill = .data$class,
      group = .data$class
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper
      ),
      alpha = 0.12,
      colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::labs(
      title = paste(
        "Association contribution by latent class:",
        association_term
      ),
      x = if (k_covariates > 0L) {
        paste(
          c(
            if (covariance_regression$k_sd > 0L) paste0("SD: ", colnames(covariance_regression$x_sd)[1L]) else NULL,
            if (covariance_regression$k_corr > 0L) paste0("correlation: ", colnames(covariance_regression$x_corr)[1L]) else NULL
          ),
          collapse = "; "
        )
      } else {
        "Reference profile"
      },
      y = "Contribution to log hazard",
      colour = "Latent class",
      fill = "Latent class"
    ) +
    theme_fn()
}

#' Posterior associations from a latent-progress mixture
#'
#' @inheritParams assoc.JoiNMeFit
#' @param trajectory If `TRUE`, return class-specific association trajectories
#'   rather than coefficient summaries.
#' @param estimand Class trajectory estimand used when `trajectory = TRUE`.
#' @param longitudinal_times,longitudinal_points Time-grid controls.
#' @param marginal_samples Within-component Monte Carlo values per posterior
#'   draw for a marginal class trajectory.
#' @param class_specific Logical. If `TRUE`, coefficient summaries also carry
#'   compact class-specific association contributions for terms affected by a
#'   fitted class-specific random-effect block.
#' @param class_points Number of reference time or covariance-covariate values
#'   used in the compact class-specific association tables. The plotting method
#'   use its requested full grid.
#'
#' @return A `PosteriorAssoc` object, or a class-trajectory ggplot when
#'   `trajectory = TRUE`. For summary-form mixture output, relevant compact
#'   class contributions are retained in the `class_association` attribute and
#'   printed beneath the common coefficient tables. 
#' @export
assoc.JoiNMeMixFit <- function(
  object,
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE,
  trajectory = FALSE,
  estimand = c("mean_per_class", "marginal_per_class"),
  longitudinal_times = NULL,
  longitudinal_points = 80L,
  marginal_samples = 32L,
  class_specific = TRUE,
  class_points = 3L,
  ...
) {
  .require_mixture_survival(object, "posterior_assoc")
  if (!isTRUE(trajectory)) {
    output <- assoc.JoiNMeFit(
      object,
      draws = draws,
      seed = seed,
      digits = digits,
      summary = summary,
      ...
    )
    attr(output, "metadata")$mixture <- object$mixture %||%
      object$config$mixture
    attr(output, "metadata")$class_association_estimand <- match.arg(
      estimand
    )
    attr(output, "metadata")$class_specific_requested <- isTRUE(
      class_specific
    )
    # Keep a plotting reference only on the mixture-facing return value.  This
    # lets `plot(posterior_assoc(fit))` recover the fitted transforms and class
    # locations without changing the established PosteriorAssoc table layout.
    attr(output, "mixture_fit") <- object
    if (isTRUE(summary) && isTRUE(class_specific)) {
      class_points <- as.integer(
        class_points
      ) # number of compact reference abscissae requested for console reporting
      if (
        length(class_points) != 1L ||
          !is.finite(class_points) ||
          class_points < 2L
      ) {
        cli::cli_abort("{.arg class_points} must be an integer of at least two.")
      }
      class_association <- .mixture_class_association_summary(
        object = object,
        association_output = output,
        estimand = match.arg(estimand),
        draws = draws %||% object$config$draws_default %||% 400L,
        seed = seed,
        digits = digits,
        longitudinal_times = longitudinal_times,
        class_points = class_points,
        marginal_samples = marginal_samples,
        ...
      ) # compact per-class log-hazard contributions for scientifically relevant terms
      attr(output, "class_association") <- class_association
    }
    return(output)
  }
  estimand <- match.arg(estimand)
  plot_object <- .plot_mixture_association_trajectory(
    object,
    estimand = estimand,
    draws = draws %||% object$config$draws_default %||% 400L,
    seed = seed,
    longitudinal_times = longitudinal_times,
    longitudinal_points = longitudinal_points,
    marginal_samples = marginal_samples,
    ci_levels = c(0.5, 0.95),
    marker = NULL,
    theme_fn = ggplot2::theme_bw,
    ...
  )
  plot_object
}

#' Determine whether a fitted class block changes an association estimand
#'
#' @param object A fitted latent-progress mixture.
#' @param association_term One fitted association term or covariance component.
#' @param estimand Either `"mean_per_class"` or `"marginal_per_class"`.
#'
#' @return `TRUE` when a fitted class-specific block enters the requested
#'   association contribution.
#' @keywords internal
#' @noRd
.mixture_class_changes_association <- function(
  object,
  association_term,
  estimand
) {
  mixture <- object$mixture %||%
    object$config$mixture # fitted latent-progress specification
  class_types <- as.character(
    mixture$class_type %||% character(0)
  ) # random-effect blocks assigned latent classes
  term_key <- sub(
    "\\[.*$",
    "",
    association_term
  ) # association channel without a covariance component suffix
  covariance_class <- any(
    class_types %in% c("corr", "vcov")
  ) # whether covariance-regression latent coordinates vary by class

  if (term_key %in% c("corr", "vcov")) {
    return(covariance_class)
  }
  source_channel <- sub(
    "^(cv|cs)_",
    "",
    term_key
  ) # longitudinal source entering the current-value or current-slope channel
  directly_relevant <- switch(
    source_channel,
    mean = "subject" %in% class_types,
    marker = "marker" %in% class_types,
    total = any(c("subject", "marker") %in% class_types),
    FALSE
  ) # class centres which alter the corresponding latent longitudinal source
  marginal_covariance_relevance <- identical(
    estimand,
    "marginal_per_class"
  ) && covariance_class && source_channel %in% c("marker", "total")
  isTRUE(directly_relevant || marginal_covariance_relevance)
}

#' Build compact class-specific posterior association tables
#'
#' @description
#' Global association coefficients are shared across classes.  Their realised
#' contribution to log hazard can nevertheless differ because a class-specific
#' random-effect block changes the longitudinal or covariance feature being
#' multiplied by that coefficient.  
#'
#' @param object A fitted latent-progress mixture.
#' @param association_output The ordinary `PosteriorAssoc` coefficient object.
#' @param estimand Class-trajectory estimand.
#' @param draws Number of posterior draws used for the derived contribution.
#' @param seed Reproducible posterior-subsetting seed.
#' @param digits Number of displayed decimal places.
#' @param longitudinal_times Optional user-supplied time values.
#' @param class_points Number of automatically selected reference values.
#' @param marginal_samples Within-class Monte Carlo values per posterior draw.
#' @param ... Optional association-term selection.
#'
#' @return A named list of posterior summary tables, one for every relevant
#'   fitted association term or covariance component.
#' @keywords internal
#' @noRd
.mixture_class_association_summary <- function(
  object,
  association_output,
  estimand,
  draws,
  seed,
  digits,
  longitudinal_times,
  class_points,
  marginal_samples,
  ...
) {
  dot_arguments <- list(
    ...
  ) # optional association component requested by the caller
  requested_term <- dot_arguments$association_term %||%
    dot_arguments$term # optional single association term or covariance component
  dot_arguments$association_term <- NULL
  dot_arguments$term <- NULL

  association_terms <- unlist(lapply(
    names(association_output),
    function(term_name) {
      if (term_name %in% c("corr", "vcov")) {
        term_table <- association_output[[term_name]] # component-labelled covariance association table
        as.character(term_table$term %||% term_name)
      } else {
        term_name
      }
    }
  ), use.names = FALSE) # every fitted scalar channel or covariance component
  if (!is.null(requested_term)) {
    association_terms <- intersect(
      association_terms,
      as.character(requested_term)
    ) # caller-selected association term, when present
  }
  association_terms <- association_terms[vapply(
    association_terms,
    function(association_term) {
      .mixture_class_changes_association(
        object,
        association_term = association_term,
        estimand = estimand
      )
    },
    logical(1)
  )] # only terms whose realised contribution can differ across fitted classes
  if (length(association_terms) == 0L) {
    return(list())
  }

  reference_times <- if (is.null(longitudinal_times)) {
    NULL
  } else {
    sort(unique(as.numeric(longitudinal_times)))
  } # optional scientifically chosen reference times supplied by the caller
  output <- list() # class-specific contribution table for each relevant term
  for (association_term in association_terms) {
    plot_arguments <- c(
      list(
        object = object,
        estimand = estimand,
        draws = draws,
        seed = seed,
        longitudinal_times = reference_times,
        longitudinal_points = class_points,
        marginal_samples = marginal_samples,
        ci_levels = 0.95,
        marker = NULL,
        theme_fn = ggplot2::theme_bw,
        association_term = association_term
      ),
      dot_arguments
    ) # exact class-trajectory calculation with a compact reference grid
    association_plot <- do.call(
      .plot_mixture_association_trajectory,
      plot_arguments
    ) # ggplot whose data contain the derived posterior contribution summaries
    plot_data <- association_plot$data # numerical values underlying the class curve
    abscissa_name <- if ("time" %in% names(plot_data)) {
      "time"
    } else {
      "x"
    } # time for longitudinal channels or covariance covariate/reference profile
    class_table <- data.frame(
      class = as.character(plot_data$class),
      reference = as.numeric(plot_data[[abscissa_name]]),
      Estimate = as.numeric(plot_data$mean),
      Est.Error = as.numeric(plot_data$sd),
      Q2.5 = as.numeric(plot_data$Q2.5),
      Q97.5 = as.numeric(plot_data$Q97.5),
      Rhat = NA_real_,
      ess_bulk = NA_real_,
      ess_tail = NA_real_,
      stringsAsFactors = FALSE
    ) # common posterior reporting schema for class-specific log-hazard contributions
    names(class_table)[2L] <- abscissa_name
    output[[association_term]] <- .round_summary_table(
      class_table,
      digits = digits
    )
  }
  output
}

#' Plot posterior association output from a latent-progress fit
#'
#' @description
#' `posterior_assoc()` and `assoc()` continue to return their familiar
#' coefficient tables or draw matrices.  When that object came from
#' `joinme_mix()`, this method can additionally display the fitted association
#' contribution along each latent-class trajectory.
#'
#' @param x A `PosteriorAssoc` object returned from a mixture fit.
#' @param trajectory Logical; request the class-specific trajectory display.
#' @param estimand Either `"mean_per_class"` or `"marginal_per_class"`.
#' @param ... Further arguments passed to [association_plot()].
#'
#' @return A ggplot or combined plot returned by [association_plot()].
#' @export
plot.PosteriorAssoc <- function(
  x,
  trajectory = TRUE,
  estimand = c("mean_per_class", "marginal_per_class"),
  ...
) {
  fitted_model <- attr(x, "mixture_fit")
  if (is.null(fitted_model) || !inherits(fitted_model, "JoiNMeMixFit")) {
    cli::cli_abort(c(
      x = "This posterior association object has no fitted mixture attached.",
      i = "Use {.fn association_plot} with the original fitted model."
    ))
  }
  if (!isTRUE(trajectory)) {
    return(association_plot(fitted_model, ...))
  }
  association_plot(
    fitted_model,
    estimand = match.arg(estimand),
    ...
  )
}
