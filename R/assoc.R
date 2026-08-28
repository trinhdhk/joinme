#' Build posterior association effects for a fitted JoiNMe model
#'
#' @description
#' Reconstructs the posterior association effects that enter the survival linear
#' predictor.
#'
#' For weighted current-value and current-slope channels 
#' (`cv_total`, `cs_total`, `cv_marker`, `cs_marker`), 
#' the returned effect is the draw-wise product of the association coefficient and the marker weight,
#' divided by the number of markers to match the scale used in the fitted hazard contribution.
#'
#' For scalar channels (`cv_mean`, `cs_mean`) the returned effect is simply the
#' posterior coefficient. For covariance-style channels (`corr`, `vcov`) the
#' returned effects are grouped by their labelled covariance component. When
#' `summary = TRUE` and an active functional transformation contains a fitted
#' affine shift, the result also contains `affine_shift`. This table begins
#' with `assoc` and `term`; `term` is `"intercept"` or `"slope"`. A row is
#' included only when that role was fitted for the corresponding active
#' association channel.
#'
#' @param object A `JoiNMeFit` object.
#' @param draws Optional number of posterior draws to retain.
#' @param seed Random seed used when subsetting posterior draws.
#' @param digits Number of digits used when `summary = TRUE`.
#' @param summary Logical. 
#' If `TRUE`, return posterior summaries. 
#' If `FALSE`, return raw MCMC sample matrices.
#' @param ... Unused.
#'
#' @return A named list with class `PosteriorAssoc`. Each association element
#'   contains either a posterior summary table or a draws-by-term matrix. In
#'   summary form, the optional `affine_shift` element contains the fitted
#'   transformation intercepts and slopes in the same posterior-summary schema.
#' @export
assoc <- function(object, ...) {
  UseMethod("assoc")
}

#' @rdname assoc
#' @export
posterior_assoc <- function(object, ...) {
  assoc(object, ...)
}

#' Describe the fitted affine shifts within association transformations
#'
#' @description
#' Functional association transformations may contain a fitted affine map,
#' \eqn{\iota_0 + \iota_1 x}, inside one or more nonlinear functions. This
#' helper identifies only the coefficients that were requested for association
#' channels present in the fitted survival model. It also translates the saved
#' Stan positions into the association and coefficient labels used in public
#' posterior tables.
#'
#' A covariance or correlation transformation is evaluated separately for
#' every applicable lower-triangular component. When a functional expression
#' contains more than one fitted nonlinear node, the node number is appended to
#' the association label so that repeated intercepts or slopes remain
#' distinguishable without introducing an internal parameter-name column.
#'
#' @param stan_data The Stan data retained in a fitted `JoiNMeFit` object.
#' @param available_variables Character vector of posterior variable names.
#'
#' @return A data frame with the posterior variable name, association label,
#'   and scientific term (`"intercept"` or `"slope"`), or `NULL` when no
#'   applicable affine shift was fitted.
#' @keywords internal
#' @noRd
.association_affine_shift_layout <- function(
  stan_data,
  available_variables
) {
  association_names <- c(
    "cv_total",
    "cv_mean",
    "cv_marker",
    "cs_total",
    "cs_mean",
    "cs_marker",
    "corr",
    "vcov"
  ) # public association order used throughout fitting and reporting
  available_variables <- as.character(
    available_variables %||% character(0)
  ) # saved posterior quantities that can genuinely be reported
  rows <- list() # one compact description for every fitted affine coefficient

  # Consider each public association channel independently. An affine flag on
  # an inactive channel is not an estimand of the fitted survival model and is
  # therefore deliberately excluded even if a similarly named quantity is
  # present in a synthetic or partially reconstructed object.
  for (association_name in association_names) {
    association_is_active <- isTRUE(
      as.integer(stan_data[[paste0("assoc_", association_name)]] %||% 0L) == 1L
    ) # whether this transformed feature contributes to the fitted log hazard
    if (!association_is_active) {
      next
    }

    association_map <- .assoc_channel_map(
      association_name
    ) # fitted affine names and counts for this association channel
    if (is.null(association_map)) {
      next
    }

    number_components <- if (association_name %in% c("corr", "vcov")) {
      .assoc_transform_component_count(
        association_name,
        stan_data$Q_idm,
        diagonal_only = identical(association_name, "vcov") &&
          isTRUE(as.integer(stan_data$indep_idmarker_cov %||% 0L) == 1L)
      )
    } else {
      1L
    } # number of separately transformed covariance features, or one scalar feature
    if (number_components < 1L) {
      next
    }

    component_labels <- if (association_name %in% c("corr", "vcov")) {
      .assoc_component_display_labels(
        association_name,
        stan_data,
        n_components = number_components
      )
    } else {
      association_name
    } # scientific association labels, including covariance-basis terms

    # Form one role at a time because an expression may estimate only its
    # intercept, only its slope, or different numbers of each. A positive count
    # is the authoritative declaration that a coefficient was fitted.
    append_role <- function(role, variable_prefix, number_nodes) {
      number_nodes <- as.integer(
        number_nodes %||% 0L
      ) # fitted nonlinear nodes carrying this coefficient role
      if (number_nodes < 1L || is.null(variable_prefix)) {
        return(NULL)
      }

      number_coefficients <- number_components * number_nodes
      expected_variables <- paste0(
        variable_prefix,
        "[",
        seq_len(number_coefficients),
        "]"
      ) # flattened component-major order used by the Stan transformation
      if (
        number_coefficients == 1L &&
          !(expected_variables[[1L]] %in% available_variables) &&
          variable_prefix %in% available_variables
      ) {
        expected_variables <- variable_prefix
      } # scalar naming admitted for lightweight fitted-object representations

      retained_positions <- which(
        expected_variables %in% available_variables
      ) # fitted and saved coefficients only
      if (length(retained_positions) == 0L) {
        return(NULL)
      }

      component_index <- ((retained_positions - 1L) %/% number_nodes) + 1L
      node_index <- ((retained_positions - 1L) %% number_nodes) + 1L
      association_label <- component_labels[component_index]
      if (number_nodes > 1L) {
        association_label <- paste0(
          association_label,
          " (node ",
          node_index,
          ")"
        )
      }

      data.frame(
        variable = expected_variables[retained_positions],
        assoc = association_label,
        term = rep(role, length(retained_positions)),
        stringsAsFactors = FALSE
      )
    }

    intercept_rows <- append_role(
      role = "intercept",
      variable_prefix = association_map$iota_intercept,
      number_nodes = stan_data[[association_map$n_iota_intercept]] %||% 0L
    ) # fitted additive shifts for this active association
    slope_rows <- append_role(
      role = "slope",
      variable_prefix = association_map$iota_slope,
      number_nodes = stan_data[[association_map$n_iota_slope]] %||% 0L
    ) # fitted multipliers within the nonlinear transformation
    channel_rows <- Filter(
      Negate(is.null),
      list(intercept_rows, slope_rows)
    )
    if (length(channel_rows) > 0L) {
      rows[[length(rows) + 1L]] <- do.call(rbind, channel_rows)
    }
  }

  if (length(rows) == 0L) {
    return(NULL)
  }
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

#' Label fitted affine-shift draws for posterior interfaces
#'
#' @description
#' Forms a unique, readable name from the association channel and affine
#' coefficient role. The same label is used by `extract()` and
#' `posterior_draws()`, including the complete renamed draw collection, so a
#' draw has one public name irrespective of the route used to obtain it.
#'
#' @param association Character vector of association labels.
#' @param term Character vector containing `"intercept"` or `"slope"`.
#'
#' @return A character vector of labels of the form
#'   `"affine_shift[association, term]"`.
#' @keywords internal
#' @noRd
.association_affine_shift_draw_label <- function(association, term) {
  paste0(
    "affine_shift[",
    as.character(association),
    ", ",
    as.character(term),
    "]"
  )
}

#' Summarise fitted affine shifts within association transformations
#'
#' @description
#' Builds the common posterior table used by `summary()`,
#' `posterior_summary()`, `assoc()`, and `posterior_assoc()`. The calculation is
#' centralised so all four interfaces apply the same association filtering,
#' labels, interval definition, rounding, and sampling diagnostics.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param available_variables Character vector of saved posterior variable
#'   names.
#' @param draws Optional number of posterior draws used in the summary.
#' @param seed Random seed used when posterior draws are reduced.
#' @param digits Number of decimal places for posterior location and interval
#'   summaries.
#'
#' @return A data frame beginning with `assoc` and `term`, followed by
#'   `Estimate`, `Est.Error`, `Q2.5`, `Q97.5`, `Rhat`, `ess_bulk`, and
#'   `ess_tail`; or `NULL` when no applicable affine coefficient was fitted.
#' @keywords internal
#' @noRd
.association_affine_shift_summary <- function(
  object,
  available_variables,
  draws = NULL,
  seed = 1,
  digits = 3
) {
  coefficient_layout <- .association_affine_shift_layout(
    stan_data = object$stan_data,
    available_variables = available_variables
  ) # exact posterior variables and their scientific association labels
  if (is.null(coefficient_layout) || nrow(coefficient_layout) == 0L) {
    return(NULL)
  }

  posterior_table <- as.data.frame(.summarise_draws_diag(
    object$fit,
    variables = coefficient_layout$variable,
    draws = draws,
    seed = seed
  )) # posterior location, uncertainty, interval, and chain diagnostics
  layout_position <- match(
    posterior_table$variable,
    coefficient_layout$variable
  ) # preserve the actual posterior ordering without joining repeated labels
  posterior_table$assoc <- coefficient_layout$assoc[layout_position]
  posterior_table$term <- coefficient_layout$term[layout_position]
  posterior_table <- posterior_table[, c(
    "assoc",
    "term",
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat",
    "ess_bulk",
    "ess_tail"
  ), drop = FALSE]

  # Use the same numerical presentation as the other model and association
  # tables whilst leaving effective sample sizes on their natural scale.
  posterior_table$Estimate <- round(posterior_table$Estimate, digits)
  posterior_table$Est.Error <- round(posterior_table$Est.Error, digits)
  posterior_table$Q2.5 <- round(posterior_table$Q2.5, digits)
  posterior_table$Q97.5 <- round(posterior_table$Q97.5, digits)
  posterior_table$Rhat <- round(posterior_table$Rhat, 3L)
  rownames(posterior_table) <- NULL
  posterior_table
}

#' @rdname assoc
#' @export
assoc.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  
  sd <- object$stan_data
  fit <- object$fit
  if (is.null(draws)) {
    draws <- object$config$draws_default
  }
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- if (isTRUE(summary)) {
    paste0("posterior_assoc_summary_", draws, "_", digits)
  } else {
    paste0("posterior_assoc_draws_", draws)
  }
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  marker_levels <- as.character(sd$marker_levels %||% paste0("marker_", seq_len(as.integer(sd$D %||% 0L))))
  n_markers <- length(marker_levels)
  weighted_divisor <- if (n_markers > 0L) n_markers else 1L

  build_scalar_term <- function(var_name, label) {
    if (is.null(var_name)) {
      return(NULL)
    }
    draw_arr <- .get_draws_array(fit, variables = var_name, draws = draws, seed = seed)
    if (isTRUE(summary)) {
      .assoc_summary_from_draw_array(draw_arr, term_labels = label, digits = digits)
    } else {
      .assoc_matrix_from_draw_array(draw_arr, term_labels = label)
    }
  }

  build_weighted_term <- function(var_name, term_key) {
    weight_array <- .association_marker_weight_array(object, term_key = term_key, draws = draws, seed = seed, all_vars = all_vars)
    if (is.null(var_name) || is.null(weight_array) || n_markers == 0L) {
      return(NULL)
    }
    alpha_arr <- .get_draws_array(fit, variables = var_name, draws = draws, seed = seed)
    effect_arr <- array(
      NA_real_,
      dim = c(dim(alpha_arr)[1], dim(alpha_arr)[2], n_markers),
      dimnames = list(
        iteration = dimnames(alpha_arr)[[1]],
        chain = dimnames(alpha_arr)[[2]],
        variable = marker_levels
      )
    )

    for (marker_index in seq_len(n_markers)) {
      effect_arr[, , marker_index] <- alpha_arr[, , 1] * weight_array[, , marker_index] / weighted_divisor
    }

    if (isTRUE(summary)) {
      .assoc_summary_from_draw_array(effect_arr, term_labels = marker_levels, digits = digits)
    } else {
      .assoc_matrix_from_draw_array(effect_arr, term_labels = marker_levels)
    }
  }

  build_vector_term <- function(var_names, term_labels, component_meta = NULL) {
    var_names <- as.character(var_names %||% character(0))
    if (!length(var_names)) {
      return(NULL)
    }
    draw_arr <- .get_draws_array(fit, variables = var_names, draws = draws, seed = seed)
    if (isTRUE(summary)) {
      out_tbl <- .assoc_summary_from_draw_array(draw_arr, term_labels = term_labels, digits = digits)
      if (!is.null(component_meta) && !is.null(out_tbl)) {
        idx <- match(out_tbl$term, component_meta$term)
        out_tbl$component <- component_meta$component[idx]
        out_tbl$row <- component_meta$row[idx]
        out_tbl$col <- component_meta$col[idx]
        out_tbl <- out_tbl[, c(
          "component",
          "row",
          "col",
          "term",
          "Estimate",
          "Est.Error",
          "Q2.5",
          "Q97.5",
          "Rhat",
          "ess_bulk",
          "ess_tail"
        ), drop = FALSE]
      }
      out_tbl
    } else {
      .assoc_matrix_from_draw_array(draw_arr, term_labels = term_labels)
    }
  }

  corr_assoc_vars <- grep("^alpha_corr\\[", all_vars, value = TRUE)
  vcov_assoc_vars <- grep("^alpha_vcov\\[", all_vars, value = TRUE)

  out <- list(
    cv_total = if (isTRUE(sd$assoc_cv_total == 1L)) {
      build_weighted_term("alpha_cv_total", term_key = "cv_total")
    } else NULL,
    cv_mean = if (isTRUE(sd$assoc_cv_mean == 1L)) {
      build_scalar_term("alpha_cv_mean", "cv_mean")
    } else NULL,
    cv_marker = if (isTRUE(sd$assoc_cv_marker == 1L)) {
      build_weighted_term("alpha_cv_marker", term_key = "cv_marker")
    } else NULL,
    cs_total = if (isTRUE(sd$assoc_cs_total == 1L)) {
      build_weighted_term("alpha_cs_total", term_key = "cs_total")
    } else NULL,
    cs_mean = if (isTRUE(sd$assoc_cs_mean == 1L)) {
      build_scalar_term("alpha_cs_mean", "cs_mean")
    } else NULL,
    cs_marker = if (isTRUE(sd$assoc_cs_marker == 1L)) {
      build_weighted_term("alpha_cs_marker", term_key = "cs_marker")
    } else NULL,
    corr = if (isTRUE(sd$assoc_corr == 1L)) {
      corr_meta <- .assoc_component_metadata("corr", sd, length(corr_assoc_vars))
      build_vector_term(
        corr_assoc_vars,
        term_labels = corr_meta$term %||% .assoc_component_display_labels("corr", sd, length(corr_assoc_vars)),
        component_meta = corr_meta
      )
    } else NULL,
    vcov = if (isTRUE(sd$assoc_vcov == 1L)) {
      vcov_meta <- .assoc_component_metadata("vcov", sd, length(vcov_assoc_vars))
      build_vector_term(
        vcov_assoc_vars,
        term_labels = vcov_meta$term %||% .assoc_component_display_labels("vcov", sd, length(vcov_assoc_vars)),
        component_meta = vcov_meta
      )
    } else NULL
  )
  out <- out[!vapply(out, is.null, logical(1))]

  if (!length(out)) {
    cli::cli_abort(c(
      x = "No association effect is available in this model.",
      i = "Fit a model with association terms such as {.val cv_total}, {.val cv_mean}, {.val corr}, or {.val vcov}."
    ))
  }

  # Add the fitted affine transformation coefficients after establishing that
  # the model contains at least one association effect. Raw association output
  # remains a collection of effect matrices; complete unsummarised affine
  # draws remain available through posterior_draws().
  if (isTRUE(summary)) {
    affine_shift <- .association_affine_shift_summary(
      object = object,
      available_variables = all_vars,
      draws = draws,
      seed = seed,
      digits = digits
    ) # only active channels and explicitly fitted intercept or slope roles
    if (!is.null(affine_shift)) {
      out$affine_shift <- affine_shift
    }
  }

  out <- structure(
    out,
    class = "PosteriorAssoc",
    metadata = list(
      summary = isTRUE(summary),
      draws = draws,
      digits = digits,
      marker_levels = marker_levels
    )
  )

  object$cache_set(cache_key, out)
  out
}


#' Print posterior association effects
#'
#' @description
#' Prints association-effect summaries or raw posterior sample availability in a
#' layout aligned with [print.summary_JoiNMeFit()].
#'
#' @param x A `PosteriorAssoc` object returned by [assoc()] or
#'   [posterior_assoc()].
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.PosteriorAssoc <- function(x, ...) {
  meta <- attr(x, "metadata") %||% list()

  .cli_summary_heading("Posterior association effects", level = 1L)
  if (!is.null(meta$draws)) {
    .cli_print_bullets(c(
      paste0("Posterior draws: ", meta$draws)
    ))
  }

  if (!length(x)) {
    return(invisible(x))
  }

  association_names <- setdiff(
    names(x),
    "affine_shift"
  ) # survival association effects, excluding the separate affine-shift table
  for (term_name in association_names) {
    .cli_summary_heading(paste0("Association term: ", term_name), level = 2L)
    if (isTRUE(meta$summary)) {
      term_tbl <- x[[term_name]]
      if (term_name %in% c("corr", "vcov") && is.data.frame(term_tbl) && all(c("row", "col") %in% names(term_tbl))) {
        .cli_print_bullets("Displayed by matrix row using fitted random-effect term labels.")
        matrix_groups <- .posterior_assoc_matrix_groups(term_tbl)
        for (row_label in names(matrix_groups)) {
          .cli_summary_heading(paste0("row = ", row_label), level = 3L)
          .cli_print_table(matrix_groups[[row_label]])
        }
      } else {
        .cli_print_table(.posterior_assoc_display_table(term_tbl))
      }
    } else {
      draw_block <- x[[term_name]]
      desc <- data.frame(
        term = colnames(draw_block) %||% character(0),
        draws = rep(nrow(draw_block), ncol(draw_block)),
        stringsAsFactors = FALSE
      )
      .cli_print_table(desc)
    }
  }

  if (isTRUE(meta$summary) && !is.null(x$affine_shift)) {
    .cli_print_table_section(
      "Association transformation affine shifts (iota)",
      .posterior_assoc_display_table(x$affine_shift),
      level = 2L
    )
  }

  class_association <- attr(x, "class_association") %||%
    list() # class-specific hazard-contribution summaries attached by joinme_mix()
  if (length(class_association) > 0L) {
    .cli_summary_heading(
      "Association contributions by latent class",
      level = 2L
    )
    .cli_print_bullets(c(
      paste0(
        "Class estimand: ",
        meta$class_association_estimand %||% "mean_per_class"
      ),
      "Only association terms affected by a fitted class-specific random-effect block are shown.",
      "The reported reference points span the fitted time or covariance-covariate range; plot() displays the full curve."
    ))
    for (association_term in names(class_association)) {
      .cli_print_table_section(
        paste0("Class association: ", association_term),
        class_association[[association_term]],
        level = 3L
      )
    }
  } else if (
    !is.null(meta$mixture) &&
      isTRUE(meta$summary) &&
      isTRUE(meta$class_specific_requested)
  ) {
    .cli_print_bullets(
      "No fitted class-specific block changes the displayed association term; the coefficient therefore remains common across classes."
    )
  }

  invisible(x)
}

#' Select the common posterior columns used for association printing
#'
#' @description
#' Association summaries use the same compact schema as the ordinary model
#' summary: posterior mean (`Estimate`), posterior standard deviation
#' (`Est.Error`), central interval and sampling diagnostics. Duplicate
#' `Mean`, `Median`, and `SD` columns are not reported.
#'
#' @param table A posterior association summary table.
#'
#' @return A data frame prepared for concise console printing.
#' @keywords internal
#' @noRd
.posterior_assoc_display_table <- function(table) {
  if (is.null(table) || !is.data.frame(table)) {
    return(table)
  }
  metric_columns <- c(
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat",
    "ess_bulk",
    "ess_tail"
  ) # explicit posterior location, spread, interval and simulation diagnostics
  identifier_columns <- setdiff(
    names(table),
    c(metric_columns, "Mean", "Median", "SD")
  ) # association labels and covariance-component identifiers
  table[, c(
    identifier_columns,
    intersect(metric_columns, names(table))
  ), drop = FALSE]
}

#' Format covariance-style posterior association summaries for printing
#'
#' @description
#' `corr` and `vcov` posterior association summaries are easier to read when the
#' lower-triangular components are displayed in row-wise groups rather than as a
#' single flat list. This helper keeps the returned object unchanged and only
#' prepares grouped tables for console output.
#'
#' @param tbl Posterior association summary table.
#'
#' @return A named list of data frames, one per matrix row label.
#' @keywords internal
#' @noRd
.posterior_assoc_matrix_groups <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(list())
  }
  if (!all(c("row", "col") %in% names(tbl))) {
    return(list())
  }

  row_levels <- unique(as.character(tbl$row))
  row_levels <- row_levels[!is.na(row_levels) & nzchar(row_levels)]
  groups <- lapply(row_levels, function(row_label) {
    block <- tbl[as.character(tbl$row) == row_label, , drop = FALSE]
    block <- block[, c(
      intersect(c(
        "col",
        "Estimate",
        "Est.Error",
        "Q2.5",
        "Q97.5",
        "Rhat",
        "ess_bulk",
        "ess_tail"
      ), names(block))
    ), drop = FALSE]
    names(block)[names(block) == "col"] <- "term"
    block
  })
  names(groups) <- row_levels
  groups
}


#' Summarise a derived posterior draw array with MCMC diagnostics
#'
#' @description
#' Derived association effects do not necessarily exist as named Stan variables.
#' This helper converts an iteration x chain x term array into the same summary
#' schema used elsewhere in JoiNMe: posterior mean (`Estimate`), posterior
#' standard deviation (`Est.Error`), central 95% interval, split-chain R-hat,
#' and bulk/tail effective sample sizes.
#'
#' @param draw_array Numeric array with dimensions iteration x chain x term.
#' @param term_labels Character vector naming the third dimension.
#' @param digits Number of decimal places used for posterior location and
#'   interval summaries.
#'
#' @return A data frame containing `Estimate`, `Est.Error`, `Q2.5`, `Q97.5`,
#'   `Rhat`, `ess_bulk`, and `ess_tail`, preceded by the term label.
#' @keywords internal
#' @noRd
.assoc_summary_from_draw_array <- function(draw_array, term_labels, digits = 3) {
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(NULL)
  }

  number_terms <- dim(
    draw_array
  )[3] # number of posterior quantities stored along the variable dimension
  term_labels <- as.character(
    term_labels %||% paste0("term_", seq_len(number_terms))
  ) # user-facing labels, which may legitimately repeat across classes or domains
  if (length(term_labels) != number_terms) {
    cli::cli_abort(
      "{.arg term_labels} must contain one label for every posterior quantity."
    )
  }

  source_variable_names <- dimnames(draw_array)[[3]]
  if (
    is.null(source_variable_names) ||
      length(source_variable_names) != number_terms ||
      any(!nzchar(source_variable_names))
  ) {
    source_variable_names <- paste0(".joinme_summary_", seq_len(number_terms))
  }
  internal_variable_names <- make.unique(
    as.character(source_variable_names),
    sep = "__"
  ) # unique diagnostic keys kept separate from possibly repeated display labels
  dimnames(draw_array) <- list(
    iteration = dimnames(draw_array)[[1]] %||% as.character(seq_len(dim(draw_array)[1])),
    chain = dimnames(draw_array)[[2]] %||% as.character(seq_len(dim(draw_array)[2])),
    variable = internal_variable_names
  )

  draws_obj <- posterior::as_draws_array(draw_array)
  draws_df <- posterior::as_draws_df(draws_obj)
  sum_df <- data.frame(
    term = term_labels,
    do.call(rbind, lapply(
      internal_variable_names,
      function(variable_name) {
        .summarize_draw_col(draws_df[[variable_name]])
      }
    )),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  diag_df <- suppressWarnings(tryCatch(
    posterior::summarise_draws(draws_obj, "rhat", "ess_bulk", "ess_tail"),
    error = function(e) NULL
  ))
  if (!is.null(diag_df)) {
    diag_df <- diag_df[, c("variable", "rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    diagnostic_position <- match(
      internal_variable_names,
      as.character(diag_df$variable)
    ) # one-to-one diagnostic row for each source variable, independent of labels
    sum_df$Rhat <- diag_df$rhat[diagnostic_position]
    sum_df$ess_bulk <- diag_df$ess_bulk[diagnostic_position]
    sum_df$ess_tail <- diag_df$ess_tail[diagnostic_position]
  } else {
    sum_df$Rhat <- NA_real_
    sum_df$ess_bulk <- NA_real_
    sum_df$ess_tail <- NA_real_
  }

  sum_df <- .round_summary_table(sum_df, digits = digits)
  sum_df[, c(
    "term",
    "Estimate",
    "Est.Error",
    "Q2.5",
    "Q97.5",
    "Rhat",
    "ess_bulk",
    "ess_tail"
  ), drop = FALSE]
}

#' Build labels for covariance-style association components
#'
#' @description
#' Converts the internal lower-triangular indexing used for covariance-style
#' association channels into term labels that refer to the underlying random
#' effect basis. This keeps downstream reports aligned with the design-matrix
#' terms seen in model summaries.
#'
#' @param term_key Association channel name. Supported values are `"corr"` and
#'   `"vcov"`.
#' @param sd Stan-data list stored in a `JoiNMeFit` object.
#' @param n_components Optional expected number of components. When supplied,
#'   the output is truncated to this length.
#'
#' @return A character vector of human-readable component labels.
#' @keywords internal
#' @noRd
.assoc_component_display_labels <- function(term_key, sd, n_components = NULL) {
  term_key <- as.character(term_key %||% "")[1]
  q_idm <- as.integer(sd$Q_idm %||% 0L)
  diagonal_only <- identical(term_key, "vcov") && as.integer(sd$indep_idmarker_cov %||% 0L) == 1L
  include_diag <- identical(term_key, "vcov")

  n_expected <- .assoc_transform_component_count(term_key, q_idm, diagonal_only = diagonal_only)
  if (is.null(n_components)) {
    n_components <- n_expected
  }
  n_components <- as.integer(n_components %||% 0L)
  if (n_components <= 0L) {
    return(character(0))
  }

  base_terms <- as.character(sd$zidm_cols %||% paste0("term_", seq_len(max(q_idm, 1L))))
  if (length(base_terms) < q_idm) {
    base_terms <- c(base_terms, paste0("term_", seq.int(length(base_terms) + 1L, q_idm)))
  }
  base_terms <- base_terms[seq_len(max(q_idm, 1L))]

  feature_map <- .assoc_cov_feature_map(
    q_idm = q_idm,
    diagonal_only = diagonal_only,
    include_diag = include_diag
  )
  if (is.null(dim(feature_map)) || nrow(feature_map) == 0L) {
    return(.assoc_transform_component_labels(term_key, n_components))
  }

  feature_map <- feature_map[seq_len(min(nrow(feature_map), n_components)), , drop = FALSE]
  row_idx <- if (!is.null(colnames(feature_map)) && "row" %in% colnames(feature_map)) {
    feature_map[, "row"]
  } else {
    feature_map[, 1]
  }
  col_idx <- if (!is.null(colnames(feature_map)) && "col" %in% colnames(feature_map)) {
    feature_map[, "col"]
  } else {
    feature_map[, 2]
  }
  paste0(
    term_key,
    "[",
    base_terms[row_idx],
    ", ",
    base_terms[col_idx],
    "]"
  )
}


#' Build explicit metadata for covariance-style association components
#'
#' @description
#' Covariance-style association channels are indexed over lower-triangular
#' coordinates of the marker-by-id random-effect basis. This helper reconstructs
#' those coordinates using the fitted random-effect term labels so posterior
#' displays can show both a compact combined term and the underlying row and
#' column terms.
#'
#' @param term_key Association channel name. Supported values are `"corr"` and
#'   `"vcov"`.
#' @param sd Stan-data list stored in the fitted object.
#' @param n_components Number of covariance-style components to report.
#'
#' @return A data frame with columns `component`, `row`, `col`, and `term`.
#' @keywords internal
#' @noRd
.assoc_component_metadata <- function(term_key, sd, n_components) {
  term_key <- as.character(term_key %||% "")[1]
  n_components <- as.integer(n_components %||% 0L)
  if (!(term_key %in% c("corr", "vcov")) || n_components <= 0L) {
    return(NULL)
  }

  q_idm <- as.integer(sd$Q_idm %||% 0L)
  diagonal_only <- identical(term_key, "vcov") && as.integer(sd$indep_idmarker_cov %||% 0L) == 1L
  feature_map <- .assoc_cov_feature_map(
    q_idm = q_idm,
    diagonal_only = diagonal_only,
    include_diag = identical(term_key, "vcov")
  )
  if (is.null(dim(feature_map)) || nrow(feature_map) == 0L) {
    return(NULL)
  }

  feature_map <- feature_map[seq_len(min(nrow(feature_map), n_components)), , drop = FALSE]
  row_idx <- if (!is.null(colnames(feature_map)) && "row" %in% colnames(feature_map)) {
    feature_map[, "row"]
  } else {
    feature_map[, 1]
  }
  col_idx <- if (!is.null(colnames(feature_map)) && "col" %in% colnames(feature_map)) {
    feature_map[, "col"]
  } else {
    feature_map[, 2]
  }
  base_terms <- as.character(sd$zidm_cols %||% paste0("term_", seq_len(max(q_idm, 1L))))
  if (length(base_terms) < q_idm) {
    base_terms <- c(base_terms, paste0("term_", seq.int(length(base_terms) + 1L, q_idm)))
  }
  base_terms <- base_terms[seq_len(max(q_idm, 1L))]

  data.frame(
    component = seq_len(nrow(feature_map)),
    row = base_terms[row_idx],
    col = base_terms[col_idx],
    term = .assoc_component_display_labels(term_key, sd, n_components = nrow(feature_map)),
    stringsAsFactors = FALSE
  )
}

#' Flatten a derived posterior draw array into an MCMC sample matrix
#'
#' @description
#' Returns posterior samples as a draws-by-term matrix so each association term
#' can be inspected directly without an additional summary step.
#'
#' @param draw_array Numeric array with dimensions iteration x chain x term.
#' @param term_labels Character vector naming the third dimension.
#'
#' @return A numeric matrix with one column per term label.
#' @keywords internal
#' @noRd
.assoc_matrix_from_draw_array <- function(draw_array, term_labels) {
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(NULL)
  }

  term_labels <- as.character(term_labels %||% paste0("term_", seq_len(dim(draw_array)[3])))
  dimnames(draw_array) <- list(
    iteration = dimnames(draw_array)[[1]] %||% as.character(seq_len(dim(draw_array)[1])),
    chain = dimnames(draw_array)[[2]] %||% as.character(seq_len(dim(draw_array)[2])),
    variable = term_labels
  )

  out <- posterior::as_draws_matrix(posterior::as_draws_array(draw_array))
  out[, term_labels, drop = FALSE]
}
