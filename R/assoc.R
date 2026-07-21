#' Build posterior association effects for a fitted JoiNMe model
#'
#' @description
#' Reconstructs the posterior association effects that enter the survival linear
#' predictor.
#'
#' For weighted current-value and current-slope channels (`cv_total`,
#' `cs_total`, `cv_marker`, `cs_marker`), the returned effect is the draw-wise
#' product of the association coefficient and the marker weight, divided by the
#' number of markers to match the scale used in the fitted hazard contribution.
#'
#' For scalar channels (`cv_mean`, `cs_mean`) the returned effect is simply the
#' posterior coefficient. For covariance-style channels (`corr`, `vcov`) the
#' returned effects are grouped by their labelled covariance component.
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
#' @return A named list with class `PosteriorAssoc`. Each list element contains
#'   either a posterior summary table or a draws-by-term matrix for one
#'   association channel.
#' @export
assoc <- function(object, ...) {
  UseMethod("assoc")
}

#' @rdname assoc
#' @export
posterior_assoc <- function(object, ...) {
  assoc(object, ...)
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
        out_tbl <- out_tbl[, c("component", "row", "col", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
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

  for (term_name in names(x)) {
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
        .cli_print_table(term_tbl)
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

  invisible(x)
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
      intersect(c("col", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), names(block))
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
#' schema used elsewhere in JoiNMe: posterior mean, posterior standard
#' deviation, central 95% interval, split-chain R-hat, and bulk/tail effective
#' sample sizes.
#'
#' @param draw_array Numeric array with dimensions iteration x chain x term.
#' @param term_labels Character vector naming the third dimension.
#' @param digits Number of decimal places used for posterior location and
#'   interval summaries.
#'
#' @return A data frame with columns `term`, `Estimate`, `Est.Error`, `Q2.5`,
#'   `Q97.5`, `Rhat`, `ess_bulk`, and `ess_tail`.
#' @keywords internal
#' @noRd
.assoc_summary_from_draw_array <- function(draw_array, term_labels, digits = 3) {
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(NULL)
  }

  term_labels <- as.character(term_labels %||% paste0("term_", seq_len(dim(draw_array)[3])))
  dimnames(draw_array) <- list(
    iteration = dimnames(draw_array)[[1]] %||% as.character(seq_len(dim(draw_array)[1])),
    chain = dimnames(draw_array)[[2]] %||% as.character(seq_len(dim(draw_array)[2])),
    variable = term_labels
  )

  draws_obj <- posterior::as_draws_array(draw_array)
  draws_df <- posterior::as_draws_df(draws_obj)
  sum_df <- data.frame(
    term = term_labels,
    do.call(rbind, lapply(term_labels, function(label) .summarize_draw_col(draws_df[[label]]))),
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
    names(diag_df) <- c("term", "Rhat", "ess_bulk", "ess_tail")
    sum_df <- merge(sum_df, diag_df, by = "term", all.x = TRUE, sort = FALSE)
  } else {
    sum_df$Rhat <- NA_real_
    sum_df$ess_bulk <- NA_real_
    sum_df$ess_tail <- NA_real_
  }

  .round_summary_table(sum_df, digits = digits)
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

