#' Marker-weighted association term keys
#'
#' @description
#' Internal helper returning the public association terms whose hazard
#' contribution is formed by a marker-specific weighted average.
#'
#' The order is used consistently across standata construction, Stan indexing,
#' simulation, summaries, and prediction. Keeping one canonical order prevents
#' accidental term swapping when weight structures are shared across terms or
#' split by term.
#'
#' @return Character vector of weighted association term keys.
#' @keywords internal
.weighted_assoc_term_keys <- function() {
  c("cv_total", "cs_total", "cv_marker", "cs_marker")
}

#' Weighted association activity flags
#'
#' @param x A standata-like list containing association flags, or a character
#'   vector of active association term keys.
#'
#' @return Named logical vector in canonical weighted-term order.
#' @keywords internal
.weighted_assoc_term_flags <- function(x) {
  keys <- .weighted_assoc_term_keys()

  if (is.character(x)) {
    active <- unique(as.character(x))
    return(stats::setNames(keys %in% active, keys))
  }

  stats::setNames(vapply(keys, function(term_key) {
    isTRUE(as.integer(x[[paste0("assoc_", term_key)]] %||% 0L) == 1L)
  }, logical(1)), keys)
}

#' Active marker-weighted association terms
#'
#' @param x A standata-like list containing association flags, or a character
#'   vector of active association term keys.
#'
#' @return Character vector of active weighted association terms.
#' @keywords internal
.active_weighted_assoc_terms <- function(x) {
  flags <- .weighted_assoc_term_flags(x)
  names(flags)[flags]
}

#' Normalise one marker-weight vector
#'
#' @param weights Numeric vector or scalar.
#' @param marker_levels Character marker labels in fitted order.
#' @param context Short context string for errors.
#' @param arg_name Argument label used in messages.
#'
#' @return Numeric vector with one entry per marker.
#' @keywords internal
.normalise_marker_weight_vector <- function(weights, marker_levels, context, arg_name = "marker_weights") {
  marker_levels <- as.character(marker_levels %||% character(0))
  n_markers <- length(marker_levels)

  if (n_markers <= 0L) {
    return(numeric(0))
  }

  if (!is.numeric(weights)) {
    cli::cli_abort(c(
      x = "{.arg {arg_name}} must be numeric in {context}.",
      i = "Provide either one scalar weight, one weight per marker, or a named list with one numeric vector per weighted association term."
    ))
  }

  if (!is.null(names(weights))) {
    weights <- weights[marker_levels]
  }

  weights <- as.numeric(weights)

  if (length(weights) == 1L) {
    weights <- rep(weights, n_markers)
  }
  if (length(weights) != n_markers) {
    cli::cli_abort(c(
      x = "{.arg {arg_name}} must have length {n_markers} in {context}.",
      i = "The expected marker order is: {.val {paste(marker_levels, collapse = ', ')}}."
    ))
  }
  if (any(!is.finite(weights))) {
    cli::cli_abort(c(
      x = "{.arg {arg_name}} must contain only finite values in {context}.",
      i = "Marker weights may be positive or negative, but they cannot be NA, NaN, or infinite."
    ))
  }

  unname(weights)
}

#' Resolve shared or term-specific marker-weight structures
#'
#' @description
#' Converts the user-supplied marker-weight specification into a term-indexed
#' internal representation. This representation is used by both fitting and
#' simulation so the data-generating process and the fitted hazard use the same
#' mapping between association terms and marker weights.
#'
#' The algorithm is:
#' 1. identify the active weighted association terms,
#' 2. normalise each supplied marker-weight vector into fitted marker order,
#' 3. either map all active terms to one shared weight set or assign one weight
#'    set per active term,
#' 4. return both per-term vectors and the compact set-index representation used
#'    by Stan.
#'
#' @param marker_weights User-supplied marker-weight specification.
#' @param marker_levels Character marker labels in fitted order.
#' @param active_terms Character vector of active weighted association terms.
#' @param shared_marker_weights Logical. If `TRUE`, all active weighted
#'   association terms share one marker-weight structure. If `FALSE`, each active
#'   weighted association term gets its own marker-weight structure.
#' @param estimate_marker_weights Logical. If `TRUE`, the default base weights
#'   are zero so the latent perturbation represents the full effective weight. If
#'   `FALSE`, the default base weights are one.
#' @param context Short context string for errors.
#'
#' @return A named list containing per-term vectors, set indices, and the compact
#'   matrix representation used by Stan.
#' @keywords internal
.resolve_marker_weight_structure <- function(marker_weights,
                                             marker_levels,
                                             active_terms,
                                             shared_marker_weights = TRUE,
                                             estimate_marker_weights = FALSE,
                                             context = "joinme") {
  keys <- .weighted_assoc_term_keys()
  marker_levels <- as.character(marker_levels %||% character(0))
  n_markers <- length(marker_levels)
  active_terms <- intersect(keys, unique(as.character(active_terms %||% character(0))))

  if (!is.logical(shared_marker_weights) || length(shared_marker_weights) != 1L || is.na(shared_marker_weights)) {
    cli::cli_abort(c(
      x = "{.arg shared_marker_weights} must be TRUE or FALSE in {context}.",
      i = "Use TRUE to share one weight structure across weighted association terms, or FALSE to estimate or fix one weight structure per active weighted term."
    ))
  }

  default_vec <- rep(if (isTRUE(estimate_marker_weights)) 0 else 1, n_markers)
  per_term <- stats::setNames(vector("list", length(keys)), keys)
  for (term_key in keys) {
    per_term[[term_key]] <- default_vec
  }

  if (!length(active_terms)) {
    return(list(
      shared_marker_weights = isTRUE(shared_marker_weights),
      active_terms = character(0),
      n_sets = 1L,
      set_index = stats::setNames(integer(length(keys)), keys),
      base_by_term = per_term,
      base_matrix = matrix(default_vec, nrow = 1L, byrow = TRUE,
        dimnames = list("inactive", marker_levels))
    ))
  }

  if (is.null(marker_weights)) {
    resolved_shared <- default_vec
    resolved_by_term <- NULL
  } else if (is.numeric(marker_weights)) {
    resolved_shared <- .normalise_marker_weight_vector(
      marker_weights,
      marker_levels = marker_levels,
      context = context,
      arg_name = "marker_weights"
    )
    resolved_by_term <- NULL
  } else if (is.matrix(marker_weights) || is.data.frame(marker_weights)) {
    resolved_by_term <- as.data.frame(marker_weights, check.names = FALSE, stringsAsFactors = FALSE)
    resolved_shared <- NULL
  } else if (is.list(marker_weights)) {
    resolved_by_term <- marker_weights
    resolved_shared <- NULL
  } else {
    cli::cli_abort(c(
      x = "{.arg marker_weights} has an unsupported type in {context}.",
      i = "Use NULL, one numeric vector, a numeric matrix/data frame with weighted-term rows, or a named list keyed by weighted association term."
    ))
  }

  if (!is.null(resolved_by_term)) {
    extract_term_weights <- function(term_key) {
      if (is.data.frame(resolved_by_term)) {
        row_names <- rownames(resolved_by_term) %||% character(0)
        if (!(term_key %in% row_names)) {
          return(NULL)
        }
        return(as.numeric(resolved_by_term[term_key, , drop = TRUE]))
      }
      resolved_by_term[[term_key]]
    }

    supplied_terms <- intersect(names(resolved_by_term) %||% rownames(resolved_by_term) %||% character(0), keys)
    if (isTRUE(shared_marker_weights)) {
      supplied_vecs <- lapply(active_terms, extract_term_weights)
      supplied_vecs <- Filter(Negate(is.null), supplied_vecs)
      if (!length(supplied_vecs)) {
        resolved_shared <- default_vec
      } else {
        normalised <- lapply(seq_along(supplied_vecs), function(i) {
          .normalise_marker_weight_vector(
            supplied_vecs[[i]],
            marker_levels = marker_levels,
            context = context,
            arg_name = paste0("marker_weights$", active_terms[min(i, length(active_terms))])
          )
        })
        ref_vec <- normalised[[1]]
        if (length(normalised) > 1L) {
          mismatch <- vapply(normalised[-1], function(x) max(abs(x - ref_vec)) > 1e-8, logical(1))
          if (any(mismatch)) {
            cli::cli_abort(c(
              x = "{.arg marker_weights} supplies different weight vectors while {.arg shared_marker_weights = TRUE} in {context}.",
              i = "Set {.arg shared_marker_weights = FALSE} to use different marker-weight structures across weighted association terms."
            ))
          }
        }
        resolved_shared <- ref_vec
      }
    } else {
      for (term_key in active_terms) {
        term_weights <- extract_term_weights(term_key)
        per_term[[term_key]] <- if (is.null(term_weights)) {
          default_vec
        } else {
          .normalise_marker_weight_vector(
            term_weights,
            marker_levels = marker_levels,
            context = context,
            arg_name = paste0("marker_weights$", term_key)
          )
        }
      }
    }

    unused_terms <- setdiff(supplied_terms, active_terms)
    if (length(unused_terms) > 0L) {
      cli::cli_warn(c(
        x = "Ignoring marker weights for inactive weighted association term{?s} in {context}: {.val {paste(unused_terms, collapse = ', ')}}.",
        i = "Only active weighted association terms contribute to the fitted or simulated hazard."
      ))
    }
  }

  if (!is.null(resolved_shared)) {
    for (term_key in active_terms) {
      per_term[[term_key]] <- resolved_shared
    }
  }

  set_index <- stats::setNames(integer(length(keys)), keys)
  if (isTRUE(shared_marker_weights)) {
    set_index[active_terms] <- 1L
    base_matrix <- matrix(per_term[[active_terms[[1L]]]], nrow = 1L, byrow = TRUE,
      dimnames = list("shared", marker_levels))
    n_sets <- 1L
  } else {
    n_sets <- length(active_terms)
    base_matrix <- matrix(0, nrow = n_sets, ncol = n_markers,
      dimnames = list(active_terms, marker_levels))
    for (set_pos in seq_along(active_terms)) {
      term_key <- active_terms[[set_pos]]
      set_index[[term_key]] <- set_pos
      base_matrix[set_pos, ] <- per_term[[term_key]]
    }
  }

  list(
    shared_marker_weights = isTRUE(shared_marker_weights),
    active_terms = active_terms,
    n_sets = as.integer(n_sets),
    set_index = set_index,
    base_by_term = per_term,
    base_matrix = base_matrix
  )
}

#' Marker-weight summary label
#'
#' @param term_key Weighted association term key.
#' @param marker_label Marker label.
#' @param shared_marker_weights Logical. If `TRUE`, omit the term key because all
#'   weighted association terms share one structure.
#'
#' @return Character label used in summaries and extraction maps.
#' @keywords internal
.marker_weight_summary_label <- function(term_key, marker_label, shared_marker_weights = TRUE) {
  if (isTRUE(shared_marker_weights)) {
    return(paste0("weight: ", marker_label))
  }
  paste0("weight[", term_key, "]: ", marker_label)
}

#' Marker-weight Stan variable prefix
#'
#' @param term_key Weighted association term key.
#' @param effective Logical; if `TRUE`, return the effective posterior variable
#'   prefix, otherwise return the base-weight prefix.
#'
#' @return Character scalar prefix.
#' @keywords internal
.marker_weight_var_prefix <- function(term_key, effective = TRUE) {
  suffix <- if (isTRUE(effective)) "eff_" else ""
  paste0("marker_weights_", suffix, term_key)
}

#' Fetch marker-weight draw array for one association term
#'
#' @param object A `JoinMeFit` object.
#' @param term_key Weighted association term key.
#' @param draws Optional posterior draw subset size.
#' @param seed Random seed used when subsetting draws.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return Iteration x chain x marker draw array, or `NULL` if unavailable.
#' @keywords internal
.association_marker_weight_array <- function(object, term_key, draws = NULL, seed = 1, all_vars = NULL) {
  sd <- object$stan_data
  n_markers <- as.integer(sd$D %||% 0L)
  if (n_markers <= 0L) {
    return(NULL)
  }

  if (!(term_key %in% .weighted_assoc_term_keys())) {
    return(NULL)
  }

  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(e) character(0))
  }

  eff_prefix <- .marker_weight_var_prefix(term_key, effective = TRUE)
  eff_vars <- paste0(eff_prefix, "[", seq_len(n_markers), "]")
  if (all(eff_vars %in% all_vars)) {
    return(.get_draws_array(object$fit, variables = eff_vars, draws = draws, seed = seed))
  }

  legacy_eff_vars <- paste0("marker_weights_eff[", seq_len(n_markers), "]")
  if (all(legacy_eff_vars %in% all_vars)) {
    return(.get_draws_array(object$fit, variables = legacy_eff_vars, draws = draws, seed = seed))
  }

  probe_var <- .first_available_draw_var(
    all_vars,
    c(
      paste0("alpha_", term_key, "_eff"),
      paste0("alpha_", term_key)
    )
  )
  if (is.null(probe_var)) {
    return(NULL)
  }

  probe_arr <- .get_draws_array(object$fit, variables = probe_var, draws = draws, seed = seed)
  base_by_term <- sd$marker_weights_by_term %||% list()
  base_weights <- base_by_term[[term_key]] %||% sd$marker_weights %||% rep(1, n_markers)
  base_weights <- as.numeric(base_weights)
  if (length(base_weights) < n_markers) {
    base_weights <- c(base_weights, rep(1, n_markers - length(base_weights)))
  }
  base_weights <- base_weights[seq_len(n_markers)]

  array(
    rep(base_weights, each = dim(probe_arr)[1] * dim(probe_arr)[2]),
    dim = c(dim(probe_arr)[1], dim(probe_arr)[2], n_markers),
    dimnames = list(
      iteration = dimnames(probe_arr)[[1]],
      chain = dimnames(probe_arr)[[2]],
      variable = paste0(eff_prefix, "[", seq_len(n_markers), "]")
    )
  )
}
