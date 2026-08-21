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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
.normalise_marker_weight_vector <- function(weights, marker_levels, context, arg_name = "marker_weights$offset") {
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
    supplied_names <- names(weights)
    named_entries <- !is.na(supplied_names) & nzchar(supplied_names)
    if (any(named_entries) && !all(named_entries)) {
      cli::cli_abort(c(
        x = "{.arg {arg_name}} must be wholly named or wholly unnamed in {context}.",
        i = "Name every marker entry, or remove every marker name and use this fitted marker order: {.val {paste(marker_levels, collapse = ', ')}}."
      ))
    }
    if (all(named_entries) && anyDuplicated(supplied_names)) {
      cli::cli_abort("{.arg {arg_name}} contains duplicated marker names in {context}.")
    }
    missing_markers <- setdiff(marker_levels, supplied_names)
    unknown_markers <- setdiff(supplied_names, marker_levels)
    if (length(missing_markers) > 0L || length(unknown_markers) > 0L) {
      cli::cli_abort(c(
        x = "Named {.arg {arg_name}} does not match the fitted markers in {context}.",
        i = "Expected marker names: {.val {paste(marker_levels, collapse = ', ')}}."
      ))
    }
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
#' Summary:
#' 1. identify the active weighted association terms,
#' 2. normalise each supplied marker-weight vector into fitted marker order,
#' 3. either map all active terms to one shared weight set or assign one weight
#'    set per active term,
#' 4. return both per-term vectors and the compact set-index representation used
#'    by Stan.
#'
#' @param marker_weight_offsets Marker-weight offset from the checked
#'   `priors$marker_weights` declaration.
#' @param marker_levels Character marker labels in fitted order.
#' @param active_terms Character vector of active weighted association terms.
#' @param marker_weight_sets_shared Logical. If `TRUE`, all active weighted
#'   association terms share one marker-weight structure. If `FALSE`, each active
#'   weighted association term gets its own marker-weight structure.
#' @param estimate_marker_weights Logical. If `TRUE`, the default offsets
#'   are zero so the latent perturbation represents the full effective weight. If
#'   `FALSE`, the default offsets are one.
#' @param context Short context string for errors.
#'
#' @return A named list containing per-term vectors, set indices, and the compact
#'   matrix representation used by Stan.
#' @keywords internal
#' @noRd
.get_marker_weight_structure <- function(marker_weight_offsets,
                                             marker_levels,
                                             active_terms,
                                             marker_weight_sets_shared = TRUE,
                                             estimate_marker_weights = FALSE,
                                             context = "JoiNMe") {
  keys <- .weighted_assoc_term_keys()
  marker_levels <- as.character(marker_levels %||% character(0))
  n_markers <- length(marker_levels)
  active_terms <- intersect(keys, unique(as.character(active_terms %||% character(0))))

  if (!is.logical(marker_weight_sets_shared) || length(marker_weight_sets_shared) != 1L || is.na(marker_weight_sets_shared)) {
    cli::cli_abort(c(
      x = "{.arg marker_weights$shared} must be TRUE or FALSE in {context}.",
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
      marker_weight_sets_shared = isTRUE(marker_weight_sets_shared),
      active_terms = character(0),
      n_sets = 1L,
      set_index = stats::setNames(integer(length(keys)), keys),
      offset_by_term = per_term,
      offset_matrix = matrix(default_vec, nrow = 1L, byrow = TRUE,
        dimnames = list("inactive", marker_levels))
    ))
  }

  if (is.null(marker_weight_offsets)) {
    resolved_shared <- default_vec
    resolved_by_term <- NULL
  } else if (is.numeric(marker_weight_offsets)) {
    resolved_shared <- .normalise_marker_weight_vector(
      marker_weight_offsets,
      marker_levels = marker_levels,
      context = context,
      arg_name = "marker_weights$offset"
    )
    resolved_by_term <- NULL
  } else if (is.matrix(marker_weight_offsets) || is.data.frame(marker_weight_offsets)) {
    resolved_by_term <- as.data.frame(marker_weight_offsets, check.names = FALSE, stringsAsFactors = FALSE)
    resolved_shared <- NULL
  } else if (is.list(marker_weight_offsets)) {
    resolved_by_term <- marker_weight_offsets
    resolved_shared <- NULL
  } else {
    cli::cli_abort(c(
      x = "{.arg marker_weights$offset} has an unsupported type in {context}.",
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
        values <- as.numeric(resolved_by_term[term_key, , drop = TRUE])
        column_names <- colnames(resolved_by_term)
        if (!is.null(column_names) && all(nzchar(column_names))) names(values) <- column_names
        return(values)
      }
      resolved_by_term[[term_key]]
    }

    declared_terms <- names(resolved_by_term) %||% rownames(resolved_by_term) %||% character(0)
    unknown_terms <- setdiff(declared_terms, keys)
    if (length(unknown_terms) > 0L) {
      cli::cli_abort(c(
        x = "{.arg marker_weights$offset} contains unknown weighted association terms in {context}: {.val {paste(unknown_terms, collapse = ', ')}}.",
        i = "Use only {.val cv_total}, {.val cs_total}, {.val cv_marker}, and {.val cs_marker}."
      ))
    }
    supplied_terms <- intersect(declared_terms, keys)
    if (isTRUE(marker_weight_sets_shared)) {
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
            arg_name = paste0("marker_weights$offset$", active_terms[min(i, length(active_terms))])
          )
        })
        ref_vec <- normalised[[1]]
        if (length(normalised) > 1L) {
          mismatch <- vapply(normalised[-1], function(x) max(abs(x - ref_vec)) > 1e-8, logical(1))
          if (any(mismatch)) {
            cli::cli_abort(c(
              x = "{.arg marker_weights$offset} supplies different weight vectors while {.arg marker_weights$shared} is TRUE in {context}.",
              i = "Set {.arg marker_weights$shared} to FALSE to use different marker-weight structures across weighted association terms."
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
            arg_name = paste0("marker_weights$offset$", term_key)
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
  if (isTRUE(marker_weight_sets_shared)) {
    set_index[active_terms] <- 1L
    offset_matrix <- matrix(per_term[[active_terms[[1L]]]], nrow = 1L, byrow = TRUE,
      dimnames = list("shared", marker_levels))
    n_sets <- 1L
  } else {
    n_sets <- length(active_terms)
    offset_matrix <- matrix(0, nrow = n_sets, ncol = n_markers,
      dimnames = list(active_terms, marker_levels))
    for (set_pos in seq_along(active_terms)) {
      term_key <- active_terms[[set_pos]]
      set_index[[term_key]] <- set_pos
      offset_matrix[set_pos, ] <- per_term[[term_key]]
    }
  }

  list(
    marker_weight_sets_shared = isTRUE(marker_weight_sets_shared),
    active_terms = active_terms,
    n_sets = as.integer(n_sets),
    set_index = set_index,
    offset_by_term = per_term,
    offset_matrix = offset_matrix
  )
}

#' Resolve the data-generating means of marker-weight sets
#'
#' @description
#' A fitted marker-weight set is decomposed into a common location and one
#' marker-specific standardised departure.  This helper gives simulation the
#' same shared-versus-term-specific indexing used by the fitted model.  A
#' scalar is recycled over all active sets; a named vector or list may instead
#' provide one value for each active weighted association term. If no value is
#' supplied, one population mean per set is drawn once and returned for storage
#' in simulation truth; the fitting intercept prior is deliberately not used as
#' its generating distribution.
#' @param marker_weight_mean Optional numeric scalar, named numeric vector, or
#'   named list specifying the common location of each simulated weight set.
#' @param marker_weight_structure Output from `.get_marker_weight_structure()`
#' @param estimate_marker_weights Logical indicating whether weights are
#'   estimated. A constant family has no fitted common location.
#' @param context Short context string used in diagnostic messages.
#'
#' @return Numeric vector with one common location per compact weight set.
#' @keywords internal
#' @noRd
.resolve_marker_weight_means <- function(marker_weight_mean,
                                         marker_weight_structure,
                                         estimate_marker_weights,
                                         context = "simulate_joinme()") {
  n_sets <- as.integer(marker_weight_structure$n_sets %||% 1L) # number of distinct means required by the shared or term-specific design
  active_terms <- marker_weight_structure$active_terms %||% character(0) # weighted association terms represented by those sets

  # Constant-family marker weights are declared offsets rather than random
  # coefficients. Their common mean and departure decomposition is therefore absent. A
  # marker_weight_mean supplied through a shared simulation call is ignored:
  # the effective weights remain exactly the declared offsets, which mirrors
  # the fitting programme's constant-family branch.
  if (!isTRUE(estimate_marker_weights) || length(active_terms) == 0L) {
    return(rep(0, n_sets))
  }

  # An omitted population location is drawn once for every fitted set, exactly
  # as omitted longitudinal and event coefficients are drawn once. These draws
  # are fixed for the remainder of data generation and are returned in truth;
  # they are not sampled from the fitting intercept prior.
  if (is.null(marker_weight_mean)) {
    return(stats::rnorm(n_sets, mean = 0, sd = 0.5))
  }

  if (is.list(marker_weight_mean)) {
    if (is.null(names(marker_weight_mean))) {
      marker_weight_mean <- unlist(marker_weight_mean, use.names = FALSE)
    } else {
      unknown_terms <- setdiff(names(marker_weight_mean), active_terms)
      if (length(unknown_terms) > 0L) {
        cli::cli_abort(c(
          x = "{.arg marker_weight_mean} contains inactive or unknown weighted association terms in {context}: {.val {paste(unknown_terms, collapse = ', ')}}.",
          i = "Use only the active terms: {.val {paste(active_terms, collapse = ', ')}}."
        ))
      }
      marker_weight_mean <- vapply(active_terms, function(term_key) {
        value <- marker_weight_mean[[term_key]]
        if (is.null(value)) NA_real_ else as.numeric(value)
      }, numeric(1))
      names(marker_weight_mean) <- active_terms
    }
  }

  if (!is.numeric(marker_weight_mean)) {
    cli::cli_abort("{.arg marker_weight_mean} must be numeric in {context}.")
  }
  supplied_names <- names(marker_weight_mean) # term labels retained before numeric coercion removes attributes
  marker_weight_mean <- as.numeric(marker_weight_mean) # candidate common locations on the effective marker-weight scale
  if (length(marker_weight_mean) == 1L) {
    marker_weight_mean <- rep(marker_weight_mean, n_sets)
  } else if (!isTRUE(marker_weight_structure$marker_weight_sets_shared) &&
             !is.null(supplied_names)) {
    marker_weight_mean <- marker_weight_mean[match(active_terms, supplied_names)]
  }
  if (length(marker_weight_mean) != n_sets || any(!is.finite(marker_weight_mean))) {
    cli::cli_abort(c(
      x = "{.arg marker_weight_mean} must contain one finite value per marker-weight set in {context}.",
      i = "The fitted structure has {n_sets} set{?s}: {.val {paste(active_terms, collapse = ', ')}}."
    ))
  }

  unname(marker_weight_mean)
}

#' Marker-weight summary label
#'
#' @param term_key Weighted association term key.
#' @param marker_label Marker label.
#' @param marker_weight_sets_shared Logical. If `TRUE`, omit the term key because all
#'   weighted association terms share one structure.
#'
#' @return Character label used in summaries and extraction maps.
#' @keywords internal
#' @noRd
.marker_weight_summary_label <- function(term_key, marker_label, marker_weight_sets_shared = TRUE) {
  if (isTRUE(marker_weight_sets_shared)) {
    return(paste0("weight: ", marker_label))
  }
  paste0("weight[", term_key, "]: ", marker_label)
}

#' Marker-weight mean summary label
#'
#' @param term_key Weighted association term key represented by the set.
#' @param marker_weight_sets_shared Logical indicating whether all active terms use
#'   one common marker-weight set.
#'
#' @return Human-readable label for the fitted common weight location.
#' @keywords internal
#' @noRd
.marker_weight_mean_summary_label <- function(term_key, marker_weight_sets_shared = TRUE) {
  if (isTRUE(marker_weight_sets_shared)) {
    return("mean weight")
  }
  paste0("mean weight[", term_key, "]")
}

#' Marker-weight spread summary label
#'
#' @param term_key Weighted association term key represented by the set.
#' @param marker_weight_sets_shared Logical indicating whether all active terms use
#'   one common marker-weight set.
#'
#' @return Human-readable label for the within-set marker-weight spread.
#' @keywords internal
#' @noRd
.marker_weight_sd_summary_label <- function(term_key, marker_weight_sets_shared = TRUE) {
  if (isTRUE(marker_weight_sets_shared)) {
    return("SD weight")
  }
  paste0("SD weight[", term_key, "]")
}

#' Marker-weight Stan variable prefix
#'
#' @param term_key Weighted association term key.
#' @param effective Logical; if `TRUE`, return the effective posterior variable
#'   prefix, otherwise return the base-weight prefix.
#'
#' @return Character scalar prefix.
#' @keywords internal
#' @noRd
.marker_weight_var_prefix <- function(term_key, effective = TRUE) {
  suffix <- if (isTRUE(effective)) "eff_" else ""
  paste0("marker_weights_", suffix, term_key)
}

#' Fetch marker-weight draw array for one association term
#'
#' @param object A `JoiNMeFit` object.
#' @param term_key Weighted association term key.
#' @param draws Optional posterior draw subset size.
#' @param seed Random seed used when subsetting draws.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return Iteration x chain x marker draw array, or `NULL` if unavailable.
#' @keywords internal
#' @noRd
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

  probe_var <- paste0("alpha_", term_key)
  if (!(probe_var %in% all_vars)) {
    return(NULL)
  }

  probe_arr <- .get_draws_array(object$fit, variables = probe_var, draws = draws, seed = seed)
  offsets_by_term <- sd$marker_weight_offsets_by_term %||% list()
  weight_offsets <- offsets_by_term[[term_key]] %||% sd$marker_weight_offsets %||% rep(1, n_markers)
  weight_offsets <- as.numeric(weight_offsets)
  if (length(weight_offsets) < n_markers) {
    weight_offsets <- c(weight_offsets, rep(1, n_markers - length(weight_offsets)))
  }
  weight_offsets <- weight_offsets[seq_len(n_markers)]

  array(
    rep(weight_offsets, each = dim(probe_arr)[1] * dim(probe_arr)[2]),
    dim = c(dim(probe_arr)[1], dim(probe_arr)[2], n_markers),
    dimnames = list(
      iteration = dimnames(probe_arr)[[1]],
      chain = dimnames(probe_arr)[[2]],
      variable = paste0(eff_prefix, "[", seq_len(n_markers), "]")
    )
  )
}

#' Fetch marker-specific marker-weight departure draws
#'
#' @description
#' The random marker-weight contribution is the effective weight after removing
#' both its declared offset and the fitted common location for the corresponding
#' weight set. This is the quantity analogous to a group-level deviation in an
#' ordinary mixed model and is consequently the quantity reported by
#' [ranef()].
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param term_key Weighted association term identifying the required set.
#' @param draws Optional posterior draw subset size.
#' @param seed Random seed used when subsetting draws.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return Iteration by chain by marker array of centred departures, or `NULL`
#'   when the marker-weight family has no fitted random component.
#' @keywords internal
#' @noRd
.marker_weight_departure_array <- function(object, term_key, draws = NULL, seed = 1, all_vars = NULL) {
  stan_data <- object$stan_data # fitted weight-set map, marker dimension and declared offsets
  number_means <- as.integer(stan_data$n_marker_weight_means %||% 0L) # number of stochastic weight sets with fitted common locations
  set_index <- as.integer(stan_data[[paste0("marker_weight_set_", term_key)]] %||% 0L) # compact fitted set used by this association term
  if (number_means < 1L || set_index < 1L || set_index > number_means) return(NULL)

  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(error) character(0)) # posterior variables available from the selected fitting engine
  }
  mean_variable <- paste0("marker_weight_mean[", set_index, "]") # common marker-weight location for this fitted set
  if (!(mean_variable %in% all_vars)) return(NULL)

  effective_array <- .association_marker_weight_array(
    object,
    term_key = term_key,
    draws = draws,
    seed = seed,
    all_vars = all_vars
  ) # effective weights equal to offset plus common location plus departure
  if (is.null(effective_array)) return(NULL)
  mean_array <- .get_draws_array(
    object$fit,
    variables = mean_variable,
    draws = draws,
    seed = seed
  ) # common-location draws aligned by iteration and chain

  number_markers <- dim(effective_array)[3L] # marker dimension represented by the effective-weight array
  offsets_by_term <- stan_data$marker_weight_offsets_by_term %||% list() # declared marker-specific offsets indexed by association term
  offsets <- as.numeric(
    offsets_by_term[[term_key]] %||%
      stan_data$marker_weight_offsets %||%
      rep(0, number_markers)
  ) # declared offset aligned to fitted marker order
  if (length(offsets) != number_markers) offsets <- rep(NA_real_, number_markers)

  departure_array <- effective_array # output array retaining iteration, chain and marker dimensions
  mean_matrix <- matrix(
    as.numeric(mean_array[, , 1L]),
    nrow = dim(mean_array)[1L],
    ncol = dim(mean_array)[2L]
  ) # ordinary numeric matrix of the common-location draws
  for (marker_index in seq_len(number_markers)) {
    effective_matrix <- matrix(
      as.numeric(effective_array[, , marker_index]),
      nrow = dim(effective_array)[1L],
      ncol = dim(effective_array)[2L]
    ) # effective draws for one marker stripped of posterior subclasses
    departure_array[, , marker_index] <- effective_matrix - mean_matrix - offsets[[marker_index]]
  }
  dimnames(departure_array)[[3L]] <- paste0(
    "marker_weight_departure_", term_key, "[", seq_len(number_markers), "]"
  ) # unique diagnostic names for the derived random-effect coordinates
  departure_array
}

#' Extract marker-weight draws on the requested mixed-model scale
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param quantity Either `"departure"` for [ranef()] or `"effective"` for
#'   [coef()].
#' @param draws Optional posterior draw subset size.
#' @param seed Random seed used when subsetting draws.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return A long data frame containing `draw`, `assoc_term`, `marker`, `term`,
#'   and `value`, or `NULL` when no corresponding quantity is fitted.
#' @keywords internal
#' @noRd
.marker_weight_long_draws <- function(object,
                                      quantity = c("departure", "effective"),
                                      draws = NULL,
                                      seed = 1,
                                      all_vars = NULL) {
  quantity <- match.arg(quantity) # requested random-effect or combined-coefficient representation
  stan_data <- object$stan_data # fitted association activation, sharing map and marker labels
  representative_terms <- .reported_marker_weight_terms(stan_data) # one term providing access to each distinct fitted weight set
  if (!length(representative_terms)) return(NULL)
  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(error) character(0)) # available posterior variables used to avoid invalid extraction requests
  }

  marker_labels <- as.character(
    stan_data$marker_levels %||% paste0("marker_", seq_len(as.integer(stan_data$D %||% 0L)))
  ) # marker labels in fitted order
  sets_shared <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # whether the representative term denotes a set shared across association terms
  output_rows <- list() # draw-level tables accumulated over distinct marker-weight sets

  for (term_key in representative_terms) {
    weight_array <- if (identical(quantity, "departure")) {
      .marker_weight_departure_array(object, term_key, draws = draws, seed = seed, all_vars = all_vars)
    } else {
      .association_marker_weight_array(object, term_key, draws = draws, seed = seed, all_vars = all_vars)
    } # posterior marker weights on the requested mixed-model scale
    if (is.null(weight_array)) next

    draw_matrix <- .assoc_matrix_from_draw_array(weight_array, term_labels = marker_labels) # flattened draws with one column per marker
    association_label <- if (sets_shared) "shared" else term_key # explicit set identity, essential when active terms have distinct weights
    metadata <- data.frame(
      assoc_term = rep(association_label, length(marker_labels)),
      marker = marker_labels,
      term = marker_labels,
      stringsAsFactors = FALSE
    ) # scientific identifiers aligned with the marker columns
    output_rows[[length(output_rows) + 1L]] <- .pivot_long_from_matrix(draw_matrix, metadata)
  }
  if (!length(output_rows)) NULL else do.call(rbind, output_rows)
}

#' Summarise marker weights on the requested mixed-model scale
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param quantity Either `"departure"` or `"effective"`.
#' @param draws Optional posterior draw subset size.
#' @param seed Random seed used when subsetting draws.
#' @param digits Number of displayed decimal places.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return A posterior summary table with association-term and marker columns,
#'   or `NULL` when no corresponding quantity is fitted.
#' @keywords internal
#' @noRd
.marker_weight_component_summary <- function(object,
                                             quantity = c("departure", "effective"),
                                             draws = NULL,
                                             seed = 1,
                                             digits = 3,
                                             all_vars = NULL) {
  quantity <- match.arg(quantity) # requested random-effect or combined-coefficient representation
  stan_data <- object$stan_data # fitted association activation, sharing map and marker labels
  representative_terms <- .reported_marker_weight_terms(stan_data) # one term providing access to each distinct fitted weight set
  if (!length(representative_terms)) return(NULL)
  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(error) character(0)) # posterior variables used for checked extraction
  }

  marker_labels <- as.character(
    stan_data$marker_levels %||% paste0("marker_", seq_len(as.integer(stan_data$D %||% 0L)))
  ) # marker labels in fitted order
  sets_shared <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # whether one common set serves every active weighted term
  output_rows <- list() # summary tables accumulated over fitted weight sets

  for (term_key in representative_terms) {
    weight_array <- if (identical(quantity, "departure")) {
      .marker_weight_departure_array(object, term_key, draws = draws, seed = seed, all_vars = all_vars)
    } else {
      .association_marker_weight_array(object, term_key, draws = draws, seed = seed, all_vars = all_vars)
    } # posterior marker weights on the requested scale
    if (is.null(weight_array)) next

    table <- .assoc_summary_from_draw_array(weight_array, term_labels = marker_labels, digits = digits) # marker-specific posterior summaries with chain diagnostics
    table$assoc_term <- if (sets_shared) "shared" else term_key
    table$marker <- table$term
    table$group <- "marker_weight"
    table <- table[, c(
      "assoc_term", "marker", "term", "Estimate", "Est.Error", "Q2.5",
      "Q97.5", "Rhat", "ess_bulk", "ess_tail", "group"
    ), drop = FALSE]
    output_rows[[length(output_rows) + 1L]] <- table
  }
  if (!length(output_rows)) NULL else do.call(rbind, output_rows)
}

#' Marker-weight sets represented in posterior reports
#'
#' @param stan_data Fitted standata and its retained presentation metadata.
#'
#' @return Character vector containing one representative association term per
#'   distinct marker-weight set.
#' @keywords internal
#' @noRd
.reported_marker_weight_terms <- function(stan_data) {
  active_terms <- .active_weighted_assoc_terms(stan_data) # weighted association channels present in the fitted event model
  shared_weights <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # whether those channels point to one common weight set
  if (shared_weights && length(active_terms) > 1L) active_terms <- active_terms[1L]
  active_terms
}

#' Map fitted marker-weight means to association terms
#'
#' @description
#' Resolves each compact marker-weight set through its explicit Stan set index.
#' This avoids assuming that posterior mean order happens to match association
#' activation order when term-specific marker weights are fitted.
#'
#' @param stan_data Fitted marker-weight set indices and sharing declaration.
#' @param all_vars Optional character vector of available posterior variables.
#'
#' @return A data frame containing posterior variable, association term, and
#'   display term for every fitted common marker-weight location.
#' @keywords internal
#' @noRd
.marker_weight_mean_map <- function(stan_data, all_vars = NULL) {
  number_means <- as.integer(stan_data$n_marker_weight_means %||% 0L) # number of fitted common marker-weight locations
  if (number_means < 1L) {
    return(data.frame(
      variable = character(0), assoc_term = character(0), term = character(0),
      stringsAsFactors = FALSE
    ))
  }

  active_terms <- .active_weighted_assoc_terms(stan_data) # weighted association terms present in the fitted event model
  term_by_set <- rep(NA_character_, number_means) # representative association term assigned through each explicit set index
  for (term_key in active_terms) {
    set_index <- as.integer(stan_data[[paste0("marker_weight_set_", term_key)]] %||% 0L) # compact set index used by this active term
    if (set_index >= 1L && set_index <= number_means && is.na(term_by_set[[set_index]])) {
      term_by_set[[set_index]] <- term_key
    }
  }
  missing_sets <- which(is.na(term_by_set)) # defensive labels for incomplete synthetic fixtures
  if (length(missing_sets)) term_by_set[missing_sets] <- paste0("set_", missing_sets)

  sets_shared <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # whether all active terms point to one common set
  variables <- paste0("marker_weight_mean[", seq_len(number_means), "]") # posterior variables in compact set order
  association_terms <- if (sets_shared) rep("shared", number_means) else term_by_set # public association identity for each fitted mean
  display_terms <- if (sets_shared) {
    rep("mean weight", number_means)
  } else {
    paste0("mean weight[", association_terms, "]")
  } # unique row and draw labels for shared or term-specific means
  output <- data.frame(
    variable = variables,
    assoc_term = association_terms,
    term = display_terms,
    stringsAsFactors = FALSE
  ) # complete mean map before optional posterior-variable filtering
  if (!is.null(all_vars)) output <- output[output$variable %in% all_vars, , drop = FALSE]
  output
}

#' Declared marker-weight offsets by fitted set
#'
#' @param stan_data Fitted standata and retained marker labels.
#'
#' @return A named list of marker-labelled numeric offset vectors, or `NULL`
#'   when no marker-weighted association is active.
#' @keywords internal
#' @noRd
.marker_weight_offset_metadata <- function(stan_data) {
  representative_terms <- .reported_marker_weight_terms(stan_data) # one term for each shared or term-specific weight set
  if (!length(representative_terms)) return(NULL)
  marker_labels <- as.character(stan_data$marker_levels %||% paste0("marker_", seq_len(as.integer(stan_data$D %||% 0L)))) # fitted marker order used by every offset vector
  offsets_by_term <- stan_data$marker_weight_offsets_by_term %||% list() # term-indexed declared offsets retained by current fits
  shared_weights <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # controls the public set label
  offsets <- lapply(representative_terms, function(term_key) {
    values <- as.numeric(offsets_by_term[[term_key]] %||% stan_data$marker_weight_offsets %||% rep(0, length(marker_labels))) # declared offsets aligned to this fitted weight set
    if (length(values) != length(marker_labels)) values <- rep(NA_real_, length(marker_labels))
    stats::setNames(values, marker_labels)
  })
  names(offsets) <- if (shared_weights) "shared" else representative_terms
  offsets
}

#' Summarise fitted marker-weight locations and within-set spread
#'
#' @description
#' For every fitted weight set, this helper reports the posterior common
#' location and the realised root-mean-square marker departure around that
#' location. The second row is a posterior descriptive quantity calculated as
#' `sqrt(mean(z_sd^2))` across the fitted markers in each draw. It is not a
#' separately fitted scale parameter: the departure family has fixed unit
#' scale, while the association slope governs the magnitude of the weighted
#' marker contribution to the log hazard.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param draws Optional posterior draw count.
#' @param seed Random seed used for posterior subsetting.
#' @param digits Number of displayed decimal places.
#' @param all_vars Optional available posterior variable names.
#'
#' @return A posterior summary table with two rows per fitted weight set, or
#'   `NULL` for fixed or inactive marker weights.
#' @keywords internal
#' @noRd
.marker_weight_set_summary <- function(object, draws = NULL, seed = 1, digits = 3, all_vars = NULL) {
  stan_data <- object$stan_data # fitted dimensions, sharing map, marker labels and declared offsets
  number_means <- as.integer(stan_data$n_marker_weight_means %||% 0L) # number of estimated common locations and therefore fitted weight sets
  representative_terms <- .reported_marker_weight_terms(stan_data) # one active association term giving access to each distinct set
  if (number_means < 1L || !length(representative_terms)) return(NULL)
  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(error) character(0))
  }
  shared_weights <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # determines compact shared or term-specific row labels
  output_rows <- list()

  for (set_index in seq_len(min(number_means, length(representative_terms)))) {
    term_key <- representative_terms[[set_index]] # representative association channel for this compact weight set
    mean_variable <- paste0("marker_weight_mean[", set_index, "]") # fitted common location variable for the set
    if (!(mean_variable %in% all_vars)) next
    mean_array <- .get_draws_array(object$fit, variables = mean_variable, draws = draws, seed = seed) # iteration-by-chain draws of the common location
    departure_array <- .marker_weight_departure_array(
      object,
      term_key = term_key,
      draws = draws,
      seed = seed,
      all_vars = all_vars
    ) # direct unit-scale departures reconstructed as effective weight minus offset and common location
    if (is.null(departure_array)) next
    departure_rms <- apply(departure_array^2, c(1L, 2L), mean) # posterior mean squared departure across the finite fitted marker collection
    departure_rms <- matrix(
      sqrt(as.numeric(departure_rms)),
      nrow = dim(departure_array)[1L],
      ncol = dim(departure_array)[2L]
    ) # root-mean-square departure for every iteration and chain, including the one-marker case
    spread_array <- array(
      departure_rms,
      dim = c(dim(departure_array)[1L], dim(departure_array)[2L], 1L),
      dimnames = list(
        iteration = dimnames(departure_array)[[1L]],
        chain = dimnames(departure_array)[[2L]],
        variable = .marker_weight_sd_summary_label(term_key, marker_weight_sets_shared = shared_weights)
      )
    ) # draws-array representation accepted by the common posterior summary calculation

    mean_label <- .marker_weight_mean_summary_label(term_key, marker_weight_sets_shared = shared_weights) # common-location row label
    spread_label <- .marker_weight_sd_summary_label(term_key, marker_weight_sets_shared = shared_weights) # realised root-mean-square departure row label
    output_rows[[length(output_rows) + 1L]] <- .assoc_summary_from_draw_array(mean_array, term_labels = mean_label, digits = digits)
    output_rows[[length(output_rows) + 1L]] <- .assoc_summary_from_draw_array(spread_array, term_labels = spread_label, digits = digits)
  }
  if (!length(output_rows)) NULL else do.call(rbind, output_rows)
}

#' Extract individual marker weights
#'
#' @description
#' Returns the effective marker weights used by marker-aggregated association
#' terms. The effective value is the declared marker-specific offset plus the
#' fitted common weight location and marker-specific departure. When
#' `marker_weights$shared = TRUE`, one shared set is returned once. Otherwise,
#' each active weighted association term receives its own set.
#'
#' Individual weights are intentionally kept out of [summary.JoiNMeFit()]. The
#' model summary reports only the common location and within-set spread. This
#' function and [coef()] return effective weights, [fixef()] returns their
#' fitted common locations, and [ranef()] returns only marker-specific
#' departures after subtracting the common location and declared offset.
#'
#' @param object A fitted JoiNMe model.
#' @param draws Optional number of posterior draws.
#' @param seed Random seed used when subsetting draws.
#' @param digits Number of decimal places used when `summary = TRUE`.
#' @param summary Logical; return posterior summaries when `TRUE`, otherwise a
#'   draws-by-weight matrix.
#' @param ... Unused.
#'
#' @return A data frame with one row per marker and weight set when
#'   `summary = TRUE`; otherwise a posterior draw matrix.
#' @export
marker_weights <- function(object, ...) {
  UseMethod("marker_weights")
}

#' @rdname marker_weights
#' @export
marker_weights.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary), msg = "summary must be TRUE or FALSE.")
  stan_data <- object$stan_data # fitted sharing map, marker labels and declared offsets
  representative_terms <- .reported_marker_weight_terms(stan_data) # one association term per distinct marker-weight set
  if (!length(representative_terms) || as.integer(stan_data$D %||% 0L) < 1L) {
    cli::cli_abort("No marker-weighted association term is available in this model.")
  }
  if (is.null(draws)) draws <- object$config$draws_default
  marker_labels <- as.character(stan_data$marker_levels %||% paste0("marker_", seq_len(as.integer(stan_data$D)))) # public marker order
  shared_weights <- isTRUE(as.integer(stan_data$marker_weight_sets_shared %||% 1L) == 1L) # whether one set is returned once
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(object$fit)), error = function(error) character(0)) # available effective-weight variables
  offset_sets <- .marker_weight_offset_metadata(stan_data) # declared marker-specific constants shown beside posterior weights
  arrays <- list()
  summaries <- list()

  for (term_key in representative_terms) {
    weight_array <- .association_marker_weight_array(object, term_key = term_key, draws = draws, seed = seed, all_vars = all_vars) # effective weights for this shared or term-specific set
    if (is.null(weight_array)) next
    set_name <- if (shared_weights) "shared" else term_key # public set identifier
    offsets <- as.numeric(offset_sets[[set_name]] %||% rep(NA_real_, length(marker_labels))) # declared offsets aligned to marker rows
    if (isTRUE(summary)) {
      table <- .assoc_summary_from_draw_array(weight_array, term_labels = marker_labels, digits = digits) # individual effective-weight posterior summaries
      table$set <- set_name
      table$marker <- table$term
      table$offset <- offsets[match(table$marker, marker_labels)]
      table <- table[, c("set", "marker", "offset", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      summaries[[length(summaries) + 1L]] <- table
    } else {
      matrix <- .assoc_matrix_from_draw_array(weight_array, term_labels = vapply(marker_labels, function(marker_label) {
        .marker_weight_summary_label(term_key, marker_label, marker_weight_sets_shared = shared_weights)
      }, character(1))) # draw matrix labelled consistently with fixef() and ranef()
      arrays[[length(arrays) + 1L]] <- matrix
    }
  }
  if (isTRUE(summary)) {
    if (!length(summaries)) return(NULL)
    return(do.call(rbind, summaries))
  }
  if (!length(arrays)) return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  do.call(cbind, arrays)
}

#' Posterior individual marker-weight alias
#'
#' @param object A fitted JoiNMe model.
#' @param ... Arguments forwarded to [marker_weights()].
#'
#' @return The draw matrix returned by `marker_weights(summary = FALSE)`.
#' @export
posterior_marker_weights <- function(object, ...) {
  marker_weights(object, summary = FALSE, ...)
}
