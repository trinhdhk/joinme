#' Posterior draws for JoiNMe objects
#'
#' @description
#' Returns posterior draws in `posterior`-compatible formats with user-facing
#' parameter names.
#'
#' Compared with [extract()], `draws()` is the convenience layer for downstream
#' posterior workflows. It is designed for tasks such as `posterior`
#' summarisation, `bayesplot` visualisation, regex-based variable selection, and
#' any workflow that expects a standard draws object.
#'
#' For fitted `JoiNMe` models, known Stan variables are relabelled from raw fit object
#' to friendly parameter names. 
#' For dynamic prediction objects, stored draw blocks are flattened into a single
#' draws object with explicit id, scale, marker, and time labels.
#'
#' Fitted-object draw arrays are cached inside the underlying R6 container after
#' the first request so later summaries, diagnostics, and plotting methods can
#' reuse the same renamed draw payload without re-reading the backend fit.
#'
#' Use `draws()` when you want:
#' - one posterior object containing renamed variables,
#' - `posterior::subset_draws()` and regex-style variable filtering,
#' - `bayesplot` directly, similar to [mcmc_plot()],
#' - a standard draws array/matrix/data frame rather than a component-specific
#'   extraction payload.
#' - `as.array` is a shorthand for `draws(format = "draws_array")`.
#'
#' Use [extract()] instead when you want:
#' - one model component at a time (`"fixef"`, `"assoc"`, `"gamma_w"`, `"basehaz"`, etc.),
#' - the explicit `term_map` telling you how user-facing labels map back to raw
#'   Stan variables,
#' - special structured payloads such as `what = "association_plot"`,
#' - prediction draw blocks separated by semantic role before flattening.
#'
#' @param object A supported JoiNMe object.
#' @param variables Optional character vector selecting variables after
#'   relabelling.
#' @param regex Logical; if `TRUE`, interpret `variables` as regular
#'   expressions.
#' @param draws Optional number of posterior draws to retain.
#' @param seed Integer seed used when subsetting draws.
#' @param what For `JoiNMeFit` objects, either `"all"` (default renamed
#'   posterior variables), `"basehaz"`, or `"baseline_hazard"`.
#' @param format Output format. Supported values are `"draws_array"`,
#'   `"draws_matrix"`, and `"draws_df"`.
#' @param ... Unused.
#'
#' @return A posterior draw object in the requested format. The result is a
#'   `posterior`-compatible object with renamed variables, suitable for
#'   `posterior::summarise_draws()`, `posterior::subset_draws()`, and
#'   `bayesplot`-style visualisation.
#' @export
draws <- function(object, ...) {
  UseMethod("draws")
}

#' @keywords internal
#' @noRd
.fit_term_map <- function(object, all_vars) {
  mapped <- lapply(
    c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "likelihood_scale"),
    function(what) .fit_component_term_map(object, what = what, all_vars = all_vars)
  )
  mapped <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, mapped)
  if (length(mapped) == 0L) {
    return(data.frame(term = all_vars, variable = all_vars, stringsAsFactors = FALSE))
  }

  map <- do.call(rbind, mapped)
  map <- map[!duplicated(map$variable), , drop = FALSE]
  leftover <- setdiff(all_vars, map$variable)
  p <- as.integer(object$stan_data$P %||% 0L)
  if (p > 0L) {
    fixed_map <- .fixed_effect_var_map(sd = object$stan_data, all_vars = all_vars)
    chosen_fixed <- as.character(fixed_map$variable)
    if (length(chosen_fixed) > 0L) {
      beta_raw <- paste0("beta[", seq_len(p), "]")
      beta_scaled <- paste0("beta_scaled[", seq_len(p), "]")
      drop_beta <- setdiff(c(beta_raw, beta_scaled), chosen_fixed)
      leftover <- setdiff(leftover, drop_beta)
    }
  }
  if (length(leftover) > 0L) {
    map <- rbind(
      map,
      data.frame(term = leftover, variable = leftover, stringsAsFactors = FALSE)
    )
  }
  map
}

#' @keywords internal
#' @noRd
.rename_draw_variables <- function(all_vars, map) {
  renamed <- all_vars
  if (is.null(map) || !nrow(map)) {
    return(make.unique(as.character(renamed)))
  }
  idx <- match(all_vars, map$variable)
  keep <- !is.na(idx)
  renamed[keep] <- map$term[idx[keep]]
  make.unique(as.character(renamed))
}

#' @keywords internal
#' @noRd
.subset_draws <- function(draws_obj, variables = NULL, regex = FALSE, draws = NULL, seed = 1) {
  out <- draws_obj
  if (!is.null(variables)) {
    out <- posterior::subset_draws(out, variable = variables, regex = regex)
    if (length(posterior::variables(out)) == 0L) {
      cli::cli_abort(c(
        x = "No posterior variables matched the requested selection.",
        i = "Check {.arg variables} after relabelling, or disable {.arg regex}."
      ))
    }
  }
  if (!is.null(draws) && is.finite(draws)) {
    nd <- posterior::ndraws(out)
    if (draws < nd) {
      set.seed(seed)
      out <- posterior::subset_draws(out, draw = sample.int(nd, size = draws))
    }
  }
  out
}

#' @keywords internal
#' @noRd
.format_draws <- function(draws_obj, format = c("draws_array", "draws_matrix", "draws_df")) {
  format <- match.arg(format)
  switch(
    format,
    draws_array = draws_obj,
    draws_matrix = posterior::as_draws_matrix(draws_obj),
    draws_df = posterior::as_draws_df(draws_obj)
  )
}

#' @keywords internal
#' @noRd
.fit_cached_draws_array <- function(object) {
  cache_key <- "renamed_draws_array"
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  draw_array <- posterior::as_draws_array(.get_draws_obj(object$fit, keep_chains = TRUE))
  all_vars <- posterior::variables(draw_array)
  term_map <- .fit_term_map(object, all_vars = all_vars)
  dimnames(draw_array)[[3]] <- .rename_draw_variables(all_vars, term_map)

  object$cache_set(cache_key, draw_array)
  draw_array
}

#' @keywords internal
#' @noRd
.flatten_dynpred_regular_block <- function(object, what) {
  ext <- tryCatch(extract.JoiNMeDynPred(object, what = what), error = function(e) NULL)
  if (is.null(ext) || is.null(ext$draws) || !length(ext$draws)) {
    return(NULL)
  }

  mats <- list()
  for (id_nm in names(ext$draws)) {
    entry <- ext$draws[[id_nm]]
    if (is.list(entry) && !is.matrix(entry)) {
      for (scale_nm in names(entry)) {
        mats[[paste0(id_nm, "::", scale_nm)]] <- entry[[scale_nm]]
      }
    } else if (is.matrix(entry)) {
      mats[[id_nm]] <- entry
    }
  }
  if (!length(mats)) {
    return(NULL)
  }
  do.call(cbind, mats)
}

#' @keywords internal
#' @noRd
.flatten_dynpred_random_effects_id <- function(object) {
  dd <- object$draws$random_effects_id
  if (is.null(dd) || !length(dd)) {
    return(NULL)
  }

  mats <- lapply(names(dd), function(id_nm) {
    entry <- dd[[id_nm]]
    mat <- entry$matrix
    if (is.null(mat) || !is.matrix(mat) || ncol(mat) == 0L) {
      return(NULL)
    }
    terms <- entry$terms %||% colnames(mat) %||% paste0("u_id[", seq_len(ncol(mat)), "]")
    colnames(mat) <- paste0("random_effects_id[id=", id_nm, "|term=", terms, "]")
    mat
  })
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) {
    return(NULL)
  }
  do.call(cbind, mats)
}

#' @keywords internal
#' @noRd
.flatten_dynpred_random_effects_marker_id <- function(object) {
  dd <- object$draws$random_effects_marker_id
  if (is.null(dd) || !length(dd)) {
    return(NULL)
  }

  mats <- list()
  for (id_nm in names(dd)) {
    entry <- dd[[id_nm]]
    mat <- entry$matrix
    if (is.matrix(mat) && ncol(mat) > 0L) {
      terms <- colnames(mat) %||% entry$terms %||% paste0("w_idm[", seq_len(ncol(mat)), "]")
      colnames(mat) <- paste0("random_effects_marker_id[id=", id_nm, "|term=", terms, "]")
      mats[[paste0(id_nm, "::matrix")]] <- mat
    }

    corr <- entry$corr
    if (!is.null(corr) && length(dim(corr)) == 3L && prod(dim(corr)[2:3]) > 0L) {
      q_dim <- dim(corr)[2]
      term_labels <- entry$terms %||% paste0("w_idm[", seq_len(q_dim), "]")
      corr_mat <- matrix(NA_real_, nrow = dim(corr)[1], ncol = q_dim * q_dim)
      corr_names <- character(q_dim * q_dim)
      pos <- 1L
      for (row_idx in seq_len(q_dim)) {
        for (col_idx in seq_len(q_dim)) {
          corr_mat[, pos] <- corr[, row_idx, col_idx]
          corr_names[pos] <- paste0(
            "corr_marker_id[id=", id_nm,
            "|row=", term_labels[row_idx],
            "|col=", term_labels[col_idx],
            "]"
          )
          pos <- pos + 1L
        }
      }
      colnames(corr_mat) <- corr_names
      mats[[paste0(id_nm, "::corr")]] <- corr_mat
    }
  }

  if (!length(mats)) {
    return(NULL)
  }
  do.call(cbind, mats)
}

#' @keywords internal
#' @noRd
.dynpred_cached_draws_array <- function(object) {
  cache_key <- "flattened_draws_array"
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  blocks <- Filter(Negate(is.null), list(
    .flatten_dynpred_regular_block(object, "longitudinal"),
    .flatten_dynpred_regular_block(object, "longitudinal_fitted"),
    .flatten_dynpred_regular_block(object, "survival"),
    .flatten_dynpred_regular_block(object, "cumhaz"),
    .flatten_dynpred_random_effects_id(object),
    .flatten_dynpred_random_effects_marker_id(object)
  ))
  if (!length(blocks)) {
    cli::cli_abort(c(
      x = "No stored posterior draws are available for this prediction object.",
      i = "Generate predictions with draw outputs enabled before calling {.fn draws}."
    ))
  }

  n_rows <- unique(vapply(blocks, nrow, integer(1)))
  if (length(n_rows) != 1L) {
    cli::cli_abort("Prediction draw blocks do not share a common number of posterior draws.")
  }

  draw_matrix <- do.call(cbind, blocks)
  colnames(draw_matrix) <- make.unique(colnames(draw_matrix))
  draw_array <- posterior::as_draws_array(array(
    draw_matrix,
    dim = c(nrow(draw_matrix), 1L, ncol(draw_matrix)),
    dimnames = list(
      iteration = as.character(seq_len(nrow(draw_matrix))),
      chain = "1",
      variable = colnames(draw_matrix)
    )
  ))

  object$cache_set(cache_key, draw_array)
  draw_array
}

#' @rdname draws
#' @export
draws.JoiNMeFit <- function(object, variables = NULL, regex = FALSE, draws = NULL, seed = 1,
                            what = c("all", "basehaz", "baseline_hazard"),
                            format = c("draws_array", "draws_matrix", "draws_df"), ...) {
  assertthat::assert_that(inherits(object, "JoiNMeFit"), msg = "Object must be a JoiNMeFit instance.")
  what <- match.arg(what)

  if (!identical(what, "all")) {
    ext_what <- switch(
      what,
      basehaz = "basehaz",
      baseline_hazard = "baseline_hazard"
    )
    ext <- extract.JoiNMeFit(
      object,
      what = ext_what,
      draws = draws,
      seed = seed,
      keep_chains = TRUE
    )
    arr <- ext$draws
    arr <- .subset_draws(arr, variables = variables, regex = regex, draws = NULL, seed = seed)
    return(.format_draws(arr, format = format))
  }

  draw_array <- .fit_cached_draws_array(object)
  draw_array <- .subset_draws(draw_array, variables = variables, regex = regex, draws = draws, seed = seed)
  .format_draws(draw_array, format = format)
}

#' @rdname draws
#' @export
draws.JoiNMeDynPred <- function(object, variables = NULL, regex = FALSE, draws = NULL, seed = 1,
                                format = c("draws_array", "draws_matrix", "draws_df"), ...) {
  assertthat::assert_that(inherits(object, "JoiNMeDynPred"), msg = "Object must be a JoiNMeDynPred instance.")
  draw_array <- .dynpred_cached_draws_array(object)
  draw_array <- .subset_draws(draw_array, variables = variables, regex = regex, draws = draws, seed = seed)
  .format_draws(draw_array, format = format)
}

#' @rdname draws
#' @export
as.array.JoiNMeFit <- function(x, ...) {
  as.array(draws.JoiNMeFit(x, format = "draws_array", ...))
}

#' @rdname draws
#' @export
as.array.JoiNMeDynPred <- function(x, ...) {
  as.array(draws.JoiNMeDynPred(x, format = "draws_array", ...))
}
