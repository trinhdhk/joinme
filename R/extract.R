#' Extract posterior draws from joinme objects
#'
#' @description
#' S3 generic to extract draw-level matrices from `JoinMeFit` and `JoinMeDynPred`
#' objects using names aligned to summary tables where feasible.
#'
#' @param object A supported joinme object.
#' @param ... Additional method-specific arguments.
#' @export
extract <- function(object, ...) {
  UseMethod("extract")
}

# File overview:
# - Extract draw matrices from JoinMeFit by summary-like components.
# - Extract draw matrices from JoinMeDynPred from stored draw payloads.

#' Extract posterior draws from a fitted joinme model
#'
#' @param object A `JoinMeFit` object.
#' @param what Character component selector. One of
#'   `"fixef"`, `"gamma_w"`, `"assoc"`, `"distributional"`,
#'   `"distributional_regression"`, or `"raw"`.
#' @param term Optional character vector of friendly term names (summary-style)
#'   to subset extracted columns.
#' @param variable Optional character vector of raw Stan variable names. This is
#'   used directly when `what = "raw"` and can also further filter mapped outputs.
#' @param draws Optional number of posterior draws to keep.
#' @param seed Integer seed used when subsetting draws.
#'
#' @return A list with fields:
#'   - `draws`: numeric matrix (rows = draws, cols = requested terms)
#'   - `term_map`: data.frame mapping `term` to Stan `variable`
#' @export
extract.JoinMeFit <- function(object,
                              what = c("fixef", "gamma_w", "assoc", "distributional", "distributional_regression", "raw"),
                              term = NULL,
                              variable = NULL,
                              draws = NULL,
                              seed = 1,
                              ...) {
  what <- match.arg(what)
  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  # Build a friendly term -> variable map following summary naming conventions.
  map <- data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE)

  if (what == "fixef") {
    beta_vars <- paste0("beta[", seq_len(sd$P), "]")
    beta_vars <- beta_vars[beta_vars %in% all_vars]
    beta_terms <- sd$x_cols %||% beta_vars
    if (length(beta_terms) != length(beta_vars)) beta_terms <- beta_vars
    map <- data.frame(term = as.character(beta_terms), variable = as.character(beta_vars), stringsAsFactors = FALSE)
  } else if (what == "gamma_w") {
    g_vars <- paste0("gamma_w[", seq_len(sd$p_w %||% 0L), "]")
    g_vars <- g_vars[g_vars %in% all_vars]
    g_terms <- sd$w_cols %||% g_vars
    if (length(g_terms) != length(g_vars)) g_terms <- g_vars
    map <- data.frame(term = as.character(g_terms), variable = as.character(g_vars), stringsAsFactors = FALSE)
  } else if (what == "assoc") {
    assoc_vars <- c(
      if (isTRUE(sd$assoc_cv_total == 1)) "alpha_cv_total",
      if (isTRUE(sd$assoc_cv_mean == 1)) "alpha_cv_mean",
      if (isTRUE(sd$assoc_cv_marker == 1)) "alpha_cv_marker",
      if (isTRUE(sd$assoc_cs_total == 1)) "alpha_cs_total",
      if (isTRUE(sd$assoc_cs_mean == 1)) "alpha_cs_mean",
      if (isTRUE(sd$assoc_cs_marker == 1)) "alpha_cs_marker",
      if (isTRUE(sd$assoc_corr == 1)) grep("^alpha_corr\\[", all_vars, value = TRUE)
    )
    assoc_vars <- assoc_vars[assoc_vars %in% all_vars]

    assoc_map <- c(
      alpha_cv_total = "cv_total",
      alpha_cv_mean = "cv_mean",
      alpha_cv_marker = "cv_marker",
      alpha_cs_total = "cs_total",
      alpha_cs_mean = "cs_mean",
      alpha_cs_marker = "cs_marker"
    )
    assoc_terms <- assoc_map[assoc_vars]
    assoc_terms[is.na(assoc_terms)] <- sub("^alpha_", "", assoc_vars[is.na(assoc_terms)])

    map <- data.frame(term = as.character(assoc_terms), variable = as.character(assoc_vars), stringsAsFactors = FALSE)

    # Marker-weight terms mirror summary naming when present.
    mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D %||% 0L), "]")
    mw_vars <- mw_vars[mw_vars %in% all_vars]
    if (length(mw_vars) > 0) {
      marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(length(mw_vars)))
      if (length(marker_terms) != length(mw_vars)) marker_terms <- paste0("marker_", seq_len(length(mw_vars)))
      mw_map <- data.frame(
        term = paste0("weight: ", marker_terms),
        variable = as.character(mw_vars),
        stringsAsFactors = FALSE
      )
      map <- rbind(map, mw_map)
    }
  } else if (what == "distributional") {
    dist_map <- .distributional_term_map(sd, cfg, all_vars)
    if (length(dist_map) > 0) {
      map <- data.frame(
        term = unname(dist_map),
        variable = names(dist_map),
        stringsAsFactors = FALSE
      )
    }
  } else if (what == "distributional_regression") {
    dist_cols <- cfg$dist$dist_cols %||% list()
    reg_specs <- list(
      sigma = list(prefix = "beta_sigma", cols = dist_cols$sigma %||% character(0)),
      nu = list(prefix = "beta_nu", cols = dist_cols$nu %||% character(0)),
      phi = list(prefix = "beta_phi", cols = dist_cols$phi %||% character(0)),
      alpha = list(prefix = "beta_alpha", cols = dist_cols$alpha %||% character(0)),
      phi_beta = list(prefix = "beta_phi_beta", cols = dist_cols$phi_beta %||% character(0)),
      tau_sde = list(prefix = "beta_tau_sde", cols = dist_cols$tau_sde %||% character(0))
    )

    map_rows <- list()
    for (nm in names(reg_specs)) {
      spec <- reg_specs[[nm]]
      vars <- grep(paste0("^", spec$prefix, "\\["), all_vars, value = TRUE)
      if (length(vars) == 0) next
      cols <- spec$cols
      if (length(cols) == length(vars)) {
        terms <- paste0(nm, ": ", cols)
      } else {
        terms <- vars
      }
      map_rows[[length(map_rows) + 1L]] <- data.frame(term = terms, variable = vars, stringsAsFactors = FALSE)
    }
    if (length(map_rows) > 0) map <- do.call(rbind, map_rows)
  } else if (what == "raw") {
    vars <- variable %||% all_vars
    vars <- vars[vars %in% all_vars]
    map <- data.frame(term = vars, variable = vars, stringsAsFactors = FALSE)
  }

  if (!is.null(variable)) {
    map <- map[map$variable %in% variable, , drop = FALSE]
  }
  if (!is.null(term)) {
    map <- map[map$term %in% term, , drop = FALSE]
  }

  if (nrow(map) == 0) {
    cli::cli_abort(c(
      x = "No matching variables found for extraction.",
      i = "Check {.arg what}, {.arg term}, and {.arg variable} filters."
    ))
  }

  # Preserve order and duplicates from the map for user-facing terms.
  var_unique <- unique(map$variable)
  dm <- .get_draws_matrix(fit, variables = var_unique, draws = draws, seed = seed)
  out <- matrix(NA_real_, nrow = nrow(dm), ncol = nrow(map))
  colnames(out) <- make.unique(map$term)
  for (j in seq_len(nrow(map))) {
    out[, j] <- dm[, map$variable[j]]
  }

  list(
    draws = out,
    term_map = map
  )
}

#' Extract draw payloads from dynamic prediction objects
#'
#' @param object A `JoinMeDynPred` object.
#' @param what Draw block selector: `"longitudinal"`, `"longitudinal_fitted"`,
#'   `"survival"`, `"cumhaz"`, `"random_effects_id"`, `"random_effects_marker_id"`.
#' @param id Optional character/integer id filter.
#' @param scale Optional scale filter for longitudinal blocks (`epred`, `linpred`, `predict`).
#'
#' @return A list with fields:
#'   - `draws`: numeric matrix or list of matrices
#'   - `meta`: extraction metadata
#' @export
extract.JoinMeDynPred <- function(object,
                                  what = c("longitudinal", "longitudinal_fitted", "survival", "cumhaz", "random_effects_id", "random_effects_marker_id"),
                                  id = NULL,
                                  scale = NULL,
                                  ...) {
  what <- match.arg(what)
  dd <- object$draws[[what]]
  if (is.null(dd) || length(dd) == 0) {
    cli::cli_abort(c(
      x = "No stored draws available for {.val {what}}.",
      i = "Generate predictions with draw outputs enabled."
    ))
  }

  ids <- names(dd)
  if (!is.null(id)) {
    ids <- intersect(ids, as.character(id))
    dd <- dd[ids]
  }

  # For random-effect payloads, return raw structured content with optional id filter.
  if (what %in% c("random_effects_id", "random_effects_marker_id")) {
    return(list(draws = dd, meta = list(what = what, ids = names(dd))))
  }

  # Flatten longitudinal/survival payloads into matrices with informative column names.
  flatten_one <- function(entry, id_label) {
    if (is.null(entry)) return(NULL)

    # Multi-scale longitudinal payload (list keyed by scale).
    if (is.list(entry) && !is.null(names(entry)) && any(names(entry) %in% c("epred", "linpred", "predict"))) {
      out_scale <- list()
      keep_scales <- if (is.null(scale)) names(entry) else intersect(names(entry), scale)
      for (sc in keep_scales) {
        e <- entry[[sc]]
        if (is.null(e$matrix)) next
        m <- e$matrix
        tt <- e$time %||% seq_len(ncol(m))
        mk <- e$marker_idx %||% rep(NA_integer_, ncol(m))
        colnames(m) <- paste0("id=", id_label, "|scale=", sc, "|marker_idx=", mk, "|time=", signif(tt, 6))
        out_scale[[sc]] <- m
      }
      return(out_scale)
    }

    if (!is.null(entry$matrix)) {
      m <- entry$matrix
      tt <- entry$time %||% seq_len(ncol(m))
      mk <- entry$marker_idx %||% rep(NA_integer_, ncol(m))
      sc <- entry$scale %||% what
      colnames(m) <- paste0("id=", id_label, "|scale=", sc, "|marker_idx=", mk, "|time=", signif(tt, 6))
      return(m)
    }

    NULL
  }

  mats <- lapply(names(dd), function(idi) flatten_one(dd[[idi]], idi))
  names(mats) <- names(dd)
  mats <- mats[!vapply(mats, is.null, logical(1))]

  list(
    draws = mats,
    meta = list(what = what, ids = names(mats), scale = scale)
  )
}
