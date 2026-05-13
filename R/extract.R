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
# - Extract draw matrices from JoinMeDynPred from stored draw components.

#' Extract posterior draws from a fitted joinme model
#'
#' @param object A `JoinMeFit` object.
#' @param what Character component selector. One of
#'   `"fixef"`, `"gamma_w"`, `"assoc"`, `"association_plot"`,
#'   `"distributional"`, `"distributional_regression"`, `"likelihood_scale"`,
#'   or `"raw"`.
#' @param term Optional character vector of friendly term names (summary-style)
#'   to subset extracted columns.
#' @param variable Optional character vector of raw Stan variable names. This is
#'   used directly when `what = "raw"` and can also further filter mapped outputs.
#' @param draws Optional number of posterior draws to keep per chain.
#' @param seed Integer seed used when subsetting draws.
#' @param keep_chains Logical; if TRUE, return draws with chains in a separate
#'   dimension (iteration x chain x term). If FALSE, return a flattened
#'   draws-by-term matrix.
#'
#' @return A list with fields:
#'   - `draws`: numeric array when `keep_chains = TRUE` (iteration x chain x term),
#'     otherwise a numeric matrix (rows = draws, cols = requested terms). For
#'     `what = "association_plot"`, this is a named list of compact draw
#'     matrices keyed by association term.
#'   - `term_map`: data.frame mapping `term` to Stan `variable`
#'   - `support`: for `what = "association_plot"`, cached model-implied raw
#'     support ranges used by association plotting.
#' @export
extract.JoinMeFit <- function(object,
                              what = c("fixef", "gamma_w", "assoc", "association_plot", "distributional", "distributional_regression", "likelihood_scale", "raw"),
                              term = NULL,
                              variable = NULL,
                              draws = NULL,
                              seed = 1,
                              keep_chains = TRUE,
                              ...) {
  what <- match.arg(what)
  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  if (what == "association_plot") {
    data <- .get_association_plot_data(object, seed = seed)
    if (is.null(data)) {
      cli::cli_abort(c(
        x = "No association plotting data is available.",
        i = "Fit a model with association terms or refit with posterior draws available."
      ))
    }

    keep_terms <- term %||% names(data$coeff_draws %||% list())
    keep_terms <- intersect(keep_terms, names(data$coeff_draws %||% list()))
    term_map <- data$term_map %||% data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE)
    if (!is.null(term)) {
      term_map <- term_map[term_map$term %in% keep_terms, , drop = FALSE]
    }
    support <- data$support %||% data.frame()
    if (!is.null(term) && nrow(support) > 0) {
      support <- support[support$term %in% keep_terms, , drop = FALSE]
    }

    return(list(
      draws = data$coeff_draws[keep_terms],
      term_map = term_map,
      support = support,
      marker_weight_draws = data$marker_weight_draws,
      transform_coeff_draws = data$transform_coeff_draws,
      transform_specs = data$transform_specs
    ))
  }

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
    corr_assoc_vars <- grep("^alpha_corr_eff\\[", all_vars, value = TRUE)
    if (length(corr_assoc_vars) == 0L) {
      corr_assoc_vars <- grep("^alpha_corr\\[", all_vars, value = TRUE)
    }
    vcov_assoc_vars <- grep("^alpha_vcov_eff\\[", all_vars, value = TRUE)
    if (length(vcov_assoc_vars) == 0L) {
      vcov_assoc_vars <- grep("^alpha_vcov\\[", all_vars, value = TRUE)
    }
    assoc_vars <- c(
      if (isTRUE(sd$assoc_cv_total == 1)) "alpha_cv_total",
      if (isTRUE(sd$assoc_cv_mean == 1)) "alpha_cv_mean",
      if (isTRUE(sd$assoc_cv_marker == 1)) "alpha_cv_marker",
      if (isTRUE(sd$assoc_cs_total == 1)) "alpha_cs_total",
      if (isTRUE(sd$assoc_cs_mean == 1)) "alpha_cs_mean",
      if (isTRUE(sd$assoc_cs_marker == 1)) "alpha_cs_marker",
      if (isTRUE(sd$assoc_corr == 1)) corr_assoc_vars,
      if (isTRUE(sd$assoc_vcov == 1)) vcov_assoc_vars
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
    assoc_terms <- sub("_eff\\[", "[", assoc_terms, perl = TRUE)

    map <- data.frame(term = as.character(assoc_terms), variable = as.character(assoc_vars), stringsAsFactors = FALSE)

    # Marker-weight terms mirror summary naming when present.
    marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D %||% 0L))
    shared_weights <- isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)
    weight_term_keys <- .active_weighted_assoc_terms(sd)
    if (shared_weights && length(weight_term_keys) > 1L) {
      weight_term_keys <- weight_term_keys[1L]
    }
    mw_rows <- list()
    for (term_key in weight_term_keys) {
      mw_vars <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(sd$D %||% 0L), "]")
      mw_vars <- mw_vars[mw_vars %in% all_vars]
      if (length(mw_vars) == 0L && shared_weights) {
        mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D %||% 0L), "]")
        mw_vars <- mw_vars[mw_vars %in% all_vars]
      }
      if (length(mw_vars) == 0L) next
      this_marker_terms <- marker_terms
      if (length(this_marker_terms) != length(mw_vars)) this_marker_terms <- paste0("marker_", seq_along(mw_vars))
      mw_rows[[length(mw_rows) + 1L]] <- data.frame(
        term = vapply(this_marker_terms, function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1)),
        variable = as.character(mw_vars),
        stringsAsFactors = FALSE
      )
    }
    if (length(mw_rows) > 0L) {
      map <- rbind(map, do.call(rbind, mw_rows))
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
  } else if (what == "likelihood_scale") {
    map_rows <- list()

    beta_vars <- paste0("beta_eff_in_likelihood[", seq_len(sd$P), "]")
    beta_vars <- beta_vars[beta_vars %in% all_vars]
    if (length(beta_vars) > 0) {
      beta_terms <- sd$x_cols %||% beta_vars
      if (length(beta_terms) != length(beta_vars)) beta_terms <- beta_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("beta_eff: ", as.character(beta_terms)),
        variable = as.character(beta_vars),
        stringsAsFactors = FALSE
      )
    }

    tau_id_vars <- paste0("tau_u_eff[", seq_len(sd$R_id %||% 0L), "]")
    tau_id_vars <- tau_id_vars[tau_id_vars %in% all_vars]
    if (length(tau_id_vars) > 0) {
      tau_id_terms <- sd$zid_cols %||% tau_id_vars
      if (length(tau_id_terms) != length(tau_id_vars)) tau_id_terms <- tau_id_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("id_sd_eff: ", as.character(tau_id_terms)),
        variable = as.character(tau_id_vars),
        stringsAsFactors = FALSE
      )
    }

    tau_marker_vars <- paste0("tau_v_eff[", seq_len(sd$R_mk %||% 0L), "]")
    tau_marker_vars <- tau_marker_vars[tau_marker_vars %in% all_vars]
    if (length(tau_marker_vars) > 0) {
      tau_marker_terms <- sd$zmk_cols %||% tau_marker_vars
      if (length(tau_marker_terms) != length(tau_marker_vars)) tau_marker_terms <- tau_marker_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("marker_sd_eff: ", as.character(tau_marker_terms)),
        variable = as.character(tau_marker_vars),
        stringsAsFactors = FALSE
      )
    }

    row_scale_vars <- paste0("marker_id_row_scale_eff[", seq_len(sd$Q_idm %||% 0L), "]")
    row_scale_vars <- row_scale_vars[row_scale_vars %in% all_vars]
    if (length(row_scale_vars) > 0) {
      row_scale_terms <- sd$zidm_cols %||% row_scale_vars
      if (length(row_scale_terms) != length(row_scale_vars)) row_scale_terms <- row_scale_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("id_marker_row_scale_eff: ", as.character(row_scale_terms)),
        variable = as.character(row_scale_vars),
        stringsAsFactors = FALSE
      )
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
  if (isTRUE(keep_chains)) {
    arr <- .get_draws_array(fit, variables = var_unique, draws = draws, seed = seed)
    var_idx <- match(map$variable, dimnames(arr)[[3]])
    if (any(is.na(var_idx))) {
      cli::cli_abort("Requested variables not found in the draw array.")
    }
    out <- array(
      NA_real_,
      dim = c(dim(arr)[1], dim(arr)[2], nrow(map)),
      dimnames = list(
        iteration = dimnames(arr)[[1]],
        chain = dimnames(arr)[[2]],
        term = make.unique(map$term)
      )
    )
    for (j in seq_len(nrow(map))) {
      out[, , j] <- arr[, , var_idx[j]]
    }
  } else {
    dm <- .get_draws_matrix(fit, variables = var_unique, draws = draws, seed = seed)
    out <- matrix(NA_real_, nrow = nrow(dm), ncol = nrow(map))
    colnames(out) <- make.unique(map$term)
    for (j in seq_len(nrow(map))) {
      out[, j] <- dm[, map$variable[j]]
    }
  }

  list(
    draws = out,
    term_map = map
  )
}

#' Extract stored posterior draw components from dynamic prediction objects
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

  # For random-effect components, return raw structured content with optional id filter.
  if (what %in% c("random_effects_id", "random_effects_marker_id")) {
    return(list(draws = dd, meta = list(what = what, ids = names(dd))))
  }

  # Flatten longitudinal/survival draw components into matrices with informative column names.
  flatten_one <- function(entry, id_label) {
    if (is.null(entry)) return(NULL)

    # Multi-scale longitudinal draw block (list keyed by scale).
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
