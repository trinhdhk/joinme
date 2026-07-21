#' Extract posterior draws from JoiNMe objects
#'
#' @description
#' S3 generic to extract component-specific posterior payloads from
#' `JoiNMeFit` and `JoiNMeDynPred` objects.
#'
#' Compared with [draws()], `extract()` provides lower-level
#' interface. It works component by component and returns the metadata needed to
#' understand how a requested summary term maps back to the stored Stan
#' variables or prediction draw blocks.
#'
#' Use `extract()` when you need a specific model component, the corresponding
#' `term_map`, or a specialised payload such as `what = "association_plot"`.
#' Use [draws()] when you want a single posterior object ready for `posterior`
#' or `bayesplot` workflows.
#'
#' @param object A supported JoiNMe object.
#' @param ... Additional method-specific arguments.
#' @seealso [draws()]
#' @export
extract <- function(object, ...) {
  UseMethod("extract")
}

#' @keywords internal
.fixed_effect_var_map <- function(sd, all_vars) {
  p <- as.integer(sd$P %||% 0L)
  if (p <= 0L) {
    return(data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE))
  }

  beta_scaled_vars <- paste0("beta_scaled[", seq_len(p), "]")
  beta_vars <- paste0("beta[", seq_len(p), "]")

  if (all(beta_vars %in% all_vars)) {
    chosen <- beta_vars
  } else if (all(beta_scaled_vars %in% all_vars)) {
    chosen <- beta_scaled_vars
  } else {
    chosen <- beta_vars[beta_vars %in% all_vars]
  }

  if (!length(chosen)) {
    return(data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE))
  }

  idx <- suppressWarnings(as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", chosen)))
  term_labels <- as.character(sd$x_cols %||% paste0("beta_", seq_len(p)))
  if (length(term_labels) < max(idx, na.rm = TRUE)) {
    term_labels <- c(term_labels, paste0("beta_", seq.int(length(term_labels) + 1L, max(idx, na.rm = TRUE))))
  }

  data.frame(
    term = term_labels[idx],
    variable = chosen,
    stringsAsFactors = FALSE
  )
}

# File overview:
# - Extract draw matrices from JoiNMeFit by summary-like components.
# - Extract draw matrices from JoiNMeDynPred from stored draw components.

#' Build a user-facing term map for extracted JoiNMeFit draws
#'
#' @param object A `JoiNMeFit` object.
#' @param what Character component selector used by [extract.JoiNMeFit()].
#' @param all_vars Optional character vector of available posterior variable
#'   names. When omitted, they are read from the fitted object.
#'
#' @return A data frame with columns `term` and `variable`.
#' @keywords internal
#' @noRd
.fit_component_term_map <- function(
  object,
  what = c(
    "fixef", 
    "gamma_w",
    "assoc",
    "distributional",
    "distributional_regression",
    "likelihood_scale",
    "raw"),
  all_vars = NULL) {
  what <- match.arg(what)

  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config
  if (is.null(all_vars)) {
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  }

  map <- data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE)

  if (what == "fixef") {
    map <- .fixed_effect_var_map(sd = sd, all_vars = all_vars)
  } else if (what == "gamma_w") {
    g_vars <- paste0("gamma_w[", seq_len(sd$p_w %||% 0L), "]")
    g_vars <- g_vars[g_vars %in% all_vars]
    if (length(g_vars) == 0L) {
      k_event <- sd$K_event %||% 1L
      g_vars <- as.vector(outer(
        seq_len(k_event),
        seq_len(sd$p_w %||% 0L),
        function(k, j) paste0("gamma_w[", k, ",", j, "]")
      ))
      g_vars <- g_vars[g_vars %in% all_vars]
    }
    if (length(g_vars) > 0L) {
      var_idx <- regmatches(g_vars, regexec("^gamma_w\\[(\\d+)(?:,(\\d+))?\\]$", g_vars))
      k_idx <- vapply(var_idx, function(x) if (length(x) >= 2L) as.integer(x[2]) else NA_integer_, integer(1))
      j_idx <- vapply(var_idx, function(x) if (length(x) >= 3L && nzchar(x[3])) as.integer(x[3]) else as.integer(x[2]), integer(1))
      k_event <- sd$K_event %||% 1L
      g_terms <- sd$w_cols %||% g_vars
      if (length(g_terms) < max(j_idx %||% 0L, 0L)) {
        g_terms <- c(g_terms, paste0("w_", seq.int(length(g_terms) + 1L, max(j_idx))))
      }
      if (length(g_terms) > 0L && any(!is.na(j_idx))) {
        g_terms <- g_terms[j_idx]
      }
      if (k_event > 1L && any(!is.na(k_idx))) {
        g_terms <- paste0("event", k_idx, ": ", g_terms)
      }
      map <- data.frame(term = as.character(g_terms), variable = as.character(g_vars), stringsAsFactors = FALSE)
    }
  } else if (what == "assoc") {
    corr_assoc_vars <- grep("^alpha_corr\\[", all_vars, value = TRUE)
    vcov_assoc_vars <- grep("^alpha_vcov\\[", all_vars, value = TRUE)
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

    map <- data.frame(term = as.character(assoc_terms), variable = as.character(assoc_vars), stringsAsFactors = FALSE)

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
      kappa = list(prefix = "beta_kappa", cols = dist_cols$kappa %||% character(0)),
      tau = list(prefix = "beta_tau", cols = dist_cols$tau %||% character(0))
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

    beta_vars <- paste0("beta_scaled[", seq_len(sd$P), "]")
    beta_vars <- beta_vars[beta_vars %in% all_vars]
    if (length(beta_vars) == 0L) {
      beta_vars <- paste0("beta[", seq_len(sd$P), "]")
    }
    beta_vars <- beta_vars[beta_vars %in% all_vars]
    if (length(beta_vars) > 0) {
      beta_terms <- sd$x_cols %||% beta_vars
      if (length(beta_terms) != length(beta_vars)) beta_terms <- beta_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("beta_scaled: ", as.character(beta_terms)),
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
    map <- data.frame(term = all_vars, variable = all_vars, stringsAsFactors = FALSE)
  }

  map
}

#' Extract posterior draws from a fitted JoiNMe model
#'
#' @description
#' Extracts one fitted-model component at a time and returns both the draw-level
#' values and the mapping that produced them.
#'
#' This differs from [draws()] in two important ways:
#' - `extract()` keeps the request scoped to one semantic component such as
#'   fixed effects, survival coefficients, association terms, distributional
#'   terms, or likelihood-scale parameters.
#' - `extract()` returns a structured list with `draws` plus `term_map`
#'   (and for `what = "association_plot"`, additional plotting support data)
#'   instead of a single `posterior` draws object.
#'
#' In short, use `extract()` when you need component-aware extraction and use
#' [draws()] when you need one renamed posterior object for general downstream
#' analysis.
#'
#' @param object A `JoiNMeFit` object.
#' @param what Character component selector. One of
#'   `"fixef"`, `"gamma_w"`, `"basehaz"`,
#'   `"baseline_hazard"`, `"assoc"`, `"association_plot"`,
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
#' @seealso [draws()] for a higher-level interface that returns a single `posterior`
#' @export
extract.JoiNMeFit <- function(object,
                              what = c("fixef", "gamma_w", "basehaz", "baseline_hazard", "assoc", "association_plot", "distributional", "distributional_regression", "likelihood_scale", "raw"),
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

  if (what %in% c("basehaz", "baseline_hazard")) {
    sd <- object$stan_data
    fit <- object$fit
    k_event <- as.integer(sd$K_event %||% 1L)
    k_bs <- as.integer(sd$Kbs %||% 0L)
    tmax <- suppressWarnings(as.numeric(sd$tmax %||% 1.0))
    if (!is.finite(tmax) || length(tmax) != 1L || tmax <= 0) {
      tmax <- 1.0
    }
    b_event <- as.matrix(sd$Bs_event_c %||% matrix(0, 0, 0))
    if (k_bs <= 0L || nrow(b_event) == 0L || ncol(b_event) != k_bs) {
      cli::cli_abort(c(
        x = "No baseline-hazard basis is available for extraction.",
        i = "Fit a survival model with baseline hazard terms before requesting basehaz extraction."
      ))
    }

    bh_vars <- as.vector(outer(seq_len(k_event), seq_len(k_bs), function(k, j) paste0("bs_gamma_c[", k, ",", j, "]")))
    all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
    bh_vars <- bh_vars[bh_vars %in% all_vars]
    if (!length(bh_vars)) {
      cli::cli_abort("No baseline-hazard coefficient draws were found in the fitted object.")
    }

    event_rows <- seq_len(nrow(b_event))
    row_labels <- if (!is.null(object$dataEvent) && nrow(object$dataEvent) == nrow(b_event)) {
      ids <- as.character(object$dataEvent$id %||% event_rows)
      tvals <- as.numeric(object$dataEvent$time %||% rep(NA_real_, nrow(b_event)))
      if (k_event > 1L && !is.null(sd$event_type) && length(sd$event_type) == nrow(b_event)) {
        paste0("etype", sd$event_type, "|id=", ids, "|time=", signif(tvals, 6))
      } else {
        paste0("id=", ids, "|time=", signif(tvals, 6))
      }
    } else {
      paste0("event_row_", event_rows)
    }

    if (isTRUE(keep_chains)) {
      arr <- .get_draws_array(fit, variables = bh_vars, draws = draws, seed = seed)
      out <- array(NA_real_, dim = c(dim(arr)[1], dim(arr)[2], k_event * nrow(b_event)),
                   dimnames = list(iteration = dimnames(arr)[[1]], chain = dimnames(arr)[[2]], term = character(k_event * nrow(b_event))))
      col_pos <- 1L
      for (k in seq_len(k_event)) {
        k_vars <- paste0("bs_gamma_c[", k, ",", seq_len(k_bs), "]")
        k_idx <- match(k_vars, dimnames(arr)[[3]])
        if (any(is.na(k_idx))) next
        for (ch in seq_len(dim(arr)[2])) {
          coef_mat <- arr[, ch, k_idx, drop = FALSE]
          out[, ch, col_pos:(col_pos + nrow(b_event) - 1L)] <- exp(as.matrix(coef_mat) %*% t(b_event)) / tmax
        }
        prefix <- if (k_event > 1L) paste0("event", k, ":") else ""
        dimnames(out)[[3]][col_pos:(col_pos + nrow(b_event) - 1L)] <- paste0(prefix, row_labels)
        col_pos <- col_pos + nrow(b_event)
      }
      map <- data.frame(term = dimnames(out)[[3]], variable = rep("basehaz", length(dimnames(out)[[3]])), stringsAsFactors = FALSE)
      if (!is.null(term)) {
        keep <- map$term %in% term
        map <- map[keep, , drop = FALSE]
        out <- out[, , keep, drop = FALSE]
      }
      return(list(draws = out, term_map = map))
    }

    dm <- .get_draws_matrix(fit, variables = bh_vars, draws = draws, seed = seed)
    out <- matrix(NA_real_, nrow = nrow(dm), ncol = k_event * nrow(b_event))
    col_names <- character(0)
    col_pos <- 1L
    for (k in seq_len(k_event)) {
      k_vars <- paste0("bs_gamma_c[", k, ",", seq_len(k_bs), "]")
      k_idx <- match(k_vars, colnames(dm))
      if (any(is.na(k_idx))) next
      out[, col_pos:(col_pos + nrow(b_event) - 1L)] <- exp(as.matrix(dm[, k_idx, drop = FALSE]) %*% t(b_event)) / tmax
      prefix <- if (k_event > 1L) paste0("event", k, ":") else ""
      col_names <- c(col_names, paste0(prefix, row_labels))
      col_pos <- col_pos + nrow(b_event)
    }
    out <- out[, seq_along(col_names), drop = FALSE]
    colnames(out) <- make.unique(col_names)
    map <- data.frame(term = colnames(out), variable = rep("basehaz", ncol(out)), stringsAsFactors = FALSE)
    if (!is.null(term)) {
      keep <- map$term %in% term
      map <- map[keep, , drop = FALSE]
      out <- out[, keep, drop = FALSE]
    }
    return(list(draws = out, term_map = map))
  }

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  map <- .fit_component_term_map(object, what = what, all_vars = all_vars)

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

  .event_time_terms <- function(sd, mapped_terms) {
    idx <- as.integer(sd$idx_time_gamma %||% integer(0))
    idx <- idx[is.finite(idx) & idx >= 1L]
    labels <- as.character(sd$w_cols %||% character(0))
    if (!length(idx) || length(labels) < max(idx)) {
      return(character(0))
    }
    base_terms <- unique(labels[idx])
    if (isTRUE(as.integer(sd$K_event %||% 1L) > 1L)) {
      base_terms <- c(
        base_terms,
        as.vector(outer(seq_len(as.integer(sd$K_event %||% 1L)), base_terms, function(k, term_label) paste0("event", k, ": ", term_label)))
      )
    }
    intersect(as.character(mapped_terms %||% character(0)), base_terms)
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
    if (identical(what, "gamma_w")) {
      tmax <- suppressWarnings(as.numeric(sd$tmax %||% 1.0))
      if (is.finite(tmax) && length(tmax) == 1L && tmax > 0) {
        time_terms <- .event_time_terms(sd, map$term)
        if (length(time_terms) > 0L) {
          keep <- map$term %in% time_terms
          out[, , keep] <- out[, , keep, drop = FALSE] / tmax
        }
      }
    }
  } else {
    dm <- .get_draws_matrix(fit, variables = var_unique, draws = draws, seed = seed)
    out <- matrix(NA_real_, nrow = nrow(dm), ncol = nrow(map))
    colnames(out) <- make.unique(map$term)
    for (j in seq_len(nrow(map))) {
      out[, j] <- dm[, map$variable[j]]
    }
    if (identical(what, "gamma_w")) {
      tmax <- suppressWarnings(as.numeric(sd$tmax %||% 1.0))
      if (is.finite(tmax) && length(tmax) == 1L && tmax > 0) {
        time_terms <- .event_time_terms(sd, map$term)
        if (length(time_terms) > 0L) {
          keep <- map$term %in% time_terms
          out[, keep] <- out[, keep, drop = FALSE] / tmax
        }
      }
    }
  }

  list(
    draws = out,
    term_map = map
  )
}

#' Extract stored posterior draw components from dynamic prediction objects
#'
#' @description
#' Extracts stored prediction draw blocks without flattening them first.
#'
#' This is the structured companion to [draws.JoiNMeDynPred()]. Use
#' `extract()` when you want to keep the original prediction block semantics
#' (`longitudinal`, `survival`, `cumhaz`, random effects, and scale/id filters).
#' Use [draws()] when you want those blocks flattened into one
#' `posterior`-compatible draw object with composite variable labels.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param what Draw block selector: `"longitudinal"`, `"longitudinal_fitted"`,
#'   `"survival"`, `"cumhaz"`, `"random_effects_id"`, `"random_effects_marker_id"`.
#' @param id Optional character/integer id filter.
#' @param scale Optional scale filter for longitudinal blocks (`epred`, `linpred`, `predict`).
#'
#' @return A list with fields:
#'   - `draws`: numeric matrix or list of matrices
#'   - `meta`: extraction metadata
#' @export
extract.JoiNMeDynPred <- function(
  object,
  what = c(
    "longitudinal",
    "longitudinal_fitted",
    "survival",
    "cumhaz",
    "random_effects_id",
    "random_effects_marker_id"
  ),
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
