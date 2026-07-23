# ---- fixef / ranef --------------------------------------------------------

#' @keywords internal
.pivot_long_from_matrix <- function(draw_matrix, meta) {
  if (is.null(draw_matrix) || !is.matrix(draw_matrix) || ncol(draw_matrix) == 0L || nrow(meta) == 0L) {
    return(data.frame())
  }
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta_rep <- meta[rep(seq_len(nrow(meta)), each = nrow(draw_matrix)), , drop = FALSE]
  data.frame(
    draw = rep(seq_len(nrow(draw_matrix)), times = ncol(draw_matrix)),
    meta_rep,
    value = as.numeric(draw_matrix),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @keywords internal
.joinme_summarize_value <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(
      Estimate = NA_real_, Est.Error = NA_real_, Q2.5 = NA_real_, Q97.5 = NA_real_,
      Rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_
    ))
  }
  qs <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
  c(
    Estimate = mean(values),
    Est.Error = stats::sd(values),
    Q2.5 = qs[1],
    Q97.5 = qs[2],
    Rhat = suppressWarnings(tryCatch(as.numeric(posterior::rhat(values)), error = function(e) NA_real_)),
    ess_bulk = suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(values)), error = function(e) NA_real_)),
    ess_tail = suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(values)), error = function(e) NA_real_))
  )
}

#' @keywords internal
.joinme_summarize_long_draws <- function(draws_df, group_cols, digits = 3, value_col = "value") {
  if (is.null(draws_df) || !nrow(draws_df)) {
    return(NULL)
  }
  group_cols <- unique(c(group_cols, value_col))
  split_key <- interaction(draws_df[group_cols[group_cols != value_col]], drop = TRUE, lex.order = TRUE)
  idx_split <- split(seq_len(nrow(draws_df)), split_key)
  out_rows <- lapply(idx_split, function(idx) {
    base_row <- draws_df[idx[1], setdiff(group_cols, value_col), drop = FALSE]
    stats_row <- .joinme_summarize_value(draws_df[[value_col]][idx])
    cbind(base_row, as.data.frame(as.list(stats_row), stringsAsFactors = FALSE), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out_rows)
  out$Estimate <- round(out$Estimate, digits)
  out$Est.Error <- round(out$Est.Error, digits)
  out$Q2.5 <- round(out$Q2.5, digits)
  out$Q97.5 <- round(out$Q97.5, digits)
  out$Rhat <- round(out$Rhat, 3)
  rownames(out) <- NULL
  out
}

#' @keywords internal
.rescale_time <- function(stan_data) {
  scale_factor <- suppressWarnings(as.numeric(stan_data$tmax %||% 1.0))
  if (!is.finite(scale_factor) || length(scale_factor) != 1L || scale_factor <= 0) {
    return(1.0)
  }
  scale_factor
}

#' @keywords internal
.rescale_public_longitudinal_draws <- function(draws_df, stan_data, idx, term_labels) {
  if (is.null(draws_df) || !nrow(draws_df)) return(draws_df)

  scale_factor <- .rescale_time(stan_data)
  if (abs(scale_factor - 1.0) < 1e-12) return(draws_df)

  idx <- as.integer(idx %||% integer(0))
  idx <- idx[is.finite(idx) & idx >= 1L]
  if (!length(idx)) return(draws_df)

  term_labels <- as.character(term_labels %||% character(0))
  if (length(term_labels) < max(idx)) return(draws_df)

  time_terms <- unique(term_labels[idx])
  keep <- draws_df$term %in% time_terms
  if (!any(keep)) return(draws_df)

  draws_df$value[keep] <- draws_df$value[keep] / scale_factor
  draws_df
}

#' @keywords internal
.rescale_public_longitudinal_summary <- function(summary_df, stan_data, idx, term_labels) {
  if (is.null(summary_df) || !nrow(summary_df)) return(summary_df)

  scale_factor <- .rescale_time(stan_data)
  if (abs(scale_factor - 1.0) < 1e-12) return(summary_df)

  idx <- as.integer(idx %||% integer(0))
  idx <- idx[is.finite(idx) & idx >= 1L]
  if (!length(idx)) return(summary_df)

  term_labels <- as.character(term_labels %||% character(0))
  if (length(term_labels) < max(idx)) return(summary_df)

  time_terms <- unique(term_labels[idx])
  keep <- summary_df$term %in% time_terms
  if (!any(keep)) return(summary_df)

  value_cols <- intersect(c("Estimate", "Est.Error", "Q2.5", "Q97.5"), names(summary_df))
  if (length(value_cols) == 0L) return(summary_df)

  summary_df[keep, value_cols] <- lapply(summary_df[keep, value_cols, drop = FALSE], function(x) x / scale_factor)
  summary_df
}

#' @keywords internal
.rescale_public_event_draws <- function(draws_df, stan_data) {
  if (is.null(draws_df) || !nrow(draws_df)) return(draws_df)
  if (!("term" %in% names(draws_df)) || !("value" %in% names(draws_df))) return(draws_df)

  scale_factor <- .rescale_time(stan_data)
  if (abs(scale_factor - 1.0) < 1e-12) return(draws_df)

  idx <- as.integer(stan_data$idx_time_gamma %||% integer(0))
  idx <- idx[is.finite(idx) & idx >= 1L]
  term_labels <- as.character(stan_data$w_cols %||% character(0))
  if (!length(idx) || length(term_labels) < max(idx)) return(draws_df)

  time_terms <- unique(term_labels[idx])
  keep <- draws_df$term %in% time_terms
  if (!any(keep)) return(draws_df)

  draws_df$value[keep] <- draws_df$value[keep] / scale_factor
  draws_df
}

#' @keywords internal
.id_labels <- function(object, n_id) {
  ids <- object$dataLong$id %||% object$dataEvent$id %||% seq_len(n_id)
  ids <- unique(as.character(ids))
  if (length(ids) < n_id) {
    ids <- c(ids, as.character(seq_len(n_id - length(ids)) + length(ids)))
  }
  ids[seq_len(n_id)]
}

#' @keywords internal
.posterior_fixef_matrix <- function(object, draws = NULL, seed = 1) {
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  beta_map <- .fixed_effect_var_map(
    sd = sd,
    all_vars = all_vars,
    term_labels = .fit_design_term_labels(object, what = "fixef", n_terms = sd$P)
  )
  beta_vars <- as.character(beta_map$variable)
  beta_terms <- as.character(beta_map$term)

  mats <- list()
  labels <- character(0)
  if (length(beta_vars) > 0L) {
    mats[[length(mats) + 1L]] <- .get_draws_matrix(fit, variables = beta_vars, draws = draws, seed = seed)
    labels <- c(labels, as.character(beta_terms))
  }

  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)
  if (isTRUE(show_marker_weights) && isTRUE((sd$D %||% 0L) > 0L)) {
    marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
    shared_weights <- isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)
    weight_term_keys <- .active_weighted_assoc_terms(sd)
    if (shared_weights && length(weight_term_keys) > 1L) weight_term_keys <- weight_term_keys[1L]
    for (term_key in weight_term_keys) {
      mw_vars <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(sd$D), "]")
      mw_vars <- mw_vars[mw_vars %in% all_vars]
      if (!length(mw_vars)) next
      mats[[length(mats) + 1L]] <- .get_draws_matrix(fit, variables = mw_vars, draws = draws, seed = seed)
      labels <- c(labels, vapply(marker_terms[seq_along(mw_vars)], function(marker_label) {
        .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
      }, character(1)))
    }
  }

  if (!length(mats)) {
    return(matrix(0, nrow = 0L, ncol = 0L))
  }
  out <- do.call(cbind, lapply(mats, as.matrix))
  colnames(out) <- labels
  out
}

#' @keywords internal
.posterior_event_coef_draws <- function(object, draws = NULL, seed = 1) {
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  if (!isTRUE((sd$p_w %||% 0L) > 0L)) return(NULL)

  g_vars <- paste0("gamma_w[", seq_len(sd$p_w), "]")
  g_vars <- g_vars[g_vars %in% all_vars]
  if (length(g_vars) == 0L) {
    k_event <- sd$K_event %||% 1L
    g_vars <- as.vector(outer(seq_len(k_event), seq_len(sd$p_w), function(k, j) paste0("gamma_w[", k, ",", j, "]")))
    g_vars <- g_vars[g_vars %in% all_vars]
  }
  if (!length(g_vars)) return(NULL)

  dmat <- .get_draws_matrix(fit, variables = g_vars, draws = draws, seed = seed)
  parse_idx <- regmatches(g_vars, regexec("^gamma_w\\[(\\d+)(?:,(\\d+))?\\]$", g_vars))
  k_idx <- vapply(parse_idx, function(x) if (length(x) >= 2L) as.integer(x[2]) else 1L, integer(1))
  j_idx <- vapply(parse_idx, function(x) if (length(x) >= 3L && nzchar(x[3])) as.integer(x[3]) else as.integer(x[2]), integer(1))
  term_labels <- sd$w_cols %||% paste0("w_", seq_len(sd$p_w))
  if (length(term_labels) < max(j_idx)) {
    term_labels <- c(term_labels, paste0("w_", seq.int(length(term_labels) + 1L, max(j_idx))))
  }
  meta <- data.frame(
    event = if ((sd$K_event %||% 1L) > 1L) paste0("event", k_idx) else rep("event", length(g_vars)),
    term = as.character(term_labels[j_idx]),
    stringsAsFactors = FALSE
  )
  out <- .pivot_long_from_matrix(dmat, meta)
  .rescale_public_event_draws(out, sd)
}

#' @keywords internal
.posterior_distreg_draws <- function(object, draws = NULL, seed = 1) {
  fit <- object$fit
  cfg <- object$config
  if (is.null(draws)) draws <- object$config$draws_default
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  dist_cols <- cfg$dist$dist_cols %||% list()
  specs <- list(
    sigma = list(prefix = "beta_sigma", cols = dist_cols$sigma %||% character(0)),
    nu = list(prefix = "beta_nu", cols = dist_cols$nu %||% character(0)),
    phi = list(prefix = "beta_phi", cols = dist_cols$phi %||% character(0)),
    alpha = list(prefix = "beta_alpha", cols = dist_cols$alpha %||% character(0)),
    kappa = list(prefix = "beta_kappa", cols = dist_cols$kappa %||% character(0)),
    tau = list(prefix = "beta_tau", cols = dist_cols$tau %||% character(0))
  )

  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    vars <- grep(paste0("^", spec$prefix, "\\["), all_vars, value = TRUE)
    if (!length(vars)) next
    dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = seed)
    term_labels <- if (length(spec$cols) == length(vars)) spec$cols else vars
    out[[nm]] <- .pivot_long_from_matrix(dmat, data.frame(term = term_labels, stringsAsFactors = FALSE))
  }
  out
}

#' @keywords internal
.posterior_fit_ranef <- function(object, draws = NULL, seed = 1) {
  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config
  if (is.null(draws)) draws <- object$config$draws_default
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  n_id <- as.integer(sd$n_id %||% 0L)
  id_labels <- .id_labels(object, n_id)
  marker_labels <- as.character(sd$marker_levels %||% paste0("marker_", seq_len(as.integer(sd$D %||% 0L))))

  extract_grouped_draws <- function(var_builder, meta_builder) {
    meta <- meta_builder()
    if (is.null(meta) || !nrow(meta)) return(NULL)
    vars <- meta$variable
    keep <- vars %in% all_vars
    meta <- meta[keep, , drop = FALSE]
    if (!nrow(meta)) return(NULL)
    dmat <- .get_draws_matrix(fit, variables = meta$variable, draws = draws, seed = seed)
    meta$variable <- NULL
    .pivot_long_from_matrix(dmat, meta)
  }

  out_long <- list(
    id = if (isTRUE((sd$R_id %||% 0L) > 0L) && isTRUE(n_id > 0L)) extract_grouped_draws(
      NULL,
      function() {
        grid <- expand.grid(id_index = seq_len(n_id), term_index = seq_len(sd$R_id), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        data.frame(
          variable = paste0("u_id[", grid$id_index, ",", grid$term_index, "]"),
          id = id_labels[grid$id_index],
          term = as.character((sd$zid_cols %||% paste0("id_re_", seq_len(sd$R_id)))[grid$term_index]),
          stringsAsFactors = FALSE
        )
      }
    ) else NULL,
    marker = if (isTRUE((sd$R_mk %||% 0L) > 0L) && isTRUE((sd$D %||% 0L) > 0L)) extract_grouped_draws(
      NULL,
      function() {
        grid <- expand.grid(marker_index = seq_len(sd$D), term_index = seq_len(sd$R_mk), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        data.frame(
          variable = paste0("v_marker[", grid$marker_index, ",", grid$term_index, "]"),
          marker = marker_labels[grid$marker_index],
          term = as.character((sd$zmk_cols %||% paste0("marker_re_", seq_len(sd$R_mk)))[grid$term_index]),
          stringsAsFactors = FALSE
        )
      }
    ) else NULL,
    marker_by_id = if (isTRUE((sd$Q_idm %||% 0L) > 0L) && isTRUE(n_id > 0L) && isTRUE((sd$D %||% 0L) > 0L)) extract_grouped_draws(
      NULL,
      function() {
        grid <- expand.grid(id_index = seq_len(n_id), marker_index = seq_len(sd$D), term_index = seq_len(sd$Q_idm), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        data.frame(
          variable = paste0("w_idm[", grid$id_index, ",", grid$marker_index, ",", grid$term_index, "]"),
          id = id_labels[grid$id_index],
          marker = marker_labels[grid$marker_index],
          term = as.character((sd$zidm_cols %||% paste0("id_marker_re_", seq_len(sd$Q_idm)))[grid$term_index]),
          stringsAsFactors = FALSE
        )
      }
    ) else NULL
  )

  out_long$id <- .rescale_public_longitudinal_draws(
    out_long$id,
    sd,
    sd$idx_time_uid,
    sd$zid_cols %||% character(0)
  )
  out_long$marker <- .rescale_public_longitudinal_draws(
    out_long$marker,
    sd,
    sd$idx_time_vmk,
    sd$zmk_cols %||% character(0)
  )
  out_long$marker_by_id <- .rescale_public_longitudinal_draws(
    out_long$marker_by_id,
    sd,
    sd$idx_time_idm,
    sd$zidm_cols %||% character(0)
  )

  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)
  if (isTRUE(show_marker_weights) && isTRUE((sd$D %||% 0L) > 0L)) {
    shared_weights <- isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)
    weight_term_keys <- .active_weighted_assoc_terms(sd)
    if (shared_weights && length(weight_term_keys) > 1L) weight_term_keys <- weight_term_keys[1L]
    weight_rows <- list()
    for (term_key in weight_term_keys) {
      vars <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(sd$D), "]")
      vars <- vars[vars %in% all_vars]
      if (!length(vars)) next
      dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = seed)
      meta <- data.frame(
        term = vapply(marker_labels[seq_along(vars)], function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      weight_rows[[length(weight_rows) + 1L]] <- .pivot_long_from_matrix(dmat, meta)
    }
    if (length(weight_rows) > 0L) out_long$assoc_weight <- do.call(rbind, weight_rows)
  }

  summarize_dist_draws <- function(param_name) {
    n_re <- as.integer(sd[[paste0("n_re_", param_name)]] %||% 0L)
    if (n_re <= 0L) return(NULL)
    K_vec <- as.integer(sd[[paste0("K_", param_name)]] %||% integer(0))
    G_vec <- as.integer(sd[[paste0("G_", param_name)]] %||% integer(0))
    if (length(K_vec) != n_re || length(G_vec) != n_re) return(NULL)

    tau_prefix <- paste0("tau_", param_name)
    z_prefix <- paste0("z_", param_name)
    term_labels <- cfg$dist$dist_re_terms[[param_name]] %||% rep(NA_character_, n_re)
    value_cols <- list()
    meta_rows <- list()
    idx <- 1L
    required_vars <- character(0)
    for (j in seq_len(n_re)) {
      required_vars <- c(required_vars, paste0(tau_prefix, "[", j, ",", seq_len(K_vec[j]), "]"))
      required_vars <- c(required_vars, as.vector(outer(seq_len(G_vec[j]), seq_len(K_vec[j]), function(g, k) {
        paste0(z_prefix, "[", j, ",", g, ",", k, "]")
      })))
    }
    required_vars <- unique(required_vars)
    required_vars <- required_vars[required_vars %in% all_vars]
    if (!length(required_vars)) return(NULL)
    dmat <- .get_draws_matrix(fit, variables = required_vars, draws = draws, seed = seed)

    for (j in seq_len(n_re)) {
      for (g in seq_len(G_vec[j])) {
        for (k in seq_len(K_vec[j])) {
          tau_nm <- paste0(tau_prefix, "[", j, ",", k, "]")
          z_nm <- paste0(z_prefix, "[", j, ",", g, ",", k, "]")
          if (!(tau_nm %in% colnames(dmat)) || !(z_nm %in% colnames(dmat))) next
          value_cols[[idx]] <- as.numeric(dmat[, tau_nm]) * as.numeric(dmat[, z_nm])
          meta_rows[[idx]] <- data.frame(
            group = as.character(g),
            term = term_labels[j] %||% paste0("re_term_", j),
            coefficient = k,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
    if (!length(value_cols)) return(NULL)
    out <- .pivot_long_from_matrix(do.call(cbind, value_cols), do.call(rbind, meta_rows))
    out$scope <- ifelse(grepl("^family=", out$term), sub("^family=([^:]+)::.*$", "\\1", out$term), "allFamilies")
    out
  }

  out_dist <- list(
    sigma = summarize_dist_draws("sigma"),
    nu = summarize_dist_draws("nu"),
    phi = summarize_dist_draws("phi"),
    alpha = summarize_dist_draws("alpha"),
    kappa = summarize_dist_draws("kappa"),
    tau = summarize_dist_draws("tau")
  )
  out_dist <- out_dist[!vapply(out_dist, is.null, logical(1))]
  if (length(out_dist) > 0L) {
    out_dist <- lapply(out_dist, function(df) split(df[, setdiff(names(df), "scope"), drop = FALSE], df$scope))
  }

  out_vcov <- NULL
  q_idm <- as.integer(sd$Q_idm %||% 0L)
  if (q_idm > 0L) {
    m_cov <- if (as.integer(sd$indep_idmarker_cov %||% 0L) == 1L) q_idm else (q_idm * (q_idm + 1L)) %/% 2L
    z_vars <- as.vector(outer(seq_len(n_id), seq_len(m_cov), function(i, m) paste0("z_L[", i, ",", m, "]")))
    z_vars <- z_vars[z_vars %in% all_vars]
    lambda_vars <- paste0("lambda_L[", seq_len(m_cov), "]")
    lambda_vars <- lambda_vars[lambda_vars %in% all_vars]
    if (length(z_vars) > 0L && length(lambda_vars) > 0L) {
      dmat <- .get_draws_matrix(fit, variables = unique(c(z_vars, lambda_vars)), draws = draws, seed = seed)
      rc_map <- matrix(NA_integer_, nrow = m_cov, ncol = 2)
      if (as.integer(sd$indep_idmarker_cov %||% 0L) == 1L) {
        for (m in seq_len(m_cov)) rc_map[m, ] <- c(m, m)
      } else {
        pos <- 1L
        for (r in seq_len(q_idm)) {
          for (c in seq_len(r)) {
            rc_map[pos, ] <- c(r, c)
            pos <- pos + 1L
          }
        }
      }
      value_cols <- list()
      meta_rows <- list()
      idx <- 1L
      for (i in seq_len(n_id)) {
        for (m in seq_len(m_cov)) {
          z_nm <- paste0("z_L[", i, ",", m, "]")
          lambda_nm <- paste0("lambda_L[", m, "]")
          if (!(z_nm %in% colnames(dmat)) || !(lambda_nm %in% colnames(dmat))) next
          value_cols[[idx]] <- as.numeric(dmat[, z_nm]) * as.numeric(dmat[, lambda_nm])
          meta_rows[[idx]] <- data.frame(
            id = id_labels[i],
            block = ifelse(rc_map[m, 1] == rc_map[m, 2], "SD[id:marker]", "K[id:marker]"),
            row = rc_map[m, 1],
            col = rc_map[m, 2],
            term = "(Intercept)",
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
      if (length(value_cols) > 0L) {
        out_vcov <- .pivot_long_from_matrix(do.call(cbind, value_cols), do.call(rbind, meta_rows))
      }
    }
  }

  list(
    formulaLong = out_long,
    formulaDist = out_dist,
    formulaVCov = out_vcov
  )
}

#' @keywords internal
.combine_fixed_random_long_draws <- function(random_draws, fixed_matrix) {
  if (is.null(random_draws) || !nrow(random_draws)) return(NULL)
  fixed_value <- numeric(nrow(random_draws))
  if (!is.null(fixed_matrix) && is.matrix(fixed_matrix) && ncol(fixed_matrix) > 0L && nrow(fixed_matrix) > 0L) {
    term_match <- match(random_draws$term, colnames(fixed_matrix))
    keep <- !is.na(term_match) & random_draws$draw <= nrow(fixed_matrix)
    if (any(keep)) {
      fixed_value[keep] <- fixed_matrix[cbind(random_draws$draw[keep], term_match[keep])]
    }
  }
  random_draws$fixed <- fixed_value
  random_draws$random <- random_draws$value
  random_draws$value <- random_draws$fixed + random_draws$random
  random_draws
}

#' @keywords internal
.posterior_fit_coef <- function(object, draws = NULL, seed = 1) {
  fixed_long <- .posterior_fixef_matrix(object, draws = draws, seed = seed)
  beta_only <- fixed_long
  if (is.matrix(beta_only) && ncol(beta_only) > 0L) {
    weight_cols <- grepl("^weight", colnames(beta_only))
    if (any(weight_cols)) beta_only <- beta_only[, !weight_cols, drop = FALSE]
  }
  ranef_draws <- .posterior_fit_ranef(object, draws = draws, seed = seed)
  dist_fixed <- .posterior_distreg_draws(object, draws = draws, seed = seed)
  event_fixed <- .posterior_event_coef_draws(object, draws = draws, seed = seed)

  out_long <- ranef_draws$formulaLong
  out_long$id <- .combine_fixed_random_long_draws(out_long$id, beta_only)
  out_long$marker <- .combine_fixed_random_long_draws(out_long$marker, beta_only)
  out_long$marker_by_id <- .combine_fixed_random_long_draws(out_long$marker_by_id, beta_only)
  out_long$population <- if (is.matrix(beta_only) && ncol(beta_only) > 0L) {
    .pivot_long_from_matrix(beta_only, data.frame(term = colnames(beta_only), stringsAsFactors = FALSE))
  } else NULL

  out_dist <- list()
  if (length(dist_fixed) > 0L) {
    out_dist$population <- dist_fixed
  }
  if (length(ranef_draws$formulaDist) > 0L) {
    out_dist$group_specific <- lapply(names(ranef_draws$formulaDist), function(param_name) {
      scope_list <- ranef_draws$formulaDist[[param_name]]
      fixed_mat <- dist_fixed[[param_name]]
      fixed_wide <- if (!is.null(fixed_mat) && nrow(fixed_mat) > 0L) {
        reshape(
          fixed_mat[, c("draw", "term", "value")],
          idvar = "draw", timevar = "term", direction = "wide"
        )
      } else NULL
      if (!is.null(fixed_wide)) {
        draw_index <- fixed_wide$draw
        fixed_wide$draw <- NULL
        fixed_wide <- as.matrix(fixed_wide)
        colnames(fixed_wide) <- sub("^value\\.", "", colnames(fixed_wide))
        fixed_wide <- fixed_wide[order(draw_index), , drop = FALSE]
      }
      lapply(scope_list, .combine_fixed_random_long_draws, fixed_matrix = fixed_wide)
    })
    names(out_dist$group_specific) <- names(ranef_draws$formulaDist)
  }

  list(
    formulaLong = out_long,
    formulaEvent = event_fixed,
    formulaDist = out_dist,
    formulaVCov = list(
      population = NULL,
      id = ranef_draws$formulaVCov
    )
  )
}

#' Extract fixed effects
#'
#' @param object A JoiNMe fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param summary Logical. If `TRUE`, return posterior summaries with the same
#'   inferential columns used throughout the package. If `FALSE`, return the
#'   posterior draw matrix.
#' @param ... Unused.
#'
#' @return When `summary = TRUE`, a data.frame of posterior summaries. When
#'   `summary = FALSE`, a draws-by-term matrix.
#' @importFrom lme4 fixef
#' @export
fixef.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoiNMeFit"), msg = "Object must be a JoiNMeFit instance.")
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!isTRUE(summary)) {
    cache_key <- paste0("posterior_fixef_", draws)
    cached <- object$cache_get(cache_key)
    if (!is.null(cached)) return(cached)
    out <- .posterior_fixef_matrix(object, draws = draws, seed = seed)
    object$cache_set(cache_key, out)
    return(out)
  }

  cache_key <- paste0("fixef_", draws, "_", digits)
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  beta_map <- .fixed_effect_var_map(sd = sd, all_vars = all_vars)
  beta_vars <- as.character(beta_map$variable)
  if (length(beta_vars) == 0L) {
    beta_vars <- paste0("beta[", seq_len(sd$P), "]")
  }

  s <- as.data.frame(.summarise_draws_diag(fit, beta_vars, draws = draws, seed = seed))
  if (nrow(beta_map) > 0L) {
    idx <- match(s$variable, beta_map$variable)
    s$term <- beta_map$term[idx]
    s$term[is.na(s$term)] <- s$variable[is.na(s$term)]
  } else {
    recovered_terms <- .fit_design_term_labels(object, what = "fixef", n_terms = nrow(s))
    s$term <- if (length(recovered_terms) == nrow(s)) recovered_terms else s$variable
  }
  out <- s[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]

  out$Estimate <- round(out$Estimate, digits)
  out$Est.Error <- round(out$Est.Error, digits)
  out$Q2.5 <- round(out$Q2.5, digits)
  out$Q97.5 <- round(out$Q97.5, digits)
  out$Rhat <- round(out$Rhat, 3)

  # Append marker-weight summaries only for marker-weighted association models
  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)

  if (sd$D > 0 && show_marker_weights) {
    marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
    shared_weights <- isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)
    weight_term_keys <- .active_weighted_assoc_terms(sd)
    if (shared_weights && length(weight_term_keys) > 1L) {
      weight_term_keys <- weight_term_keys[1L]
    }
    weight_rows <- list()
    for (term_key in weight_term_keys) {
      mw_vars <- paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(sd$D), "]")
      mw_vars <- mw_vars[mw_vars %in% all_vars]
      if (length(mw_vars) == 0L) next
      mw <- as.data.frame(.summarise_draws_diag(fit, mw_vars, draws = draws, seed = seed))
      mw$term <- vapply(marker_terms[seq_len(nrow(mw))], function(marker_label) {
        .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
      }, character(1))
      mw <- mw[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      mw$Estimate <- round(mw$Estimate, digits)
      mw$Est.Error <- round(mw$Est.Error, digits)
      mw$Q2.5 <- round(mw$Q2.5, digits)
      mw$Q97.5 <- round(mw$Q97.5, digits)
      mw$Rhat <- round(mw$Rhat, 3)
      weight_rows[[length(weight_rows) + 1L]] <- mw
    }
    if (length(weight_rows) > 0L) {
      out <- rbind(out, do.call(rbind, weight_rows))
    }
  }
  object$cache_set(cache_key, out)
  out
}

#' Posterior fixed-effect alias for fitted JoiNMe models
#'
#' @param object A `JoiNMeFit` object.
#' @param ... Additional arguments forwarded to [fixef()].
#'
#' @return The same object returned by `fixef(object, summary = FALSE, ...)`.
#' @export
posterior_fixef <- function(object, ...) {
  fixef(object, summary = FALSE, ...)
}

#' Extract random effects
#'
#' @param object A JoiNMeFit fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param summary Logical. If `TRUE`, return posterior summaries. If `FALSE`,
#'   return the posterior extraction on the coefficient scale.
#' @param ... Unused.
#'
#' @return A nested list with top-level entries `formulaLong` and `formulaDist`.
#'   `formulaLong` contains random-effect summaries for longitudinal model
#'   components (`id`, `marker`, `marker_by_id_latent`, when present).
#'   `formulaDist` contains distributional random-effect summaries organised by
#'   parameter and family scope (e.g., `sigma$student_t`, `nu$allFamilies`).
#' @importFrom lme4 ranef
#' @export
ranef.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoiNMeFit"), msg = "Object must be a JoiNMeFit instance.")
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!isTRUE(summary)) {
    cache_key <- paste0("posterior_ranef_", draws)
    cached <- object$cache_get(cache_key)
    if (!is.null(cached)) return(cached)
    out <- .posterior_fit_ranef(object, draws = draws, seed = seed)
    object$cache_set(cache_key, out)
    return(out)
  }

  cache_key <- paste0("ranef_", draws, "_", digits)
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  summarize_vector <- function(values) {
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      return(c(Estimate = NA_real_, Est.Error = NA_real_, Q2.5 = NA_real_, Q97.5 = NA_real_, Rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_))
    }
    qs <- stats::quantile(values, probs = c(0.025, 0.975), names = FALSE)
    c(
      Estimate = mean(values),
      Est.Error = stats::sd(values),
      Q2.5 = qs[1],
      Q97.5 = qs[2],
      Rhat = suppressWarnings(tryCatch(as.numeric(posterior::rhat(values)), error = function(e) NA_real_)),
      ess_bulk = suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(values)), error = function(e) NA_real_)),
      ess_tail = suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(values)), error = function(e) NA_real_))
    )
  }

  summarize_block <- function(varnames, group, term_labels = NULL) {
    if (length(varnames) == 0) return(NULL)
    s <- as.data.frame(.summarise_draws_diag(fit, varnames, draws = draws, seed = seed))
    out <- s[, c("variable", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    names(out)[1] <- "variable"
    term_out <- out$variable
    term_labels <- as.character(term_labels %||% character(0))
    if (length(term_labels) > 0L) {
      term_index <- suppressWarnings(as.integer(sub("^.*,(\\d+)\\]$", "\\1", out$variable)))
      keep <- is.finite(term_index) & term_index >= 1L & term_index <= length(term_labels)
      term_out[keep] <- term_labels[term_index[keep]]
    }
    out$term <- term_out
    out$group <- group
    out$variable <- NULL
    out <- out[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail", "group"), drop = FALSE]

    out$Estimate <- round(out$Estimate, digits)
    out$Est.Error <- round(out$Est.Error, digits)
    out$Q2.5 <- round(out$Q2.5, digits)
    out$Q97.5 <- round(out$Q97.5, digits)
    out$Rhat <- round(out$Rhat, 3)
    out
  }

  u_vars <- vars[grepl("^u_id\\[", vars)]
  v_vars <- vars[grepl("^v_marker\\[", vars)]
  zw_vars <- vars[grepl("^z_w\\[", vars)]

  out_long <- list(
    id = summarize_block(u_vars, "id", term_labels = sd$zid_cols %||% character(0)),
    marker = if (sd$R_mk > 0) summarize_block(v_vars, "marker", term_labels = sd$zmk_cols %||% character(0)) else NULL,
    marker_by_id_latent = summarize_block(zw_vars, "marker_by_id_latent", term_labels = sd$zidm_cols %||% character(0))
  )

  out_long$id <- .rescale_public_longitudinal_summary(
    out_long$id,
    sd,
    sd$idx_time_uid,
    sd$zid_cols %||% character(0)
  )
  out_long$marker <- .rescale_public_longitudinal_summary(
    out_long$marker,
    sd,
    sd$idx_time_vmk,
    sd$zmk_cols %||% character(0)
  )

  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)

  weight_term_keys <- .active_weighted_assoc_terms(sd)
  shared_weights <- isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)
  if (shared_weights && length(weight_term_keys) > 1L) {
    weight_term_keys <- weight_term_keys[1L]
  }
  mw_vars <- unlist(lapply(weight_term_keys, function(term_key) {
    paste0(.marker_weight_var_prefix(term_key, effective = TRUE), "[", seq_len(sd$D), "]")
  }), use.names = FALSE)
  mw_vars <- mw_vars[mw_vars %in% vars]
  if (show_marker_weights && length(mw_vars) > 0) {
    mw <- summarize_block(mw_vars, "assoc_weight")
    if (!is.null(mw) && !is.null(sd$marker_levels)) {
      marker_terms <- unlist(lapply(weight_term_keys, function(term_key) {
        vapply(sd$marker_levels, function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1))
      }), use.names = FALSE)
      if (length(marker_terms) == nrow(mw)) {
        mw$term <- marker_terms
      }
    }
    out_long$assoc_weight <- mw
  } else if (show_marker_weights && isTRUE(as.logical(sd$fixed_marker_weights)) && !is.null(sd$marker_weights)) {
    marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
    mw <- data.frame(
      term = paste0("weight: ", marker_terms),
      Estimate = as.numeric(sd$marker_weights),
      Est.Error = NA_real_,
      Q2.5 = NA_real_,
      Q97.5 = NA_real_,
      Rhat = NA_real_,
      ess_bulk = NA_real_,
      ess_tail = NA_real_,
      group = "assoc_weight",
      stringsAsFactors = FALSE
    )
    out_long$assoc_weight <- mw
  }

  # Distributional random effects are represented as latent z and scale tau.
  # For each distributional parameter and RE term j, compute b = z * tau per
  # (group, coefficient) cell and summarize across draws.
  summarize_dist_ranef <- function(param_name) {
    n_re <- as.integer(sd[[paste0("n_re_", param_name)]] %||% 0L)
    if (n_re <= 0L) return(NULL)

    K_vec <- as.integer(sd[[paste0("K_", param_name)]] %||% integer(0))
    G_vec <- as.integer(sd[[paste0("G_", param_name)]] %||% integer(0))
    if (length(K_vec) != n_re || length(G_vec) != n_re) return(NULL)

    tau_prefix <- paste0("tau_", param_name)
    z_prefix <- paste0("z_", param_name)
    term_labels <- object$config$dist$dist_re_terms[[param_name]] %||% rep(NA_character_, n_re)

    required_vars <- character(0)
    for (j in seq_len(n_re)) {
      required_vars <- c(required_vars, paste0(tau_prefix, "[", j, ",", seq_len(K_vec[j]), "]"))
      required_vars <- c(required_vars, as.vector(outer(seq_len(G_vec[j]), seq_len(K_vec[j]), function(g, k) {
        paste0(z_prefix, "[", j, ",", g, ",", k, "]")
      })))
    }
    required_vars <- unique(required_vars)
    required_vars <- required_vars[required_vars %in% vars]
    if (length(required_vars) == 0) return(NULL)

    dmat <- .get_draws_matrix(fit, variables = required_vars, draws = draws, seed = seed)
    rows <- list()
    idx <- 1L
    for (j in seq_len(n_re)) {
      for (g in seq_len(G_vec[j])) {
        for (k in seq_len(K_vec[j])) {
          tau_nm <- paste0(tau_prefix, "[", j, ",", k, "]")
          z_nm <- paste0(z_prefix, "[", j, ",", g, ",", k, "]")
          if (!(tau_nm %in% colnames(dmat)) || !(z_nm %in% colnames(dmat))) next
          vals <- as.numeric(dmat[, tau_nm]) * as.numeric(dmat[, z_nm])
          ss <- summarize_vector(vals)
          rows[[idx]] <- data.frame(
            term_index = j,
            group = g,
            coefficient = k,
            term = term_labels[j] %||% paste0("re_term_", j),
            Estimate = ss[["Estimate"]],
            Est.Error = ss[["Est.Error"]],
            Q2.5 = ss[["Q2.5"]],
            Q97.5 = ss[["Q97.5"]],
            Rhat = ss[["Rhat"]],
            ess_bulk = ss[["ess_bulk"]],
            ess_tail = ss[["ess_tail"]],
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
    if (length(rows) == 0) return(NULL)
    out_df <- do.call(rbind, rows)
    out_df$Estimate <- round(out_df$Estimate, digits)
    out_df$Est.Error <- round(out_df$Est.Error, digits)
    out_df$Q2.5 <- round(out_df$Q2.5, digits)
    out_df$Q97.5 <- round(out_df$Q97.5, digits)
    out_df$Rhat <- round(out_df$Rhat, 3)
    out_df
  }

  split_dist_scopes <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    scope_names <- ifelse(grepl("^family=", df$term),
                          sub("^family=([^:]+)::.*$", "\\1", df$term),
                          "allFamilies")
    out <- split(df, scope_names)
    out <- out[order(names(out))]
    out
  }

  out_dist <- list(
    sigma = split_dist_scopes(summarize_dist_ranef("sigma")),
    nu = split_dist_scopes(summarize_dist_ranef("nu")),
    phi = split_dist_scopes(summarize_dist_ranef("phi")),
    alpha = split_dist_scopes(summarize_dist_ranef("alpha")),
    kappa = split_dist_scopes(summarize_dist_ranef("kappa")),
    tau = split_dist_scopes(summarize_dist_ranef("tau"))
  )
  out_dist <- out_dist[!vapply(out_dist, is.null, logical(1))]

  out <- list(
    formulaLong = out_long,
    formulaDist = out_dist
  )
  object$cache_set(cache_key, out)
  out
}

#' Posterior random-effect alias for fitted JoiNMe models
#'
#' @param object A `JoiNMeFit` object.
#' @param ... Additional arguments forwarded to [ranef()].
#'
#' @return The same object returned by `ranef(object, summary = FALSE, ...)`.
#' @export
posterior_ranef <- function(object, ...) {
  ranef(object, summary = FALSE, ...)
}

#' Combined posterior coefficients for fitted JoiNMe models
#'
#' @description
#' Returns posterior coefficients on the scale used by each model component.
#'
#' The guiding rule is simple:
#' for every coefficient carried by a group-specific model matrix, the returned
#' value is the sum of the population-level contribution and the matching
#' group-level deviation. When no group-level deviation exists, the returned
#' coefficient is the population-level coefficient itself.
#'
#' This mirrors the interpretation used in multilevel modelling:
#' a subject-specific or marker-specific coefficient is the coefficient that
#' would multiply the corresponding column of the model matrix for that unit.
#'
#' @param object A `JoiNMeFit` object.
#' @param draws Optional number of posterior draws to retain.
#' @param seed Integer seed used when subsetting posterior draws.
#' @param digits Number of digits used when `summary = TRUE`.
#' @param summary Logical. If `TRUE`, return posterior summaries. If `FALSE`,
#'   return posterior draw-level extractions.
#' @param ... Unused.
#'
#' @return When `summary = FALSE`, a nested list of draw-level data frames for
#'   the longitudinal, event, distributional, and covariance-regression parts of
#'   the model. When `summary = TRUE`, the same structure is returned after
#'   summarising each coefficient with posterior means, posterior uncertainty,
#'   interval estimates, and MCMC diagnostics.
#' @method coef JoiNMeFit
#' @export
coef.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoiNMeFit"), msg = "Object must be a JoiNMeFit instance.")
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  if (is.null(draws)) draws <- object$config$draws_default

  cache_key <- if (isTRUE(summary)) {
    paste0("coef_summary_", draws, "_", digits)
  } else {
    paste0("coef_draws_", draws)
  }
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  out <- .posterior_fit_coef(object, draws = draws, seed = seed)
  if (isTRUE(summary)) {
    summarize_component <- function(x) {
      if (is.null(x)) return(NULL)
      if (is.data.frame(x) && all(c("draw", "value") %in% names(x))) {
        keep_cols <- setdiff(names(x), c("draw", "value", "fixed", "random"))
        tbl <- .joinme_summarize_long_draws(x, group_cols = c(keep_cols, "value"), digits = digits)
        if ("event" %in% names(tbl)) {
          tbl$Hazard.Ratio <- round(exp(tbl$Estimate), digits)
          tbl$HR.Q2.5 <- round(exp(tbl$Q2.5), digits)
          tbl$HR.Q97.5 <- round(exp(tbl$Q97.5), digits)
        }
        return(tbl)
      }
      if (is.list(x)) return(lapply(x, summarize_component))
      x
    }
    out <- summarize_component(out)
  }

  object$cache_set(cache_key, out)
  out
}

#' Posterior coefficient alias for fitted JoiNMe models
#'
#' @param object A `JoiNMeFit` object.
#' @param ... Additional arguments forwarded to [coef()].
#'
#' @return The same object returned by `coef(object, summary = FALSE, ...)`.
#' @export
posterior_coef <- function(object, ...) {
  coef(object, summary = FALSE, ...)
}

#' Extract predicted random effects from dynamic predictions
#'
#' @description
#' Returns predicted random effects from a `JoiNMeDynPred` object. Marker-by-id
#' random effects are available only when marker covariance is configured to be
#' subject-dependent.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param ... Unused.
#'
#' @return A named list containing random-effects summary tables.
#' @importFrom lme4 ranef
#' @export
ranef.JoiNMeDynPred <- function(object, ...) {
    assertthat::assert_that(inherits(object, "JoiNMeDynPred"), msg = "Object must be a JoiNMeDynPred instance.")

    if (!isTRUE(object$metadata$marker_corr_depends_on_id)) {
        cli::cli_abort(c(
            x = "Predicted marker-by-id random effects are only available when marker covariance depends on id.",
            i = "Refit with subject-dependent covariance structure in {.arg formulaVCov} and marker-by-id random effects (Q_idm > 0)."
        ))
    }

    sum_obj <- summary(object)
    list(
        formulaLong = list(
            marker_by_id = sum_obj$tables$random_effects_marker_id
        )
    )
}
