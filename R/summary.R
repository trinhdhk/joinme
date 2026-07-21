#' Summary of a JoiNMe object
#'
#' @description
#' Builds posterior summaries for the longitudinal process, survival process,
#' association terms, and optional covariance blocks. For survival, the summary
#' reports `survival_process` when the event model
#' contains non-intercept covariates beyond association features.
#'
#' @param object A JoiNMe fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param include_corr Logical; include covariance summaries.
#' @param ... Unused.
#'
#' @return A `summary_JoiNMeFit` object 
#' @export
summary.JoiNMeFit <- function(object, draws = NULL, seed = 1, digits = 3,
                           include_corr = TRUE, ...) {

  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  if (is.null(draws)) draws <- cfg$draws_default
  assertthat::assert_that(is.null(draws) || (is.numeric(draws) && draws > 0),
                          msg = "draws must be NULL or a positive number.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- paste0("summary_draws=", draws, "_digits=", digits, "_corr=", as.integer(include_corr))
  # Cache by draw count + digits to avoid repeat summaries
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  diag <- .joinme_sampler_diagnostics(fit)
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  s_beta <- .extract_fit_summary(object, what = "fixef", draws = draws, seed = seed, digits = digits)

  g_ext <- tryCatch(
    extract.JoiNMeFit(object, what = "gamma_w", draws = draws, seed = seed, keep_chains = TRUE),
    error = function(e) NULL
  )
  s_g <- .summarise_named_draws(g_ext$draws %||% NULL, digits = digits)

  s_basehaz <- NULL
  k_event <- as.integer(sd$K_event %||% 1L)
  k_bs <- as.integer(sd$Kbs %||% 0L)
  if (k_bs > 0L) {
    bh_vars <- as.vector(outer(seq_len(k_event), seq_len(k_bs), function(k, j) paste0("bs_gamma_c[", k, ",", j, "]")))
    bh_vars <- bh_vars[bh_vars %in% all_vars]
    if (length(bh_vars) > 0L) {
      s_basehaz <- as.data.frame(.summarise_draws_diag(fit, bh_vars, draws = draws, seed = seed))
      var_idx <- regmatches(as.character(s_basehaz$variable), regexec("^bs_gamma_c\\[(\\d+),(\\d+)\\]$", as.character(s_basehaz$variable)))
      k_idx <- vapply(var_idx, function(x) as.integer(x[2]), integer(1))
      j_idx <- vapply(var_idx, function(x) as.integer(x[3]), integer(1))
      bh_terms <- as.character(sd$basehaz_cols %||% paste0("basis_", seq_len(k_bs)))
      if (length(bh_terms) < max(j_idx)) {
        bh_terms <- c(bh_terms, paste0("basis_", seq.int(length(bh_terms) + 1L, max(j_idx))))
      }
      term_labels <- bh_terms[j_idx]
      if (k_event > 1L) {
        term_labels <- paste0("event", k_idx, ": ", term_labels)
      }
      s_basehaz$term <- as.character(term_labels)
      s_basehaz <- s_basehaz[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      s_basehaz$Estimate <- round(s_basehaz$Estimate, digits)
      s_basehaz$Est.Error <- round(s_basehaz$Est.Error, digits)
      s_basehaz$Q2.5 <- round(s_basehaz$Q2.5, digits)
      s_basehaz$Q97.5 <- round(s_basehaz$Q97.5, digits)
      s_basehaz$Rhat <- round(s_basehaz$Rhat, 3)
    }
  }
  if (!is.null(s_basehaz) && nrow(s_basehaz) > 0L) {
    tmax <- suppressWarnings(as.numeric(sd$tmax %||% 1.0))
    if (!is.finite(tmax) || length(tmax) != 1L || tmax <= 0) {
      tmax <- 1.0
    }
    # browser()
    term_base <- sub("^event[0-9]+:\\s*", "", as.character(s_basehaz$term))
    # is_intercept <- tolower(trimws(term_base)) %in% c("(intercept)", "intercept", "1")
    is_intercept <- grepl(",1]", trimws(bh_vars), fixed=TRUE)
    if (any(is_intercept) && tmax != 1.0) {
      log_tmax <- log(tmax)
      s_basehaz$Estimate[is_intercept] <- s_basehaz$Estimate[is_intercept] - log_tmax
      s_basehaz$Q2.5[is_intercept] <- s_basehaz$Q2.5[is_intercept] - log_tmax
      s_basehaz$Q97.5[is_intercept] <- s_basehaz$Q97.5[is_intercept] - log_tmax
    }
    s_basehaz$Hazard.Ratio <- round(exp(s_basehaz$Estimate), digits)
    s_basehaz$HR.Q2.5 <- round(exp(s_basehaz$Q2.5), digits)
    s_basehaz$HR.Q97.5 <- round(exp(s_basehaz$Q97.5), digits)
    s_basehaz <- s_basehaz[, c(
      "term", "Estimate", "Hazard.Ratio", "Est.Error", "Q2.5", "Q97.5",
      "HR.Q2.5", "HR.Q97.5", "Rhat", "ess_bulk", "ess_tail"
    ), drop = FALSE]
  }

  # Survival-process report (non-association baseline covariates only)
  #
  # - If and only if baseline event covariates exist in the survival model,
  #   `survival_process` must be reported.
  # - These terms must match the event design columns (e.g., x1, x2) and must
  #   not be silently dropped.
  #
  # 1) Start from gamma_w summaries (the survival linear predictor terms).
  # 2) Re-label terms from standata event columns when available (`w_cols`).
  # 3) Keep only non-intercept baseline covariates.
  # 4) Add hazard-ratio summaries exp(beta) for direct interpretation.
  s_surv <- NULL
  if (!is.null(s_g) && nrow(s_g) > 0 && isTRUE((sd$p_w %||% 0L) > 0L)) {
    raw_terms <- as.character(s_g$term)
    is_intercept_term <- function(term) {
      term_chr <- trimws(as.character(term))
      tolower(term_chr) %in% c("(intercept)", "intercept", "1")
    }

    baseline_terms <- sub("^event[0-9]+:\\s*", "", raw_terms)

    keep_idx <- !vapply(baseline_terms, is_intercept_term, logical(1))
    if (!any(keep_idx) && length(baseline_terms) > 0) {
      # if intercept filtering drops everything, keep all
      # baseline terms rather than hiding survival covariates.
      keep_idx <- rep(TRUE, length(baseline_terms))
    }
    if (any(keep_idx)) {
      s_surv <- s_g[keep_idx, , drop = FALSE]
      s_surv$term <- baseline_terms[keep_idx]
      s_surv$Hazard.Ratio <- round(exp(s_surv$Estimate), digits)
      s_surv$HR.Q2.5 <- round(exp(s_surv$Q2.5), digits)
      s_surv$HR.Q97.5 <- round(exp(s_surv$Q97.5), digits)
      s_surv <- s_surv[, c(
        "term", "Estimate", "Hazard.Ratio", "Est.Error", "Q2.5", "Q97.5",
        "HR.Q2.5", "HR.Q97.5", "Rhat", "ess_bulk", "ess_tail"
      ), drop = FALSE]
    }
  }

  a_ext <- tryCatch(
    extract.JoiNMeFit(object, what = "assoc", draws = draws, seed = seed, keep_chains = TRUE),
    error = function(e) NULL
  )
  s_a <- .summarise_named_draws(a_ext$draws %||% NULL, digits = digits)

  # Marker-weight association summaries are shown only when marker-weighted
  # association terms are active (cv_total/cv_marker/cs_total/cs_marker).
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

    weight_tables <- lapply(weight_term_keys, function(term_key) {
      weight_arr <- .association_marker_weight_array(object, term_key = term_key, draws = draws, seed = seed, all_vars = all_vars)
      if (!is.null(weight_arr)) {
        labels <- vapply(marker_terms, function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1))
        return(.assoc_summary_from_draw_array(weight_arr, term_labels = labels, digits = digits))
      }

      base_by_term <- sd$marker_weights_by_term %||% list()
      base_weights <- as.numeric(base_by_term[[term_key]] %||% sd$marker_weights %||% rep(1, sd$D))
      if (length(base_weights) != sd$D) {
        return(NULL)
      }

      data.frame(
        term = vapply(marker_terms, function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1)),
        Estimate = round(base_weights, digits),
        Est.Error = NA_real_,
        Q2.5 = NA_real_,
        Q97.5 = NA_real_,
        Rhat = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    weight_tables <- Filter(Negate(is.null), weight_tables)
    if (length(weight_tables) > 0L) {
      s_mw <- do.call(rbind, weight_tables)
      if (is.null(s_a)) {
        s_a <- s_mw
      } else {
        all_cols <- union(names(s_a), names(s_mw))
        add_missing_cols <- function(tbl, cols) {
          miss <- setdiff(cols, names(tbl))
          if (length(miss) > 0L) {
            for (nm in miss) tbl[[nm]] <- NA_real_
          }
          tbl[, cols, drop = FALSE]
        }
        s_a <- add_missing_cols(s_a, all_cols)
        s_mw <- add_missing_cols(s_mw, all_cols)
        s_a <- rbind(s_a, s_mw)
      }
    }
  }

  # Distributional parameter summaries (family-aware, marker-labeled)
  #
  # - Include only parameters required by each marker family.
  # - Replace numeric marker indices with marker names in term labels.
  s_d <- .extract_fit_summary(object, what = "distributional", draws = draws, seed = seed, digits = digits)

  s_dr <- .extract_fit_summary(object, what = "distributional_regression", draws = draws, seed = seed, digits = digits)
  if (!is.null(s_dr) && nrow(s_dr) > 0L) {
    s_dr$parameter <- sub(":.*$", "", s_dr$term)
    s_dr <- s_dr[, c("parameter", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
  }

  # Optional variance/covariance summaries (costly)
  corr_tables <- NULL
  id_marker_cov_tables <- NULL
  if (isTRUE(include_corr)) {
    any_re_indep <- any(as.integer(c(
      sd$indep_id_re %||% 0L,
      sd$indep_marker_re %||% 0L,
      sd$indep_idmarker_cov %||% 0L
    )) == 1L)
    .filter_diag_rows <- function(tbl) {
      if (is.null(tbl) || !all(c("row", "col") %in% names(tbl))) {
        return(tbl)
      }
      tbl[tbl$row == tbl$col, , drop = FALSE]
    }

    q_idm <- as.integer(sd$Q_idm %||% 0L)

    if (q_idm > 0) {
      m_cov <- if (as.integer(sd$indep_idmarker_cov %||% 0L) == 1L) q_idm else (q_idm * (q_idm + 1L)) %/% 2L

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

      reg_rows <- list()

      .summarize_block_parameters <- function(var_names,
                                             block_labels,
                                             term_labels,
                                             row_labels = rep(NA_integer_, length(var_names)),
                                             col_labels = rep(NA_integer_, length(var_names))) {
        keep <- var_names %in% all_vars
        var_names <- var_names[keep]
        block_labels <- block_labels[keep]
        term_labels <- term_labels[keep]
        row_labels <- row_labels[keep]
        col_labels <- col_labels[keep]
        if (length(var_names) == 0) return(NULL)
        out <- as.data.frame(.summarise_draws_diag(fit, var_names, draws = draws, seed = seed))
        idx <- match(out$variable, var_names)
        out$block <- block_labels[idx]
        out$row <- row_labels[idx]
        out$col <- col_labels[idx]
        out$term <- term_labels[idx]
        out <- out[, c("block", "row", "col", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
        out$Estimate <- round(out$Estimate, digits)
        out$Est.Error <- round(out$Est.Error, digits)
        out$Q2.5 <- round(out$Q2.5, digits)
        out$Q97.5 <- round(out$Q97.5, digits)
        out$Rhat <- round(out$Rhat, 3)
        out
      }

      alpha_vars <- paste0("alpha_L[", seq_len(m_cov), "]")
      alpha_blocks <- ifelse(rc_map[, 1] == rc_map[, 2], "SD[id:marker]", "K[id:marker]")
      alpha_tbl <- .summarize_block_parameters(
        alpha_vars,
        block_labels = alpha_blocks,
        term_labels = rep("(Intercept)", m_cov),
        row_labels = rc_map[, 1],
        col_labels = rc_map[, 2]
      )
      if (!is.null(alpha_tbl)) reg_rows[[length(reg_rows) + 1L]] <- alpha_tbl

      k_cov <- as.integer(sd$K_cov %||% 0L)
      if (m_cov > 0 && k_cov > 0) {
        cov_labels <- colnames(sd$Xcov)
        if (is.null(cov_labels) || length(cov_labels) != k_cov) {
          cov_labels <- paste0("k", seq_len(k_cov))
        }
        beta_vars <- as.vector(outer(seq_len(m_cov), seq_len(k_cov), function(m, k) paste0("beta_L[", m, ",", k, "]")))
        beta_blocks <- rep(alpha_blocks, each = k_cov)
        beta_rows <- as.vector(outer(seq_len(m_cov), seq_len(k_cov), function(m, k) rc_map[m, 1]))
        beta_cols <- as.vector(outer(seq_len(m_cov), seq_len(k_cov), function(m, k) rc_map[m, 2]))
        beta_terms <- as.vector(outer(seq_len(m_cov), seq_len(k_cov), function(m, k) {
          cov_labels[k]
        }))
        beta_tbl <- .summarize_block_parameters(
          beta_vars,
          block_labels = beta_blocks,
          term_labels = beta_terms,
          row_labels = beta_rows,
          col_labels = beta_cols
        )
        if (!is.null(beta_tbl)) reg_rows[[length(reg_rows) + 1L]] <- beta_tbl
      }

      lambda_vars <- paste0("lambda_L[", seq_len(m_cov), "]")
      lambda_blocks <- alpha_blocks
      lambda_tbl <- .summarize_block_parameters(
        lambda_vars,
        block_labels = lambda_blocks,
        term_labels = rep("lambda", m_cov),
        row_labels = rc_map[, 1],
        col_labels = rc_map[, 2]
      )
      if (!is.null(lambda_tbl)) reg_rows[[length(reg_rows) + 1L]] <- lambda_tbl

      reg_rows <- Filter(Negate(is.null), reg_rows)
      regression_tbl <- if (length(reg_rows) > 0) do.call(rbind, reg_rows) else NULL
      if (!is.null(regression_tbl)) {
        regression_tbl <- regression_tbl[, c("block", "row", "col", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      }
      id_marker_cov_tables <- list(
        regression = regression_tbl,
        hyperparameters = NULL
      )
      id_marker_cov_tables <- id_marker_cov_tables[!vapply(id_marker_cov_tables, is.null, logical(1))]
      if (length(id_marker_cov_tables) == 0) id_marker_cov_tables <- NULL
    }

    corr_tables <- list(
      id = vcov(object, what = "id", draws = draws),
      marker = if (sd$R_mk > 0) vcov(object, what = "marker", draws = draws) else NULL
    )
    if (isTRUE(any_re_indep)) {
      corr_tables <- lapply(corr_tables, .filter_diag_rows)
    }
    corr_tables$id <- .label_covariance_summary_table(
      corr_tables$id,
      term_labels = as.character(sd$zid_cols %||% paste0("id_re_", seq_len(as.integer(sd$R_id %||% 0L))))
    )
    corr_tables$marker <- .label_covariance_summary_table(
      corr_tables$marker,
      term_labels = as.character(sd$zmk_cols %||% paste0("marker_re_", seq_len(as.integer(sd$R_mk %||% 0L))))
    )
  }

  transform_param_specs <- list(
    cv_total = list(prefix = "coeff_cv_eff", n = sd$n_coeff_cv %||% 0L),
    cs_total = list(prefix = "coeff_cs_eff", n = sd$n_coeff_cs %||% 0L),
    corr = list(prefix = "coeff_corr_eff", n = sd$n_coeff_corr %||% 0L),
    vcov = list(prefix = "coeff_vcov_eff", n = sd$n_coeff_vcov %||% 0L),
    cv_mean = list(prefix = "coeff_cv_mean_eff", n = sd$n_coeff_cv_mean %||% 0L),
    cv_marker = list(prefix = "coeff_cv_marker_eff", n = sd$n_coeff_cv_marker %||% 0L),
    cs_mean = list(prefix = "coeff_cs_mean_eff", n = sd$n_coeff_cs_mean %||% 0L),
    cs_marker = list(prefix = "coeff_cs_marker_eff", n = sd$n_coeff_cs_marker %||% 0L)
  )
  transform_param_tables <- list()
  for (channel in names(transform_param_specs)) {
    spec <- transform_param_specs[[channel]]
    if (channel %in% c("corr", "vcov")) {
      n_components <- .assoc_transform_component_count(
        channel,
        sd$Q_idm,
        diagonal_only = identical(channel, "vcov") && isTRUE(as.integer(sd$indep_idmarker_cov %||% 0L) == 1L)
      )
      if (n_components < 1L || as.integer(spec$n) < 1L) {
        next
      }
      component_labels <- .assoc_transform_component_labels(channel, n_components)
      var_names <- as.vector(outer(seq_len(n_components), seq_len(as.integer(spec$n)), function(m, j) paste0(spec$prefix, "[", m, ",", j, "]")))
      var_names <- var_names[var_names %in% all_vars]
      if (length(var_names) == 0L && as.integer(spec$n) > 0L) {
        var_names <- paste0(spec$prefix, "[", seq_len(as.integer(spec$n)), "]")
        var_names <- var_names[var_names %in% all_vars]
      }
    } else {
      component_labels <- channel
      var_names <- paste0(spec$prefix, "[", seq_len(as.integer(spec$n)), "]")
      var_names <- var_names[var_names %in% all_vars]
    }
    if (length(var_names) == 0) next
    tmp <- as.data.frame(.summarise_draws_diag(fit, var_names, draws = draws, seed = seed))
    if (channel %in% c("corr", "vcov") && any(grepl(paste0("^", spec$prefix, "\\[\\d+,\\d+\\]$"), tmp$variable))) {
      comp_idx <- as.integer(sub(paste0("^", spec$prefix, "\\[(\\d+),\\d+\\]$"), "\\1", tmp$variable))
      basis_idx <- as.integer(sub(paste0("^", spec$prefix, "\\[\\d+,(\\d+)\\]$"), "\\1", tmp$variable))
      tmp$channel <- component_labels[comp_idx]
      tmp$term <- paste0("coeff_", basis_idx)
    } else {
      idx <- match(tmp$variable, var_names)
      tmp$channel <- channel
      tmp$term <- paste0("coeff_", idx)
    }
    tmp <- tmp[, c("channel", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    tmp <- .omit_fixed_transform_endpoint_rows(tmp, channel = channel, spec = spec, sd = sd)
    if (nrow(tmp) == 0L) {
      next
    }
    tmp$Estimate <- round(tmp$Estimate, digits)
    tmp$Est.Error <- round(tmp$Est.Error, digits)
    tmp$Q2.5 <- round(tmp$Q2.5, digits)
    tmp$Q97.5 <- round(tmp$Q97.5, digits)
    tmp$Rhat <- round(tmp$Rhat, 3)
    ord_tmp <- order(tmp$channel, tmp$term)
    transform_param_tables[[channel]] <- tmp[ord_tmp, , drop = FALSE]
  }
  iota_param_specs <- list(
    cv_total = list(intercept = "iota_intercept_cv_eff", slope = "iota_slope_cv_eff", n_intercept = sd$estimate_iota_intercept_cv %||% 0L, n_slope = sd$estimate_iota_slope_cv %||% 0L),
    cs_total = list(intercept = "iota_intercept_cs_eff", slope = "iota_slope_cs_eff", n_intercept = sd$estimate_iota_intercept_cs %||% 0L, n_slope = sd$estimate_iota_slope_cs %||% 0L),
    corr = list(intercept = "iota_intercept_corr_eff", slope = "iota_slope_corr_eff", n_intercept = sd$estimate_iota_intercept_corr %||% 0L, n_slope = sd$estimate_iota_slope_corr %||% 0L),
    vcov = list(intercept = "iota_intercept_vcov_eff", slope = "iota_slope_vcov_eff", n_intercept = sd$estimate_iota_intercept_vcov %||% 0L, n_slope = sd$estimate_iota_slope_vcov %||% 0L),
    cv_mean = list(intercept = "iota_intercept_cv_mean_eff", slope = "iota_slope_cv_mean_eff", n_intercept = sd$estimate_iota_intercept_cv_mean %||% 0L, n_slope = sd$estimate_iota_slope_cv_mean %||% 0L),
    cv_marker = list(intercept = "iota_intercept_cv_marker_eff", slope = "iota_slope_cv_marker_eff", n_intercept = sd$estimate_iota_intercept_cv_marker %||% 0L, n_slope = sd$estimate_iota_slope_cv_marker %||% 0L),
    cs_mean = list(intercept = "iota_intercept_cs_mean_eff", slope = "iota_slope_cs_mean_eff", n_intercept = sd$estimate_iota_intercept_cs_mean %||% 0L, n_slope = sd$estimate_iota_slope_cs_mean %||% 0L),
    cs_marker = list(intercept = "iota_intercept_cs_marker_eff", slope = "iota_slope_cs_marker_eff", n_intercept = sd$estimate_iota_intercept_cs_marker %||% 0L, n_slope = sd$estimate_iota_slope_cs_marker %||% 0L)
  )
  for (channel in names(iota_param_specs)) {
    spec <- iota_param_specs[[channel]]
    if (channel %in% c("corr", "vcov")) {
      n_components <- .assoc_transform_component_count(
        channel,
        sd$Q_idm,
        diagonal_only = identical(channel, "vcov") && isTRUE(as.integer(sd$indep_idmarker_cov %||% 0L) == 1L)
      )
      if (n_components < 1L) {
        next
      }
      component_labels <- .assoc_transform_component_labels(channel, n_components)
      n_intercept <- as.integer(spec$n_intercept %||% 0L)
      n_slope <- as.integer(spec$n_slope %||% 0L)
      var_map <- list(
        intercept = if (n_intercept > 0L) paste0(spec$intercept, "[", seq_len(n_components * n_intercept), "]") else character(0),
        slope = if (n_slope > 0L) paste0(spec$slope, "[", seq_len(n_components * n_slope), "]") else character(0)
      )
      var_names <- unlist(var_map, use.names = FALSE)
      var_names <- var_names[var_names %in% all_vars]
      if (!length(var_names)) {
        if (n_intercept <= 1L) {
          var_names <- c(var_names, paste0(spec$intercept, "[", seq_len(n_components), "]"))
        }
        if (n_slope <= 1L) {
          var_names <- c(var_names, paste0(spec$slope, "[", seq_len(n_components), "]"))
        }
        var_names <- var_names[var_names %in% all_vars]
      }
      if (!length(var_names)) {
        next
      }
      tmp <- as.data.frame(.summarise_draws_diag(fit, var_names, draws = draws, seed = seed))
      raw_idx <- as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", tmp$variable))
      is_intercept <- grepl(paste0("^", spec$intercept, "\\["), tmp$variable)
      count_use <- ifelse(is_intercept, max(1L, n_intercept), max(1L, n_slope))
      tmp$channel <- component_labels[((raw_idx - 1L) %/% count_use) + 1L]
      tmp$term <- paste0(ifelse(is_intercept, "iota_1", "iota_2"), "[", ((raw_idx - 1L) %% count_use) + 1L, "]")
    } else {
      var_map <- c(
        if (as.integer(spec$n_intercept %||% 0L) > 0L) paste0(spec$intercept, "[", seq_len(as.integer(spec$n_intercept)), "]") else character(0),
        if (as.integer(spec$n_slope %||% 0L) > 0L) paste0(spec$slope, "[", seq_len(as.integer(spec$n_slope)), "]") else character(0)
      )
      if (!length(var_map)) {
        var_map <- c(spec$intercept, spec$slope)
      }
      var_names <- unname(var_map)[unname(var_map) %in% all_vars]
      if (!length(var_names)) {
        next
      }
      tmp <- as.data.frame(.summarise_draws_diag(fit, var_names, draws = draws, seed = seed))
      tmp$channel <- channel
      tmp$term <- ifelse(
        grepl(paste0("^", spec$intercept, "(\\[|$)"), tmp$variable),
        paste0("iota_1[", ifelse(grepl("\\[", tmp$variable), sub(paste0("^", spec$intercept, "\\[(\\d+)\\]$"), "\\1", tmp$variable), "1"), "]"),
        paste0("iota_2[", ifelse(grepl("\\[", tmp$variable), sub(paste0("^", spec$slope, "\\[(\\d+)\\]$"), "\\1", tmp$variable), "1"), "]")
      )
    }
    tmp <- tmp[, c("channel", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    tmp$Estimate <- round(tmp$Estimate, digits)
    tmp$Est.Error <- round(tmp$Est.Error, digits)
    tmp$Q2.5 <- round(tmp$Q2.5, digits)
    tmp$Q97.5 <- round(tmp$Q97.5, digits)
    tmp$Rhat <- round(tmp$Rhat, 3)
    ord_tmp <- order(tmp$channel, tmp$term)
    transform_param_tables[[paste0(channel, "_iota")]] <- tmp[ord_tmp, , drop = FALSE]
  }
  transform_params <- if (length(transform_param_tables) > 0) {
    do.call(rbind, unname(transform_param_tables))
  } else {
    NULL
  }

  transform_specs <- cfg$transforms_spec %||% object$call$transforms
  transform_formulas <- .transform_formulas_from_specs(transform_specs, sd = sd)

  if (!is.null(id_marker_cov_tables) && !is.null(id_marker_cov_tables$regression)) {
    id_marker_cov_tables$regression <- .label_covariance_summary_table(
      id_marker_cov_tables$regression,
      term_labels = as.character(sd$zidm_cols %||% paste0("id_marker_re_", seq_len(as.integer(sd$Q_idm %||% 0L))))
    )
  }

  term_diag <- .term_diagnostics_from_tables(list(s_beta, s_basehaz, s_surv, s_a, transform_params, s_d, s_dr, corr_tables, id_marker_cov_tables))
  diag_table <- .build_common_diagnostics_table(
    draws = as.numeric(diag$draws %||% NA_real_),
    divergences = as.numeric(diag$divergences %||% NA_real_),
    treedepth_hits = as.numeric(diag$treedepth_hits %||% NA_real_),
    ebfmi_min = as.numeric(diag$ebfmi_min %||% NA_real_),
    max_rhat = as.numeric(term_diag$max_rhat %||% diag$max_rhat %||% NA_real_),
    min_ess_bulk = as.numeric(term_diag$min_ess_bulk %||% diag$min_ess_bulk %||% diag$min_ess %||% NA_real_),
    min_ess_tail = as.numeric(term_diag$min_ess_tail %||% diag$min_ess_tail %||% diag$min_ess %||% NA_real_),
    n_terms_total = as.numeric(term_diag$n_terms_total %||% 0),
    n_terms_bad_rhat = as.numeric(term_diag$n_terms_bad_rhat %||% 0),
    n_terms_low_ess_bulk = as.numeric(term_diag$n_terms_low_ess_bulk %||% 0),
    n_terms_low_ess_tail = as.numeric(term_diag$n_terms_low_ess_tail %||% 0)
  )

  summary_obj <- SummaryJoiNMeFit$new(
    tables = list(
      diagnostics = diag_table,
      fixef = s_beta,
      baseline_hazard = s_basehaz,
      survival_process = s_surv,
      assoc = s_a,
      transform_parameters = transform_params,
      distributional = s_d,
      distributional_regression = s_dr,
      corr = corr_tables,
      id_marker_cov = id_marker_cov_tables
    ),
    diagnostics = diag,
    metadata = list(
      call = if (!is.null(object$call)) paste(deparse(object$call, width.cutoff = 500L), collapse = " ") else NULL,
      family = sd$family_names %||% .family_code_to_name(cfg$family_long),
      basehaz = sd$basehaz %||% NULL,
      tmax = sd$tmax %||% 1.0,
      draws = draws,
      transforms = cfg$transforms,
      transform_formulas = transform_formulas
    )
  )

  object$cache_set(cache_key, summary_obj)
  summary_obj
}

#' @importFrom brms posterior_summary
#' @rdname summary.JoiNMeFit
#' @export
posterior_summary.JoiNMeFit <- function(object, ...) {
  summary(object, ...)
}

#' @rdname summary.JoiNMeDynPred
#' @importFrom brms posterior_summary
#' @export
posterior_summary.JoiNMeDynPred <- function(object, ...) {
  summary(object, ...)
}

#' Summary of JoiNMe dynamic prediction
#'
#' @description
#' Summarises a `JoiNMeDynPred` object
#'
#' @details
#' The summary includes:
#' - an overview table of row/subject counts by process,
#' - median survival time per subject (draw-wise crossing of `S(t)=0.5`,
#'   summarised with estimate, uncertainty, and diagnostics),
#' - predicted id-level random effects per subject (estimate, uncertainty,
#'   interval, and diagnostics),
#' - predicted marker-by-id random effects and covariance per subject when
#'   marker covariance depends on id,
#' - a compact diagnostics table counting potential convergence/ESS issues.
#'
#' @param object A `JoiNMeDynPred` object produced by [predict].
#' @param ... Unused.
#'
#' @return A `summary_JoiNMeDynPred` object containing tabular summaries.
#' @method summary JoiNMeDynPred
#' @export
summary.JoiNMeDynPred <- function(object, ...) {
    
    cached <- object$cache_get("summary")
    if (!is.null(cached)) return(cached)

    n_long_rows <- if (!is.null(object$predictions$longitudinal)) nrow(object$predictions$longitudinal) else 0
    n_surv_rows <- if (!is.null(object$predictions$survival)) nrow(object$predictions$survival) else 0
    n_cumhaz_rows <- if (!is.null(object$predictions$cumhaz)) nrow(object$predictions$cumhaz) else 0
    n_long_subjects <- if (!is.null(object$predictions$longitudinal)) length(unique(object$predictions$longitudinal$id)) else 0
    n_surv_subjects <- if (!is.null(object$predictions$survival)) length(unique(object$predictions$survival$id)) else 0
    n_cumhaz_subjects <- if (!is.null(object$predictions$cumhaz)) length(unique(object$predictions$cumhaz$id)) else 0

    overview_table <- data.frame(
        metric = c(
            "longitudinal_rows",
            "longitudinal_subjects",
            "survival_rows",
            "survival_subjects",
            "cumhaz_rows",
            "cumhaz_subjects",
            "prediction_scale",
            "posterior_draws"
        ),
        value = c(
            n_long_rows,
            n_long_subjects,
            n_surv_rows,
            n_surv_subjects,
            n_cumhaz_rows,
            n_cumhaz_subjects,
            if (!is.null(object$metadata$scales)) paste(object$metadata$scales, collapse = ",") else object$metadata$scale %||% NA_character_,
            object$n_samples %||% object$metadata$n_samples %||% NA_real_
        ),
        stringsAsFactors = FALSE
    )

    median_survival_table <- NULL
    if (!is.null(object$draws$survival) && length(object$draws$survival) > 0) {
        median_survival_rows <- lapply(names(object$draws$survival), function(subject_id) {
            subject_survival <- object$draws$survival[[subject_id]]
            draw_median_times <- .median_survival_time_by_draw(
                survival_draw_matrix = subject_survival$matrix,
                time_grid = subject_survival$time,
                threshold = 0.5
            )
            draw_summary <- .summarize_draw_vector_with_diagnostics(draw_median_times)
            cbind(
                id = subject_id,
                term = "median_survival_time",
                n_reached = sum(is.finite(draw_median_times)),
                n_total = length(draw_median_times),
                draw_summary,
                stringsAsFactors = FALSE
            )
        })
        median_survival_table <- do.call(rbind, median_survival_rows)
    }

    random_effects_id_table <- NULL
    if (!is.null(object$draws$random_effects_id) && length(object$draws$random_effects_id) > 0) {
        random_effect_rows <- lapply(names(object$draws$random_effects_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_id[[subject_id]]
            u_id_matrix <- subject_effects$matrix
            u_id_terms <- subject_effects$terms %||% paste0("u_id[", seq_len(ncol(u_id_matrix)), "]")

            per_term <- lapply(seq_len(ncol(u_id_matrix)), function(term_index) {
                draw_summary <- .summarize_draw_vector_with_diagnostics(u_id_matrix[, term_index])
                cbind(
                    id = subject_id,
                    term = u_id_terms[term_index],
                    draw_summary,
                    stringsAsFactors = FALSE
                )
            })
            do.call(rbind, per_term)
        })
        random_effects_id_table <- do.call(rbind, random_effect_rows)
    }

    marker_corr_depends_on_id <- isTRUE(object$metadata$marker_corr_depends_on_id)
    any_re_indep <- any(as.integer(c(
        object$metadata$indep_id_re %||% 0L,
        object$metadata$indep_marker_re %||% 0L,
        object$metadata$indep_idmarker_cov %||% 0L
    )) == 1L)

    random_effects_marker_id_table <- NULL
    if (marker_corr_depends_on_id && !is.null(object$draws$random_effects_marker_id) && length(object$draws$random_effects_marker_id) > 0) {
        marker_id_rows <- lapply(names(object$draws$random_effects_marker_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_marker_id[[subject_id]]
            w_idm_matrix <- subject_effects$matrix
            if (is.null(w_idm_matrix) || ncol(w_idm_matrix) == 0) return(NULL)

            marker_terms <- colnames(w_idm_matrix)
            if (is.null(marker_terms)) {
                marker_terms <- paste0("marker_id[", seq_len(ncol(w_idm_matrix)), "]")
            }

            per_term <- lapply(seq_len(ncol(w_idm_matrix)), function(term_index) {
                draw_summary <- .summarize_draw_vector_with_diagnostics(w_idm_matrix[, term_index])
                term_label <- as.character(marker_terms[term_index])
                marker_label <- if (grepl("::", term_label, fixed = TRUE)) {
                    sub("::.*$", "", term_label)
                } else {
                    NA_character_
                }
                basis_term <- if (grepl("::", term_label, fixed = TRUE)) {
                    sub("^.*::", "", term_label)
                } else {
                    term_label
                }
                cbind(
                    id = subject_id,
                    marker = marker_label,
                    term = basis_term,
                    draw_summary,
                    stringsAsFactors = FALSE
                )
            })
            do.call(rbind, per_term)
        })
        marker_id_rows <- Filter(Negate(is.null), marker_id_rows)
        if (length(marker_id_rows) > 0) {
            random_effects_marker_id_table <- do.call(rbind, marker_id_rows)
        }
    }

    corr_marker_id_table <- NULL
    if (marker_corr_depends_on_id && !is.null(object$draws$random_effects_marker_id) && length(object$draws$random_effects_marker_id) > 0) {
        corr_rows <- lapply(names(object$draws$random_effects_marker_id), function(subject_id) {
            subject_effects <- object$draws$random_effects_marker_id[[subject_id]]
            corr_draws <- subject_effects$corr
            if (is.null(corr_draws) || length(dim(corr_draws)) != 3) return(NULL)

            terms <- as.character(subject_effects$terms %||% paste0("w_idm[", seq_len(dim(corr_draws)[2]), "]"))
            q_dim <- dim(corr_draws)[2]

            per_cell <- lapply(seq_len(q_dim), function(r_idx) {
                lapply(seq_len(q_dim), function(c_idx) {
                    draw_summary <- .summarize_draw_vector_with_diagnostics(corr_draws[, r_idx, c_idx])
                    cbind(
                        id = subject_id,
                        row = terms[r_idx],
                        col = terms[c_idx],
                        draw_summary,
                        stringsAsFactors = FALSE
                    )
                })
            })
            do.call(rbind, unlist(per_cell, recursive = FALSE))
        })
        corr_rows <- Filter(Negate(is.null), corr_rows)
        if (length(corr_rows) > 0) {
            corr_marker_id_table <- do.call(rbind, corr_rows)
            if (isTRUE(any_re_indep)) {
                corr_marker_id_table <- corr_marker_id_table[corr_marker_id_table$row == corr_marker_id_table$col, , drop = FALSE]
            }
        }
    }

    diagnostics_table <- NULL
    bind_rows_safe <- function(...) {
        row_list <- Filter(Negate(is.null), list(...))
        if (length(row_list) == 0) return(NULL)
        all_columns <- unique(unlist(lapply(row_list, names), use.names = FALSE))
        row_list <- lapply(row_list, function(df) {
            missing_columns <- setdiff(all_columns, names(df))
            if (length(missing_columns) > 0) {
                for (column_name in missing_columns) df[[column_name]] <- NA
            }
            df[, all_columns, drop = FALSE]
        })
        do.call(rbind, row_list)
    }
    diagnostics_source <- bind_rows_safe(
        if (!is.null(random_effects_id_table)) transform(random_effects_id_table, section = "random_effects_id") else NULL,
        if (!is.null(random_effects_marker_id_table)) transform(random_effects_marker_id_table, section = "random_effects_marker_id") else NULL,
        if (!is.null(median_survival_table)) transform(median_survival_table, section = "median_survival_time") else NULL
    )
    term_diag <- .term_diagnostics_from_tables(diagnostics_source)
    sampler_diag <- object$metadata$sampler_diagnostics %||% list()

    diagnostics_table <- .build_common_diagnostics_table(
        draws = as.numeric(object$n_samples %||% object$metadata$n_samples %||% NA_real_),
        divergences = as.numeric(sampler_diag$divergences %||% 0),
        treedepth_hits = as.numeric(sampler_diag$treedepth_hits %||% 0),
        ebfmi_min = as.numeric(sampler_diag$ebfmi_min %||% NA_real_),
        max_rhat = as.numeric(term_diag$max_rhat %||% sampler_diag$max_rhat %||% NA_real_),
        min_ess_bulk = as.numeric(term_diag$min_ess_bulk %||% sampler_diag$min_ess_bulk %||% NA_real_),
        min_ess_tail = as.numeric(term_diag$min_ess_tail %||% sampler_diag$min_ess_tail %||% NA_real_),
        n_terms_total = as.numeric(term_diag$n_terms_total %||% 0),
        n_terms_bad_rhat = as.numeric(term_diag$n_terms_bad_rhat %||% 0),
        n_terms_low_ess_bulk = as.numeric(term_diag$n_terms_low_ess_bulk %||% 0),
        n_terms_low_ess_tail = as.numeric(term_diag$n_terms_low_ess_tail %||% 0)
    )

    tables <- list(
        diagnostics = diagnostics_table,
        overview = overview_table,
        median_survival_time = median_survival_table,
        random_effects_id = random_effects_id_table,
        random_effects_marker_id = random_effects_marker_id_table,
        corr_marker_id = corr_marker_id_table
    )

    summary_obj <- SummaryJoiNMeDynPred$new(tables = tables, metadata = object$metadata)
    object$cache_set("summary", summary_obj)
    summary_obj
}


#' @export
print.summary_JoiNMeDynPred <- function(x, ...) {
    .cli_summary_heading("Prediction summary", level = 1L)
    meta_lines <- character(0)
    if (!is.null(x$metadata$pred_type)) {
        meta_lines <- c(meta_lines, paste0("Prediction type: ", x$metadata$pred_type))
    }
    if (!is.null(x$metadata$n_subjects)) {
        meta_lines <- c(meta_lines, paste0("Subjects: ", x$metadata$n_subjects))
    }
    .cli_print_bullets(meta_lines)
    .cli_print_table_section("Diagnostics", x$tables$diagnostics, level = 2L, formatter = .format_common_diagnostics_for_print)
    .cli_print_table_section("Overview", x$tables$overview, level = 2L)
    .cli_print_table_section("Median survival time", x$tables$median_survival_time, level = 2L)
    .cli_print_table_section("Predicted id random effects", x$tables$random_effects_id, level = 2L)
    .cli_print_table_section("Predicted marker-by-id random effects", x$tables$random_effects_marker_id, level = 2L)
    .cli_print_table_section("Predicted marker-by-id covariance", x$tables$corr_marker_id, level = 2L)
    invisible(x)
}
