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
#' @return A `summary_JoiNMeFit` object. Its `tables` element includes
#'   posterior association coefficients, a compact marker-weight table with
#'   common locations and within-set spreads, transform parameters, and, when an
#'   ordered piecewise-linear association is active, `piecewise_ordinates`
#'   containing the relative log-hazard and hazard-ratio contribution at every
#'   knot. Declared marker-weight offsets are retained in
#'   `metadata$marker_weight_offsets` and printed above the posterior tables.
#' @export
summary.JoiNMeFit <- function(object, draws = NULL, seed = .Random.seed[[1]], digits = 3,
                           include_corr = TRUE, ...) {

  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  if (is.null(draws)) draws <- cfg$draws_default
  assertthat::assert_that(is.null(draws) || (is.numeric(draws) && draws > 0),
                          msg = "draws must be NULL or a positive number.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- paste0(
    "summary_draws=", draws,
    "_seed=", seed,
    "_digits=", digits,
    "_corr=", as.integer(include_corr)
  )
  # Cache by every argument that changes the selected draws or reported table.
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  diag <- .joinme_sampler_diagnostics(fit)
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  s_beta <- .extract_fit_summary(object, what = "fixef", draws = draws, seed = seed, digits = digits)

  g_ext <- tryCatch(
    extract.JoiNMeFit(object, what = "gamma_w", draws = draws, seed = seed, keep_chains = TRUE),
    error = function(e) NULL
  )
  s_g <- .summarise_named_draws(g_ext$posterior_draws %||% NULL, digits = digits)

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
      bh_terms <- .basehaz_term_labels(
        basehaz = sd$basehaz %||% cfg$basehaz %||% "bs",
        n_terms = k_bs,
        supplied_names = sd$basehaz_cols
      ) # also repairs bare integer spline labels retained by previously fitted objects
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
  s_a <- .summarise_named_draws(a_ext$posterior_draws %||% NULL, digits = digits)

  # Marker weights have their own compact model summary. Individual effective
  # weights are deliberately excluded here: marker_weights() and coef() return
  # the effective values, fixef() returns fitted set means, and ranef() returns
  # marker-specific departures after removing the means and declared offsets.
  # For each fitted shared or term-specific set, report only its common location
  # and the posterior root-mean-square realised departure around that location.
  s_marker_weights <- .marker_weight_set_summary(
    object,
    draws = draws,
    seed = seed,
    digits = digits,
    all_vars = all_vars
  )

  # Distributional parameter summaries (family-aware, marker-labeled)
  #
  # - Include only parameters required by each marker family.
  # - Replace numeric marker indices with marker names in term labels.
  s_d <- .extract_fit_summary(object, what = "distributional", draws = draws, seed = seed, digits = digits)
  s_tau_fixed <- .fixed_tau_summary_table(object, digits = digits)
  if (!is.null(s_tau_fixed)) {
    s_d <- if (is.null(s_d)) {
      s_tau_fixed
    } else {
      rbind(s_d, s_tau_fixed)
    }
  }

  s_dr <- .extract_fit_summary(object, what = "distributional_regression", draws = draws, seed = seed, digits = digits)
  if (!is.null(s_dr) && nrow(s_dr) > 0L) {
    s_dr$parameter <- sub(":.*$", "", s_dr$term)
    s_dr <- s_dr[, c("parameter", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
  }

  # Optional variance/covariance summaries (costly)
  corr_tables <- NULL
  id_marker_cov_tables <- NULL
  if (isTRUE(include_corr)) {
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

      # The two formulaVCov components may have unrelated model matrices, so
      # their posterior slope names and scientific labels are assembled
      # independently. Earlier fitted objects stored one shared beta_L array;
      # that representation remains readable without changing new fits.
      covariance_design <- .stored_vcov_design(sd) # current split design or a read-only earlier shared-design view
      k_cov_sd <- covariance_design$k_sd
      k_cov_corr <- covariance_design$k_corr
      if (isTRUE(covariance_design$shared_format) && m_cov > 0L && k_cov_sd > 0L) {
        shared_labels <- colnames(covariance_design$x_sd) %||% paste0("k", seq_len(k_cov_sd))
        beta_shared_vars <- as.vector(outer(
          seq_len(m_cov), seq_len(k_cov_sd),
          function(m, k) paste0("beta_L[", m, ",", k, "]")
        ))
        beta_shared_tbl <- .summarize_block_parameters(
          beta_shared_vars,
          block_labels = rep(alpha_blocks, times = k_cov_sd),
          term_labels = rep(shared_labels, each = m_cov),
          row_labels = rep(rc_map[, 1], times = k_cov_sd),
          col_labels = rep(rc_map[, 2], times = k_cov_sd)
        )
        if (!is.null(beta_shared_tbl)) reg_rows[[length(reg_rows) + 1L]] <- beta_shared_tbl
      } else if (q_idm > 0L && k_cov_sd > 0L) {
        sd_labels <- colnames(covariance_design$x_sd) %||% paste0("k", seq_len(k_cov_sd))
        beta_sd_vars <- as.vector(outer(
          seq_len(q_idm), seq_len(k_cov_sd),
          function(r, k) paste0("beta_L_sd[", r, ",", k, "]")
        ))
        beta_sd_tbl <- .summarize_block_parameters(
          beta_sd_vars,
          block_labels = rep("SD[id:marker]", q_idm * k_cov_sd),
          term_labels = rep(sd_labels, each = q_idm),
          row_labels = rep(seq_len(q_idm), times = k_cov_sd),
          col_labels = rep(seq_len(q_idm), times = k_cov_sd)
        )
        if (!is.null(beta_sd_tbl)) reg_rows[[length(reg_rows) + 1L]] <- beta_sd_tbl
      }

      correlation_map <- rc_map[rc_map[, 1] != rc_map[, 2], , drop = FALSE]
      if (!isTRUE(covariance_design$shared_format) && nrow(correlation_map) > 0L && k_cov_corr > 0L) {
        corr_labels <- colnames(covariance_design$x_corr) %||% paste0("k", seq_len(k_cov_corr))
        beta_corr_vars <- as.vector(outer(
          seq_len(nrow(correlation_map)), seq_len(k_cov_corr),
          function(m, k) paste0("beta_L_corr[", m, ",", k, "]")
        ))
        beta_corr_tbl <- .summarize_block_parameters(
          beta_corr_vars,
          block_labels = rep("K[id:marker]", nrow(correlation_map) * k_cov_corr),
          term_labels = rep(corr_labels, each = nrow(correlation_map)),
          row_labels = rep(correlation_map[, 1], times = k_cov_corr),
          col_labels = rep(correlation_map[, 2], times = k_cov_corr)
        )
        if (!is.null(beta_corr_tbl)) reg_rows[[length(reg_rows) + 1L]] <- beta_corr_tbl
      }

      lambda_vars <- paste0("lambda_L[", seq_len(m_cov), "]")
      lambda_blocks <- alpha_blocks
      lambda_tbl <- .summarize_block_parameters(
        lambda_vars,
        block_labels = lambda_blocks,
        term_labels = rep("latent SD (lambda)", m_cov),
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
    # Retain the complete symmetric covariance summaries even when a block was
    # fitted as independent. In that case the off-diagonal entries are exact
    # zeros, which is meaningful model information rather than redundant
    # output and prevents an independent block in one domain from hiding
    # covariance entries belonging to another domain.
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

#   browser()
  transform_specs <- cfg$transform_spec %||% object$call$transforms
  transform_formulas <- .transform_formulas_from_specs(transform_specs, sd = sd)
  piecewise_ordinates <- .summarise_pwlin_ordinates(
    object = object,
    transform_specs = transform_specs,
    draws = draws,
    seed = seed,
    digits = digits
  )

  if (!is.null(id_marker_cov_tables) && !is.null(id_marker_cov_tables$regression)) {
    id_marker_cov_tables$regression <- .label_covariance_summary_table(
      id_marker_cov_tables$regression,
      term_labels = as.character(sd$zidm_cols %||% paste0("id_marker_re_", seq_len(as.integer(sd$Q_idm %||% 0L))))
    )
  }

  has_survival_process <- .fit_includes_survival(
    object
  ) # whether event observations contributed information to this fit
  if (!has_survival_process) {
    # The common Stan programme still declares event-process parameters for a
    # longitudinal-only analysis. Their draws come solely from their priors
    # and are computational scaffolding, not estimands supported by observed
    # survival data. Excluding them also keeps convergence counts focused on
    # the nested longitudinal mixed model that was actually fitted.
    s_basehaz <- NULL
    s_surv <- NULL
    s_a <- NULL
    s_marker_weights <- NULL
    transform_params <- NULL
    piecewise_ordinates <- NULL
  }

  reported_tables <- list(
    s_beta,
    s_basehaz,
    s_surv,
    s_a,
    s_marker_weights,
    transform_params,
    s_d,
    s_dr,
    corr_tables,
    id_marker_cov_tables
  ) # posterior tables whose parameter-level convergence counts are reported
  diag_table <- .summary_diagnostics_table(
    sampler_diagnostics = diag,
    reported_tables = reported_tables
  ) # sampler-wide extrema combined with counts from the displayed estimands

  summary_obj <- SummaryJoiNMeFit$new(
    tables = list(
      diagnostics = diag_table,
      fixef = s_beta,
      baseline_hazard = s_basehaz,
      survival_process = s_surv,
      assoc = s_a,
      marker_weights = s_marker_weights,
      transform_parameters = transform_params,
      piecewise_ordinates = piecewise_ordinates,
      distributional = s_d,
      distributional_regression = s_dr,
      corr = corr_tables,
      id_marker_cov = id_marker_cov_tables
    ),
    diagnostics = diag,
    metadata = list(
      call = if (!is.null(object$call)) paste(deparse(object$call, width.cutoff = 500L), collapse = " ") else NULL,
      family = sd$family_names %||% .family_code_to_name(cfg$family_long),
      basehaz = if (has_survival_process) {
        sd$basehaz %||% NULL
      } else {
        NULL
      },
      tmax = sd$tmax %||% 1.0,
      draws = draws,
      event_process = if (has_survival_process) {
        "joint longitudinal-survival"
      } else {
        "not fitted"
      },
      marker_weight_offsets = if (has_survival_process) {
        .marker_weight_offset_metadata(sd)
      } else {
        NULL
      },
      transforms = if (has_survival_process) {
        cfg$transforms
      } else {
        NULL
      },
      transform_formulas = if (has_survival_process) {
        transform_formulas
      } else {
        NULL
      }
    )
  )

  object$cache_set(cache_key, summary_obj)
  summary_obj
}

#' Build diagnostics for the parameters represented in a model summary
#'
#' @description
#' The sampler knows the most extreme R-hat and effective sample size over all
#' saved variables, whereas the summary tables identify how many displayed
#' scientific estimands breach the reporting thresholds. Both views matter.
#' This helper retains the sampler-wide extrema and obtains the counts from the
#' tables that the reader can inspect. In particular, a poorly mixed latent
#' class parameter must not disappear from the headline diagnostics merely
#' because the ordinary joint-model summary was assembled first.
#'
#' @param sampler_diagnostics Named diagnostics returned by .joinme_sampler_diagnostics().
#' @param reported_tables A possibly nested list of posterior summary tables.
#'
#' @return A data frame in the common diagnostics schema.
#' @keywords internal
#' @noRd
.summary_diagnostics_table <- function(
  sampler_diagnostics,
  reported_tables
) {
  term_diagnostics <- .term_diagnostics_from_tables(
    reported_tables
  ) # convergence extrema and threshold counts among displayed parameters
  finite_maximum <- function(values) {
    finite_values <- values[
      is.finite(values)
    ] # available finite diagnostics from the sampler and displayed tables
    if (length(finite_values) == 0L) NA_real_ else max(finite_values)
  }
  finite_minimum <- function(values) {
    finite_values <- values[
      is.finite(values)
    ] # available finite diagnostics from the sampler and displayed tables
    if (length(finite_values) == 0L) NA_real_ else min(finite_values)
  }

  sampler_bulk_ess <- sampler_diagnostics$min_ess_bulk %||%
    sampler_diagnostics$min_ess %||%
    NA_real_ # smallest sampler-wide bulk ESS, with support for older fits
  sampler_tail_ess <- sampler_diagnostics$min_ess_tail %||%
    sampler_diagnostics$min_ess %||%
    NA_real_ # smallest sampler-wide tail ESS, with support for older fits

  .build_common_diagnostics_table(
    draws = as.numeric(sampler_diagnostics$draws %||% NA_real_),
    divergences = as.numeric(
      sampler_diagnostics$divergences %||% NA_real_
    ),
    treedepth_hits = as.numeric(
      sampler_diagnostics$treedepth_hits %||% NA_real_
    ),
    ebfmi_min = as.numeric(sampler_diagnostics$ebfmi_min %||% NA_real_),
    max_rhat = finite_maximum(c(
      term_diagnostics$max_rhat,
      sampler_diagnostics$max_rhat
    )),
    min_ess_bulk = finite_minimum(c(
      term_diagnostics$min_ess_bulk,
      sampler_bulk_ess
    )),
    min_ess_tail = finite_minimum(c(
      term_diagnostics$min_ess_tail,
      sampler_tail_ess
    )),
    n_terms_total = as.numeric(term_diagnostics$n_terms_total %||% 0),
    n_terms_bad_rhat = as.numeric(
      term_diagnostics$n_terms_bad_rhat %||% 0
    ),
    n_terms_low_ess_bulk = as.numeric(
      term_diagnostics$n_terms_low_ess_bulk %||% 0
    ),
    n_terms_low_ess_tail = as.numeric(
      term_diagnostics$n_terms_low_ess_tail %||% 0
    )
  )
}

#' Summarise fitted piecewise-linear log-hazard ordinates
#'
#' @description
#' Reconstruct each ordered piecewise-linear transform at its knots and combine
#' it draw by draw with the corresponding association coefficient. The result
#' is the posterior relative log-hazard contribution at every knot, together
#' with its hazard-ratio representation. The first ordinate is the zero
#' reference imposed for identifiability with the baseline hazard.
#'
#' @param object A fitted JoiNMeFit object.
#' @param transform_specs Named association transform specifications.
#' @param draws Optional posterior draw count.
#' @param seed Integer seed used when draws are subsampled.
#' @param digits Number of decimal places used in the returned table.
#'
#' @return A data frame with channel/component, knot number and location,
#'   direction, and posterior summaries on relative log-hazard and hazard-ratio
#'   scales; NULL when no active fitted piecewise-linear transform is present.
#' @keywords internal
#' @noRd
.summarise_pwlin_ordinates <- function(object, transform_specs, draws = NULL, seed = 1, digits = 3) {
  if (is.null(transform_specs) || !length(transform_specs)) {
    return(NULL)
  }
  pwlin_channels <- names(transform_specs)[vapply(transform_specs, function(spec) {
    identical(.canonicalise_transform_type(spec$type %||% "identity"), "pwlin")
  }, logical(1))]
  if (!length(pwlin_channels)) {
    return(NULL)
  }

  plot_data <- tryCatch(.get_association_plot_data(object, seed = seed), error = function(e) NULL)
  if (is.null(plot_data)) {
    return(NULL)
  }
  available_terms <- unique(as.character(plot_data$term_map$term %||% character(0)))
  available_terms <- available_terms[!grepl("^weight:\\s", available_terms)]
  rows <- list()

  for (channel in pwlin_channels) {
    channel_terms <- available_terms[
      available_terms == channel | startsWith(available_terms, paste0(channel, "["))
    ]
    if (!length(channel_terms)) {
      next
    }
    map <- .assoc_channel_map(channel)
    knots <- as.numeric(
      object$stan_data[[map$knots]] %||%
        transform_specs[[channel]]$knots %||%
        transform_specs[[channel]]$cutpoints %||%
        transform_specs[[channel]]$x
    )
    if (length(knots) < 2L) {
      next
    }

    mode_name <- paste0("tf_mode_", .get_transform_term(channel)$mode_suffix)
    direction <- .monotone_direction_label(
      transform_specs[[channel]]$direction %||%
        if (as.integer(object$stan_data[[mode_name]] %||% 3L) == 7L) -1L else 1L
    )
    for (term in channel_terms) {
      alpha <- tryCatch(
        .assoc_coeff_draws(object, term = term, data = plot_data, seed = seed),
        error = function(e) NULL
      )
      if (is.null(alpha) || !length(alpha)) {
        next
      }
      transform_at_knots <- tryCatch(
        .association_transform_matrix(
          object,
          term_key = channel,
          term = term,
          x_grid = knots,
          n_draws = length(alpha),
          seed = seed,
          data = plot_data
        ),
        error = function(e) NULL
      )
      if (is.null(transform_at_knots) || nrow(transform_at_knots) != length(alpha)) {
        next
      }

      draw_index <- seq_along(alpha)
      if (!is.null(draws) && length(draw_index) > draws) {
        set.seed(seed)
        draw_index <- sample(draw_index, draws)
      }
      alpha <- alpha[draw_index]
      transform_at_knots <- transform_at_knots[draw_index, , drop = FALSE]
      log_hazard <- transform_at_knots * alpha
      q <- t(apply(
        log_hazard,
        2L,
        stats::quantile,
        probs = c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      ))
      estimate <- colMeans(log_hazard)
      rows[[length(rows) + 1L]] <- data.frame(
        channel = term,
        ordinate = seq_along(knots),
        knot = knots,
        direction = direction,
        Estimate = round(estimate, digits),
        Est.Error = round(apply(log_hazard, 2L, stats::sd), digits),
        Q2.5 = round(q[, 1L], digits),
        Q97.5 = round(q[, 2L], digits),
        Hazard.Ratio = round(exp(estimate), digits),
        HR.Q2.5 = round(exp(q[, 1L]), digits),
        HR.Q97.5 = round(exp(q[, 2L]), digits),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

#' @importFrom brms posterior_summary
#' @rdname summary.JoiNMeFit
#' @param what Posterior view to return. `"model"` gives the complete fitted
#'   model summary. `"fixef"`, `"ranef"`, and `"coef"` provide the common
#'   reporting layer used by the corresponding high-level methods.
#' @param summary Logical. For a coefficient view, `TRUE` returns posterior
#'   summaries and `FALSE` returns its structured draw-level extraction.
#' @details
#' For coefficient views, this is the reporting layer between [extract()] and
#' the conventional [fixef()], [ranef()], and [coef()] methods. The high-level
#' methods delegate here; requests with `summary = FALSE` continue to the
#' component-aware draw representation owned by `extract()`.
#' @export
posterior_summary.JoiNMeFit <- function(
  object,
  what = c("model", "fixef", "ranef", "coef"),
  draws = NULL,
  seed = 1,
  digits = 3,
  summary = TRUE,
  ...
) {
  what <- match.arg(what) # requested reporting layer within the fitted model
  if (identical(what, "model")) {
    return(summary.JoiNMeFit(
      object,
      draws = draws,
      seed = seed,
      digits = digits,
      ...
    ))
  }

  # Draw-level coefficient requests belong to the extraction layer. Routing
  # them here, before any tabulation, keeps the public hierarchy strictly
  # one-directional: coefficient method -> posterior_summary() -> extract().
  if (!isTRUE(summary)) {
    extraction_view <- switch(
      what,
      fixef = "fixed_effects",
      ranef = "random_effects",
      coef = "coefficients"
    ) # component-aware selector owned by extract.JoiNMeFit()
    return(extract(
      object,
      what = extraction_view,
      draws = draws,
      seed = seed,
      keep_chains = FALSE
    )$posterior_draws)
  }

  component_function <- switch(
    what,
    fixef = .summarise_fixed_effect_posterior,
    ranef = .summarise_random_effect_posterior,
    coef = .summarise_combined_coefficient_posterior
  ) # one implementation owner for every coefficient view
  component_function(
    object,
    draws = draws,
    seed = seed,
    digits = digits,
    summary = summary,
    ...
  )
}

#' @rdname summary.JoiNMeDynPred
#' @importFrom brms posterior_summary
#' @export
posterior_summary.JoiNMeDynPred <- function(object, ...) {
  summary(object, ...)
}

#' Central posterior credible intervals for JoiNMe models
#'
#' @description
#' Computes central posterior credible intervals through
#' [rstantools::posterior_interval()]. The method uses the same scientific
#' coefficient views as [posterior_summary()]: the complete fitted posterior,
#' population-level coefficients, group-specific deviations, or their combined
#' coefficients.
#'
#' @details
#' For `what = "model"`, the returned matrix has one row per friendly posterior
#' parameter name. The `variables` argument selects those names after JoiNMe has
#' translated the Stan coordinates into their reported statistical terms.
#'
#' The `fixef` view is also rectangular and therefore returns an interval
#' matrix. The `ranef` and `coef` views contain several statistically distinct
#' coefficient tables. Their nested list structure is retained, while each
#' draw-level `value` column is replaced by the two interval limits. This keeps
#' subject, marker, event, distributional-family, association-term, and
#' covariance identities explicit.
#'
#' Every interval is calculated by the default matrix method of
#' [rstantools::posterior_interval()]. Consequently, `prob` is the total
#' posterior probability contained between the two central quantiles.
#'
#' @param object A fitted `JoiNMeFit` object. Latent-class fits inherit this
#'   method through `JoiNMeMixFit`.
#' @param prob A single number strictly between zero and one giving the
#'   posterior probability contained in the interval. The default is `0.9`,
#'   following `rstantools`.
#' @param what Posterior view to interval: `"model"`, `"fixef"`, `"ranef"`,
#'   or `"coef"`.
#' @param variables Optional character vector selecting friendly parameter
#'   names when `what = "model"`.
#' @param regex Logical; when `TRUE`, interpret `variables` as regular
#'   expressions. This argument applies only to `what = "model"`.
#' @param draws Optional number of posterior draws to retain before computing
#'   the intervals.
#' @param seed Integer seed used when posterior draws are subsampled.
#' @param ... Additional arguments passed to
#'   [rstantools::posterior_interval()].
#'
#' @return For `what = "model"` or `what = "fixef"`, a numeric matrix with one
#'   row per term and two probability-labelled columns. For `what = "ranef"`
#'   or `what = "coef"`, a nested list of data frames retaining the coefficient
#'   identifiers and containing the same two probability-labelled columns.
#'
#' @importFrom rstantools posterior_interval
#' @export
posterior_interval.JoiNMeFit <- function(
  object,
  prob = 0.9,
  what = c("model", "fixef", "ranef", "coef"),
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  ...
) {
  what <- match.arg(what) # scientific posterior view whose uncertainty is requested

  # The complete posterior is already represented by one friendly-named draws
  # matrix. Passing it directly to rstantools preserves its established row and
  # probability-column convention without duplicating interval calculations.
  if (identical(what, "model")) {
    posterior_matrix <- posterior_draws(
      object,
      variables = variables,
      regex = regex,
      draws = draws,
      seed = seed,
      format = "draws_matrix"
    ) # retained posterior draws with reported parameter names
    return(rstantools::posterior_interval(
      as.matrix(posterior_matrix),
      prob = prob,
      ...
    ))
  }

  if (!is.null(variables)) {
    cli::cli_abort(c(
      x = "{.arg variables} is available only when {.code what = 'model'}.",
      i = "The coefficient views retain several identifier columns; select the desired table from the returned structure."
    ))
  }

  # Reuse the draw-level layer underlying posterior_summary(), fixef(), ranef(),
  # and coef(). This ensures intervals and conventional summaries refer to
  # precisely the same coefficients and scientific scales.
  coefficient_draws <- posterior_summary(
    object,
    what = what,
    draws = draws,
    seed = seed,
    summary = FALSE
  ) # matrix for fixed effects; nested draw tables for random or combined effects

  .joinme_posterior_interval_structure(
    coefficient_draws,
    prob = prob,
    ...
  )
}

#' @keywords internal
#' @noRd
.joinme_posterior_interval_structure <- function(object, prob, ...) {
  # Rectangular components can be delegated immediately to the authoritative
  # rstantools calculation. Coercion removes posterior-specific matrix classes
  # while retaining the term names required for interval row labels.
  if (is.matrix(object)) {
    return(rstantools::posterior_interval(as.matrix(object), prob = prob, ...))
  }

  # Draw-level coefficient tables use one row per posterior draw and coefficient
  # identity. All columns other than the numerical decomposition are therefore
  # grouping variables that must remain in the interval result.
  if (is.data.frame(object) && all(c("draw", "value") %in% names(object))) {
    identity_columns <- setdiff(
      names(object),
      c("draw", "value", "fixed", "random")
    ) # columns uniquely identifying a scientific coefficient

    if (length(identity_columns) == 0L) {
      interval <- rstantools::posterior_interval(
        matrix(object$value, ncol = 1L, dimnames = list(NULL, "value")),
        prob = prob,
        ...
      )
      return(as.data.frame(interval, check.names = FALSE))
    }

    # Interaction creates a stable group index without converting the original
    # identifiers themselves; their types and displayed values are copied from
    # the first posterior row in each group below.
    grouping_factors <- lapply(
      object[identity_columns],
      factor,
      exclude = NULL
    ) # factors retaining missing identifiers as an explicit level
    coefficient_group <- do.call(
      interaction,
      c(grouping_factors, list(drop = TRUE, lex.order = TRUE))
    ) # draw-table row membership for each distinct coefficient
    group_rows <- split(seq_len(nrow(object)), coefficient_group)

    interval_rows <- lapply(group_rows, function(row_index) {
      coefficient_values <- as.numeric(object$value[row_index]) # posterior sample for one coefficient
      coefficient_name <- paste(
        vapply(object[identity_columns][row_index[1L], , drop = FALSE], as.character, character(1L)),
        collapse = ": "
      ) # temporary matrix label used only during the rstantools calculation
      interval <- rstantools::posterior_interval(
        matrix(coefficient_values, ncol = 1L, dimnames = list(NULL, coefficient_name)),
        prob = prob,
        ...
      )
      cbind(
        object[row_index[1L], identity_columns, drop = FALSE],
        as.data.frame(interval, check.names = FALSE),
        row.names = NULL
      )
    })
    return(do.call(rbind, interval_rows))
  }

  # Nested coefficient views are traversed without changing their names or
  # statistical hierarchy. Null and empty components remain absent or empty,
  # matching the corresponding posterior_summary() result.
  if (is.list(object)) {
    return(lapply(
      object,
      .joinme_posterior_interval_structure,
      prob = prob,
      ...
    ))
  }

  object
}

#' Central posterior credible intervals for dynamic predictions
#'
#' @description
#' Flattens the stored dynamic-prediction draws using [posterior_draws()] and
#' computes central credible intervals with
#' [rstantools::posterior_interval()]. Composite parameter names retain the
#' subject, marker, prediction scale, and evaluation-time identities.
#'
#' @inheritParams posterior_interval.JoiNMeFit
#' @param object A `JoiNMeDynPred` dynamic-prediction object. Latent-class
#'   predictions inherit this method through `JoiNMeMixDynPred`.
#'
#' @return A numeric matrix with one row per selected prediction quantity and
#'   two probability-labelled interval columns.
#'
#' @importFrom rstantools posterior_interval
#' @export
posterior_interval.JoiNMeDynPred <- function(
  object,
  prob = 0.9,
  variables = NULL,
  regex = FALSE,
  draws = NULL,
  seed = 1,
  ...
) {
  posterior_matrix <- posterior_draws(
    object,
    variables = variables,
    regex = regex,
    draws = draws,
    seed = seed,
    format = "draws_matrix"
  ) # flattened dynamic-prediction draws with explicit scientific identities
  rstantools::posterior_interval(
    as.matrix(posterior_matrix),
    prob = prob,
    ...
  )
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
