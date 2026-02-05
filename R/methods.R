#' @name joinme_methods
#' @title joinme S3 Methods
#'
#' @description
#' S3 methods for joinme R6 objects with cached summaries
#' @keywords internal
#' @author Trinh Dong
NULL

# ---- internal helpers -----------------------------------------------------

#' @keywords internal
.summarise_draws_diag <- function(fit, variables, draws = NULL, seed = 1) {
  if (length(variables) == 0) {
    return(data.frame(variable = character(0)))
  }

  draws_obj <- .get_draws_obj(fit, variables = variables, draws = draws, seed = seed)

  ddf <- posterior::as_draws_df(draws_obj)
  sum_list <- lapply(variables, function(v) .summarize_draw_col(ddf[[v]]))
  sum_df <- data.frame(
    variable = variables,
    do.call(rbind, sum_list),
    row.names = NULL,
    check.names = FALSE
  )

  diag_df <- suppressWarnings(tryCatch(
    posterior::summarise_draws(draws_obj, "rhat", "ess_bulk"),
    error = function(e) NULL
  ))
  if (!is.null(diag_df)) {
    diag_df <- diag_df[, c("variable", "rhat", "ess_bulk"), drop = FALSE]
    names(diag_df) <- c("variable", "Rhat", "ESS")
    sum_df <- merge(sum_df, diag_df, by = "variable", all.x = TRUE, sort = FALSE)
  } else {
    sum_df$Rhat <- NA_real_
    sum_df$ESS <- NA_real_
  }
  sum_df
}

#' @keywords internal
.joinme_sampler_diagnostics <- function(fit) {
  out <- list(
    draws = NA_integer_,
    divergences = NA_integer_,
    treedepth_hits = NA_integer_,
    ebfmi_min = NA_real_,
    max_rhat = NA_real_,
    min_ess = NA_real_
  )

  out$draws <- tryCatch(posterior::ndraws(.get_draws_obj(fit)), error = function(e) NA_integer_)

  if (.is_cmdstanr_fit(fit)) {
    diag <- tryCatch(fit$diagnostic_summary(), error = function(e) NULL)
    if (!is.null(diag)) {
      if ("num_divergent" %in% names(diag)) {
        out$divergences <- sum(diag$num_divergent, na.rm = TRUE)
      }
      if ("num_max_treedepth" %in% names(diag)) {
        out$treedepth_hits <- sum(diag$num_max_treedepth, na.rm = TRUE)
      }
      if ("ebfmi" %in% names(diag)) {
        out$ebfmi_min <- suppressWarnings(min(diag$ebfmi, na.rm = TRUE))
        if (!is.finite(out$ebfmi_min)) out$ebfmi_min <- NA_real_
      }
    }

    sumdf <- tryCatch(fit$summary(), error = function(e) NULL)
    if (!is.null(sumdf) && all(c("rhat", "ess_bulk") %in% names(sumdf))) {
      out$max_rhat <- suppressWarnings(max(sumdf$rhat, na.rm = TRUE))
      out$min_ess <- suppressWarnings(min(sumdf$ess_bulk, na.rm = TRUE))
      if (!is.finite(out$max_rhat)) out$max_rhat <- NA_real_
      if (!is.finite(out$min_ess)) out$min_ess <- NA_real_
    }
  } else if (.is_rstan_fit(fit)) {
    params <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    if (!is.null(params)) {
      out$divergences <- sum(vapply(params, function(x) sum(x[, "divergent__"], na.rm = TRUE), numeric(1)))
    }
    sumdf <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
    if (!is.null(sumdf)) {
      if ("Rhat" %in% colnames(sumdf)) {
        out$max_rhat <- suppressWarnings(max(sumdf[, "Rhat"], na.rm = TRUE))
        if (!is.finite(out$max_rhat)) out$max_rhat <- NA_real_
      }
      if ("n_eff" %in% colnames(sumdf)) {
        out$min_ess <- suppressWarnings(min(sumdf[, "n_eff"], na.rm = TRUE))
        if (!is.finite(out$min_ess)) out$min_ess <- NA_real_
      }
    }
  }

  out
}

# ---- print/summary --------------------------------------------------------

#' Print a joinme object
#'
#' @param x A joinme fit object.
#' @param ... Unused.
#'
#' @return Invisibly returns the object.
#' @export
print.JoinMeFit <- function(x, ...) {
  assertthat::assert_that(inherits(x, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  cat("JoinMe model fit\n")
  cat("===============\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  family <- x$stan_data$family_names %||% .family_code_to_name(x$config$family_long)
  if (!is.null(family)) {
    cat("Family: ", paste(family, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$tmax)) {
    cat("tmax: ", x$tmax, "\n", sep = "")
  }
  cat("Use summary() for parameter summaries.\n")
  invisible(x)
}

#' Summarize a joinme object
#'
#' @param object A joinme fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param include_vcov Logical; include covariance summaries.
#' @param ... Unused.
#'
#' @return A summary_JoinMeFit object.
#' @export
summary.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3,
                           include_vcov = FALSE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  if (is.null(draws)) draws <- cfg$draws_default
  assertthat::assert_that(is.null(draws) || (is.numeric(draws) && draws > 0),
                          msg = "draws must be NULL or a positive number.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- paste0("summary_", draws, "_", digits, "_", include_vcov)
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  diag <- .joinme_sampler_diagnostics(fit)
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  beta_vars <- paste0("beta[", seq_len(sd$P), "]")
  s_beta <- as.data.frame(.summarise_draws_diag(fit, beta_vars, draws = draws, seed = seed))
  s_beta$term <- if (!is.null(sd$x_cols) && length(sd$x_cols) == nrow(s_beta)) sd$x_cols else s_beta$variable
  s_beta <- s_beta[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
  s_beta$Estimate <- round(s_beta$Estimate, digits)
  s_beta$Est.Error <- round(s_beta$Est.Error, digits)
  s_beta$Q2.5 <- round(s_beta$Q2.5, digits)
  s_beta$Q97.5 <- round(s_beta$Q97.5, digits)
  s_beta$Rhat <- round(s_beta$Rhat, 3)
  s_beta$ESS <- round(s_beta$ESS, 1)

  s_g <- NULL
  if (sd$p_w > 0) {
    g_vars <- paste0("gamma_w[", seq_len(sd$p_w), "]")
    g_vars <- g_vars[g_vars %in% all_vars]
    if (length(g_vars) > 0) {
      s_g <- as.data.frame(.summarise_draws_diag(fit, g_vars, draws = draws, seed = seed))
      s_g$term <- if (!is.null(sd$w_cols) && length(sd$w_cols) == nrow(s_g)) sd$w_cols else s_g$variable
      s_g <- s_g[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
      s_g$Estimate <- round(s_g$Estimate, digits)
      s_g$Est.Error <- round(s_g$Est.Error, digits)
      s_g$Q2.5 <- round(s_g$Q2.5, digits)
      s_g$Q97.5 <- round(s_g$Q97.5, digits)
      s_g$Rhat <- round(s_g$Rhat, 3)
      s_g$ESS <- round(s_g$ESS, 1)
    }
  }

  a_vars <- c("alpha_cv_total", "alpha_cv_mean", "alpha_cv_marker", "alpha_cs_total", "alpha_cs_mean", "alpha_cs_marker")
  a_vars <- c(a_vars, grep("^alpha_vcov_var\\[", all_vars, value = TRUE))
  a_vars <- a_vars[a_vars %in% all_vars]
  active_vars <- c(
    if (isTRUE(sd$assoc_cv_total == 1)) "alpha_cv_total",
    if (isTRUE(sd$assoc_cv_mean == 1)) "alpha_cv_mean",
    if (isTRUE(sd$assoc_cv_marker == 1)) "alpha_cv_marker",
    if (isTRUE(sd$assoc_cs_total == 1)) "alpha_cs_total",
    if (isTRUE(sd$assoc_cs_mean == 1)) "alpha_cs_mean",
    if (isTRUE(sd$assoc_cs_marker == 1)) "alpha_cs_marker",
    if (isTRUE(sd$assoc_vcov == 1)) grep("^alpha_vcov_var\\[", a_vars, value = TRUE)
  )
  if (length(active_vars) > 0) {
    a_vars <- a_vars[a_vars %in% active_vars]
  }
  s_a <- NULL
  if (length(a_vars) > 0) {
    s_a <- as.data.frame(.summarise_draws_diag(fit, a_vars, draws = draws, seed = seed))
    assoc_map <- c(
      alpha_cv_total = "cv_total",
      alpha_cv_mean = "cv_mean",
      alpha_cv_marker = "cv_marker",
      alpha_cs_total = "cs_total",
      alpha_cs_mean = "cs_mean",
      alpha_cs_marker = "cs_marker"
    )
    s_a$term <- assoc_map[s_a$variable]
    s_a$term[is.na(s_a$term)] <- sub("^alpha_", "", s_a$variable[is.na(s_a$term)])
    s_a <- s_a[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
    s_a$Estimate <- round(s_a$Estimate, digits)
    s_a$Est.Error <- round(s_a$Est.Error, digits)
    s_a$Q2.5 <- round(s_a$Q2.5, digits)
    s_a$Q97.5 <- round(s_a$Q97.5, digits)
    s_a$Rhat <- round(s_a$Rhat, 3)
    s_a$ESS <- round(s_a$ESS, 1)
  }

  # Marker-weight association summaries (when estimated)
  if (sd$D > 0) {
    mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D), "]")
    mw_vars <- mw_vars[mw_vars %in% all_vars]
    if (length(mw_vars) > 0) {
      s_mw <- as.data.frame(.summarise_draws_diag(fit, mw_vars, draws = draws, seed = seed))
      marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
      if (length(marker_terms) == nrow(s_mw)) {
        s_mw$term <- paste0("weight: ", marker_terms)
      } else {
        s_mw$term <- s_mw$variable
      }
      s_mw <- s_mw[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
      s_mw$Estimate <- round(s_mw$Estimate, digits)
      s_mw$Est.Error <- round(s_mw$Est.Error, digits)
      s_mw$Q2.5 <- round(s_mw$Q2.5, digits)
      s_mw$Q97.5 <- round(s_mw$Q97.5, digits)
      s_mw$Rhat <- round(s_mw$Rhat, 3)
      s_mw$ESS <- round(s_mw$ESS, 1)
      s_a <- if (is.null(s_a)) s_mw else rbind(s_a, s_mw)
    } else if (isTRUE(as.logical(sd$estimate_marker_weights)) && !is.null(sd$marker_weights)) {
      marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
      s_mw <- data.frame(
        term = paste0("weight: ", marker_terms),
        Estimate = round(as.numeric(sd$marker_weights), digits),
        Est.Error = NA_real_,
        Q2.5 = NA_real_,
        Q97.5 = NA_real_,
        Rhat = NA_real_,
        ESS = NA_real_,
        stringsAsFactors = FALSE
      )
      s_a <- if (is.null(s_a)) s_mw else rbind(s_a, s_mw)
    }
  }

  dist_candidates <- c(
    "sigma_y",
    grep("^sigma_marker\\[", all_vars, value = TRUE),
    grep("^nu_marker\\[", all_vars, value = TRUE),
    grep("^phi_nb_marker\\[", all_vars, value = TRUE),
    grep("^alpha_skew_marker\\[", all_vars, value = TRUE)
  )
  dist_vars <- dist_candidates[dist_candidates %in% all_vars]
  s_d <- NULL
  if (length(dist_vars) > 0) {
    s_d <- as.data.frame(.summarise_draws_diag(fit, dist_vars, draws = draws, seed = seed))
    s_d$term <- s_d$variable
    s_d <- s_d[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
    s_d$Estimate <- round(s_d$Estimate, digits)
    s_d$Est.Error <- round(s_d$Est.Error, digits)
    s_d$Q2.5 <- round(s_d$Q2.5, digits)
    s_d$Q97.5 <- round(s_d$Q97.5, digits)
    s_d$Rhat <- round(s_d$Rhat, 3)
    s_d$ESS <- round(s_d$ESS, 1)
  }

  s_dr <- NULL
  dist_cols <- cfg$dist$dist_cols %||% list()
  dist_reg_specs <- list(
    sigma = list(prefix = "beta_sigma", cols = dist_cols$sigma %||% character(0)),
    nu = list(prefix = "beta_nu", cols = dist_cols$nu %||% character(0)),
    phi = list(prefix = "beta_phi", cols = dist_cols$phi %||% character(0)),
    alpha = list(prefix = "beta_alpha", cols = dist_cols$alpha %||% character(0))
  )
  dist_reg_tables <- list()
  for (nm in names(dist_reg_specs)) {
    spec <- dist_reg_specs[[nm]]
    vars <- grep(paste0("^", spec$prefix, "\\["), all_vars, value = TRUE)
    if (length(vars) > 0) {
      tmp <- as.data.frame(.summarise_draws_diag(fit, vars, draws = draws, seed = seed))
      if (length(spec$cols) == nrow(tmp)) {
        tmp$term <- paste0(nm, ": ", spec$cols)
      } else {
        tmp$term <- tmp$variable
      }
      tmp$parameter <- nm
      tmp <- tmp[, c("parameter", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
      tmp$Estimate <- round(tmp$Estimate, digits)
      tmp$Est.Error <- round(tmp$Est.Error, digits)
      tmp$Q2.5 <- round(tmp$Q2.5, digits)
      tmp$Q97.5 <- round(tmp$Q97.5, digits)
      tmp$Rhat <- round(tmp$Rhat, 3)
      tmp$ESS <- round(tmp$ESS, 1)
      dist_reg_tables[[nm]] <- tmp
    }
  }
  if (length(dist_reg_tables) > 0) {
    s_dr <- do.call(rbind, dist_reg_tables)
  }

  vcov_tables <- NULL
  if (isTRUE(include_vcov)) {
    vcov_tables <- list(
      id = vcov(object, what = "id", draws = draws),
      marker = if (sd$R_mk > 0) vcov(object, what = "marker", draws = draws) else NULL,
      marker_by_id_latent = vcov(object, what = "marker_by_id_latent", draws = draws)
    )
  }

  summary_obj <- SummaryJoinMeFit$new(
    tables = list(
      fixef = s_beta,
      gamma_w = s_g,
      assoc = s_a,
      distributional = s_d,
      distributional_regression = s_dr,
      vcov = vcov_tables
    ),
    diagnostics = diag,
    metadata = list(
      family = sd$family_names %||% .family_code_to_name(cfg$family_long),
      tmax = sd$tmax %||% cfg$tmax %||% 1.0,
      draws = draws,
      transforms = cfg$transforms
    )
  )

  object$cache_set(cache_key, summary_obj)
  summary_obj
}

#' Update a joinme fit
#'
#' @description
#' Refits a joinme model using the stored call, with optional updates to formulas,
#' data, and control arguments. This mirrors the pattern used by brms::update.
#'
#' @param object A joinme fit object.
#' @param formulaLong Optional updated longitudinal formula. Use an update formula
#'   (e.g., `~ . + x`) to modify the existing model.
#' @param dataLong Optional updated longitudinal dataset.
#' @param formulaEvent Optional updated survival formula (full or update form).
#' @param dataEvent Optional updated event dataset.
#' @param formulaVcov Optional updated covariance formula (full or update form).
#' @param formulaDist Optional distributional regression formulas with parameter
#'   names on the LHS (e.g., `sigma ~ 1 + time`).
#' @param control Optional updated control list.
#' @param draws Optional draws override.
#' @param families Optional updated families specification.
#' @param transforms Optional updated transforms specification.
#' @param priors Optional updated priors list.
#' @param ... Additional arguments passed to `joinme()`.
#'
#' @return A refitted joinme object.
#' @export
update.JoinMeFit <- function(
  object,
  formulaLong = NULL,
  dataLong = NULL,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaVcov = NULL,
  formulaDist = NULL,
  control = NULL,
  draws = NULL,
  families = NULL,
  transforms = NULL,
  priors = NULL,
  ...
) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  call_obj <- object$call
  if (is.null(call_obj) || !is.call(call_obj)) {
    call_obj <- call("joinme")
  }
  call_obj[[1]] <- quote(joinme)

  update_formula <- function(current, updated, name) {
    if (is.null(updated)) return(current)
    if (!inherits(updated, "formula")) {
      cli::cli_abort("{.arg {name}} must be a formula.")
    }
    if (any(all.vars(updated) == ".")) {
      return(stats::update(current, updated))
    }
    updated
  }

  call_obj$formulaLong <- update_formula(object$formulaLong, formulaLong, "formulaLong")
  call_obj$formulaEvent <- update_formula(object$formulaEvent, formulaEvent, "formulaEvent")
  call_obj$formulaVcov <- update_formula(object$formulaVcov, formulaVcov, "formulaVcov")
  if (!is.null(formulaDist)) call_obj$formulaDist <- formulaDist

  if (!is.null(dataLong)) call_obj$dataLong <- dataLong
  if (!is.null(dataEvent)) call_obj$dataEvent <- dataEvent
  if (!is.null(control)) call_obj$control <- control
  if (!is.null(draws)) call_obj$draws <- draws
  if (!is.null(families)) call_obj$families <- families
  if (!is.null(transforms)) call_obj$transforms <- transforms
  if (!is.null(priors)) call_obj$priors <- priors

  dots <- list(...)
  if (length(dots) > 0) {
    if (is.null(names(dots)) || any(names(dots) == "")) {
      cli::cli_abort("All additional arguments must be named.")
    }
    for (nm in names(dots)) {
      call_obj[[nm]] <- dots[[nm]]
    }
  }

  if (is.null(call_obj$dataLong)) {
    cli::cli_abort("{.arg dataLong} must be supplied when it is not available in the stored call.")
  }
  if (is.null(call_obj$dataEvent)) {
    cli::cli_abort("{.arg dataEvent} must be supplied when it is not available in the stored call.")
  }

  eval(call_obj, parent.frame())
}

#' Print a summary_JoinMeFit object
#'
#' @param x A summary_JoinMeFit object.
#' @param ... Unused.
#'
#' @return Invisibly returns the summary object.
#' @export
print.summary_JoinMeFit <- function(x, ...) {
  cat("Joint mixed effects model summary\n")
  cat("=================================\n")
  if (!is.null(x$metadata$family)) {
    family <- x$metadata$family
    family <- if (length(family) == 1) {
      unique(family)
    } else {
      fm_collapsed <- paste(head(family, 5), collapse = ", ")
      if (length(family) > 5) {
        fm_collapsed <- paste0(fm_collapsed, ", ... (total ", length(family), " families)")
      }
      fm_collapsed
    }

    cat("Family: ", family, "\n", sep = "")
  }
  if (!is.null(x$metadata$tmax)) {
    cat("tmax: ", x$metadata$tmax, "\n", sep = "")
  }
  cat("\n")
  if (!is.null(x$diagnostics)) {
    cat("Sampler diagnostics\n")
    cat("-------------------\n")
    cat(sprintf("Draws: %s\n", x$diagnostics$draws %||% NA))
    cat(sprintf("Divergences: %s\n", x$diagnostics$divergences %||% NA))
    cat(sprintf("Max treedepth hits: %s\n", x$diagnostics$treedepth_hits %||% NA))
    cat(sprintf("Min E-BFMI: %s\n", ifelse(is.finite(x$diagnostics$ebfmi_min), sprintf("%.3f", x$diagnostics$ebfmi_min), "NA")))
    cat(sprintf("Max R-hat: %s\n", ifelse(is.finite(x$diagnostics$max_rhat), sprintf("%.3f", x$diagnostics$max_rhat), "NA")))
    cat(sprintf("Min bulk ESS: %s\n", ifelse(is.finite(x$diagnostics$min_ess), sprintf("%.1f", x$diagnostics$min_ess), "NA")))
  }
  if (!is.null(x$tables$fixef)) {
    cat("\nFixed effects (beta)\n")
    cat("---------------------\n")
    print(x$tables$fixef, row.names = FALSE)
  }
  if (!is.null(x$tables$gamma_w)) {
    cat("\nSurvival effects (gamma_w)\n")
    cat("---------------------------\n")
    print(x$tables$gamma_w, row.names = FALSE)
  }
  if (!is.null(x$tables$assoc)) {
    cat("\nAssociation parameters\n")
    cat("-----------------------\n")
    print(x$tables$assoc, row.names = FALSE)
  }
  if (!is.null(x$tables$distributional)) {
    cat("\nDistributional parameters\n")
    cat("---------------------------\n")
    print(x$tables$distributional, row.names = FALSE)
  }
  if (!is.null(x$tables$distributional_regression)) {
    cat("\nDistributional regression\n")
    cat("---------------------------\n")
    print(x$tables$distributional_regression, row.names = FALSE)
  }
  invisible(x)
}

# ---- fixef / ranef --------------------------------------------------------

#' Extract fixed effects
#'
#' @param object A joinme fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param ... Unused.
#'
#' @return A data.frame of fixed effects summaries.
#' @importFrom lme4 fixef
#' @export
fixef.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  cache_key <- paste0("fixef_", draws, "_", digits)
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  beta_vars <- all_vars[grepl("^beta\\[", all_vars)]
  if (length(beta_vars) == 0) beta_vars <- paste0("beta[", seq_len(sd$P), "]")

  s <- as.data.frame(.summarise_draws_diag(fit, beta_vars, draws = draws, seed = seed))
  s$term <- if (!is.null(sd$x_cols) && length(sd$x_cols) == nrow(s)) sd$x_cols else s$variable
  out <- s[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]

  out$Estimate <- round(out$Estimate, digits)
  out$Est.Error <- round(out$Est.Error, digits)
  out$Q2.5 <- round(out$Q2.5, digits)
  out$Q97.5 <- round(out$Q97.5, digits)
  out$Rhat <- round(out$Rhat, 3)
  out$ESS <- round(out$ESS, 1)

  # Append marker-weight summaries as association-like terms
  if (sd$D > 0) {
    mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D), "]")
    mw_vars <- mw_vars[mw_vars %in% all_vars]
    if (length(mw_vars) > 0) {
      mw <- as.data.frame(.summarise_draws_diag(fit, mw_vars, draws = draws, seed = seed))
      marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
      if (length(marker_terms) == nrow(mw)) {
        mw$term <- paste0("weight: ", marker_terms)
      } else {
        mw$term <- mw$variable
      }
      mw <- mw[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
      mw$Estimate <- round(mw$Estimate, digits)
      mw$Est.Error <- round(mw$Est.Error, digits)
      mw$Q2.5 <- round(mw$Q2.5, digits)
      mw$Q97.5 <- round(mw$Q97.5, digits)
      mw$Rhat <- round(mw$Rhat, 3)
      mw$ESS <- round(mw$ESS, 1)
      out <- rbind(out, mw)
    }
  }
  object$cache_set(cache_key, out)
  out
}

#' Extract random effects
#'
#' @param object A joinme fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param ... Unused.
#'
#' @return A list of data.frames for random effects.
#' @importFrom lme4 ranef
#' @export
ranef.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  cache_key <- paste0("ranef_", draws, "_", digits)
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  summarize_block <- function(varnames, group) {
    if (length(varnames) == 0) return(NULL)
    s <- as.data.frame(.summarise_draws_diag(fit, varnames, draws = draws, seed = seed))
    out <- s[, c("variable", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ESS"), drop = FALSE]
    names(out)[1] <- "term"
    out$group <- group

    out$Estimate <- round(out$Estimate, digits)
    out$Est.Error <- round(out$Est.Error, digits)
    out$Q2.5 <- round(out$Q2.5, digits)
    out$Q97.5 <- round(out$Q97.5, digits)
    out$Rhat <- round(out$Rhat, 3)
    out$ESS <- round(out$ESS, 1)
    out
  }

  u_vars <- vars[grepl("^u_id\\[", vars)]
  v_vars <- vars[grepl("^v_marker\\[", vars)]
  zw_vars <- vars[grepl("^z_w\\[", vars)]

  out <- list(
    id = summarize_block(u_vars, "id"),
    marker = if (sd$R_mk > 0) summarize_block(v_vars, "marker") else NULL,
    marker_byid_latent = summarize_block(zw_vars, "marker_byid_latent")
  )

  mw_vars <- vars[grepl("^marker_weights_eff\\[", vars)]
  if (length(mw_vars) > 0) {
    mw <- summarize_block(mw_vars, "assoc_weight")
    if (!is.null(mw) && !is.null(sd$marker_levels)) {
      marker_terms <- sd$marker_levels
      if (length(marker_terms) == nrow(mw)) {
        mw$term <- paste0("weight: ", marker_terms)
      }
    }
    out$assoc_weight <- mw
  } else if (isTRUE(as.logical(sd$estimate_marker_weights)) && !is.null(sd$marker_weights)) {
    marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
    mw <- data.frame(
      term = paste0("weight: ", marker_terms),
      Estimate = as.numeric(sd$marker_weights),
      Est.Error = NA_real_,
      Q2.5 = NA_real_,
      Q97.5 = NA_real_,
      Rhat = NA_real_,
      ESS = NA_real_,
      group = "assoc_weight",
      stringsAsFactors = FALSE
    )
    out$assoc_weight <- mw
  }
  object$cache_set(cache_key, out)
  out
}

# ---- vcov -----------------------------------------------------------------

#' Extract covariance summaries
#'
#' @param object A joinme fit object.
#' @param what Covariance block: "id", "marker", or "marker_by_id_latent".
#' @param draws Number of draws to use for summaries.
#' @param ... Unused.
#'
#' @return A data.frame with covariance summaries.
#' @export
vcov.JoinMeFit <- function(object, what = c("id", "marker", "marker_by_id_latent"), draws = NULL, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  what <- match.arg(what)
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  summarize_cov_matrix <- function(tau_prefix, Lcorr_prefix, dim, label) {
    if (dim <= 0) {
      return(data.frame(
        block = label, row = integer(0), col = integer(0),
        Estimate = numeric(0), Est.Error = numeric(0), Q2.5 = numeric(0), Q97.5 = numeric(0)
      ))
    }
    tau_vars <- paste0(tau_prefix, "[", seq_len(dim), "]")
    L_vars <- as.vector(outer(seq_len(dim), seq_len(dim),
                              function(r, c) paste0(Lcorr_prefix, "[", r, ",", c, "]")))
    vars <- c(tau_vars, L_vars)

    ddf <- .get_draws_df(fit, variables = vars, draws = draws, seed = 1)
    tau_draws <- as.matrix(ddf[, tau_vars, drop = FALSE])
    L_draws <- as.matrix(ddf[, L_vars, drop = FALSE])
    nD <- nrow(tau_draws)

    Sigma_list <- vector("list", nD)
    for (i in seq_len(nD)) {
      tau <- as.numeric(tau_draws[i, ])
      L <- matrix(as.numeric(L_draws[i, ]), nrow = dim, ncol = dim, byrow = FALSE)
      L[upper.tri(L)] <- 0
      Corr <- L %*% t(L)
      Sigma_list[[i]] <- diag(tau, dim) %*% Corr %*% diag(tau, dim)
    }

    out <- list()
    for (r in seq_len(dim)) for (c in seq_len(dim)) {
      vals <- vapply(Sigma_list, function(S) S[r, c], numeric(1))
      out[[length(out) + 1]] <- cbind(block = label, row = r, col = c, t(.summarize_draw_col(vals)))
    }
    df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
    df$row <- as.integer(df$row)
    df$col <- as.integer(df$col)
    df$Estimate <- as.numeric(df$Estimate)
    df$Est.Error <- as.numeric(df$Est.Error)
    df$Q2.5 <- as.numeric(df$Q2.5)
    df$Q97.5 <- as.numeric(df$Q97.5)
    df
  }

  if (what == "id") return(summarize_cov_matrix("tau_u", "Lcorr_u", sd$R_id, "Sigma_u"))
  if (what == "marker") {
    if (sd$R_mk <= 0) {
      cli::cli_abort("Marker vcov requested but R_mk = 0.")
    }
    return(summarize_cov_matrix("tau_v", "Lcorr_v", sd$R_mk, "Sigma_v"))
  }
  summarize_cov_matrix("tau_w", "Lcorr_w", sd$Q_idm, "Sigma_w")
}
