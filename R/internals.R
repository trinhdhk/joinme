
# ---- internal helpers -----------------------------------------------------

#' @keywords internal
.summarise_draws_diag <- function(fit, variables, draws = NULL, seed = 1) {
  # Summarise draws with point estimates + diagnostics in one table
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
    posterior::summarise_draws(draws_obj, "rhat", "ess_bulk", "ess_tail"),
    error = function(e) NULL
  ))
  if (!is.null(diag_df)) {
    diag_df <- diag_df[, c("variable", "rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    names(diag_df) <- c("variable", "Rhat", "ess_bulk", "ess_tail")
    sum_df <- merge(sum_df, diag_df, by = "variable", all.x = TRUE, sort = FALSE)
  } else {
    sum_df$Rhat <- NA_real_
    sum_df$ess_bulk <- NA_real_
    sum_df$ess_tail <- NA_real_
  }
  sum_df
}

#' @keywords internal
#' @noRd
.format_transform_expr <- function(expr) {
  if (is.null(expr)) return("identity")
  if (rlang::is_quosure(expr)) expr <- rlang::get_expr(expr)
  if (inherits(expr, "formula")) expr <- expr[[2]]
  if (is.character(expr)) return(paste0("~ ", expr))
  if (is.call(expr) || is.name(expr) || is.numeric(expr)) {
    return(paste0("~ ", paste(deparse(expr), collapse = "")))
  }
  "functional"
}

#' @keywords internal
#' @noRd
.transform_expr_with_affine_shift <- function(expr, iota_nodes = list()) {
  if (is.null(iota_nodes) || !length(iota_nodes)) {
    return(expr)
  }

  if (rlang::is_quosure(expr)) expr <- rlang::get_expr(expr)
  if (inherits(expr, "formula")) expr <- expr[[2]]

  iota_map <- if (length(iota_nodes)) {
    stats::setNames(iota_nodes, vapply(iota_nodes, function(node) node$path %||% "", character(1)))
  } else {
    list()
  }
  build_iota_symbol <- function(prefix, idx) {
    parse(text = paste0(prefix, "[", idx, "]"))[[1]]
  }

  rewrite_node <- function(node, path = "root") {
    if (!is.call(node)) {
      return(node)
    }
    parts <- as.list(node)
    if (length(parts) > 1L) {
      for (idx in seq_along(parts[-1])) {
        parts[[idx + 1L]] <- rewrite_node(parts[[idx + 1L]], path = paste0(path, "/", idx))
      }
    }
    info <- iota_map[[path]]
    if (!is.null(info) && length(parts) >= 2L) {
      base_expr <- parts[[2L]]
      intercept_expr <- if (isTRUE((info$intercept_index %||% 0L) > 0L)) build_iota_symbol("iota_1", info$intercept_index) else NULL
      slope_expr <- if (isTRUE((info$slope_index %||% 0L) > 0L)) build_iota_symbol("iota_2", info$slope_index) else NULL
      shifted <- if (!is.null(intercept_expr) && !is.null(slope_expr)) {
        bquote(.(intercept_expr) + .(slope_expr) * .(base_expr))
      } else if (!is.null(intercept_expr)) {
        bquote(.(intercept_expr) + .(base_expr))
      } else if (!is.null(slope_expr)) {
        bquote(.(slope_expr) * .(base_expr))
      } else {
        base_expr
      }
      parts[[2L]] <- shifted
    }
    as.call(parts)
  }

  rewrite_node(expr)
}

#' @keywords internal
#' @noRd
.format_transform_spec <- function(spec) {
  if (is.null(spec) || is.null(spec$type) || spec$type == "identity") {
    return("identity")
  }
  if (spec$type == "functional") {
    return(.format_transform_expr(.transform_expr_with_affine_shift(
      spec$expr,
      iota_nodes = .transform_iota_nodes(spec)
    )))
  }
  spec$type <- .canonicalise_transform_type(spec$type)
  if (spec$type %in% c("ispline", "ispline_penalised", "pmonospline", "pmono", "ispline_expit", "ispline_expit_penalised")) {
    knots <- spec$raw_knots %||% spec$knots %||% spec$x
    degree <- spec$degree %||% 3L
    label <- if (spec$type %in% c("ispline_expit", "ispline_expit_penalised")) "ispline_expit" else "ispline"
    direction <- spec$direction %||% if (!is.null(spec$spline_direction)) .monotone_direction_label(spec$spline_direction) else NULL
    if (!is.null(knots)) {
      knot_text <- paste(format(knots, digits = 3, trim = TRUE), collapse = ", ")
      direction_text <- if (!is.null(direction)) paste0(", direction = ", direction) else ""
      return(paste0(label, "(knots = c(", knot_text, "), degree = ", degree, direction_text, ")"))
    }
    direction_text <- if (!is.null(direction)) paste0(", direction = ", direction) else ""
    return(paste0(label, "(degree = ", degree, direction_text, ")"))
  }
  if (spec$type == "pwlin") {
    knots <- spec$knots %||% spec$cutpoints %||% spec$x
    knot_text <- if (!is.null(knots)) paste(format(knots, digits = 3, trim = TRUE), collapse = ", ") else ""
    direction <- .monotone_direction_label(spec$direction %||% if (!is.null(spec$y)) .infer_monotone_direction(knots, spec$y) else 1L)
    return(paste0("pwlin(knots = c(", knot_text, "), direction = ", direction, ")"))
  }
  spec$type
}

#' @keywords internal
#' @noRd
.omit_fixed_transform_endpoint_rows <- function(tbl, channel, spec, sd) {
  if (is.null(tbl) || nrow(tbl) == 0L) {
    return(tbl)
  }

  channel_map <- list(
    cv_total = list(estimate = "estimate_spline_cv", n_free = "n_free_spline_cv"),
    cs_total = list(estimate = "estimate_spline_cs", n_free = "n_free_spline_cs"),
    corr = list(estimate = "estimate_spline_corr", n_free = "n_free_spline_corr"),
    vcov = list(estimate = "estimate_spline_vcov", n_free = "n_free_spline_vcov"),
    cv_mean = list(estimate = "estimate_spline_cv_mean", n_free = "n_free_spline_cv_mean"),
    cv_marker = list(estimate = "estimate_spline_cv_marker", n_free = "n_free_spline_cv_marker"),
    cs_mean = list(estimate = "estimate_spline_cs_mean", n_free = "n_free_spline_cs_mean"),
    cs_marker = list(estimate = "estimate_spline_cs_marker", n_free = "n_free_spline_cs_marker")
  )

  map <- channel_map[[channel]]
  if (is.null(map)) {
    return(tbl)
  }

  estimate_flag <- as.integer(sd[[map$estimate]] %||% 0L)
  n_free <- as.integer(sd[[map$n_free]] %||% 0L)
  n_coeff <- as.integer(spec$n %||% 0L)

  if (!isTRUE(estimate_flag == 1L) || n_coeff < 2L || n_free != (n_coeff - 1L)) {
    return(tbl)
  }

  anchored_terms <- c("coeff_1", paste0("coeff_", n_coeff))
  tbl[!(tbl$term %in% anchored_terms), , drop = FALSE]
}

#' @keywords internal
.transform_formulas_from_specs <- function(transforms, sd = NULL) {
  if (is.null(transforms) || length(transforms) == 0) return(NULL)
  terms <- names(transforms)
  if (is.null(terms)) return(NULL)

  rows <- lapply(terms, function(term_name) {
    labels <- term_name
    if (!is.null(sd) && term_name %in% c("corr", "vcov")) {
      labels <- .assoc_transform_component_labels(
        term_name,
        .assoc_transform_component_count(
          term_name,
          sd$Q_idm,
          diagonal_only = identical(term_name, "vcov") && as.integer(sd$indep_idmarker_cov %||% 0L) == 1L
        )
      )
    }
    data.frame(
      term = labels,
      formula = rep(.format_transform_spec(transforms[[term_name]]), length(labels)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' @keywords internal
.distributional_term_map <- function(sd, cfg, all_vars) {
  # Build a family-aware map from Stan variable names to human-readable
  # distributional summary terms.
  #
  # Key calculation rule:
  # - when no distributional regression is used for a parameter, report one
  #   family-level parameter per required family (shared across markers of the
  #   same family),
  # - when regression is used, fixed/random regression summaries are reported in
  #   distributional_regression and baseline family constants are omitted there.
  family_long <- sd$family_long %||% cfg$family_long %||% integer(0)
  D <- sd$D %||% length(family_long)
  if (D <= 0) {
    return(character(0))
  }

  if (length(family_long) < D) {
    family_long <- c(as.integer(family_long), rep(NA_integer_, D - length(family_long)))
  } else {
    family_long <- as.integer(family_long[seq_len(D)])
  }

  marker_levels <- sd$marker_levels %||% paste0("marker_", seq_len(D))
  marker_levels <- as.character(marker_levels)
  if (length(marker_levels) < D) {
    marker_levels <- c(marker_levels, paste0("marker_", seq.int(length(marker_levels) + 1L, D)))
  }
  marker_levels <- marker_levels[seq_len(D)]
  fixed_tau <- .normalise_fixed_tau_by_marker(
    use_tau_fixed = sd$use_tau_fixed,
    tau_fixed = sd$tau_fixed,
    family_codes = family_long
  )

  # Family -> parameter requirements via .family_distrib_params().
  req_by_marker <- lapply(seq_len(D), function(d) {
    fam <- family_long[d]
    if (is.na(fam)) return(character(0))
    tryCatch(.family_distrib_params(fam), error = function(e) character(0))
  })

  term_map <- character(0)

  # Resolve which parameters are regression-driven.
  P_sigma <- as.integer(sd$P_sigma %||% 0L)
  P_nu <- as.integer(sd$P_nu %||% 0L)
  P_phi <- as.integer(sd$P_phi %||% 0L)
  P_alpha <- as.integer(sd$P_alpha %||% 0L)
  P_kappa <- as.integer(sd$P_kappa %||% 0L)
  P_tau <- as.integer(sd$P_tau %||% 0L)
  n_re_sigma <- as.integer(sd$n_re_sigma %||% 0L)
  n_re_nu <- as.integer(sd$n_re_nu %||% 0L)
  n_re_phi <- as.integer(sd$n_re_phi %||% 0L)
  n_re_alpha <- as.integer(sd$n_re_alpha %||% 0L)
  n_re_kappa <- as.integer(sd$n_re_kappa %||% 0L)
  n_re_tau <- as.integer(sd$n_re_tau %||% 0L)

  has_reg <- list(
    sigma = (P_sigma > 0 || n_re_sigma > 0),
    nu = (P_nu > 0 || n_re_nu > 0),
    phi = (P_phi > 0 || n_re_phi > 0),
    alpha = (P_alpha > 0 || n_re_alpha > 0),
    kappa = (P_kappa > 0 || n_re_kappa > 0),
    tau = (P_tau > 0 || n_re_tau > 0)
  )

  # Family names for readable labels.
  fam_name <- function(code) {
    nm <- tryCatch(.family_name(code), error = function(e) as.character(code))
    as.character(nm)
  }

  fam_sets <- list(
    sigma = sort(unique(family_long[vapply(req_by_marker, function(x) "sigma" %in% x, logical(1))])),
    nu = sort(unique(family_long[vapply(req_by_marker, function(x) "nu" %in% x, logical(1))])),
    phi = sort(unique(family_long[vapply(req_by_marker, function(x) "phi" %in% x, logical(1))])),
    alpha = sort(unique(family_long[vapply(req_by_marker, function(x) "alpha" %in% x, logical(1))])),
    kappa = sort(unique(family_long[vapply(req_by_marker, function(x) "kappa" %in% x, logical(1))])),
    tau = sort(unique(family_long[
      vapply(req_by_marker, function(x) "tau" %in% x, logical(1)) &
        fixed_tau$use_tau_fixed == 0L
    ]))
  )

  marker_label_prefix <- list(
    sigma = "sigma_marker",
    nu = "nu_marker",
    phi = "phi_nb_marker",
    alpha = "alpha_skew_marker",
    kappa = "kappa_marker",
    tau = "tau_marker"
  )
  family_label_prefix <- list(
    sigma = "sigma",
    nu = "nu",
    phi = "phi_nb",
    alpha = "alpha_skew",
    kappa = "kappa",
    tau = "tau"
  )

  markers_for_param <- function(param, fam_code) {
    marker_index <- which(
      family_long == fam_code &
        vapply(req_by_marker, function(x) param %in% x, logical(1))
    )
    if (identical(param, "tau")) {
      marker_index <- marker_index[fixed_tau$use_tau_fixed[marker_index] == 0L]
    }
    marker_index
  }

  add_family_terms <- function(param, stan_prefix, fam_codes, use_regression) {
    if (length(fam_codes) == 0 || isTRUE(use_regression)) return(invisible(NULL))
    for (idx in seq_along(fam_codes)) {
      var_nm <- paste0(stan_prefix, "[", idx, "]")
      if (var_nm %in% all_vars) {
        fam_code <- fam_codes[idx]
        if (identical(param, "sigma")) {
          mk_idx <- markers_for_param(param, fam_code)
          if (length(mk_idx) == 1) {
            mk <- marker_levels[mk_idx]
            term_map[var_nm] <<- paste0("sigma[", mk, "]")
          } else {
            term_map[var_nm] <<- paste0("sigma[family=", fam_name(fam_code), "]")
          }
        } else {
          mk_idx <- markers_for_param(param, fam_code)
          if (length(mk_idx) == 1) {
            mk <- marker_levels[mk_idx]
            term_map[var_nm] <<- paste0(marker_label_prefix[[param]], "[", mk, "]")
          } else {
            term_map[var_nm] <<- paste0(family_label_prefix[[param]], "[family=", fam_name(fam_code), "]")
          }
        }
      }
    }
  }

  add_family_terms("sigma", "sigma_family", fam_sets$sigma, has_reg$sigma)
  add_family_terms("nu", "nu_family", fam_sets$nu, has_reg$nu)
  add_family_terms("phi", "phi_family", fam_sets$phi, has_reg$phi)
  add_family_terms("alpha", "alpha_family", fam_sets$alpha, has_reg$alpha)
  add_family_terms("kappa", "kappa_family", fam_sets$kappa, has_reg$kappa)
  add_family_terms("tau", "tau_family", fam_sets$tau, has_reg$tau)

  # Backward-compatible fallback for older fits that still expose marker-level
  # distributional constants.
  for (d in seq_len(D)) {
    req <- req_by_marker[[d]]
    mk <- marker_levels[d]

    if ("sigma" %in% req && !isTRUE(has_reg$sigma)) {
      var_nm <- paste0("sigma_marker[", d, "]")
      term_map[var_nm] <- paste0("sigma[", mk, "]")
    }
    if ("nu" %in% req && !isTRUE(has_reg$nu)) {
      var_nm <- paste0("nu_marker[", d, "]")
      term_map[var_nm] <- paste0("nu_marker[", mk, "]")
    }
    if ("phi" %in% req && !isTRUE(has_reg$phi)) {
      var_nm <- paste0("phi_nb_marker[", d, "]")
      term_map[var_nm] <- paste0("phi_nb_marker[", mk, "]")
    }
    if ("alpha" %in% req && !isTRUE(has_reg$alpha)) {
      var_nm <- paste0("alpha_skew_marker[", d, "]")
      term_map[var_nm] <- paste0("alpha_skew_marker[", mk, "]")
    }
    if ("kappa" %in% req && !isTRUE(has_reg$kappa)) {
      var_nm <- paste0("kappa_marker[", d, "]")
      term_map[var_nm] <- paste0("kappa_marker[", mk, "]")
    }
    if ("tau" %in% req && !isTRUE(has_reg$tau)) {
      var_nm <- paste0("tau_marker[", d, "]")
      term_map[var_nm] <- paste0("tau_marker[", mk, "]")
    }
  }

  # Keep only variables that are present in posterior draws.
  term_map[names(term_map) %in% all_vars]
}

#' Summarise marker-specific fixed skew-Laplace quantiles
#'
#' @description
#' Constructs deterministic summary rows for skew-Laplace `tau` values supplied
#' through [jm_family()]. These constants are part of the statistical model but
#' are not posterior variables; reporting them alongside estimated
#' distributional parameters makes that distinction visible in
#' [summary.JoiNMeFit()].
#'
#' The posterior standard deviation is zero because the value is fixed.
#' Convergence and effective-sample-size diagnostics are undefined and are
#' therefore reported as missing.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param digits Number of decimal places used for reported values.
#'
#' @return A standard JoiNMe posterior-summary data frame, or `NULL` when no
#'   marker has a fixed skew-Laplace quantile.
#' @keywords internal
#' @noRd
.fixed_tau_summary_table <- function(object, digits = 3) {
  sd <- object$stan_data
  family_codes <- as.integer(
    sd$family_long %||%
      object$config$family_long %||%
      integer(0)
  )
  if (length(family_codes) == 0L) {
    return(NULL)
  }

  fixed_tau <- .normalise_fixed_tau_by_marker(
    use_tau_fixed = sd$use_tau_fixed,
    tau_fixed = sd$tau_fixed,
    family_codes = family_codes
  )
  marker_index <- which(fixed_tau$use_tau_fixed == 1L)
  if (length(marker_index) == 0L) {
    return(NULL)
  }

  marker_levels <- as.character(
    sd$marker_levels %||% paste0("marker_", seq_along(family_codes))
  )
  if (length(marker_levels) != length(family_codes)) {
    marker_levels <- paste0("marker_", seq_along(family_codes))
  }
  value <- fixed_tau$tau_fixed[marker_index]

  data.frame(
    term = paste0("tau_fixed[", marker_levels[marker_index], "]"),
    Estimate = round(value, digits),
    Est.Error = 0,
    Q2.5 = round(value, digits),
    Q97.5 = round(value, digits),
    Rhat = NA_real_,
    ess_bulk = NA_real_,
    ess_tail = NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.sampler_diagnostic_draws <- function(fit) {
  # Read the sampler-state variables in a common iteration-by-chain layout.
  # These variables are deliberately kept separate from model parameters by
  # both Stan backends, although the two backends expose them through different
  # methods.
  if (.is_cmdstanr_fit(fit)) {
    sampler_draws <- tryCatch(
      fit$sampler_diagnostics(format = "draws_array"),
      error = function(error) NULL
    ) # CmdStanR sampler states retained in memory by .import_cmdstanr_fit()
    if (is.null(sampler_draws)) {
      return(NULL)
    }
    return(posterior::as_draws_array(sampler_draws))
  }

  if (.is_rstan_fit(fit)) {
    sampler_parameters <- tryCatch(
      rstan::get_sampler_params(fit, inc_warmup = FALSE),
      error = function(error) NULL
    ) # one sampler-state matrix for each rstan chain
    if (is.null(sampler_parameters) || length(sampler_parameters) == 0L) {
      return(NULL)
    }
    common_variables <- Reduce(
      intersect,
      lapply(sampler_parameters, colnames)
    ) # sampler variables available in every chain
    if (length(common_variables) == 0L) {
      return(NULL)
    }
    number_iterations <- min(vapply(
      sampler_parameters,
      nrow,
      integer(1)
    )) # common post-warm-up iteration count across chains
    sampler_array <- array(
      NA_real_,
      dim = c(
        number_iterations,
        length(sampler_parameters),
        length(common_variables)
      ),
      dimnames = list(
        iteration = as.character(seq_len(number_iterations)),
        chain = as.character(seq_along(sampler_parameters)),
        variable = common_variables
      )
    ) # common sampler-state array matching posterior::draws_array conventions
    for (chain_index in seq_along(sampler_parameters)) {
      sampler_array[, chain_index, ] <- sampler_parameters[[chain_index]][
        seq_len(number_iterations),
        common_variables,
        drop = FALSE
      ]
    }
    return(posterior::as_draws_array(sampler_array))
  }

  NULL
}

#' Summarise Hamiltonian energy exploration by chain
#'
#' @description
#' Computes E-BFMI directly from the retained energy sequence and reports the
#' path quantities that help distinguish an isolated warning from persistent
#' slow movement through the posterior energy distribution.  The backend's
#' own E-BFMI is substituted later when it is available, so displayed values
#' agree with CmdStanR or rstan exactly.
#'
#' @param fit A CmdStanR or rstan fit object.
#'
#' @return A data frame with one row per chain, or an empty data frame when the
#'   sampler states are unavailable.
#' @keywords internal
#' @noRd
.sampler_energy_by_chain <- function(fit) {
  sampler_draws <- .sampler_diagnostic_draws(
    fit
  ) # iteration-by-chain sampler states supplied by either supported backend
  if (is.null(sampler_draws)) {
    return(data.frame())
  }
  sampler_array <- as.array(
    sampler_draws
  ) # numeric array used for direct sequential energy calculations
  sampler_variables <- dimnames(sampler_array)[[3L]] %||%
    character(0) # names of the retained sampler-state variables
  if (!"energy__" %in% sampler_variables) {
    return(data.frame())
  }

  value_or_missing <- function(chain_index, variable_name) {
    if (!variable_name %in% sampler_variables) {
      return(rep(NA_real_, dim(sampler_array)[1L]))
    }
    as.numeric(sampler_array[, chain_index, variable_name])
  }

  chain_rows <- vector(
    "list",
    dim(sampler_array)[2L]
  ) # one transparent energy/path summary for each Markov chain
  for (chain_index in seq_len(dim(sampler_array)[2L])) {
    energy <- value_or_missing(
      chain_index,
      "energy__"
    ) # Hamiltonian energy in retained sampling order
    finite_energy <- is.finite(
      energy
    ) # iterations that can contribute to sequential energy diagnostics
    energy <- energy[finite_energy]
    energy_variance <- stats::var(
      energy
    ) # marginal energy variance in this chain
    energy_change <- diff(
      energy
    ) # successive energy movement used in the E-BFMI numerator
    ebfmi <- if (
      length(energy_change) > 0L &&
        is.finite(energy_variance) &&
        energy_variance > 0
    ) {
      mean(energy_change^2) / energy_variance
    } else {
      NA_real_
    } # direct E-BFMI calculation following the Stan diagnostic definition
    energy_lag_one <- if (
      length(energy) > 2L && stats::sd(energy) > 0
    ) {
      stats::cor(energy[-length(energy)], energy[-1L])
    } else {
      NA_real_
    } # persistence of energy from one retained transition to the next
    treedepth <- value_or_missing(
      chain_index,
      "treedepth__"
    ) # binary-tree depth used by the No-U-Turn sampler at each transition
    leapfrog <- value_or_missing(
      chain_index,
      "n_leapfrog__"
    ) # number of numerical integration steps at each transition
    acceptance <- value_or_missing(
      chain_index,
      "accept_stat__"
    ) # realised Metropolis acceptance statistic at each transition
    step_size <- value_or_missing(
      chain_index,
      "stepsize__"
    ) # adapted numerical integration step size recorded per transition
    divergences <- value_or_missing(
      chain_index,
      "divergent__"
    ) # indicator that the numerical Hamiltonian path diverged

    chain_rows[[chain_index]] <- data.frame(
      chain = as.integer(chain_index),
      draws = length(energy),
      E_BFMI = ebfmi,
      energy_SD = stats::sd(energy),
      energy_lag_one = energy_lag_one,
      divergences = if (all(is.na(divergences))) {
        NA_integer_
      } else {
        sum(divergences, na.rm = TRUE)
      },
      median_treedepth = stats::median(treedepth, na.rm = TRUE),
      maximum_treedepth = suppressWarnings(max(treedepth, na.rm = TRUE)),
      median_leapfrog = stats::median(leapfrog, na.rm = TRUE),
      mean_acceptance = mean(acceptance, na.rm = TRUE),
      step_size = stats::median(step_size, na.rm = TRUE),
      stringsAsFactors = FALSE
    ) # chain-specific evidence about energy movement and integration effort
  }
  output <- do.call(
    rbind,
    chain_rows
  ) # combined chain-level diagnostic table
  numeric_columns <- setdiff(
    names(output),
    c("chain", "draws", "divergences")
  ) # continuous path summaries which may be entirely unavailable
  for (column_name in numeric_columns) {
    output[[column_name]][!is.finite(output[[column_name]])] <- NA_real_
  }
  output
}

#' @keywords internal
.joinme_sampler_diagnostics <- function(fit) {
  # Collect sampler diagnostics across cmdstanr or rstan backends
  out <- list(
    draws = NA_integer_,
    divergences = NA_integer_,
    treedepth_hits = NA_integer_,
    ebfmi_min = NA_real_,
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    min_ess_tail = NA_real_,
    min_ess = NA_real_,
    energy_by_chain = data.frame()
  )

  out$energy_by_chain <- .sampler_energy_by_chain(
    fit
  ) # detailed energy and numerical-path assessment retained for diagnosis()

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
        if (
          nrow(out$energy_by_chain) == length(diag$ebfmi)
        ) {
          out$energy_by_chain$E_BFMI <- as.numeric(
            diag$ebfmi
          ) # use the backend value verbatim when CmdStanR supplies it
        }
      }
      if (
        "num_max_treedepth" %in% names(diag) &&
          nrow(out$energy_by_chain) == length(diag$num_max_treedepth)
      ) {
        out$energy_by_chain$treedepth_hits <- as.integer(
          diag$num_max_treedepth
        ) # transitions which reached the configured rather than observed limit
      }
    }

    sumdf <- tryCatch(fit$summary(), error = function(e) NULL)
    if (!is.null(sumdf) && all(c("rhat", "ess_bulk") %in% names(sumdf))) {
      out$max_rhat <- suppressWarnings(max(sumdf$rhat, na.rm = TRUE))
      out$min_ess_bulk <- suppressWarnings(min(sumdf$ess_bulk, na.rm = TRUE))
      if ("ess_tail" %in% names(sumdf)) {
        out$min_ess_tail <- suppressWarnings(min(sumdf$ess_tail, na.rm = TRUE))
      }
      if (!is.finite(out$max_rhat)) out$max_rhat <- NA_real_
      if (!is.finite(out$min_ess_bulk)) out$min_ess_bulk <- NA_real_
      if (!is.finite(out$min_ess_tail)) out$min_ess_tail <- NA_real_
    }
  } else if (.is_rstan_fit(fit)) {
    params <- tryCatch(rstan::get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    if (!is.null(params)) {
      out$divergences <- sum(vapply(params, function(x) sum(x[, "divergent__"], na.rm = TRUE), numeric(1)))
      out$treedepth_hits <- rstan::get_num_max_treedepth(fit)
      out$ebfmi_min <- min(rstan::get_bfmi(fit))
      backend_ebfmi <- as.numeric(
        rstan::get_bfmi(fit)
      ) # rstan's own chain-specific E-BFMI calculation
      if (nrow(out$energy_by_chain) == length(backend_ebfmi)) {
        out$energy_by_chain$E_BFMI <- backend_ebfmi
      }
    }
    sumdf <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
    if (!is.null(sumdf)) {
      if ("Rhat" %in% colnames(sumdf)) {
        out$max_rhat <- suppressWarnings(max(sumdf[, "Rhat"], na.rm = TRUE))
        if (!is.finite(out$max_rhat)) out$max_rhat <- NA_real_
      }
      if ("n_eff" %in% colnames(sumdf)) {
        out$min_ess_bulk <- suppressWarnings(min(sumdf[, "n_eff"], na.rm = TRUE))
        if (!is.finite(out$min_ess_bulk)) out$min_ess_bulk <- NA_real_
      }
    }
  }

  if (is.na(out$min_ess) || !is.finite(out$min_ess)) {
    out$min_ess <- out$min_ess_bulk
  }
  if (!is.finite(out$min_ess_tail)) {
    out$min_ess_tail <- out$min_ess_bulk
  }

  out
}

#' @keywords internal
.build_common_diagnostics_table <- function(draws = NA_real_,
                                            divergences = NA_real_,
                                            treedepth_hits = NA_real_,
                                            ebfmi_min = NA_real_,
                                            max_rhat = NA_real_,
                                            min_ess_bulk = NA_real_,
                                            min_ess_tail = NA_real_,
                                            n_terms_total = NA_real_,
                                            n_terms_bad_rhat = NA_real_,
                                            n_terms_low_ess_bulk = NA_real_,
                                            n_terms_low_ess_tail = NA_real_) {
  data.frame(
    metric = c(
      "draws",
      "divergences",
      "treedepth_hits",
      "ebfmi_min",
      "max_rhat",
      "min_ess_bulk",
      "min_ess_tail",
      "n_terms_total",
      "n_terms_bad_rhat",
      "n_terms_low_ess_bulk",
      "n_terms_low_ess_tail"
    ),
    value = as.numeric(c(
      draws,
      divergences,
      treedepth_hits,
      ebfmi_min,
      max_rhat,
      min_ess_bulk,
      min_ess_tail,
      n_terms_total,
      n_terms_bad_rhat,
      n_terms_low_ess_bulk,
      n_terms_low_ess_tail
    )),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.format_common_diagnostics_for_print <- function(diag_tbl) {
  # Prepare the diagnostics table for console printing.
  #
  # R stores the `value` column as numeric so downstream code can keep using
  # arithmetic on the table. For printing, a subset of metrics are counts and
  # should be displayed as integers rather than floating-point values.
  if (is.null(diag_tbl) || !is.data.frame(diag_tbl) || !all(c("metric", "value") %in% names(diag_tbl))) {
    return(diag_tbl)
  }

  int_metrics <- c(
    "draws",
    "divergences",
    "treedepth_hits",
    "n_terms_total",
    "n_terms_bad_rhat",
    "n_terms_low_ess_bulk",
    "n_terms_low_ess_tail"
  )

  out <- diag_tbl[, c("metric", "value"), drop = FALSE]
  out$value <- vapply(seq_len(nrow(out)), function(i) {
    metric_i <- as.character(out$metric[i])
    value_i <- out$value[i]
    if (!is.finite(value_i)) return(NA_character_)
    if (metric_i %in% int_metrics) {
      return(as.character(as.integer(round(value_i))))
    }
    trimws(format(round(value_i, 3), nsmall = 0, scientific = FALSE))
  }, character(1))
  out
}

#' @keywords internal
.cli_summary_heading <- function(text, level = 1L) {
  level <- as.integer(level %||% 1L)
  captured_output <- sink.number(type = "output") > 0L
  if (!captured_output) {
    if (level <= 1L) {
      cli::cli_h1(text)
    } else if (level == 2L) {
      cli::cli_h2(text)
    } else {
      cli::cli_h3(text)
    }
  } else {
    cat("\n", text, "\n", sep = "")
  }
  invisible(NULL)
}

#' @keywords internal
.cli_print_table <- function(tbl, formatter = identity) {
  if (is.null(tbl)) {
    return(invisible(NULL))
  }
  formatted_tbl <- tibble::as_tibble(formatter(tbl))
  print(formatted_tbl, n = nrow(formatted_tbl), width = Inf)
  invisible(formatted_tbl)
}

#' @keywords internal
.cli_print_bullets <- function(lines) {
  lines <- lines[!is.na(lines) & nzchar(lines)]
  if (!length(lines)) {
    return(invisible(NULL))
  }
  bullet <- cli::symbol$bullet %||% "-"
  cat(paste0(bullet, " ", lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

#' @keywords internal
.cli_print_table_section <- function(title, tbl, level = 2L, formatter = identity) {
  if (is.null(tbl)) {
    return(invisible(NULL))
  }
  .cli_summary_heading(title, level = level)
  .cli_print_table(tbl, formatter = formatter)
  invisible(tbl)
}

#' @keywords internal
.diagnostics_table_from_sampler <- function(diag) {
  .build_common_diagnostics_table(
    draws = as.numeric(diag$draws %||% NA_real_),
    divergences = as.numeric(diag$divergences %||% NA_real_),
    treedepth_hits = as.numeric(diag$treedepth_hits %||% NA_real_),
    ebfmi_min = as.numeric(diag$ebfmi_min %||% NA_real_),
    max_rhat = as.numeric(diag$max_rhat %||% NA_real_),
    min_ess_bulk = as.numeric(diag$min_ess_bulk %||% diag$min_ess %||% NA_real_),
    min_ess_tail = as.numeric(diag$min_ess_tail %||% diag$min_ess %||% NA_real_),
    n_terms_total = NA_real_,
    n_terms_bad_rhat = NA_real_,
    n_terms_low_ess_bulk = NA_real_,
    n_terms_low_ess_tail = NA_real_
  )
}

#' @keywords internal
.term_diagnostics_from_tables <- function(x) {
  extract_df <- function(obj) {
    if (is.null(obj)) return(list())
    if (is.data.frame(obj)) {
      if (any(c("Rhat", "ess_bulk", "ess_tail") %in% names(obj))) {
        return(list(obj))
      }
      return(list())
    }
    if (is.list(obj)) {
      return(unlist(lapply(obj, extract_df), recursive = FALSE))
    }
    list()
  }

  dfs <- extract_df(x)
  if (length(dfs) == 0) {
    return(list(
      n_terms_total = 0,
      n_terms_bad_rhat = 0,
      n_terms_low_ess_bulk = 0,
      n_terms_low_ess_tail = 0,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss) > 0) {
      for (m in miss) df[[m]] <- rep(NA_real_, nrow(df))
    }
    df[, all_cols, drop = FALSE]
  })
  d <- do.call(rbind, dfs)
  if (!"Rhat" %in% names(d)) d$Rhat <- NA_real_
  if (!"ess_bulk" %in% names(d)) d$ess_bulk <- NA_real_
  if (!"ess_tail" %in% names(d)) d$ess_tail <- NA_real_

  n_terms_total <- nrow(d)
  n_terms_bad_rhat <- sum(is.finite(d$Rhat) & d$Rhat > 1.01, na.rm = TRUE)
  n_terms_low_ess_bulk <- sum(is.finite(d$ess_bulk) & d$ess_bulk < 100, na.rm = TRUE)
  n_terms_low_ess_tail <- sum(is.finite(d$ess_tail) & d$ess_tail < 100, na.rm = TRUE)

  max_rhat <- suppressWarnings(max(d$Rhat, na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(d$ess_bulk, na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(d$ess_tail, na.rm = TRUE))
  if (!is.finite(max_rhat)) max_rhat <- NA_real_
  if (!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
  if (!is.finite(min_ess_tail)) min_ess_tail <- NA_real_

  list(
    n_terms_total = as.numeric(n_terms_total),
    n_terms_bad_rhat = as.numeric(n_terms_bad_rhat),
    n_terms_low_ess_bulk = as.numeric(n_terms_low_ess_bulk),
    n_terms_low_ess_tail = as.numeric(n_terms_low_ess_tail),
    max_rhat = as.numeric(max_rhat),
    min_ess_bulk = as.numeric(min_ess_bulk),
    min_ess_tail = as.numeric(min_ess_tail)
  )
}

#' @keywords internal
.aggregate_sampler_diagnostics <- function(diags) {
  if (is.null(diags) || length(diags) == 0) {
    return(list(
      draws = NA_real_,
      divergences = 0,
      treedepth_hits = 0,
      ebfmi_min = NA_real_,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_
    ))
  }

  get_num <- function(name) as.numeric(vapply(diags, function(d) d[[name]] %||% NA_real_, numeric(1)))
  draws <- get_num("draws")
  divergences <- get_num("divergences")
  treedepth_hits <- get_num("treedepth_hits")
  ebfmi_min <- get_num("ebfmi_min")
  max_rhat <- get_num("max_rhat")
  min_ess_bulk <- get_num("min_ess_bulk")
  min_ess_tail <- get_num("min_ess_tail")

  out <- list(
    draws = if (any(is.finite(draws))) max(draws, na.rm = TRUE) else NA_real_,
    divergences = if (any(is.finite(divergences))) sum(divergences[is.finite(divergences)]) else 0,
    treedepth_hits = if (any(is.finite(treedepth_hits))) sum(treedepth_hits[is.finite(treedepth_hits)]) else 0,
    ebfmi_min = if (any(is.finite(ebfmi_min))) min(ebfmi_min, na.rm = TRUE) else NA_real_,
    max_rhat = if (any(is.finite(max_rhat))) max(max_rhat, na.rm = TRUE) else NA_real_,
    min_ess_bulk = if (any(is.finite(min_ess_bulk))) min(min_ess_bulk, na.rm = TRUE) else NA_real_,
    min_ess_tail = if (any(is.finite(min_ess_tail))) min(min_ess_tail, na.rm = TRUE) else NA_real_
  )
  out
}

# ---- print/summary --------------------------------------------------------

#' Print a JoiNMe object
#'
#' @param x A JoiNMe fit object.
#' @param ... Unused.
#'
#' @return Invisibly returns the object.
#' @export
print.JoiNMeFit <- function(x, ...) {
  # Delegate to R6 print.
  NextMethod('print')
}


#' Print a summary_JoiNMeFit object
#'
#' @param x A summary_JoiNMeFit object.
#' @param ... Unused.
#'
#' @return Invisibly returns the summary object.
#' @export
print.summary_JoiNMeFit <- function(x, ...) {
  longitudinal_only_fit <- identical(
    x$metadata$event_process,
    "not fitted"
  ) # whether the summary contains only the nested longitudinal process
  .cli_summary_heading(
    if (longitudinal_only_fit) {
      "Nested longitudinal mixed effects model summary"
    } else {
      "Joint mixed effects model summary"
    },
    level = 1L
  )
  meta_lines <- character(0)
  if (!is.null(x$metadata$call) && nzchar(x$metadata$call)) {
    meta_lines <- c(meta_lines, paste0("Call: ", x$metadata$call))
  }
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

    meta_lines <- c(meta_lines, paste0("Family: ", family))
  }
  if (!is.null(x$metadata$event_process)) {
    meta_lines <- c(
      meta_lines,
      paste0("Event process: ", x$metadata$event_process)
    )
  }
  if (!is.null(x$metadata$basehaz)) {
    meta_lines <- c(meta_lines, paste0("Baseline hazard type: ", x$metadata$basehaz))
  }
  if (!is.null(x$metadata$tmax)) {
    meta_lines <- c(meta_lines, paste0("tmax: ", x$metadata$tmax))
  }
  marker_weight_offsets <- x$metadata$marker_weight_offsets %||% list() # declared marker-specific constants displayed before posterior tables
  if (length(marker_weight_offsets) > 0L) {
    for (set_name in names(marker_weight_offsets)) {
      offset_values <- marker_weight_offsets[[set_name]] # named marker offsets for one shared or term-specific fitted weight set
      offset_text <- paste0(names(offset_values), "=", format(signif(as.numeric(offset_values), 6), trim = TRUE), collapse = ", ") # compact marker-by-offset declaration
      meta_lines <- c(meta_lines, paste0("Marker-weight offsets [", set_name, "]: ", offset_text))
    }
  }
  .cli_print_bullets(meta_lines)
  diag_tbl <- x$tables$diagnostics
  if (is.null(diag_tbl) && !is.null(x$diagnostics)) {
    diag_tbl <- .diagnostics_table_from_sampler(x$diagnostics)
  }
  .cli_print_table_section("Sampler diagnostics", diag_tbl, level = 2L, formatter = .format_common_diagnostics_for_print)
  .cli_print_table_section("Transformations", x$metadata$transform_formulas, level = 2L)
  .cli_print_table_section("Fixed effects (beta)", x$tables$fixef, level = 2L)
  .cli_print_table_section("Baseline hazard coefficients", x$tables$baseline_hazard, level = 2L)
  .cli_print_table_section("Survival process (non-association covariates)", x$tables$survival_process, level = 2L)
  .cli_print_table_section("Association parameters", x$tables$assoc, level = 2L)
  .cli_print_table_section("Marker-weight distribution", x$tables$marker_weights, level = 2L)
  .cli_print_table_section("Transform parameters", x$tables$transform_parameters, level = 2L)
  .cli_print_table_section("Piecewise-linear relative log-hazard ordinates", x$tables$piecewise_ordinates, level = 2L)
  .cli_print_table_section("Distributional parameters", x$tables$distributional, level = 2L)
  .cli_print_table_section("Distributional regression", x$tables$distributional_regression, level = 2L)
  if (!is.null(x$tables$corr)) {
    .cli_summary_heading("Covariance summaries (diagonal entries are variances)", level = 2L)
    if (!is.null(x$tables$corr$id)) {
      .cli_summary_heading("id", level = 3L)
      .cli_print_table(x$tables$corr$id)
    }
    if (!is.null(x$tables$corr$marker)) {
      .cli_summary_heading("marker", level = 3L)
      .cli_print_table(x$tables$corr$marker)
    }
  }
  if (!is.null(x$tables$id_marker_cov)) {
    .cli_summary_heading("id:marker covariance parameters", level = 2L)
    if (!is.null(x$tables$id_marker_cov$regression)) {
      .cli_print_table_section("covariance regression coefficients", x$tables$id_marker_cov$regression, level = 3L)
    }
    if (!is.null(x$tables$id_marker_cov$hyperparameters)) {
      .cli_print_table_section("covariance regression hyperparameters", x$tables$id_marker_cov$hyperparameters, level = 3L)
    }
  } else if (!is.null(x$tables$corr) && !is.null(x$tables$corr$id_marker)) {
    .cli_print_table_section("id:marker covariance parameters", x$tables$corr$id_marker, level = 2L)
  }
  invisible(x)
}


#' Print dynamic prediction results
#'
#' @param x A dynamic prediction object.
#' @param ... Unused.
#'
#' @return Invisibly returns the object.
#' @export
print.JoiNMeDynPred <- function(x, ...) {
    
    .cli_summary_heading("Joint mixed effects dynamic prediction", level = 1L)
    if (!is.null(x$call)) {
        cat("Call:\n")
        print(x$call)
    }
    if (!is.null(x$metadata$pred_type)) {
        cat("Prediction type: ", x$metadata$pred_type, "\n", sep = "")
    }
    scales <- x$metadata$scales %||% x$metadata$scale
    if (!is.null(scales)) {
        if (length(scales) > 1L) {
            cat("Scales: ", paste(scales, collapse = ", "), "\n", sep = "")
        } else {
            cat("Scale: ", scales, "\n", sep = "")
        }
    }
    if (!is.null(x$metadata$n_subjects)) {
        cat("Subjects: ", x$metadata$n_subjects, "\n", sep = "")
    }
    if (!is.null(x$n_samples)) {
        cat("Posterior draws: ", x$n_samples, "\n", sep = "")
    }
    cat("Use summary() for prediction summaries.\n")
    cat("Use plot() for trajectory and interval visualisation.\n")
    invisible(x)
}

#' Summarise renamed posterior draws
#'
#' @param draws_obj Posterior draw matrix or array with renamed variables.
#' @param digits Number of decimal places used for posterior summaries.
#'
#' @return A summary table with the common JoiNMe schema.
#' @keywords internal
#' @noRd
.summarise_named_draws <- function(draws_obj, digits = 3) {
  draw_array <- .make_named_draw_array(draws_obj)
  if (is.null(draw_array) || !length(dim(draw_array)) || dim(draw_array)[3] == 0L) {
    return(NULL)
  }
  term_labels <- dimnames(draw_array)[[3]] %||% paste0("term_", seq_len(dim(draw_array)[3]))
  .assoc_summary_from_draw_array(draw_array, term_labels = term_labels, digits = digits)
}

#' Round summary tables using the standard JoiNMe summary schema
#'
#' @description
#' Applies the same rounding rules used throughout JoiNMe summaries
#'
#' The rounded columns are the inferential columns that users typically inspect
#' first: posterior mean, posterior standard deviation, interval bounds, and
#' convergence diagnostics. Columns not listed in the schema are left unchanged.
#'
#' @param tbl A data frame built from posterior draws.
#' @param digits Number of decimal places used for posterior location and
#'   interval summaries.
#'
#' @return The same data frame with rounded summary columns.
#' @keywords internal
#' @noRd
.round_summary_table <- function(tbl, digits = 3) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(tbl)
  }

  rounded <- tbl
  for (col_name in intersect(
    c(
      "Estimate",
      "Mean",
      "Median",
      "Est.Error",
      "SD",
      "Q2.5",
      "Q97.5"
    ),
    names(rounded)
  )) {
    rounded[[col_name]] <- round(rounded[[col_name]], digits)
  }
  if ("Rhat" %in% names(rounded)) {
    rounded$Rhat <- round(rounded$Rhat, 3)
  }
  rounded
}

#' Convert renamed draws into a chains-aware array
#'
#' @param x Posterior draws returned by [extract()] or [posterior_draws()].
#'
#' @return A three-dimensional array with dimensions iteration x chain x term.
#' @keywords internal
#' @noRd
.make_named_draw_array <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (length(dim(x)) == 3L) {
    return(x)
  }
  if (is.matrix(x)) {
    return(array(
      x,
      dim = c(nrow(x), 1L, ncol(x)),
      dimnames = list(
        iteration = rownames(x) %||% as.character(seq_len(nrow(x))),
        chain = "1",
        variable = colnames(x) %||% as.character(seq_len(ncol(x)))
      )
    ))
  }
  NULL
}

#' Summarise extracted draws for one JoiNMeFit component
#'
#' @param object A `JoiNMeFit` object.
#' @param what Component selector passed to [extract.JoiNMeFit()].
#' @param draws Optional number of posterior draws to keep.
#' @param seed Integer seed used when subsetting draws.
#' @param digits Number of decimal places used for posterior summaries.
#' @param term Optional term filter passed through to [extract.JoiNMeFit()].
#'
#' @return A posterior summary table, or `NULL` when no draws match.
#' @keywords internal
#' @noRd
.extract_fit_summary <- function(object, what, draws = NULL, seed = 1, digits = 3, term = NULL) {
  ext <- tryCatch(
    extract.JoiNMeFit(object, what = what, term = term, draws = draws, seed = seed, keep_chains = TRUE),
    error = function(e) NULL
  )
  if (is.null(ext) || is.null(ext$posterior_draws)) {
    return(NULL)
  }
  .summarise_named_draws(ext$posterior_draws, digits = digits)
}

#' Convert integer covariance indices into model-term labels
#'
#' @description
#' The low-level covariance extractors return matrix coordinates as integer row
#' and column positions. For console summaries, those positions are harder to
#' interpret than the corresponding random-effect terms. This helper replaces the
#' raw indices with term labels recovered from the same model matrix basis.
#'
#' @param tbl Covariance summary table containing `row` and `col` columns.
#' @param term_labels Character vector of basis labels.
#'
#' @return The same data frame with labelled `row` and `col` columns.
#' @keywords internal
#' @noRd
.label_covariance_summary_table <- function(tbl, term_labels) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(tbl)
  }
  if (!all(c("row", "col") %in% names(tbl))) {
    return(tbl)
  }

  term_labels <- as.character(term_labels %||% character(0))
  if (!length(term_labels)) {
    return(tbl)
  }

  out <- tbl
  row_index <- suppressWarnings(as.integer(out$row))
  col_index <- suppressWarnings(as.integer(out$col))
  row_default <- as.character(out$row)
  col_default <- as.character(out$col)

  valid_row <- !is.na(row_index) & row_index >= 1L & row_index <= length(term_labels)
  valid_col <- !is.na(col_index) & col_index >= 1L & col_index <= length(term_labels)

  row_default[valid_row] <- term_labels[row_index[valid_row]]
  col_default[valid_col] <- term_labels[col_index[valid_col]]

  out$row <- row_default
  out$col <- col_default
  out
}
