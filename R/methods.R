#' @name joinme_methods
#' @title joinme S3 Methods
#'
#' @description
#' S3 methods for joinme R6 objects with cached summaries
#' @keywords internal
#' @author Trinh Dong
NULL

# File overview:
# - S3 methods for JoinMeFit and summary output.
# - Internal helpers for diagnostics and draw summaries.

#' Extract covariance summaries
#'
#' Generic for extracting covariance summaries from joinme objects.
#'
#' @param object A supported joinme object.
#' @param ... Additional arguments passed to methods.
#' @export
corr <- function(object, ...) {
  UseMethod("corr")
}

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
.format_transform_spec <- function(spec) {
  if (is.null(spec) || is.null(spec$type) || spec$type == "identity") {
    return("identity")
  }
  if (spec$type == "functional") {
    return(.format_transform_expr(spec$expr))
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
    x_vals <- spec$x
    y_vals <- spec$y
    x_text <- if (!is.null(x_vals)) paste(format(x_vals, digits = 3, trim = TRUE), collapse = ", ") else ""
    y_text <- if (!is.null(y_vals)) paste(format(y_vals, digits = 3, trim = TRUE), collapse = ", ") else ""
    return(paste0("pwlin(x = c(", x_text, "), y = c(", y_text, "))"))
  }
  spec$type
}

#' @keywords internal
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
  # Key contract:
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
  P_phi_beta <- as.integer(sd$P_phi_beta %||% 0L)
  P_tau_sde <- as.integer(sd$P_tau_sde %||% 0L)
  n_re_sigma <- as.integer(sd$n_re_sigma %||% 0L)
  n_re_nu <- as.integer(sd$n_re_nu %||% 0L)
  n_re_phi <- as.integer(sd$n_re_phi %||% 0L)
  n_re_alpha <- as.integer(sd$n_re_alpha %||% 0L)
  n_re_phi_beta <- as.integer(sd$n_re_phi_beta %||% 0L)
  n_re_tau_sde <- as.integer(sd$n_re_tau_sde %||% 0L)

  has_reg <- list(
    sigma = (P_sigma > 0 || n_re_sigma > 0),
    nu = (P_nu > 0 || n_re_nu > 0),
    phi = (P_phi > 0 || n_re_phi > 0),
    alpha = (P_alpha > 0 || n_re_alpha > 0),
    phi_beta = (P_phi_beta > 0 || n_re_phi_beta > 0),
    tau_sde = (P_tau_sde > 0 || n_re_tau_sde > 0)
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
    phi_beta = sort(unique(family_long[vapply(req_by_marker, function(x) "phi_beta" %in% x, logical(1))])),
    tau_sde = sort(unique(family_long[vapply(req_by_marker, function(x) "tau_sde" %in% x, logical(1))]))
  )

  marker_label_prefix <- list(
    sigma = "sigma_marker",
    nu = "nu_marker",
    phi = "phi_nb_marker",
    alpha = "alpha_skew_marker",
    phi_beta = "phi_beta_marker",
    tau_sde = "tau_sde_marker"
  )
  family_label_prefix <- list(
    sigma = "sigma",
    nu = "nu",
    phi = "phi_nb",
    alpha = "alpha_skew",
    phi_beta = "phi_beta",
    tau_sde = "tau_sde"
  )

  markers_for_param <- function(param, fam_code) {
    which(family_long == fam_code & vapply(req_by_marker, function(x) param %in% x, logical(1)))
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
  add_family_terms("phi_beta", "phi_beta_family", fam_sets$phi_beta, has_reg$phi_beta)
  add_family_terms("tau_sde", "tau_sde_family", fam_sets$tau_sde, has_reg$tau_sde)

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
    if ("phi_beta" %in% req && !isTRUE(has_reg$phi_beta)) {
      var_nm <- paste0("phi_beta_marker[", d, "]")
      term_map[var_nm] <- paste0("phi_beta_marker[", mk, "]")
    }
    if ("tau_sde" %in% req && !isTRUE(has_reg$tau_sde)) {
      var_nm <- paste0("tau_sde_marker[", d, "]")
      term_map[var_nm] <- paste0("tau_sde_marker[", mk, "]")
    }
  }

  # Keep only variables that are present in posterior draws.
  term_map[names(term_map) %in% all_vars]
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
      for (m in miss) df[[m]] <- NA_real_
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

#' Print a joinme object
#'
#' @param x A joinme fit object.
#' @param ... Unused.
#'
#' @return Invisibly returns the object.
#' @export
print.JoinMeFit <- function(x, ...) {
  # User-facing summary header for fits
  assertthat::assert_that(inherits(x, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  cat("JoinMe model fit\n")
  cat("===============\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  family <- x$stan_data$family_names %||% .family_code_to_name(x$config$family_long)
  if (!is.null(family)) {
    fml <-
      if (length(unique(family)) == 1) family[1] else
      if (length(family) > 5) paste(paste(head(family, 5), collapse = ", "), "...") else
      paste(family, collapse = ", ")
    cat("Family: ", fml, "\n", sep = "")
  }
  if (!is.null(x$tmax)) {
    cat("tmax: ", x$tmax, "\n", sep = "")
  }
  cat("Use summary() for parameter summaries.\n")
  invisible(x)
}

#' Summarise a joinme object
#'
#' @description
#' Builds posterior summaries for the longitudinal process, survival process,
#' association terms, and optional covariance blocks. For survival, the summary
#' now includes a dedicated `survival_process` report whenever the event model
#' contains non-intercept covariates beyond association features.
#'
#' @param object A joinme fit object.
#' @param draws Number of draws to use for summaries.
#' @param seed Random seed for subsetting draws.
#' @param digits Number of digits to round summary values.
#' @param include_corr Logical; include covariance summaries.
#' @param ... Unused.
#'
#' @return A `summary_JoinMeFit` object containing tables such as `fixef`,
#'   `survival_process` (when applicable), `assoc`, covariance
#'   summaries (`id`, `marker`), dedicated `id:marker` covariance-parameter
#'   summaries (latent + covariance-regression blocks when `Q_idm > 0`), and diagnostics.
#' @export
summary.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3,
                           include_corr = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config

  if (is.null(draws)) draws <- cfg$draws_default
  assertthat::assert_that(is.null(draws) || (is.numeric(draws) && draws > 0),
                          msg = "draws must be NULL or a positive number.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- paste0("summary_", draws, "_", digits, "_", include_corr)
  # Cache by draw count + digits to avoid repeat summaries
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) return(cached)

  diag <- .joinme_sampler_diagnostics(fit)
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))

  beta_vars <- paste0("beta[", seq_len(sd$P), "]")
  s_beta <- as.data.frame(.summarise_draws_diag(fit, beta_vars, draws = draws, seed = seed))
  s_beta$term <- if (!is.null(sd$x_cols) && length(sd$x_cols) == nrow(s_beta)) sd$x_cols else s_beta$variable
  s_beta <- s_beta[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
  s_beta$Estimate <- round(s_beta$Estimate, digits)
  s_beta$Est.Error <- round(s_beta$Est.Error, digits)
  s_beta$Q2.5 <- round(s_beta$Q2.5, digits)
  s_beta$Q97.5 <- round(s_beta$Q97.5, digits)
  s_beta$Rhat <- round(s_beta$Rhat, 3)

  s_g <- NULL
  if (sd$p_w > 0) {
    g_vars <- paste0("gamma_w[", seq_len(sd$p_w), "]")
    g_vars <- g_vars[g_vars %in% all_vars]

    if (length(g_vars) == 0) {
      k_event <- sd$K_event %||% 1L
      g_vars_mat <- outer(
        seq_len(k_event),
        seq_len(sd$p_w),
        function(k, j) paste0("gamma_w[", k, ",", j, "]")
      )
      g_vars <- as.vector(g_vars_mat)
      g_vars <- g_vars[g_vars %in% all_vars]
    }

    if (length(g_vars) > 0) {
      s_g <- as.data.frame(.summarise_draws_diag(fit, g_vars, draws = draws, seed = seed))
      if (!is.null(sd$w_cols)) {
        var_idx <- regmatches(s_g$variable, regexec("^gamma_w\\[(\\d+)(?:,(\\d+))?\\]$", s_g$variable))
        k_idx <- vapply(var_idx, function(x) if (length(x) >= 2) as.integer(x[2]) else NA_integer_, integer(1))
        j_idx <- vapply(var_idx, function(x) if (length(x) >= 3 && nzchar(x[3])) as.integer(x[3]) else NA_integer_, integer(1))
        k_event <- sd$K_event %||% 1L
        term_base <- s_g$variable
        if (any(!is.na(j_idx))) term_base[!is.na(j_idx)] <- sd$w_cols[j_idx[!is.na(j_idx)]]
        if (k_event > 1 && !all(is.na(k_idx))) {
          s_g$term <- paste0("event", k_idx, ": ", term_base)
        } else {
          s_g$term <- term_base
        }
      } else {
        s_g$term <- s_g$variable
      }
      s_g <- s_g[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      s_g$Estimate <- round(s_g$Estimate, digits)
      s_g$Est.Error <- round(s_g$Est.Error, digits)
      s_g$Q2.5 <- round(s_g$Q2.5, digits)
      s_g$Q97.5 <- round(s_g$Q97.5, digits)
      s_g$Rhat <- round(s_g$Rhat, 3)
    }
  }

  # Survival-process report (non-association baseline covariates only)
  #
  # IMPORTANT:
  # - If and only if baseline event covariates exist in the survival model,
  #   `survival_process` must be reported.
  # - These terms must match the event design columns (e.g., x1, x2) and must
  #   not be silently dropped.
  #
  # Algorithm:
  # 1) Start from gamma_w summaries (the survival linear predictor terms).
  # 2) Re-label terms from standata event columns when available (`w_cols`).
  # 3) Keep only non-intercept baseline covariates.
  # 4) Add hazard-ratio summaries exp(beta) for direct interpretation.
  s_surv <- NULL
  if (!is.null(s_g) && nrow(s_g) > 0 && isTRUE((sd$p_w %||% 0L) > 0L)) {
    is_intercept_term <- function(term) {
      term_chr <- trimws(as.character(term))
      tolower(term_chr) %in% c("(intercept)", "intercept", "1")
    }

    baseline_terms <- if (!is.null(sd$w_cols) && length(sd$w_cols) == nrow(s_g)) {
      as.character(sd$w_cols)
    } else {
      as.character(s_g$term)
    }

    keep_idx <- !vapply(baseline_terms, is_intercept_term, logical(1))
    if (!any(keep_idx) && length(baseline_terms) > 0) {
      # Defensive fallback: if intercept filtering drops everything, keep all
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

  a_vars <- c("alpha_cv_total", "alpha_cv_mean", "alpha_cv_marker", "alpha_cs_total", "alpha_cs_mean", "alpha_cs_marker")
  a_vars <- c(a_vars, grep("^alpha_corr\\[", all_vars, value = TRUE))
  a_vars <- c(a_vars, grep("^alpha_vcov\\[", all_vars, value = TRUE))
  a_vars <- a_vars[a_vars %in% all_vars]
  active_vars <- c(
    if (isTRUE(sd$assoc_cv_total == 1)) "alpha_cv_total",
    if (isTRUE(sd$assoc_cv_mean == 1)) "alpha_cv_mean",
    if (isTRUE(sd$assoc_cv_marker == 1)) "alpha_cv_marker",
    if (isTRUE(sd$assoc_cs_total == 1)) "alpha_cs_total",
    if (isTRUE(sd$assoc_cs_mean == 1)) "alpha_cs_mean",
    if (isTRUE(sd$assoc_cs_marker == 1)) "alpha_cs_marker",
    if (isTRUE(sd$assoc_corr == 1)) grep("^alpha_corr\\[", a_vars, value = TRUE),
    if (isTRUE(sd$assoc_vcov == 1)) grep("^alpha_vcov\\[", a_vars, value = TRUE)
  )
  if (length(active_vars) > 0) {
    a_vars <- a_vars[a_vars %in% active_vars]
  }
  # Association coefficient summaries (respect assoc flags)
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
    s_a <- s_a[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    s_a$Estimate <- round(s_a$Estimate, digits)
    s_a$Est.Error <- round(s_a$Est.Error, digits)
    s_a$Q2.5 <- round(s_a$Q2.5, digits)
    s_a$Q97.5 <- round(s_a$Q97.5, digits)
    s_a$Rhat <- round(s_a$Rhat, 3)
  }

  # Marker-weight association summaries are shown only when marker-weighted
  # association terms are active (cv_total/cv_marker/cs_total/cs_marker).
  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)

  if (sd$D > 0 && show_marker_weights) {
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
      s_mw <- s_mw[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      s_mw$Estimate <- round(s_mw$Estimate, digits)
      s_mw$Est.Error <- round(s_mw$Est.Error, digits)
      s_mw$Q2.5 <- round(s_mw$Q2.5, digits)
      s_mw$Q97.5 <- round(s_mw$Q97.5, digits)
      s_mw$Rhat <- round(s_mw$Rhat, 3)
      s_a <- if (is.null(s_a)) s_mw else rbind(s_a, s_mw)
    } else if (isTRUE(as.logical(sd$fixed_marker_weights)) && !is.null(sd$marker_weights)) {
      marker_terms <- sd$marker_levels %||% paste0("marker_", seq_len(sd$D))
      s_mw <- data.frame(
        term = paste0("weight: ", marker_terms),
        Estimate = round(as.numeric(sd$marker_weights), digits),
        Est.Error = NA_real_,
        Q2.5 = NA_real_,
        Q97.5 = NA_real_,
        Rhat = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        stringsAsFactors = FALSE
      )
      s_a <- if (is.null(s_a)) s_mw else rbind(s_a, s_mw)
    }
  }

  # Distributional parameter summaries (family-aware, marker-labeled)
  #
  # Reporting rules:
  # - Include only parameters required by each marker family.
  # - Replace numeric marker indices with marker names in term labels.
  dist_term_map <- .distributional_term_map(sd, cfg, all_vars)
  dist_vars <- names(dist_term_map)
  s_d <- NULL
  if (length(dist_vars) > 0) {
    s_d <- as.data.frame(.summarise_draws_diag(fit, dist_vars, draws = draws, seed = seed))
    s_d$term <- unname(dist_term_map[s_d$variable])
    s_d <- s_d[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    s_d$Estimate <- round(s_d$Estimate, digits)
    s_d$Est.Error <- round(s_d$Est.Error, digits)
    s_d$Q2.5 <- round(s_d$Q2.5, digits)
    s_d$Q97.5 <- round(s_d$Q97.5, digits)
    s_d$Rhat <- round(s_d$Rhat, 3)
  }

  s_dr <- NULL
  dist_cols <- cfg$dist$dist_cols %||% list()
  # Distributional regression summaries (fixed effects)
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
      tmp <- tmp[, c("parameter", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      tmp$Estimate <- round(tmp$Estimate, digits)
      tmp$Est.Error <- round(tmp$Est.Error, digits)
      tmp$Q2.5 <- round(tmp$Q2.5, digits)
      tmp$Q97.5 <- round(tmp$Q97.5, digits)
      tmp$Rhat <- round(tmp$Rhat, 3)
      dist_reg_tables[[nm]] <- tmp
    }
  }
  if (length(dist_reg_tables) > 0) {
    s_dr <- do.call(rbind, dist_reg_tables)
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

      summarize_latent_cov_matrix <- function() {
        tau_vars <- paste0("tau_w[", seq_len(q_idm), "]")
        L_vars <- as.vector(outer(seq_len(q_idm), seq_len(q_idm), function(r, c) paste0("Lcorr_w[", r, ",", c, "]")))
        vars_needed <- c(tau_vars, L_vars)
        vars_needed <- vars_needed[vars_needed %in% all_vars]
        if (length(tau_vars[tau_vars %in% all_vars]) != q_idm) return(NULL)
        if (length(L_vars[L_vars %in% all_vars]) != q_idm * q_idm) return(NULL)
        dmat <- .get_draws_matrix(fit, variables = c(tau_vars, L_vars), draws = draws, seed = seed)
        out <- list()
        idx <- 1L
        for (r in seq_len(q_idm)) {
          for (c in seq_len(q_idm)) {
            vals <- vapply(seq_len(nrow(dmat)), function(i) {
              tau <- as.numeric(dmat[i, tau_vars])
              L <- matrix(as.numeric(dmat[i, L_vars]), nrow = q_idm, ncol = q_idm, byrow = FALSE)
              L[upper.tri(L)] <- 0
              Corr <- L %*% t(L)
              Sigma <- diag(tau, q_idm, q_idm) %*% Corr %*% diag(tau, q_idm, q_idm)
              as.numeric(Sigma[r, c])
            }, numeric(1))
            ss <- .summarize_draw_col(vals)
            rhat <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(vals)), error = function(e) NA_real_))
            ess_bulk <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(vals)), error = function(e) NA_real_))
            ess_tail <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(vals)), error = function(e) NA_real_))
            out[[idx]] <- data.frame(
              block = "sigma_latent",
              row = r,
              col = c,
              Estimate = round(as.numeric(ss[["Estimate"]]), digits),
              Est.Error = round(as.numeric(ss[["Est.Error"]]), digits),
              Q2.5 = round(as.numeric(ss[["Q2.5"]]), digits),
              Q97.5 = round(as.numeric(ss[["Q97.5"]]), digits),
              Rhat = round(rhat, 3),
              ess_bulk = as.numeric(ess_bulk),
              ess_tail = as.numeric(ess_tail),
              stringsAsFactors = FALSE
            )
            idx <- idx + 1L
          }
        }
        if (length(out) == 0) return(NULL)
        do.call(rbind, out)
      }

      latent_tbl <- summarize_latent_cov_matrix()

      alpha_vars <- paste0("alpha_L[", seq_len(m_cov), "]")
      alpha_blocks <- rep("L[id:marker]", m_cov)
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
        beta_blocks <- rep("L[id:marker]", length(beta_vars))
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
      lambda_blocks <- rep("L[id:marker]", m_cov)
      lambda_tbl <- .summarize_block_parameters(
        lambda_vars,
        block_labels = lambda_blocks,
        term_labels = rep("lambda", m_cov),
        row_labels = rc_map[, 1],
        col_labels = rc_map[, 2]
      )
      if (!is.null(lambda_tbl)) reg_rows[[length(reg_rows) + 1L]] <- lambda_tbl

      hyperparameter_tbl <- .summarize_block_parameters("tau_L", "id:marker", "sd_u")
      if (!is.null(hyperparameter_tbl)) {
        hyperparameter_tbl <- hyperparameter_tbl[, c("block", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      }

      reg_rows <- Filter(Negate(is.null), reg_rows)
      regression_tbl <- if (length(reg_rows) > 0) do.call(rbind, reg_rows) else NULL
      if (!is.null(regression_tbl)) {
        regression_tbl <- regression_tbl[, c("block", "row", "col", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      }
      if (isTRUE(any_re_indep)) {
        latent_tbl <- .filter_diag_rows(latent_tbl)
        regression_tbl <- NULL
      }
      id_marker_cov_tables <- list(
        latent = latent_tbl,
        regression = regression_tbl,
        hyperparameters = hyperparameter_tbl
      )
      id_marker_cov_tables <- id_marker_cov_tables[!vapply(id_marker_cov_tables, is.null, logical(1))]
      if (length(id_marker_cov_tables) == 0) id_marker_cov_tables <- NULL
    }

    corr_tables <- list(
      id = corr(object, what = "id", draws = draws),
      marker = if (sd$R_mk > 0) corr(object, what = "marker", draws = draws) else NULL
    )
    if (isTRUE(any_re_indep)) {
      corr_tables <- lapply(corr_tables, .filter_diag_rows)
    }
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
    transform_param_tables[[channel]] <- dplyr::arrange(tmp, channel, term)
  }
  transform_params <- if (length(transform_param_tables) > 0) {
    do.call(rbind, unname(transform_param_tables))
  } else {
    NULL
  }

  transform_specs <- cfg$transforms_spec %||% object$call$transforms
  transform_formulas <- .transform_formulas_from_specs(transform_specs, sd = sd)

  term_diag <- .term_diagnostics_from_tables(list(s_beta, s_surv, s_a, transform_params, s_d, s_dr, corr_tables, id_marker_cov_tables))
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

  summary_obj <- SummaryJoinMeFit$new(
    tables = list(
      diagnostics = diag_table,
      fixef = s_beta,
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
      tmax = sd$tmax %||% cfg$tmax %||% 1.0,
      draws = draws,
      transforms = cfg$transforms,
      transform_formulas = transform_formulas
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
#' @param formulaCorr Optional updated covariance formula (full or update form).
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
#' @method update JoinMeFit
#' @export
update.JoinMeFit <- function(
  object,
  formulaLong = NULL,
  dataLong = NULL,
  formulaEvent = NULL,
  dataEvent = NULL,
  formulaCorr = NULL,
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
    # Support update formulas ("~ . + x") while preserving original
    if (is.null(updated)) return(current)
    if (!inherits(updated, "formula")) {
      cli::cli_abort("{.arg {name}} must be a formula.")
    }

    # Robust dot detection without regex; the dot symbol is not reported by all.vars().
    expr_has_dot <- function(expr) {
      if (is.null(expr)) return(FALSE)
      if (is.name(expr) && identical(as.character(expr), ".")) return(TRUE)
      if (!is.call(expr)) return(FALSE)
      any(vapply(as.list(expr)[-1], expr_has_dot, logical(1)))
    }

    replace_dot <- function(expr, replacement) {
      if (is.null(expr)) return(expr)
      if (is.name(expr) && identical(as.character(expr), ".")) return(replacement)
      if (!is.call(expr)) return(expr)
      as.call(c(expr[[1]], lapply(as.list(expr)[-1], replace_dot, replacement = replacement)))
    }

    has_lhs <- length(updated) == 3
    rhs <- if (has_lhs) updated[[3]] else updated[[2]]
    if (expr_has_dot(rhs) || !has_lhs) {
      current_exp <- reformulas::expandDoubleVerts(current)
      current_bars <- reformulas::findbars(current_exp)
      current_fix <- reformulas::nobars(current_exp)
      current_fix_rhs <- if (length(current_fix) == 3) current_fix[[3]] else current_fix[[2]]
      new_rhs <- if (expr_has_dot(rhs)) replace_dot(rhs, current_fix_rhs) else rhs

      tmp_form <- rlang::new_formula(NULL, new_rhs, env = rlang::`%||%`(environment(current), environment(updated)))
      has_bars <- length(reformulas::findbars(reformulas::expandDoubleVerts(tmp_form))) > 0
      if (!has_bars && length(current_bars) > 0) {
        new_rhs <- Reduce(function(acc, bar) call("+", acc, bar), current_bars, init = new_rhs)
      }

      lhs <- if (has_lhs) updated[[2]] else if (length(current) == 3) current[[2]] else NULL
      return(rlang::new_formula(lhs, new_rhs, env = rlang::`%||%`(environment(current), environment(updated))))
    }
    updated
  }

  call_obj$formulaLong <- update_formula(object$formulaLong, formulaLong, "formulaLong")
  call_obj$formulaEvent <- update_formula(object$formulaEvent, formulaEvent, "formulaEvent")
  call_obj$formulaCorr <- update_formula(object$formulaCorr, formulaCorr, "formulaCorr")
  if (!is.null(formulaDist)) call_obj$formulaDist <- formulaDist

  if (!is.null(call_obj$formulaLong) && inherits(call_obj$formulaLong, "formula")) {
    f_exp <- reformulas::expandDoubleVerts(call_obj$formulaLong)
    if (length(reformulas::findbars(f_exp)) == 0 && !is.null(object$formulaLong)) {
      current_exp <- reformulas::expandDoubleVerts(object$formulaLong)
      current_bars <- reformulas::findbars(current_exp)
      if (length(current_bars) > 0) {
        fe_form <- reformulas::nobars(f_exp)
        fe_rhs <- if (length(fe_form) == 3) fe_form[[3]] else fe_form[[2]]
        rhs <- Reduce(function(acc, bar) call("+", acc, bar), current_bars, init = fe_rhs)
        lhs <- if (length(call_obj$formulaLong) == 3) call_obj$formulaLong[[2]] else object$formulaLong[[2]]
        call_obj$formulaLong <- rlang::new_formula(lhs, rhs, env = rlang::`%||%`(environment(call_obj$formulaLong), environment(object$formulaLong)))
      }
    }
  }

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
  .cli_summary_heading("Joint mixed effects model summary", level = 1L)
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
  if (!is.null(x$metadata$tmax)) {
    meta_lines <- c(meta_lines, paste0("tmax: ", x$metadata$tmax))
  }
  .cli_print_bullets(meta_lines)
  diag_tbl <- x$tables$diagnostics
  if (is.null(diag_tbl) && !is.null(x$diagnostics)) {
    diag_tbl <- .diagnostics_table_from_sampler(x$diagnostics)
  }
  .cli_print_table_section("Sampler diagnostics", diag_tbl, level = 2L, formatter = .format_common_diagnostics_for_print)
  .cli_print_table_section("Transformations", x$metadata$transform_formulas, level = 2L)
  .cli_print_table_section("Fixed effects (beta)", x$tables$fixef, level = 2L)
  .cli_print_table_section("Survival process (non-association covariates)", x$tables$survival_process, level = 2L)
  .cli_print_table_section("Association parameters", x$tables$assoc, level = 2L)
  .cli_print_table_section("Transform parameters", x$tables$transform_parameters, level = 2L)
  .cli_print_table_section("Distributional parameters", x$tables$distributional, level = 2L)
  .cli_print_table_section("Distributional regression", x$tables$distributional_regression, level = 2L)
  if (!is.null(x$tables$corr)) {
    .cli_summary_heading("Covariance summaries", level = 2L)
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
    if (!is.null(x$tables$id_marker_cov$latent)) {
      .cli_print_table_section("latent covariance matrix", x$tables$id_marker_cov$latent, level = 3L)
    }
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
      mw <- mw[, c("term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      mw$Estimate <- round(mw$Estimate, digits)
      mw$Est.Error <- round(mw$Est.Error, digits)
      mw$Q2.5 <- round(mw$Q2.5, digits)
      mw$Q97.5 <- round(mw$Q97.5, digits)
      mw$Rhat <- round(mw$Rhat, 3)
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
#' @return A nested list with top-level entries `formulaLong` and `formulaDist`.
#'   `formulaLong` contains random-effect summaries for longitudinal model
#'   components (`id`, `marker`, `marker_by_id_latent`, when present).
#'   `formulaDist` contains distributional random-effect summaries organised by
#'   parameter and family scope (e.g., `sigma$student_t`, `nu$allFamilies`).
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

  summarize_block <- function(varnames, group) {
    if (length(varnames) == 0) return(NULL)
    s <- as.data.frame(.summarise_draws_diag(fit, varnames, draws = draws, seed = seed))
    out <- s[, c("variable", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    names(out)[1] <- "term"
    out$group <- group

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
    id = summarize_block(u_vars, "id"),
    marker = if (sd$R_mk > 0) summarize_block(v_vars, "marker") else NULL,
    marker_by_id_latent = summarize_block(zw_vars, "marker_by_id_latent")
  )

  show_marker_weights <- isTRUE(sd$assoc_cv_total == 1) ||
    isTRUE(sd$assoc_cv_marker == 1) ||
    isTRUE(sd$assoc_cs_total == 1) ||
    isTRUE(sd$assoc_cs_marker == 1)

  mw_vars <- vars[grepl("^marker_weights_eff\\[", vars)]
  if (show_marker_weights && length(mw_vars) > 0) {
    mw <- summarize_block(mw_vars, "assoc_weight")
    if (!is.null(mw) && !is.null(sd$marker_levels)) {
      marker_terms <- sd$marker_levels
      if (length(marker_terms) == nrow(mw)) {
        mw$term <- paste0("weight: ", marker_terms)
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
    phi_beta = split_dist_scopes(summarize_dist_ranef("phi_beta")),
    tau_sde = split_dist_scopes(summarize_dist_ranef("tau_sde"))
  )
  out_dist <- out_dist[!vapply(out_dist, is.null, logical(1))]

  out <- list(
    formulaLong = out_long,
    formulaDist = out_dist
  )
  object$cache_set(cache_key, out)
  out
}

# ---- corr -----------------------------------------------------------------

#' Extract covariance summaries
#'
#' @param object A joinme fit object.
#' @param what Optional covariance block selector (`"id"`, `"marker"`).
#'   When `NULL` (default), returns a nested list for
#'   all covariance components in `formulaLong` and `formulaDist`.
#' @param draws Number of draws to use for summaries.
#' @param ... Unused.
#'
#' @return If `what` is supplied, a data frame for the requested covariance
#'   block. Otherwise, a nested list with `formulaLong` and `formulaDist`
#'   covariance summaries.
#' @export
corr.JoinMeFit <- function(object, what = NULL, draws = NULL, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!is.null(what)) {
    what <- match.arg(what, choices = c("id", "marker"))
  }

  summarize_cov_matrix <- function(tau_prefix, Lcorr_prefix, dim, label) {
    if (dim <= 0) {
      return(data.frame(
        block = label, row = integer(0), col = integer(0),
        Estimate = numeric(0), Est.Error = numeric(0), Q2.5 = numeric(0), Q97.5 = numeric(0),
        Rhat = numeric(0), ess_bulk = numeric(0), ess_tail = numeric(0)
      ))
    }
    tau_vars <- paste0(tau_prefix, "[", seq_len(dim), "]")
    L_vars <- as.vector(outer(seq_len(dim), seq_len(dim),
                              function(r, c) paste0(Lcorr_prefix, "[", r, ",", c, "]")))
    vars <- c(tau_vars, L_vars)

    dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = 1)
    tau_draws <- dmat[, tau_vars, drop = FALSE]
    L_draws <- dmat[, L_vars, drop = FALSE]
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
      ss <- .summarize_draw_col(vals)
      rhat <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(vals)), error = function(e) NA_real_))
      ess_bulk <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(vals)), error = function(e) NA_real_))
      ess_tail <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(vals)), error = function(e) NA_real_))
      out[[length(out) + 1]] <- cbind(block = label, row = r, col = c, t(ss), Rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail)
    }
    df <- as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE)
    df$row <- as.integer(df$row)
    df$col <- as.integer(df$col)
    df$Estimate <- as.numeric(df$Estimate)
    df$Est.Error <- as.numeric(df$Est.Error)
    df$Q2.5 <- as.numeric(df$Q2.5)
    df$Q97.5 <- as.numeric(df$Q97.5)
    df$Rhat <- as.numeric(df$Rhat)
    df$ess_bulk <- as.numeric(df$ess_bulk)
    df$ess_tail <- as.numeric(df$ess_tail)
    df
  }

  if (!is.null(what)) {
    if (what == "id") return(summarize_cov_matrix("tau_u", "Lcorr_u", sd$R_id, "id"))
    if (what == "marker") {
      if (sd$R_mk <= 0) {
        cli::cli_abort("Marker corr requested but R_mk = 0.")
      }
      return(summarize_cov_matrix("tau_v", "Lcorr_v", sd$R_mk, "marker"))
    }
  }

  summarize_dist_corr <- function(param_name) {
    n_re <- as.integer(sd[[paste0("n_re_", param_name)]] %||% 0L)
    if (n_re <= 0L) return(NULL)

    K_vec <- as.integer(sd[[paste0("K_", param_name)]] %||% integer(0))
    if (length(K_vec) != n_re) return(NULL)

    tau_prefix <- paste0("tau_", param_name)
    term_labels <- object$config$dist$dist_re_terms[[param_name]] %||% rep(NA_character_, n_re)

    out_terms <- vector("list", n_re)
    for (j in seq_len(n_re)) {
      tau_vars <- paste0(tau_prefix, "[", j, ",", seq_len(K_vec[j]), "]")
      dmat <- .get_draws_matrix(fit, variables = tau_vars, draws = draws, seed = 1)
      if (ncol(dmat) == 0) next

      rows <- list()
      idx <- 1L
      for (r in seq_len(K_vec[j])) {
        for (c in seq_len(K_vec[j])) {
          if (r == c) {
            vals <- as.numeric(dmat[, tau_vars[r]])^2
          } else {
            vals <- rep(0, nrow(dmat))
          }
          ss <- .summarize_draw_col(vals)
          rhat <- suppressWarnings(tryCatch(as.numeric(posterior::rhat(vals)), error = function(e) NA_real_))
          ess_bulk <- suppressWarnings(tryCatch(as.numeric(posterior::ess_basic(vals)), error = function(e) NA_real_))
          ess_tail <- suppressWarnings(tryCatch(as.numeric(posterior::ess_tail(vals)), error = function(e) NA_real_))
          rows[[idx]] <- data.frame(
            block = term_labels[j] %||% paste0("re_term_", j),
            row = r,
            col = c,
            Estimate = as.numeric(ss[["Estimate"]]),
            Est.Error = as.numeric(ss[["Est.Error"]]),
            Q2.5 = as.numeric(ss[["Q2.5"]]),
            Q97.5 = as.numeric(ss[["Q97.5"]]),
            Rhat = rhat,
            ess_bulk = ess_bulk,
            ess_tail = ess_tail,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
      out_terms[[j]] <- do.call(rbind, rows)
    }
    out_terms <- Filter(Negate(is.null), out_terms)
    if (length(out_terms) == 0) return(NULL)
    out_df <- do.call(rbind, out_terms)
    scope_names <- ifelse(grepl("^family=", out_df$block),
                          sub("^family=([^:]+)::.*$", "\\1", out_df$block),
                          "allFamilies")
    split(out_df, scope_names)
  }

  out <- list(
    formulaLong = list(
      id = summarize_cov_matrix("tau_u", "Lcorr_u", sd$R_id, "id"),
      marker = if (sd$R_mk > 0) summarize_cov_matrix("tau_v", "Lcorr_v", sd$R_mk, "marker") else NULL
    ),
    formulaDist = list(
      sigma = summarize_dist_corr("sigma"),
      nu = summarize_dist_corr("nu"),
      phi = summarize_dist_corr("phi"),
      alpha = summarize_dist_corr("alpha"),
      phi_beta = summarize_dist_corr("phi_beta"),
      tau_sde = summarize_dist_corr("tau_sde")
    )
  )
  out$formulaLong <- out$formulaLong[!vapply(out$formulaLong, is.null, logical(1))]
  out$formulaDist <- out$formulaDist[!vapply(out$formulaDist, is.null, logical(1))]
  out
}

#' Extract posterior covariance summaries
#'
#' @description
#' Returns posterior covariance summaries for fitted joinme models.
#' This method delegates to [corr()] and supports selecting specific blocks
#' through `...` (for example, `what = "id"`).
#'
#' @param object A joinme fit object.
#' @param ... Additional arguments passed to [corr()].
#'
#' @return Posterior covariance summary table(s).
#' @method vcov JoinMeFit
#' @export
vcov.JoinMeFit <- function(object, ...) {
  corr(object, ...)
}
