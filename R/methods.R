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

#' Round summary tables using the standard joinme summary schema
#'
#' @description
#' Applies the same rounding rules used throughout joinme summaries so new
#' posterior reports remain directly comparable to [summary.JoinMeFit()].
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
.round_joinme_summary_table <- function(tbl, digits = 3) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(tbl)
  }

  rounded <- tbl
  for (col_name in intersect(c("Estimate", "Est.Error", "Q2.5", "Q97.5"), names(rounded))) {
    rounded[[col_name]] <- round(rounded[[col_name]], digits)
  }
  if ("Rhat" %in% names(rounded)) {
    rounded$Rhat <- round(rounded$Rhat, 3)
  }
  rounded
}

#' Find the first posterior variable available in a candidate set
#'
#' @description
#' Many fitted objects retain both raw and effective-scale parameter names. This
#' helper selects the first available variable name from a priority-ordered
#' character vector.
#'
#' @param all_vars Character vector of posterior variable names available in the
#'   fitted object.
#' @param candidates Character vector ordered from preferred to fallback names.
#'
#' @return A single character scalar, or `NULL` when none of the candidates are
#'   present.
#' @keywords internal
.first_available_draw_var <- function(all_vars, candidates) {
  candidates <- as.character(candidates %||% character(0))
  hit <- candidates[candidates %in% all_vars]
  if (!length(hit)) {
    return(NULL)
  }
  hit[[1]]
}

#' Build readable labels for covariance-style association components
#'
#' @description
#' Converts the internal lower-triangular indexing used for covariance-style
#' association channels into term labels that refer to the underlying random
#' effect basis. This keeps downstream reports aligned with the design-matrix
#' terms seen in model summaries.
#'
#' @param term_key Association channel name. Supported values are `"corr"` and
#'   `"vcov"`.
#' @param sd Stan-data list stored in a `JoinMeFit` object.
#' @param n_components Optional expected number of components. When supplied,
#'   the output is truncated to this length.
#'
#' @return A character vector of human-readable component labels.
#' @keywords internal
.assoc_component_display_labels <- function(term_key, sd, n_components = NULL) {
  term_key <- as.character(term_key %||% "")[1]
  q_idm <- as.integer(sd$Q_idm %||% 0L)
  diagonal_only <- identical(term_key, "vcov") && as.integer(sd$indep_idmarker_cov %||% 0L) == 1L
  include_diag <- identical(term_key, "vcov")

  n_expected <- .assoc_transform_component_count(term_key, q_idm, diagonal_only = diagonal_only)
  if (is.null(n_components)) {
    n_components <- n_expected
  }
  n_components <- as.integer(n_components %||% 0L)
  if (n_components <= 0L) {
    return(character(0))
  }

  base_terms <- as.character(sd$zidm_cols %||% paste0("term_", seq_len(max(q_idm, 1L))))
  if (length(base_terms) < q_idm) {
    base_terms <- c(base_terms, paste0("term_", seq.int(length(base_terms) + 1L, q_idm)))
  }
  base_terms <- base_terms[seq_len(max(q_idm, 1L))]

  feature_map <- .assoc_cov_feature_map(
    q_idm = q_idm,
    diagonal_only = diagonal_only,
    include_diag = include_diag
  )
  if (is.null(dim(feature_map)) || nrow(feature_map) == 0L) {
    return(.assoc_transform_component_labels(term_key, n_components))
  }

  feature_map <- feature_map[seq_len(min(nrow(feature_map), n_components)), , drop = FALSE]
  row_idx <- if (!is.null(colnames(feature_map)) && "row" %in% colnames(feature_map)) {
    feature_map[, "row"]
  } else {
    feature_map[, 1]
  }
  col_idx <- if (!is.null(colnames(feature_map)) && "col" %in% colnames(feature_map)) {
    feature_map[, "col"]
  } else {
    feature_map[, 2]
  }
  paste0(
    term_key,
    "[",
    base_terms[row_idx],
    ", ",
    base_terms[col_idx],
    "]"
  )
}

#' Build explicit metadata for covariance-style association components
#'
#' @description
#' Covariance-style association channels are indexed over lower-triangular
#' coordinates of the marker-by-id random-effect basis. This helper reconstructs
#' those coordinates using the fitted random-effect term labels so posterior
#' displays can show both a compact combined term and the underlying row and
#' column terms.
#'
#' @param term_key Association channel name. Supported values are `"corr"` and
#'   `"vcov"`.
#' @param sd Stan-data list stored in the fitted object.
#' @param n_components Number of covariance-style components to report.
#'
#' @return A data frame with columns `component`, `row`, `col`, and `term`.
#' @keywords internal
.assoc_component_metadata <- function(term_key, sd, n_components) {
  term_key <- as.character(term_key %||% "")[1]
  n_components <- as.integer(n_components %||% 0L)
  if (!(term_key %in% c("corr", "vcov")) || n_components <= 0L) {
    return(NULL)
  }

  q_idm <- as.integer(sd$Q_idm %||% 0L)
  diagonal_only <- identical(term_key, "vcov") && as.integer(sd$indep_idmarker_cov %||% 0L) == 1L
  feature_map <- .assoc_cov_feature_map(
    q_idm = q_idm,
    diagonal_only = diagonal_only,
    include_diag = identical(term_key, "vcov")
  )
  if (is.null(dim(feature_map)) || nrow(feature_map) == 0L) {
    return(NULL)
  }

  feature_map <- feature_map[seq_len(min(nrow(feature_map), n_components)), , drop = FALSE]
  row_idx <- if (!is.null(colnames(feature_map)) && "row" %in% colnames(feature_map)) {
    feature_map[, "row"]
  } else {
    feature_map[, 1]
  }
  col_idx <- if (!is.null(colnames(feature_map)) && "col" %in% colnames(feature_map)) {
    feature_map[, "col"]
  } else {
    feature_map[, 2]
  }
  base_terms <- as.character(sd$zidm_cols %||% paste0("term_", seq_len(max(q_idm, 1L))))
  if (length(base_terms) < q_idm) {
    base_terms <- c(base_terms, paste0("term_", seq.int(length(base_terms) + 1L, q_idm)))
  }
  base_terms <- base_terms[seq_len(max(q_idm, 1L))]

  data.frame(
    component = seq_len(nrow(feature_map)),
    row = base_terms[row_idx],
    col = base_terms[col_idx],
    term = .assoc_component_display_labels(term_key, sd, n_components = nrow(feature_map)),
    stringsAsFactors = FALSE
  )
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

#' Weighted association terms constrained positive in summary output
#'
#' @param sd Standata list stored in a fitted object.
#'
#' @return Character vector of weighted association term keys whose raw `alpha`
#'   parameter is constrained positive in the fitted model.
#' @keywords internal
.positive_weighted_assoc_terms <- function(sd) {
  active_terms <- .active_weighted_assoc_terms(sd)
  if (!length(active_terms)) {
    return(character(0))
  }

  if (isTRUE(as.integer(sd$shared_marker_weights %||% 1L) == 1L)) {
    active_terms[[1L]]
  } else {
    active_terms
  }
}

#' Summarise a derived posterior draw array with MCMC diagnostics
#'
#' @description
#' Derived association effects do not necessarily exist as named Stan variables.
#' This helper converts an iteration x chain x term array into the same summary
#' schema used elsewhere in joinme: posterior mean, posterior standard
#' deviation, central 95% interval, split-chain R-hat, and bulk/tail effective
#' sample sizes.
#'
#' @param draw_array Numeric array with dimensions iteration x chain x term.
#' @param term_labels Character vector naming the third dimension.
#' @param digits Number of decimal places used for posterior location and
#'   interval summaries.
#'
#' @return A data frame with columns `term`, `Estimate`, `Est.Error`, `Q2.5`,
#'   `Q97.5`, `Rhat`, `ess_bulk`, and `ess_tail`.
#' @keywords internal
.assoc_summary_from_draw_array <- function(draw_array, term_labels, digits = 3) {
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(NULL)
  }

  term_labels <- as.character(term_labels %||% paste0("term_", seq_len(dim(draw_array)[3])))
  dimnames(draw_array) <- list(
    iteration = dimnames(draw_array)[[1]] %||% as.character(seq_len(dim(draw_array)[1])),
    chain = dimnames(draw_array)[[2]] %||% as.character(seq_len(dim(draw_array)[2])),
    variable = term_labels
  )

  draws_obj <- posterior::as_draws_array(draw_array)
  draws_df <- posterior::as_draws_df(draws_obj)
  sum_df <- data.frame(
    term = term_labels,
    do.call(rbind, lapply(term_labels, function(label) .summarize_draw_col(draws_df[[label]]))),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  diag_df <- suppressWarnings(tryCatch(
    posterior::summarise_draws(draws_obj, "rhat", "ess_bulk", "ess_tail"),
    error = function(e) NULL
  ))
  if (!is.null(diag_df)) {
    diag_df <- diag_df[, c("variable", "rhat", "ess_bulk", "ess_tail"), drop = FALSE]
    names(diag_df) <- c("term", "Rhat", "ess_bulk", "ess_tail")
    sum_df <- merge(sum_df, diag_df, by = "term", all.x = TRUE, sort = FALSE)
  } else {
    sum_df$Rhat <- NA_real_
    sum_df$ess_bulk <- NA_real_
    sum_df$ess_tail <- NA_real_
  }

  .round_joinme_summary_table(sum_df, digits = digits)
}

#' Flatten a derived posterior draw array into an MCMC sample matrix
#'
#' @description
#' Returns posterior samples as a draws-by-term matrix so each association term
#' can be inspected directly without an additional summary step.
#'
#' @param draw_array Numeric array with dimensions iteration x chain x term.
#' @param term_labels Character vector naming the third dimension.
#'
#' @return A numeric matrix with one column per term label.
#' @keywords internal
.assoc_matrix_from_draw_array <- function(draw_array, term_labels) {
  if (is.null(draw_array) || length(dim(draw_array)) != 3L || dim(draw_array)[3] == 0L) {
    return(NULL)
  }

  term_labels <- as.character(term_labels %||% paste0("term_", seq_len(dim(draw_array)[3])))
  dimnames(draw_array) <- list(
    iteration = dimnames(draw_array)[[1]] %||% as.character(seq_len(dim(draw_array)[1])),
    chain = dimnames(draw_array)[[2]] %||% as.character(seq_len(dim(draw_array)[2])),
    variable = term_labels
  )

  out <- posterior::as_draws_matrix(posterior::as_draws_array(draw_array))
  out[, term_labels, drop = FALSE]
}

#' Build posterior association effects for a fitted joinme model
#'
#' @description
#' Reconstructs the posterior association effects that enter the survival linear
#' predictor.
#'
#' For weighted current-value and current-slope channels (`cv_total`,
#' `cs_total`, `cv_marker`, `cs_marker`), the returned effect is the draw-wise
#' product of the association coefficient and the marker weight, divided by the
#' number of markers to match the scale used in the fitted hazard contribution.
#'
#' For scalar channels (`cv_mean`, `cs_mean`) the returned effect is simply the
#' posterior coefficient. For covariance-style channels (`corr`, `vcov`) the
#' returned effects are grouped by their labelled covariance component.
#'
#' @param object A `JoinMeFit` object.
#' @param draws Optional number of posterior draws to retain.
#' @param seed Random seed used when subsetting posterior draws.
#' @param digits Number of digits used when `summary = TRUE`.
#' @param summary Logical. If `TRUE`, return posterior summaries with the same
#'   inferential columns as [summary.JoinMeFit()]. If `FALSE`, return raw MCMC
#'   sample matrices.
#' @param ... Unused.
#'
#' @return A named list with class `PosteriorAssoc`. Each list element contains
#'   either a posterior summary table or a draws-by-term matrix for one
#'   association channel.
#' @export
assoc <- function(object, ...) {
  UseMethod("assoc")
}

#' @rdname assoc
#' @export
posterior_assoc <- function(object, ...) {
  assoc(object, ...)
}

#' @rdname assoc
#' @export
assoc.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")

  sd <- object$stan_data
  fit <- object$fit
  if (is.null(draws)) {
    draws <- object$config$draws_default
  }
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  assertthat::assert_that(is.numeric(digits) && digits >= 0, msg = "digits must be non-negative.")

  cache_key <- if (isTRUE(summary)) {
    paste0("posterior_assoc_summary_", draws, "_", digits)
  } else {
    paste0("posterior_assoc_draws_", draws)
  }
  cached <- object$cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  marker_levels <- as.character(sd$marker_levels %||% paste0("marker_", seq_len(as.integer(sd$D %||% 0L))))
  n_markers <- length(marker_levels)
  weighted_divisor <- if (n_markers > 0L) n_markers else 1L

  build_scalar_term <- function(var_name, label) {
    if (is.null(var_name)) {
      return(NULL)
    }
    draw_arr <- .get_draws_array(fit, variables = var_name, draws = draws, seed = seed)
    if (isTRUE(summary)) {
      .assoc_summary_from_draw_array(draw_arr, term_labels = label, digits = digits)
    } else {
      .assoc_matrix_from_draw_array(draw_arr, term_labels = label)
    }
  }

  build_weighted_term <- function(var_name, term_key) {
    weight_array <- .association_marker_weight_array(object, term_key = term_key, draws = draws, seed = seed, all_vars = all_vars)
    if (is.null(var_name) || is.null(weight_array) || n_markers == 0L) {
      return(NULL)
    }
    alpha_arr <- .get_draws_array(fit, variables = var_name, draws = draws, seed = seed)
    effect_arr <- array(
      NA_real_,
      dim = c(dim(alpha_arr)[1], dim(alpha_arr)[2], n_markers),
      dimnames = list(
        iteration = dimnames(alpha_arr)[[1]],
        chain = dimnames(alpha_arr)[[2]],
        variable = marker_levels
      )
    )

    for (marker_index in seq_len(n_markers)) {
      effect_arr[, , marker_index] <- alpha_arr[, , 1] * weight_array[, , marker_index] / weighted_divisor
    }

    if (isTRUE(summary)) {
      .assoc_summary_from_draw_array(effect_arr, term_labels = marker_levels, digits = digits)
    } else {
      .assoc_matrix_from_draw_array(effect_arr, term_labels = marker_levels)
    }
  }

  build_vector_term <- function(var_names, term_labels, component_meta = NULL) {
    var_names <- as.character(var_names %||% character(0))
    if (!length(var_names)) {
      return(NULL)
    }
    draw_arr <- .get_draws_array(fit, variables = var_names, draws = draws, seed = seed)
    if (isTRUE(summary)) {
      out_tbl <- .assoc_summary_from_draw_array(draw_arr, term_labels = term_labels, digits = digits)
      if (!is.null(component_meta) && !is.null(out_tbl)) {
        idx <- match(out_tbl$term, component_meta$term)
        out_tbl$component <- component_meta$component[idx]
        out_tbl$row <- component_meta$row[idx]
        out_tbl$col <- component_meta$col[idx]
        out_tbl <- out_tbl[, c("component", "row", "col", "term", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), drop = FALSE]
      }
      out_tbl
    } else {
      .assoc_matrix_from_draw_array(draw_arr, term_labels = term_labels)
    }
  }

  corr_assoc_vars <- grep("^alpha_corr_eff\\[", all_vars, value = TRUE)
  if (length(corr_assoc_vars) == 0L) {
    corr_assoc_vars <- grep("^alpha_corr\\[", all_vars, value = TRUE)
  }
  vcov_assoc_vars <- grep("^alpha_vcov_eff\\[", all_vars, value = TRUE)
  if (length(vcov_assoc_vars) == 0L) {
    vcov_assoc_vars <- grep("^alpha_vcov\\[", all_vars, value = TRUE)
  }

  out <- list(
    cv_total = if (isTRUE(sd$assoc_cv_total == 1L)) {
      build_weighted_term(.first_available_draw_var(all_vars, c("alpha_cv_total_eff", "alpha_cv_total")), term_key = "cv_total")
    } else NULL,
    cv_mean = if (isTRUE(sd$assoc_cv_mean == 1L)) {
      build_scalar_term(.first_available_draw_var(all_vars, c("alpha_cv_mean_eff", "alpha_cv_mean")), "cv_mean")
    } else NULL,
    cv_marker = if (isTRUE(sd$assoc_cv_marker == 1L)) {
      build_weighted_term(.first_available_draw_var(all_vars, c("alpha_cv_marker_eff", "alpha_cv_marker")), term_key = "cv_marker")
    } else NULL,
    cs_total = if (isTRUE(sd$assoc_cs_total == 1L)) {
      build_weighted_term(.first_available_draw_var(all_vars, c("alpha_cs_total_eff", "alpha_cs_total")), term_key = "cs_total")
    } else NULL,
    cs_mean = if (isTRUE(sd$assoc_cs_mean == 1L)) {
      build_scalar_term(.first_available_draw_var(all_vars, c("alpha_cs_mean_eff", "alpha_cs_mean")), "cs_mean")
    } else NULL,
    cs_marker = if (isTRUE(sd$assoc_cs_marker == 1L)) {
      build_weighted_term(.first_available_draw_var(all_vars, c("alpha_cs_marker_eff", "alpha_cs_marker")), term_key = "cs_marker")
    } else NULL,
    corr = if (isTRUE(sd$assoc_corr == 1L)) {
      corr_meta <- .assoc_component_metadata("corr", sd, length(corr_assoc_vars))
      build_vector_term(
        corr_assoc_vars,
        term_labels = corr_meta$term %||% .assoc_component_display_labels("corr", sd, length(corr_assoc_vars)),
        component_meta = corr_meta
      )
    } else NULL,
    vcov = if (isTRUE(sd$assoc_vcov == 1L)) {
      vcov_meta <- .assoc_component_metadata("vcov", sd, length(vcov_assoc_vars))
      build_vector_term(
        vcov_assoc_vars,
        term_labels = vcov_meta$term %||% .assoc_component_display_labels("vcov", sd, length(vcov_assoc_vars)),
        component_meta = vcov_meta
      )
    } else NULL
  )
  out <- out[!vapply(out, is.null, logical(1))]

  if (!length(out)) {
    cli::cli_abort(c(
      x = "No association effects are available in this fitted object.",
      i = "Fit a model with association terms such as {.val cv_total}, {.val cv_mean}, {.val corr}, or {.val vcov}."
    ))
  }

  out <- structure(
    out,
    class = "PosteriorAssoc",
    metadata = list(
      summary = isTRUE(summary),
      draws = draws,
      digits = digits,
      marker_levels = marker_levels
    )
  )

  object$cache_set(cache_key, out)
  out
}

#' Posterior summary alias for joinme objects
#'
#' @description
#' Provides a user-facing alias to [summary()] so posterior summaries can be
#' requested with terminology that emphasizes Bayesian output.
#'
#' @param object A joinme object.
#' @param ... Additional arguments forwarded to [summary()].
#'
#' @return The same object that [summary()] would return for the supplied class.
#' @export
posterior_summary <- function(object, ...) {
  UseMethod("posterior_summary")
}

#' @rdname posterior_summary
#' @export
posterior_summary.JoinMeFit <- function(object, ...) {
  summary(object, ...)
}

#' @rdname posterior_summary
#' @export
posterior_summary.JoinMeDynPred <- function(object, ...) {
  summary(object, ...)
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

#' Format covariance-style posterior association summaries for printing
#'
#' @description
#' `corr` and `vcov` posterior association summaries are easier to read when the
#' lower-triangular components are displayed in row-wise groups rather than as a
#' single flat list. This helper keeps the returned object unchanged and only
#' prepares grouped tables for console output.
#'
#' @param tbl Posterior association summary table.
#'
#' @return A named list of data frames, one per matrix row label.
#' @keywords internal
.posterior_assoc_matrix_groups <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(list())
  }
  if (!all(c("row", "col") %in% names(tbl))) {
    return(list())
  }

  row_levels <- unique(as.character(tbl$row))
  row_levels <- row_levels[!is.na(row_levels) & nzchar(row_levels)]
  groups <- lapply(row_levels, function(row_label) {
    block <- tbl[as.character(tbl$row) == row_label, , drop = FALSE]
    block <- block[, c(
      intersect(c("col", "Estimate", "Est.Error", "Q2.5", "Q97.5", "Rhat", "ess_bulk", "ess_tail"), names(block))
    ), drop = FALSE]
    names(block)[names(block) == "col"] <- "term"
    block
  })
  names(groups) <- row_levels
  groups
}

#' Print posterior association effects
#'
#' @description
#' Prints association-effect summaries or raw posterior sample availability in a
#' layout aligned with [print.summary_JoinMeFit()].
#'
#' @param x A `PosteriorAssoc` object returned by [assoc()] or
#'   [posterior_assoc()].
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.PosteriorAssoc <- function(x, ...) {
  meta <- attr(x, "metadata") %||% list()

  .cli_summary_heading("Posterior association effects", level = 1L)
  .cli_print_bullets(c(
    if (!is.null(meta$draws)) paste0("Posterior draws: ", meta$draws) else NULL,
    paste0("Returned scale: ", if (isTRUE(meta$summary)) "posterior summary" else "MCMC samples")
  ))

  if (!length(x)) {
    return(invisible(x))
  }

  for (term_name in names(x)) {
    .cli_summary_heading(paste0("Association term: ", term_name), level = 2L)
    if (isTRUE(meta$summary)) {
      term_tbl <- x[[term_name]]
      if (term_name %in% c("corr", "vcov") && is.data.frame(term_tbl) && all(c("row", "col") %in% names(term_tbl))) {
        .cli_print_bullets("Displayed by matrix row using fitted random-effect term labels.")
        matrix_groups <- .posterior_assoc_matrix_groups(term_tbl)
        for (row_label in names(matrix_groups)) {
          .cli_summary_heading(paste0("row = ", row_label), level = 3L)
          .cli_print_table(matrix_groups[[row_label]])
        }
      } else {
        .cli_print_table(term_tbl)
      }
    } else {
      draw_block <- x[[term_name]]
      desc <- data.frame(
        term = colnames(draw_block) %||% character(0),
        draws = rep(nrow(draw_block), ncol(draw_block)),
        stringsAsFactors = FALSE
      )
      .cli_print_table(desc)
    }
  }

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

  corr_assoc_vars <- grep("^alpha_corr_eff\\[", all_vars, value = TRUE)
  if (length(corr_assoc_vars) == 0L) corr_assoc_vars <- grep("^alpha_corr\\[", all_vars, value = TRUE)
  vcov_assoc_vars <- grep("^alpha_vcov_eff\\[", all_vars, value = TRUE)
  if (length(vcov_assoc_vars) == 0L) vcov_assoc_vars <- grep("^alpha_vcov\\[", all_vars, value = TRUE)
  a_vars <- c("alpha_cv_total", "alpha_cv_mean", "alpha_cv_marker", "alpha_cs_total", "alpha_cs_mean", "alpha_cs_marker")
  a_vars <- c(a_vars, corr_assoc_vars)
  a_vars <- c(a_vars, vcov_assoc_vars)
  a_vars <- a_vars[a_vars %in% all_vars]
  active_vars <- c(
    if (isTRUE(sd$assoc_cv_total == 1)) "alpha_cv_total",
    if (isTRUE(sd$assoc_cv_mean == 1)) "alpha_cv_mean",
    if (isTRUE(sd$assoc_cv_marker == 1)) "alpha_cv_marker",
    if (isTRUE(sd$assoc_cs_total == 1)) "alpha_cs_total",
    if (isTRUE(sd$assoc_cs_mean == 1)) "alpha_cs_mean",
    if (isTRUE(sd$assoc_cs_marker == 1)) "alpha_cs_marker",
    if (isTRUE(sd$assoc_corr == 1)) corr_assoc_vars,
    if (isTRUE(sd$assoc_vcov == 1)) vcov_assoc_vars
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
    s_a$term <- sub("_eff\\[", "[", s_a$term, perl = TRUE)
    positive_terms <- .positive_weighted_assoc_terms(sd)
    if (length(positive_terms) > 0L) {
      positive_idx <- s_a$term %in% positive_terms
      s_a$term[positive_idx] <- paste0(s_a$term[positive_idx], " (+)")
    }
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
      tmax = sd$tmax_internal %||% 1.0,
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
#' @param formulaVCov Optional updated covariance formula (full or update form).
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
  formulaVCov = NULL,
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
  call_obj$formulaVCov <- update_formula(object$formulaVCov, formulaVCov, "formulaVCov")
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

# ---- fixef / ranef --------------------------------------------------------

#' @keywords internal
.joinme_draw_long_from_matrix <- function(draw_matrix, meta) {
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
.joinme_public_time_scale <- function(stan_data) {
  scale_factor <- suppressWarnings(as.numeric(stan_data$tmax_internal %||% 1.0))
  if (!is.finite(scale_factor) || length(scale_factor) != 1L || scale_factor <= 0) {
    return(1.0)
  }
  scale_factor
}

#' @keywords internal
.rescale_public_longitudinal_draws <- function(draws_df, stan_data, idx, term_labels) {
  if (is.null(draws_df) || !nrow(draws_df)) return(draws_df)

  scale_factor <- .joinme_public_time_scale(stan_data)
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

  scale_factor <- .joinme_public_time_scale(stan_data)
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
.joinme_id_labels <- function(object, n_id) {
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

  beta_vars <- paste0("beta[", seq_len(sd$P %||% 0L), "]")
  beta_vars <- beta_vars[beta_vars %in% all_vars]
  beta_terms <- sd$x_cols %||% beta_vars
  if (length(beta_terms) != length(beta_vars)) beta_terms <- beta_vars

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
      if (length(mw_vars) == 0L && shared_weights) {
        mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D), "]")
        mw_vars <- mw_vars[mw_vars %in% all_vars]
      }
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
  .joinme_draw_long_from_matrix(dmat, meta)
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
    phi_beta = list(prefix = "beta_phi_beta", cols = dist_cols$phi_beta %||% character(0)),
    tau_sde = list(prefix = "beta_tau_sde", cols = dist_cols$tau_sde %||% character(0))
  )

  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    vars <- grep(paste0("^", spec$prefix, "\\["), all_vars, value = TRUE)
    if (!length(vars)) next
    dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = seed)
    term_labels <- if (length(spec$cols) == length(vars)) spec$cols else vars
    out[[nm]] <- .joinme_draw_long_from_matrix(dmat, data.frame(term = term_labels, stringsAsFactors = FALSE))
  }
  out
}

#' @keywords internal
.posterior_ranef_joinmefit <- function(object, draws = NULL, seed = 1) {
  fit <- object$fit
  sd <- object$stan_data
  cfg <- object$config
  if (is.null(draws)) draws <- object$config$draws_default
  all_vars <- tryCatch(posterior::variables(.get_draws_obj(fit)), error = function(e) character(0))
  n_id <- as.integer(sd$n_id %||% 0L)
  id_labels <- .joinme_id_labels(object, n_id)
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
    .joinme_draw_long_from_matrix(dmat, meta)
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
      if (length(vars) == 0L && shared_weights) {
        vars <- paste0("marker_weights_eff[", seq_len(sd$D), "]")
        vars <- vars[vars %in% all_vars]
      }
      if (!length(vars)) next
      dmat <- .get_draws_matrix(fit, variables = vars, draws = draws, seed = seed)
      meta <- data.frame(
        term = vapply(marker_labels[seq_along(vars)], function(marker_label) {
          .marker_weight_summary_label(term_key, marker_label, shared_marker_weights = shared_weights)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      weight_rows[[length(weight_rows) + 1L]] <- .joinme_draw_long_from_matrix(dmat, meta)
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
    out <- .joinme_draw_long_from_matrix(do.call(cbind, value_cols), do.call(rbind, meta_rows))
    out$scope <- ifelse(grepl("^family=", out$term), sub("^family=([^:]+)::.*$", "\\1", out$term), "allFamilies")
    out
  }

  out_dist <- list(
    sigma = summarize_dist_draws("sigma"),
    nu = summarize_dist_draws("nu"),
    phi = summarize_dist_draws("phi"),
    alpha = summarize_dist_draws("alpha"),
    phi_beta = summarize_dist_draws("phi_beta"),
    tau_sde = summarize_dist_draws("tau_sde")
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
        out_vcov <- .joinme_draw_long_from_matrix(do.call(cbind, value_cols), do.call(rbind, meta_rows))
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
.posterior_coef_joinmefit <- function(object, draws = NULL, seed = 1) {
  fixed_long <- .posterior_fixef_matrix(object, draws = draws, seed = seed)
  beta_only <- fixed_long
  if (is.matrix(beta_only) && ncol(beta_only) > 0L) {
    weight_cols <- grepl("^weight", colnames(beta_only))
    if (any(weight_cols)) beta_only <- beta_only[, !weight_cols, drop = FALSE]
  }
  ranef_draws <- .posterior_ranef_joinmefit(object, draws = draws, seed = seed)
  dist_fixed <- .posterior_distreg_draws(object, draws = draws, seed = seed)
  event_fixed <- .posterior_event_coef_draws(object, draws = draws, seed = seed)

  out_long <- ranef_draws$formulaLong
  out_long$id <- .combine_fixed_random_long_draws(out_long$id, beta_only)
  out_long$marker <- .combine_fixed_random_long_draws(out_long$marker, beta_only)
  out_long$marker_by_id <- .combine_fixed_random_long_draws(out_long$marker_by_id, beta_only)
  out_long$population <- if (is.matrix(beta_only) && ncol(beta_only) > 0L) {
    .joinme_draw_long_from_matrix(beta_only, data.frame(term = colnames(beta_only), stringsAsFactors = FALSE))
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
#' @param object A joinme fit object.
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
fixef.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
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
      if (length(mw_vars) == 0L && shared_weights) {
        mw_vars <- paste0("marker_weights_eff[", seq_len(sd$D), "]")
        mw_vars <- mw_vars[mw_vars %in% all_vars]
      }
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

#' Posterior fixed-effect alias for fitted joinme models
#'
#' @param object A `JoinMeFit` object.
#' @param ... Additional arguments forwarded to [fixef()].
#'
#' @return The same object returned by `fixef(object, summary = FALSE, ...)`.
#' @export
posterior_fixef <- function(object, ...) {
  fixef(object, summary = FALSE, ...)
}

#' Extract random effects
#'
#' @param object A joinme fit object.
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
ranef.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  assertthat::assert_that(is.logical(summary) && length(summary) == 1L && !is.na(summary),
                          msg = "summary must be TRUE or FALSE.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!isTRUE(summary)) {
    cache_key <- paste0("posterior_ranef_", draws)
    cached <- object$cache_get(cache_key)
    if (!is.null(cached)) return(cached)
    out <- .posterior_ranef_joinmefit(object, draws = draws, seed = seed)
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

#' Posterior random-effect alias for fitted joinme models
#'
#' @param object A `JoinMeFit` object.
#' @param ... Additional arguments forwarded to [ranef()].
#'
#' @return The same object returned by `ranef(object, summary = FALSE, ...)`.
#' @export
posterior_ranef <- function(object, ...) {
  ranef(object, summary = FALSE, ...)
}

#' Combined posterior coefficients for fitted joinme models
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
#' @param object A `JoinMeFit` object.
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
#' @method coef JoinMeFit
#' @export
coef.JoinMeFit <- function(object, draws = NULL, seed = 1, digits = 3, summary = TRUE, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
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

  out <- .posterior_coef_joinmefit(object, draws = draws, seed = seed)
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

#' Posterior coefficient alias for fitted joinme models
#'
#' @param object A `JoinMeFit` object.
#' @param ... Additional arguments forwarded to [coef()].
#'
#' @return The same object returned by `coef(object, summary = FALSE, ...)`.
#' @export
posterior_coef <- function(object, ...) {
  coef(object, summary = FALSE, ...)
}

# ---- corr -----------------------------------------------------------------

#' Extract posterior covariance summaries
#'
#' @description
#' Returns posterior covariance summaries for fitted joinme models. Supports
#' selecting specific blocks through `what` (for example, `what = "id"`).
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
#' @method vcov JoinMeFit
#' @export
vcov.JoinMeFit <- function(object, what = NULL, draws = NULL, ...) {
  assertthat::assert_that(inherits(object, "JoinMeFit"), msg = "Object must be a JoinMeFit instance.")
  fit <- object$fit
  sd <- object$stan_data
  if (is.null(draws)) draws <- object$config$draws_default

  if (!is.null(what)) {
    what <- match.arg(what, choices = c("id", "marker"))
  }

  summarize_cov_matrix <- function(tau_prefix, Lcorr_prefix, dim, label, diagonal_only = FALSE) {
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
      if (isTRUE(diagonal_only)) {
        Corr <- diag(dim)
      } else {
        L <- matrix(as.numeric(L_draws[i, ]), nrow = dim, ncol = dim, byrow = FALSE)
        L[upper.tri(L)] <- 0
        Corr <- L %*% t(L)
      }
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
    if (what == "id") {
      return(summarize_cov_matrix(
        "tau_u",
        "Lcorr_u",
        sd$R_id,
        "id",
        diagonal_only = as.integer(sd$indep_id_re %||% 0L) == 1L
      ))
    }
    if (what == "marker") {
      if (sd$R_mk <= 0) {
        cli::cli_abort("Marker corr requested but R_mk = 0.")
      }
      return(summarize_cov_matrix(
        "tau_v",
        "Lcorr_v",
        sd$R_mk,
        "marker",
        diagonal_only = as.integer(sd$indep_marker_re %||% 0L) == 1L
      ))
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
      id = summarize_cov_matrix(
        "tau_u",
        "Lcorr_u",
        sd$R_id,
        "id",
        diagonal_only = as.integer(sd$indep_id_re %||% 0L) == 1L
      ),
      marker = if (sd$R_mk > 0) summarize_cov_matrix(
        "tau_v",
        "Lcorr_v",
        sd$R_mk,
        "marker",
        diagonal_only = as.integer(sd$indep_marker_re %||% 0L) == 1L
      ) else NULL
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
