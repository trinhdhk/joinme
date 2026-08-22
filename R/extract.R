#' Extract posterior draws from JoiNMe objects
#'
#' @description
#' S3 generic to extract component-specific posterior results from
#' `JoiNMeFit` and `JoiNMeDynPred` objects.
#'
#' Together with [posterior_draws()], `extract()` is the lowest public posterior
#' interface. It works component by component and returns the metadata needed to
#' understand how a requested summary term maps back to the stored Stan
#' variables or prediction draw blocks.
#'
#' The coefficient hierarchy is deliberately one-directional. `extract()` owns
#' the draw-level representations selected by `what = "fixed_effects"`,
#' `"random_effects"`, and `"coefficients"`; [posterior_summary()] summarises
#' those representations; and [fixef()], [ranef()], and [coef()] are the
#' conventional high-level entry points. This single path avoids parallel
#' coefficient APIs with competing semantics.
#'
#' Use `extract()` when you need a specific model component, the corresponding
#' `term_map`, or a specialised result such as `what = "association_plot"`.
#' Use [posterior_draws()] when you want a single posterior object ready for `posterior`
#' or `bayesplot` workflows.
#'
#' @param object A supported JoiNMe object.
#' @param ... Additional method-specific arguments.
#' @seealso [posterior_draws()]
#' @export
extract <- function(object, ...) {
  UseMethod("extract")
}

#' Recover model-matrix column labels for posterior reporting
#'
#' @description
#' Returns the statistical term labels associated with one fitted design
#' matrix. Current fitted objects store these labels in `stan_data`. Objects
#' created by an earlier fitting route may contain only the numeric Stan data,
#' because character metadata were removed before sampling. For those objects,
#' this helper reconstructs the same fixed-effect or event-process model matrix
#' from the recorded formula and observed data.
#'
#' Reconstruction is deliberately restricted to column names. Posterior values
#' are never recalculated, and the observed data are not modified. Consequently,
#' the helper can repair presentation of an existing fit without changing its
#' posterior distribution.
#'
#' @param object A fitted `JoiNMeFit` object.
#' @param what Design matrix to label: `"fixef"` for the longitudinal
#'   population-level design or `"gamma_w"` for the event-process design.
#' @param n_terms Expected number of columns in the fitted design.
#'
#' @return A character vector of length `n_terms`. Stored model-matrix labels
#'   are preferred; formula-derived labels are used when stored labels are
#'   absent; stable indexed labels are the final fallback.
#' @keywords internal
#' @noRd
.fit_design_term_labels <- function(object, what = c("fixef", "gamma_w"), n_terms) {
  what <- match.arg(what)
  n_terms <- as.integer(n_terms %||% 0L)
  if (n_terms < 1L) {
    return(character(0))
  }

  # Labels saved when the design matrix was constructed are authoritative:
  # they include factor contrasts, interactions, and basis-function columns in
  # precisely the order used by Stan.
  stored_name <- if (identical(what, "fixef")) "x_cols" else "w_cols"
  stored <- as.character(object$stan_data[[stored_name]] %||% character(0))
  if (length(stored) == n_terms && all(!is.na(stored)) && all(nzchar(stored))) {
    return(stored)
  }

  # Older fitted objects may lack character metadata. Reconstruct only the
  # design-column vocabulary from the formula and original data. Scaling the
  # time variable changes values but not ordinary model-matrix column names.
  recovered <- tryCatch({
    if (identical(what, "fixef")) {
      expanded <- reformulas::expandDoubleVerts(object$formulaLong)
      fixed_formula <- reformulas::nobars(expanded)
      fixed_rhs <- stats::update(fixed_formula, . ~ .)
      fixed_rhs[[2L]] <- NULL
      colnames(stats::model.matrix(fixed_rhs, data = object$dataLong))
    } else {
      colnames(.mm_event(object$formulaEvent, object$dataEvent))
    }
  }, error = function(e) character(0))

  recovered <- as.character(recovered %||% character(0))
  if (length(recovered) == n_terms && all(!is.na(recovered)) && all(nzchar(recovered))) {
    return(recovered)
  }

  prefix <- if (identical(what, "fixef")) "beta_" else "w_"
  paste0(prefix, seq_len(n_terms))
}

#' Map longitudinal fixed-effect variables to statistical terms
#'
#' @description
#' Selects the posterior coefficient representation on the original time scale
#' and pairs every selected Stan variable with its model-matrix term.
#'
#' @param sd Stan data and fitted-model metadata.
#' @param all_vars Character vector of available posterior variables.
#' @param term_labels Optional model-matrix labels aligned with the `P`
#'   longitudinal fixed-effect columns.
#'
#' @return A two-column data frame containing friendly `term` labels and raw
#'   posterior `variable` names.
#' @keywords internal
#' @noRd
.fixed_effect_var_map <- function(sd, all_vars, term_labels = NULL) {
  p <- as.integer(sd$P %||% 0L)
  if (p <= 0L) {
    return(data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE))
  }

  beta_vars <- paste0("beta[", seq_len(p), "]")
  chosen <- beta_vars[beta_vars %in% all_vars]

  if (!length(chosen)) {
    return(data.frame(term = character(0), variable = character(0), stringsAsFactors = FALSE))
  }

  idx <- suppressWarnings(as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", chosen)))
  term_labels <- as.character(term_labels %||% sd$x_cols %||% paste0("beta_", seq_len(p)))
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
    fixed_terms <- .fit_design_term_labels(object, what = "fixef", n_terms = sd$P)
    map <- .fixed_effect_var_map(
      sd = sd,
      all_vars = all_vars,
      term_labels = fixed_terms
    )
  } else if (what == "gamma_w") {
    p_w <- as.integer(sd$p_w %||% 0L)
    if (p_w < 1L) {
      return(map)
    }
    g_vars <- paste0("gamma_w[", seq_len(p_w), "]")
    g_vars <- g_vars[g_vars %in% all_vars]
    if (length(g_vars) == 0L) {
      k_event <- sd$K_event %||% 1L
      g_vars <- as.vector(outer(
        seq_len(k_event),
        seq_len(p_w),
        function(k, j) paste0("gamma_w[", k, ",", j, "]")
      ))
      g_vars <- g_vars[g_vars %in% all_vars]
    }
    if (length(g_vars) > 0L) {
      var_idx <- regmatches(g_vars, regexec("^gamma_w\\[(\\d+)(?:,(\\d+))?\\]$", g_vars))
      k_idx <- vapply(var_idx, function(x) if (length(x) >= 2L) as.integer(x[2]) else NA_integer_, integer(1))
      j_idx <- vapply(var_idx, function(x) if (length(x) >= 3L && nzchar(x[3])) as.integer(x[3]) else as.integer(x[2]), integer(1))
      k_event <- sd$K_event %||% 1L
      g_terms <- .fit_design_term_labels(object, what = "gamma_w", n_terms = sd$p_w)
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

    beta_vars <- paste0("beta[", seq_len(sd$P), "]")
    beta_vars <- beta_vars[beta_vars %in% all_vars]
    if (length(beta_vars) > 0) {
      beta_terms <- .fit_design_term_labels(object, what = "fixef", n_terms = sd$P)
      if (length(beta_terms) != length(beta_vars)) beta_terms <- beta_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("beta: ", as.character(beta_terms)),
        variable = as.character(beta_vars),
        stringsAsFactors = FALSE
      )
    }

    tau_id_vars <- paste0("tau_u[", seq_len(sd$R_id %||% 0L), "]")
    tau_id_vars <- tau_id_vars[tau_id_vars %in% all_vars]
    if (length(tau_id_vars) > 0) {
      tau_id_terms <- sd$zid_cols %||% tau_id_vars
      if (length(tau_id_terms) != length(tau_id_vars)) tau_id_terms <- tau_id_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("id_sd: ", as.character(tau_id_terms)),
        variable = as.character(tau_id_vars),
        stringsAsFactors = FALSE
      )
    }

    tau_marker_vars <- paste0("tau_v[", seq_len(sd$R_mk %||% 0L), "]")
    tau_marker_vars <- tau_marker_vars[tau_marker_vars %in% all_vars]
    if (length(tau_marker_vars) > 0) {
      tau_marker_terms <- sd$zmk_cols %||% tau_marker_vars
      if (length(tau_marker_terms) != length(tau_marker_vars)) tau_marker_terms <- tau_marker_vars
      map_rows[[length(map_rows) + 1L]] <- data.frame(
        term = paste0("marker_sd: ", as.character(tau_marker_terms)),
        variable = as.character(tau_marker_vars),
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
#' This differs from [posterior_draws()] in two important ways:
#' - `extract()` keeps the request scoped to one semantic component such as
#'   fixed effects, survival coefficients, association terms, distributional
#'   terms, or likelihood-scale parameters.
#' - `extract()` returns a structured list with `posterior_draws` plus `term_map`
#'   (and for `what = "association_plot"`, additional plotting support data)
#'   instead of a single `posterior` draws object.
#'
#' In short, use `extract()` when you need component-aware extraction and use
#' [posterior_draws()] when you need one renamed posterior object for general downstream
#' analysis.
#'
#' @param object A `JoiNMeFit` object.
#' @param what Character component selector. One of
#'   `"fixed_effects"`, `"random_effects"`, `"coefficients"`,
#'   `"marker_weights"`, `"fixef"`,
#'   `"gamma_w"`, `"basehaz"`,
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
#'   - `posterior_draws`: numeric array when `keep_chains = TRUE` (iteration x chain x term),
#'     otherwise a numeric matrix (rows = draws, cols = requested terms). For
#'     `what = "association_plot"`, this is a named list of compact draw
#'     matrices keyed by association term.
#'   - `term_map`: data.frame mapping `term` to Stan `variable`
#'   - `support`: for `what = "association_plot"`, cached model-implied raw
#'     support ranges used by association plotting.
#'   - `transform_coeff_draws`: for `what = "association_plot"`,
#'     draw-specific fitted I-spline or ordered piecewise-linear ordinates,
#'     including one matrix per covariance component when applicable.
#' @seealso [posterior_draws()] for a higher-level interface that returns a single `posterior`
#' @export
extract.JoiNMeFit <- function(object,
                              what = c("fixed_effects", "random_effects", "coefficients", "marker_weights", "fixef", "gamma_w", "basehaz", "baseline_hazard", "assoc", "association_plot", "distributional", "distributional_regression", "likelihood_scale", "raw"),
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

  # These three selectors own the complete draw-level inputs used by the
  # posterior-summary layer and corresponding high-level coefficient methods. They deliberately return
  # scientific structures rather than Stan storage coordinates. The shorter
  # `what = "fixef"` selector below remains the longitudinal design-matrix
  # component and is useful when a rectangular posterior object is required.
  if (what %in% c("fixed_effects", "random_effects", "coefficients")) {
    structured_draws <- switch(
      what,
      fixed_effects = .extract_fixed_effect_draws(object, draws = draws, seed = seed),
      random_effects = .extract_random_effect_draws(object, draws = draws, seed = seed),
      coefficients = .extract_combined_coefficient_draws(object, draws = draws, seed = seed)
    )
    return(list(
      posterior_draws = structured_draws,
      term_map = NULL,
      metadata = list(
        what = what,
        scale = "scientific",
        keep_chains = FALSE
      )
    ))
  }

  # Marker weights are a derived posterior quantity rather than one rectangular
  # Stan variable. This selector makes their chain-preserving reconstruction
  # available without requiring callers to reach into package internals.
  if (identical(what, "marker_weights")) {
    active_terms <- .reported_marker_weight_terms(sd) # one representative term per fitted shared or term-specific set
    if (!is.null(term)) active_terms <- intersect(active_terms, as.character(term))
    if (!length(active_terms)) {
      cli::cli_abort("No matching marker-weight set is available for extraction.")
    }
    all_variables <- tryCatch(
      posterior::variables(.get_draws_obj(fit)),
      error = function(error) character(0)
    ) # posterior vocabulary used to select fitted rather than constant weights
    weight_arrays <- lapply(active_terms, function(term_key) {
      .association_marker_weight_array(
        object,
        term_key = term_key,
        draws = draws,
        seed = seed,
        all_vars = all_variables
      )
    })
    names(weight_arrays) <- active_terms
    if (!isTRUE(keep_chains)) {
      weight_arrays <- lapply(weight_arrays, function(weight_array) {
        posterior::as_draws_matrix(posterior::as_draws_array(weight_array))
      })
    }
    return(list(
      posterior_draws = if (length(weight_arrays) == 1L) weight_arrays[[1L]] else weight_arrays,
      term_map = data.frame(
        term = active_terms,
        variable = "effective_marker_weight",
        stringsAsFactors = FALSE
      ),
      metadata = list(what = what, keep_chains = isTRUE(keep_chains))
    ))
  }

  # These three selectors own the complete draw-level inputs used by the
  # posterior-summary layer and corresponding high-level coefficient methods. They deliberately return
  # scientific structures rather than Stan storage coordinates. The shorter
  # `what = "fixef"` selector remains the longitudinal design-matrix
  # component and is useful when a rectangular posterior object is required.
  if (what %in% c("fixed_effects", "random_effects", "coefficients")) {
    structured_draws <- switch(
      what,
      fixed_effects = .extract_fixed_effect_draws(object, draws = draws, seed = seed),
      random_effects = .extract_random_effect_draws(object, draws = draws, seed = seed),
      coefficients = .extract_combined_coefficient_draws(object, draws = draws, seed = seed)
    )
    return(list(
      posterior_draws = structured_draws,
      term_map = NULL,
      metadata = list(
        what = what,
        scale = "scientific",
        keep_chains = FALSE
      )
    ))
  }

  # Marker weights are a derived posterior quantity rather than one rectangular
  # Stan variable. This selector makes their chain-preserving reconstruction
  # available without requiring callers to reach into package internals.
  if (identical(what, "marker_weights")) {
    active_terms <- .reported_marker_weight_terms(sd) # one representative term per fitted shared or term-specific set
    if (!is.null(term)) active_terms <- intersect(active_terms, as.character(term))
    if (!length(active_terms)) {
      cli::cli_abort("No matching marker-weight set is available for extraction.")
    }
    all_variables <- tryCatch(
      posterior::variables(.get_draws_obj(fit)),
      error = function(error) character(0)
    ) # posterior vocabulary used to select fitted rather than constant weights
    weight_arrays <- lapply(active_terms, function(term_key) {
      .association_marker_weight_array(
        object,
        term_key = term_key,
        draws = draws,
        seed = seed,
        all_vars = all_variables
      )
    })
    names(weight_arrays) <- active_terms
    if (!isTRUE(keep_chains)) {
      weight_arrays <- lapply(weight_arrays, function(weight_array) {
        posterior::as_draws_matrix(posterior::as_draws_array(weight_array))
      })
    }
    return(list(
      posterior_draws = if (length(weight_arrays) == 1L) weight_arrays[[1L]] else weight_arrays,
      term_map = data.frame(
        term = active_terms,
        variable = "effective_marker_weight",
        stringsAsFactors = FALSE
      ),
      metadata = list(what = what, keep_chains = isTRUE(keep_chains))
    ))
  }

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
      posterior_draws = data$coeff_draws[keep_terms],
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
      tvals <- .event_ordinate_to_original_time(
        sd$S_event,
        tmax
      ) # endpoint labels in original study-time units, independent of the Surv variable's name
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
          coef_mat <- matrix(
            arr[, ch, k_idx, drop = FALSE],
            nrow = dim(arr)[1],
            ncol = length(k_idx)
          ) # iteration-by-basis coefficient matrix after removing the singleton chain dimension
          out[, ch, col_pos:(col_pos + nrow(b_event) - 1L)] <- exp(coef_mat %*% t(b_event)) / tmax
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
      return(list(posterior_draws = out, term_map = map))
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
    return(list(posterior_draws = out, term_map = map))
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
    posterior_draws = out,
    term_map = map
  )
}

#' Extract stored posterior draw components from dynamic prediction objects
#'
#' @description
#' Extracts stored prediction draw blocks without flattening them first.
#'
#' This is the structured companion to [posterior_draws.JoiNMeDynPred()]. Use
#' `extract()` when you want to keep the original prediction block semantics
#' (`longitudinal`, `survival`, `cumhaz`, random effects, and scale/id filters).
#' Use [posterior_draws()] when you want those blocks flattened into one
#' `posterior`-compatible draw object with composite variable labels.
#'
#' @param object A `JoiNMeDynPred` object.
#' @param what Draw block selector: `"longitudinal"`, `"longitudinal_fitted"`,
#'   `"survival"`, `"cumhaz"`, `"random_effects_id"`, `"random_effects_marker_id"`.
#' @param id Optional character/integer id filter.
#' @param scale Optional scale filter for longitudinal blocks (`epred`, `linpred`, `predict`).
#'
#' @return A list with fields:
#'   - `posterior_draws`: numeric matrix or list of matrices
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
    return(list(posterior_draws = dd, meta = list(what = what, ids = names(dd))))
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
    posterior_draws = mats,
    meta = list(what = what, ids = names(mats), scale = scale)
  )
}
