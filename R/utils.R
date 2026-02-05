#' @name joinme_utils
#' @title joinme Utility Functions
#'
#' @importFrom stats terms
#'
#' @description
#' Helper functions shared by standata builder, simulation, and methods.
#'
#' @keywords internal
NULL

suppressPackageStartupMessages({
  library(lme4)
  library(reformulas)
  library(splines)
  library(posterior)
})

#' @keywords internal
.gk15_nodes <- function() {
  c(
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
    -0.7415311855993945, -0.5860872354676911, -0.4058451513773972,
    -0.2077849550078985, 0.0, 0.2077849550078985,
    0.4058451513773972, 0.5860872354676911, 0.7415311855993945,
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126
  )
}

#' Safe model.matrix wrapper
#' @keywords internal
.mm <- function(formula, data) {
  X <- stats::model.matrix(formula, data = data)
  storage.mode(X) <- "double"
  X
}

#' Normalize distributional formula input
#' @keywords internal
.normalize_formula_dist <- function(formulaDist) {
  if (is.null(formulaDist)) return(list())
  if (inherits(formulaDist, "formula") || is.character(formulaDist)) {
    formulaDist <- list(formulaDist)
  }
  if (!is.list(formulaDist)) {
    cli::cli_abort(c(
      x = "{.arg formulaDist} must be a formula, character, or list of formulas.",
      i = "Example: list(sigma = ~ 1 + time, nu = nu ~ marker)."
    ))
  }
  if (is.null(names(formulaDist))) names(formulaDist) <- rep("", length(formulaDist))

  out <- list()
  for (i in seq_along(formulaDist)) {
    f <- formulaDist[[i]]
    if (!inherits(f, "formula")) f <- stats::as.formula(f)

    nm <- names(formulaDist)[i]
    has_name <- !is.null(nm) && nzchar(nm)
    has_lhs <- length(f) >= 3

    if (!has_lhs && has_name) {
      rhs_txt <- deparse(f[[2]])
      f <- stats::as.formula(paste(nm, "~", rhs_txt))
      has_lhs <- TRUE
    }
    if (!has_lhs) {
      cli::cli_abort(c(
        x = "Each distributional regression must have a parameter on the LHS unless the list entry is named.",
        i = "Example: list(sigma = ~ 1 + time) or list(sigma ~ 1 + time)."
      ))
    }

    lhs <- all.vars(f[[2]])
    if (length(lhs) != 1) {
      cli::cli_abort(c(
        x = "Distributional regression LHS must be a single parameter name.",
        i = "Allowed: sigma, nu, phi, alpha (or skew)."
      ))
    }
    param <- tolower(lhs)
    if (param == "skew") param <- "alpha"
    if (param %in% c("phi_beta", "beta_phi", "precision")) param <- "phi_beta"
    if (param %in% c("tau", "tau_sde", "skew_sde")) param <- "tau_sde"
    if (!param %in% c("sigma", "nu", "phi", "alpha", "phi_beta", "tau_sde")) {
      cli::cli_abort(c(
        x = "Unknown distributional parameter: {lhs}.",
        i = "Use one of: sigma, nu, phi, alpha (or skew), phi_beta, tau_sde."
      ))
    }

    if (has_name && tolower(nm) != param) {
      cli::cli_abort(c(
        x = "Distributional formula name '{nm}' does not match LHS '{param}'.",
        i = "Use names matching the LHS, or leave list entries unnamed."
      ))
    }
    out[[param]] <- f
  }

  out
}

#' Build distributional design matrix
#' @keywords internal
.build_dist_matrix <- function(formula, data) {
  if (is.null(formula)) {
    return(list(P = 0L, X = matrix(0.0, nrow(data), 0), cols = character(0)))
  }
  if (!inherits(formula, "formula")) formula <- stats::as.formula(formula)
  rhs <- stats::update(formula, . ~ .)
  rhs[[2]] <- NULL
  rhs <- reformulas::nobars(rhs)
  X <- .mm(rhs, data)
  list(P = ncol(X), X = X, cols = colnames(X))
}

#' Parse mixed-effects terms for distributional regression
#' @keywords internal
.build_dist_re_terms <- function(formula, data) {
  if (is.null(formula)) {
    return(list(n_re = 0L, K = integer(0), G = integer(0), Z = list(), J = list(), terms = character(0)))
  }
  if (!inherits(formula, "formula")) formula <- stats::as.formula(formula)

  f_exp <- reformulas::expandDoubleVerts(formula)
  bars <- reformulas::findbars(f_exp)
  if (length(bars) == 0) {
    return(list(n_re = 0L, K = integer(0), G = integer(0), Z = list(), J = list(), terms = character(0)))
  }

  Z_list <- list()
  J_list <- list()
  K_list <- integer(0)
  G_list <- integer(0)
  term_names <- character(0)

  for (bt in bars) {
    rhs_expr <- bt[[2]]
    grp_expr <- bt[[3]]

    if (any(grepl("\\|", deparse(rhs_expr)))) {
      cli::cli_abort(c(
        x = "Nested random effects are not supported in {.arg formulaDist}.",
        i = "Use flat terms like (1 + t | id) + (1 | region)."
      ))
    }

    rhs_formula <- stats::as.formula(paste0("~ ", paste(deparse(rhs_expr, width.cutoff = 500L), collapse = " ")))
    Z <- .mm(rhs_formula, data)

    grp_df <- stats::model.frame(
      stats::as.formula(paste0("~ ", paste(deparse(grp_expr, width.cutoff = 500L), collapse = " "))),
      data = data,
      drop.unused.levels = TRUE
    )
    if (ncol(grp_df) == 1) {
      grp <- grp_df[[1]]
    } else {
      grp <- interaction(grp_df, drop = TRUE)
    }
    grp <- factor(grp)
    J <- as.integer(grp)

    Z_list[[length(Z_list) + 1L]] <- Z
    J_list[[length(J_list) + 1L]] <- J
    K_list <- c(K_list, ncol(Z))
    G_list <- c(G_list, nlevels(grp))
    term_names <- c(term_names, paste0("(", paste(deparse(rhs_expr, width.cutoff = 500L), collapse = " "), "|", paste(deparse(grp_expr, width.cutoff = 500L), collapse = " "), ")"))
  }

  list(
    n_re = length(Z_list),
    K = K_list,
    G = G_list,
    Z = Z_list,
    J = J_list,
    terms = term_names
  )
}

#' Pad random-effects term matrices to max dimensions
#' @keywords internal
.pad_re_terms <- function(re_terms, n_rows) {
  if (re_terms$n_re == 0) {
    return(list(
      n_re = 0L,
      K = integer(0),
      G = integer(0),
      K_max = 0L,
      G_max = 0L,
      Z = list(),
      J = list(),
      J_mat = matrix(0L, nrow = 0, ncol = n_rows),
      terms = re_terms$terms
    ))
  }
  K_max <- max(re_terms$K)
  G_max <- max(re_terms$G)
  Z_padded <- lapply(seq_len(re_terms$n_re), function(i) {
    Z <- re_terms$Z[[i]]
    if (ncol(Z) < K_max) {
      Z <- cbind(Z, matrix(0.0, nrow(Z), K_max - ncol(Z)))
    }
    if (nrow(Z) != n_rows) {
      cli::cli_abort(c(
        x = "Random-effects design matrix row count mismatch.",
        i = "Check distributional regression design construction."
      ))
    }
    Z
  })
  J_mat <- do.call(rbind, lapply(re_terms$J, function(j) as.integer(j)))
  list(
    n_re = re_terms$n_re,
    K = re_terms$K,
    G = re_terms$G,
    K_max = K_max,
    G_max = G_max,
    Z = Z_padded,
    J = re_terms$J,
    J_mat = J_mat,
    terms = re_terms$terms
  )
}

#' Parse association include flags
#' @keywords internal
.parse_assoc <- function(assoc) {
  assoc <- unique(assoc)
  list(
    assoc_cv_total  = as.integer("cv_total" %in% assoc),
    assoc_cv_mean   = as.integer("cv_mean" %in% assoc),
    assoc_cv_marker = as.integer("cv_marker" %in% assoc),
    assoc_cs_total  = as.integer("cs_total" %in% assoc),
    assoc_cs_mean   = as.integer("cs_mean" %in% assoc),
    assoc_cs_marker = as.integer("cs_marker" %in% assoc),
    assoc_vcov      = as.integer("vcov" %in% assoc)
  )
}

#' Validate transform flags
#'
#' Transform codes: 0 id, 1 exp, 2 log, 3 inv_logit, 4 probit, 5 sqrt, 6 cbrt.
#' Restrictions: log and sqrt are restricted to cv_mean (legacy constraint).
#'
#' @keywords internal
.validate_tf <- function(tf) {
  nm <- c("cv_mean", "cv_marker", "cs_mean", "cs_marker", "vcov")
  tf_map <- c(
    identity = 0L,
    id = 0L,
    exp = 1L,
    log = 2L,
    inv_logit = 3L,
    probit = 4L,
    sqrt = 5L,
    cbrt = 6L
  )
  for (k in nm) if (is.null(tf[[k]])) tf[[k]] <- 0L
  for (k in nm) {
    if (is.character(tf[[k]])) {
      key <- tolower(tf[[k]])
      if (!key %in% names(tf_map)) {
        cli::cli_abort(c(
          x = "Unknown transform name for {k}: {tf[[k]]}.",
          i = "Use one of: {paste(names(tf_map), collapse = ', ')}."
        ))
      }
      tf[[k]] <- tf_map[[key]]
    }
    tf[[k]] <- as.integer(tf[[k]])
    if (!(tf[[k]] %in% 0:6)) {
      cli::cli_abort(c(
        x = "Transform flag {k} must be in 0..6.",
        i = "Use a supported name or integer code."
      ))
    }
  }
  bad <- c(tf$cv_marker, tf$cs_mean, tf$cs_marker, tf$vcov)
  if (any(bad %in% c(2L, 5L))) {
    cli::cli_abort(c(
      x = "log/sqrt transforms are restricted to cv_mean only.",
      i = "Use identity or other transforms for cs_mean/cs_marker/vcov."
    ))
  }
  tf
}

#' Detect whether RHS expression includes a time variable
#' @keywords internal
.expr_has_time <- function(expr, time_var) {
  f <- if (is.character(expr)) {
    stats::as.formula(paste0("~ ", expr))
  } else {
    stats::as.formula(paste0("~ ", paste(deparse(expr, width.cutoff = 500L), collapse = " ")))
  }
  time_var %in% all.vars(f)
}

#' Disable an association component with a warning
#' @keywords internal
.disable_assoc <- function(assoc, what, reason) {
  if (what %in% assoc) {
    warning(sprintf("Association '%s' disabled: %s", what, reason), call. = FALSE)
    assoc <- setdiff(assoc, what)
  }
  assoc
}

#' Convert lme4 bar terms into RHS-only formula list
#' @keywords internal
.bar_terms_to_rhs_list <- function(bar_terms) {
  lapply(bar_terms, function(bt) {
    expr <- bt[[2]]
    stats::as.formula(paste0("~ ", paste(deparse(expr, width.cutoff = 500L), collapse = " ")))
  })
}

#' Evaluate RHS list on an n_id x K matrix of times
#' @keywords internal
.eval_rhs_list_on_times <- function(rhs_list, dataEvent, time_var, times_mat) {
  n_id <- nrow(dataEvent)
  K <- ncol(times_mat)
  dd <- dataEvent[rep(seq_len(n_id), each = K), , drop = FALSE]
  dd[[time_var]] <- as.vector(t(times_mat))
  mats <- lapply(rhs_list, function(rhs) .mm(rhs, dd))
  if (length(mats) == 0) {
    return(array(0.0, dim = c(n_id, K, 0)))
  }
  Xbig <- do.call(cbind, mats)
  array(Xbig, dim = c(n_id, K, ncol(Xbig)))
}

#' Evaluate RHS list at event times
#' @keywords internal
.eval_rhs_list_at_event <- function(rhs_list, dataEvent, time_var, t_vec) {
  dd <- dataEvent
  dd[[time_var]] <- t_vec
  mats <- lapply(rhs_list, function(rhs) .mm(rhs, dd))
  if (length(mats) == 0) {
    return(matrix(0.0, nrow(dd), 0))
  }
  do.call(cbind, mats)
}

#' Center baseline hazard basis columns
#' @keywords internal
.center_baseline <- function(Bs_event_raw, Bs_gk_raw) {
  n_id <- nrow(Bs_event_raw)
  Kbs <- ncol(Bs_event_raw)
  big <- rbind(
    Bs_event_raw,
    matrix(aperm(Bs_gk_raw, c(1, 3, 2)), nrow = n_id * 15, ncol = Kbs)
  )
  colm <- colMeans(big)
  Bs_event_c <- sweep(Bs_event_raw, 2, colm, "-")
  Bs_gk_c <- sweep(Bs_gk_raw, MARGIN = 3, STATS = colm, FUN = "-")
  list(Bs_event_c = Bs_event_c, Bs_gk_c = Bs_gk_c, col_means = colm)
}

#' Construct a zero-dimension marker-only block
#' @keywords internal
.zero_marker_block <- function(N, n_id) {
  list(
    R_mk = 0L,
    Z_mk_obs = matrix(0.0, N, 0),
    Z_mk_gk_now = array(0.0, dim = c(n_id, 15, 0)),
    Z_mk_gk_fwd = array(0.0, dim = c(n_id, 15, 0)),
    Z_mk_event_now = matrix(0.0, n_id, 0),
    Z_mk_event_fwd = matrix(0.0, n_id, 0)
  )
}

#' Baseline hazard basis (bs/ns)
#' @keywords internal
.make_basehaz_basis <- function(x, basis = c("bs", "ns"), knots, degree = 3, boundary = c(0, 1)) {
  basis <- match.arg(basis)
  # Suppress warnings about x values beyond boundary knots
  # (This is expected in prediction when using Gauss-Kronrod quadrature near boundaries)
  if (basis == "bs") {
    suppressWarnings(splines::bs(x, knots = knots, Boundary.knots = boundary, degree = degree, intercept = FALSE))
  } else {
    suppressWarnings(splines::ns(x, knots = knots, Boundary.knots = boundary, intercept = FALSE))
  }
}

#' Extract nested marker syntax terms
#'
#' Outer: ( ... | marker )
#' Inner: ( ... | id ) terms inside outer define marker-by-id basis (optional).
#' These terms create id-specific random intercepts/slopes per marker when present.
#' Remaining no-bars part defines marker-only terms (optional).
#'
#' @keywords internal
.extract_nested_marker_terms <- function(formulaLong, marker_var, id_var) {
  f_exp <- reformulas::expandDoubleVerts(formulaLong)
  bars <- reformulas::findbars(f_exp)
  if (length(bars) == 0) {
    return(list(mk_rhs_list = list(), idm_rhs_list = list(), marker_terms = list()))
  }
  
  # Extract grouping variable name - handle both symbols and compound expressions
  grp <- character(length(bars))
  for (i in seq_along(bars)) {
    b <- bars[[i]]
    grp_expr <- b[[3]]
    if (is.call(grp_expr) && as.character(grp_expr)[[1]] == ":") {
      # Handle marker:id or marker*id patterns - extract first part
      grp[i] <- deparse(grp_expr[[2]], width.cutoff = 500L)[1]
    } else if (is.call(grp_expr) && as.character(grp_expr)[[1]] == "*") {
      # Handle marker*id patterns - extract first part
      grp[i] <- deparse(grp_expr[[2]], width.cutoff = 500L)[1]
    } else {
      # Simple symbol
      grp[i] <- as.character(grp_expr)[[1]]
    }
  }
  
  mk_idx <- which(grp == marker_var)
  if (length(mk_idx) == 0) {
    return(list(mk_rhs_list = list(), idm_rhs_list = list(), marker_terms = list()))
  }

  mk_rhs_list <- list()
  idm_rhs_list <- list()
  marker_terms <- list()

  for (j in mk_idx) {
    bt <- bars[[j]]
    outer_expr <- bt[[2]]
    marker_terms <- c(marker_terms, list(outer_expr))

    f_inner <- stats::as.formula(paste0("~ ", paste(deparse(outer_expr, width.cutoff = 500L), collapse = " ")))
    f_inner <- reformulas::expandDoubleVerts(f_inner)

    inner_bars <- reformulas::findbars(f_inner)
    if (length(inner_bars) > 0) {
      inner_grp <- vapply(inner_bars, function(b) as.character(b[[3]]), character(1))
      id_inner_idx <- which(inner_grp == id_var)
      if (length(id_inner_idx) > 0) {
        idm_rhs_list <- c(idm_rhs_list, .bar_terms_to_rhs_list(inner_bars[id_inner_idx]))
      }
    }

    f_mk <- reformulas::nobars(f_inner)
    rhs_terms <- terms(f_mk)
    term_labels <- attr(rhs_terms, "term.labels")
    has_only_intercept <- (length(term_labels) == 0) && (attr(rhs_terms, "intercept") == 1)

    outer_txt <- paste(deparse(outer_expr, width.cutoff = 500L), collapse = " ")
    user_explicit_intercept <- grepl("(\\+\\s*1\\s*(\\+|$))", outer_txt)

    if (has_only_intercept && !user_explicit_intercept) next

    rhs_txt <- paste(deparse(f_mk[[2]], width.cutoff = 500L), collapse = " ")
    rhs_txt <- trimws(rhs_txt)
    if (rhs_txt == "0") next

    mk_rhs_list <- c(mk_rhs_list, list(f_mk))
  }

  list(mk_rhs_list = mk_rhs_list, idm_rhs_list = idm_rhs_list, marker_terms = marker_terms)
}

# ---- NEW: time-index metadata for internal scaling in Stan -----------------

#' Detect indices of time-related columns in a design matrix
#'
#' @param colnames_vec Character vector of column names.
#' @param time_var Time variable name.
#' @return Integer vector of indices.
#' @keywords internal
.detect_time_cols <- function(colnames_vec, time_var) {
  if (is.null(colnames_vec) || length(colnames_vec) == 0) {
    return(integer(0))
  }
  idx <- which(colnames_vec == time_var)
  if (length(idx) == 0) {
    pat <- paste0("(^", time_var, "$)|(^", time_var, "\\b)|\\b", time_var, "\\b")
    idx <- grep(pat, colnames_vec)
  }
  as.integer(idx)
}

#' Build time-index metadata list for Stan
#'
#' @return named list with n_time_* and idx_time_* values.
#' @keywords internal
.make_time_index_metadata <- function(x_cols, zid_cols, zmk_cols, zidm_cols, time_var) {
  idx_beta <- .detect_time_cols(x_cols, time_var)
  idx_uid <- .detect_time_cols(zid_cols, time_var)
  idx_vmk <- .detect_time_cols(zmk_cols, time_var)
  idx_widm <- .detect_time_cols(zidm_cols, time_var)

  list(
    n_time_beta = length(idx_beta),
    idx_time_beta = idx_beta,
    n_time_uid = length(idx_uid),
    idx_time_uid = idx_uid,
    n_time_vmk = length(idx_vmk),
    idx_time_vmk = idx_vmk,
    n_time_widm = length(idx_widm),
    idx_time_widm = idx_widm
  )
}

# ---- Stan engine helpers ---------------------------------------------------

#' Resolve preferred Stan engine
#' @keywords internal
.resolve_stan_engine <- function(engine = NULL) {
  resolved <- engine %||% getOption("stan_preferred_engine", "cmdstanr")
  resolved <- tolower(as.character(resolved))
  if (!resolved %in% c("cmdstanr", "rstan")) {
    cli::cli_abort(c(
      x = "Unknown Stan engine: {resolved}.",
      i = "Use 'cmdstanr' or 'rstan' via options(stan_preferred_engine=...)."
    ))
  }
  resolved
}

#' Detect CmdStanR fit
#' @keywords internal
.is_cmdstanr_fit <- function(fit) {
  inherits(fit, "CmdStanMCMC")
}

#' Detect rstan fit
#' @keywords internal
.is_rstan_fit <- function(fit) {
  inherits(fit, "stanfit")
}

#' Unified draws object helper
#' @keywords internal
.get_draws_obj <- function(fit, variables = NULL, draws = NULL, seed = 1) {
  if (.is_cmdstanr_fit(fit)) {
    d <- fit$draws(variables = variables)
  } else if (.is_rstan_fit(fit)) {
    d <- posterior::as_draws_array(fit)
    if (!is.null(variables)) {
      vars_avail <- intersect(variables, posterior::variables(d))
      d <- posterior::subset_draws(d, variable = vars_avail)
    }
  } else {
    cli::cli_abort("Unsupported Stan fit object; expected CmdStanR or rstan.")
  }
  if (!is.null(draws) && is.finite(draws)) {
    nd <- posterior::ndraws(d)
    if (draws < nd) {
      set.seed(seed)
      idx <- sample.int(nd, size = draws)
      d <- posterior::subset_draws(d, draw = idx)
    }
  }
  d
}

#' Unified draws matrix helper
#' @keywords internal
.get_draws_matrix <- function(fit, variables = NULL, draws = NULL, seed = 1) {
  d <- .get_draws_obj(fit, variables = variables, draws = draws, seed = seed)
  posterior::as_draws_matrix(d)
}

# ---- Posterior draw helpers for methods -----------------------------------

#' Subset draws safely before conversion to draws_df
#' @keywords internal
.get_draws_df <- function(fit, variables, draws = NULL, seed = 1) {
  d <- .get_draws_obj(fit, variables = variables, draws = draws, seed = seed)
  posterior::as_draws_df(d)
}

#' Summarize a numeric vector of draws
#' @keywords internal
.summarize_draw_col <- function(x) {
  c(
    Estimate = mean(x),
    Est.Error = stats::sd(x),
    Q2.5 = stats::quantile(x, 0.025, names = FALSE),
    Q97.5 = stats::quantile(x, 0.975, names = FALSE)
  )
}

#' Clean summary wrapper
#' @keywords internal
.summarize_draws_cleaned <- function(fit, variables, draws = NULL, seed = 1) {
  ddf <- .get_draws_df(fit, variables = variables, draws = draws, seed = seed)
  out <- lapply(variables, function(v) .summarize_draw_col(ddf[[v]]))
  mat <- do.call(rbind, out)
  df <- data.frame(variable = variables, mat, row.names = NULL, check.names = FALSE)
  df
}

#' TF name helper
#' @keywords internal
.tf_name <- function(code) {
  map <- c("identity", "exp", "log", "inv_logit", "probit", "sqrt", "cbrt")
  if (code >= 0 && code <= 6) map[code + 1] else "unknown"
}

# ---- Stan data parsing (rstan compatibility) ------------------------------

#' Expand Stan includes
#' @keywords internal
.read_stan_with_includes <- function(file, visited = character()) {
  if (!file.exists(file) || file %in% visited) return(character(0))
  visited <- c(visited, file)
  lines <- readLines(file, warn = FALSE)
  out <- character(0)
  for (line in lines) {
    if (grepl("^\\s*#include\\s+", line)) {
      inc <- sub("^\\s*#include\\s+\"?([^\"\\s]+)\"?.*$", "\\1", line)
      inc_path <- file.path(dirname(file), inc)
      if (file.exists(inc_path)) {
        out <- c(out, .read_stan_with_includes(inc_path, visited = visited))
      }
    } else {
      out <- c(out, line)
    }
  }
  out
}

#' Extract data block variable names from Stan file
#' @keywords internal
.stan_data_names <- function(stan_file) {
  lines <- .read_stan_with_includes(stan_file)
  if (length(lines) == 0) return(character(0))
  data_start <- grep("^\\s*data\\s*\\{", lines)
  if (length(data_start) == 0) return(character(0))
  start_idx <- data_start[1]
  depth <- 1L
  vars <- character(0)
  count_char <- function(x, ch) {
    hits <- gregexpr(ch, x, fixed = TRUE)
    if (length(hits) == 1 && hits[[1]][1] == -1) return(0L)
    sum(vapply(hits, length, integer(1)))
  }
  for (i in seq.int(start_idx + 1L, length(lines))) {
    line <- lines[i]
    line_clean <- gsub("//.*$", "", line)
    if (nchar(trimws(line_clean)) == 0) {
      depth <- depth + count_char(line, "{") - count_char(line, "}")
      if (depth <= 0) break
      next
    }
    if (grepl(";", line_clean)) {
      decl <- sub(";.*$", "", line_clean)
      decl <- trimws(decl)
      if (nzchar(decl)) {
        var <- sub(".*\\b([A-Za-z_][A-Za-z0-9_]*)\\s*(\\[.*\\])?$", "\\1", decl)
        if (nzchar(var)) vars <- c(vars, var)
      }
    }
    depth <- depth + count_char(line, "{") - count_char(line, "}")
    if (depth <= 0) break
  }
  unique(vars)
}

#' Coerce distributional regression arrays for rstan
#' @keywords internal
.coerce_rstan_dist_arrays <- function(sd) {
  make_num_array <- function(x, dims) {
    if (length(dims) == 0) return(x)
    if (is.null(x) || (is.list(x) && length(x) == 0)) {
      return(array(0, dim = dims))
    }
    if (is.list(x)) {
      arr <- array(unlist(x), dim = dims)
      storage.mode(arr) <- "double"
      return(arr)
    }
    if (is.array(x)) {
      storage.mode(x) <- "double"
      return(x)
    }
    arr <- array(x, dim = dims)
    storage.mode(arr) <- "double"
    arr
  }
  make_int_array <- function(x, dims) {
    if (length(dims) == 0) return(x)
    if (is.null(x) || (is.list(x) && length(x) == 0)) {
      return(array(1L, dim = dims))
    }
    if (is.list(x)) {
      arr <- array(unlist(x), dim = dims)
      storage.mode(arr) <- "integer"
      return(arr)
    }
    if (is.array(x)) {
      storage.mode(x) <- "integer"
      return(x)
    }
    arr <- array(x, dim = dims)
    storage.mode(arr) <- "integer"
    arr
  }

  N <- sd$N %||% 0L
  sd$Z_sigma <- make_num_array(sd$Z_sigma, c(sd$n_re_sigma %||% 0L, N, sd$K_sigma_max %||% 0L))
  sd$J_sigma <- make_int_array(sd$J_sigma, c(sd$n_re_sigma %||% 0L, N))
  sd$Z_nu <- make_num_array(sd$Z_nu, c(sd$n_re_nu %||% 0L, N, sd$K_nu_max %||% 0L))
  sd$J_nu <- make_int_array(sd$J_nu, c(sd$n_re_nu %||% 0L, N))
  sd$Z_phi <- make_num_array(sd$Z_phi, c(sd$n_re_phi %||% 0L, N, sd$K_phi_max %||% 0L))
  sd$J_phi <- make_int_array(sd$J_phi, c(sd$n_re_phi %||% 0L, N))
  sd$Z_alpha <- make_num_array(sd$Z_alpha, c(sd$n_re_alpha %||% 0L, N, sd$K_alpha_max %||% 0L))
  sd$J_alpha <- make_int_array(sd$J_alpha, c(sd$n_re_alpha %||% 0L, N))
  sd$Z_phi_beta <- make_num_array(sd$Z_phi_beta, c(sd$n_re_phi_beta %||% 0L, N, sd$K_phi_beta_max %||% 0L))
  sd$J_phi_beta <- make_int_array(sd$J_phi_beta, c(sd$n_re_phi_beta %||% 0L, N))
  sd$Z_tau_sde <- make_num_array(sd$Z_tau_sde, c(sd$n_re_tau_sde %||% 0L, N, sd$K_tau_sde_max %||% 0L))
  sd$J_tau_sde <- make_int_array(sd$J_tau_sde, c(sd$n_re_tau_sde %||% 0L, N))
  sd
}

#' Coerce time index fields for rstan
#' @keywords internal
.coerce_rstan_time_indices <- function(sd) {
  idx_names <- c("idx_time_beta", "idx_time_uid", "idx_time_vmk", "idx_time_widm")
  n_names <- c("n_time_beta", "n_time_uid", "n_time_vmk", "n_time_widm")
  for (i in seq_along(idx_names)) {
    idx <- sd[[idx_names[i]]] %||% integer(0)
    idx <- as.integer(idx)
    sd[[n_names[i]]] <- length(idx)
    sd[[idx_names[i]]] <- array(idx, dim = c(length(idx)))
  }
  sd
}

#' Coerce vector fields to 1D arrays for rstan
#' @keywords internal
.coerce_rstan_vectors <- function(sd, fields) {
  for (nm in fields) {
    x <- sd[[nm]]
    if (is.null(x)) next
    if (is.atomic(x) && is.null(dim(x))) {
      sd[[nm]] <- array(x, dim = c(length(x)))
    }
  }
  sd
}

#' Family name helper
#' @keywords internal
.family_name <- function(code) {
  map <- c("gaussian", "student", "bernoulli", "binomial", "poisson", "negbin")
  if (code > 0 && code <= 6) map[code] else "unknown"
}
