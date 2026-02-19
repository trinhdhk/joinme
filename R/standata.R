#' Build Stan data for joinme (time-internal scaling)
#'
#' @importFrom stats predict setNames
#' @importFrom dplyr %>% .data all_of
#'
#' @description
#' Standata builder for the joinme joint model (multivariate longitudinal + survival).
#'
#' This function follows the notebook `joinme_standata()` design, with these updates:
#'
#' 1. Time is still scaled in the design matrices by tmax for numerical stability.
#' 2. The Stan program internally rescales time-related coefficients using `tmax` and
#'    index vectors indicating which columns correspond to time.
#' 3. Therefore, this builder computes and supplies:
#'    - n_time_beta, idx_time_beta
#'    - n_time_uid,  idx_time_uid
#'    - n_time_vmk,  idx_time_vmk
#'    - n_time_widm, idx_time_widm
#'
#' @details
#' This builder performs three key steps:
#' 1. Constructs fixed-effect and random-effect design matrices on a scaled time axis.
#' 2. Creates distributional regression matrices (including optional random effects).
#' 3. Assembles spline bases and Gauss–Kronrod nodes for survival integration.
#'
#' Time is scaled internally as `t_scaled = t/t_max` for numerical stability; indices
#' are recorded to rescale time-associated coefficients back to the original units
#' inside Stan.
#'
#' @param formulaLong Longitudinal formula defining fixed effects, id-level effects,
#'   and the marker block. The marker block may optionally include an inner
#'   `( ... | id )` term for marker-by-id random effects. When the inner term is
#'   omitted, marker-by-id random effects are disabled (Q_idm = 0). Use `||` to
#'   enforce diagonal random-effect covariance (e.g., `(1 + time || id)`).
#' @param dataLong Long-format longitudinal data with columns for id, marker, time,
#'   outcome, and covariates referenced in `formulaLong`.
#' @param formulaEvent Survival formula for baseline covariates and event model.
#' @param dataEvent One row per id event data with event time, event indicator, and
#'   covariates referenced in `formulaEvent`.
#' @param formulaVcov Covariance regression formula for id-specific marker-by-id effects.
#'   If the marker block omits the inner `( ... | id )`, then marker-by-id effects
#'   are absent and `vcov` associations are not allowed.
#' @param formulaDist Optional list of formulas for distributional regression. Two
#'   forms are supported:
#'   1. Named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#'   2. Unnamed list with LHS parameter names, e.g. `list(sigma ~ 1 + time)`.
#'   Random-effects terms with `|` are supported; nested random-effects formulas are not.
#' @param id_var Column name for subject id in both longitudinal and event data.
#' @param marker_var Column name for marker/biomarker id in longitudinal data.
#' @param time_var Column name for longitudinal time in `dataLong`.
#' @param y_var Column name for longitudinal outcome in `dataLong`.
#' @param event_time_var Column name for event/censoring time in `dataEvent`.
#' @param event_var Column name for event indicator in `dataEvent`.
#' @param event_type_var Optional column name for competing-risk event type in `dataEvent`.
#' @param stage_from_var Optional column name for multi-stage "from" state in `dataEvent`.
#' @param stage_to_var Optional column name for multi-stage "to" state in `dataEvent`.
#' @param eps_fd Positive finite-difference step for association derivatives.
#' @param assoc Character vector specifying association components
#'   (e.g., "cv_mean", "cs_total", "vcov").
#' @param families Optional marker-specific family names. If NULL, all markers
#'   use Gaussian responses.
#' @param transforms Optional list specifying transformations for association terms.
#'   Each element (cv_total, cs_total, vcov) is a list with a `type` and fields
#'   required by that type (see `build_standata_transforms()`).
#' @param beta_prior Prior specification for longitudinal fixed effects.
#' @param alpha_prior Prior specification for association parameters.
#' @param lkj_prior Prior specification for correlation structures.
#' @param allow_marker_crosscorr Integer flag; 1 allows cross-marker correlation in marker RE.
#' @param shrinkage Integer flag controlling shrinkage behavior for marker-by-id effects.
#' @param marker_weights Optional numeric vector of length D giving base weights
#'   for marker-specific association components. These weights are applied to both
#'   current value (CV) and current slope (CS) marker summaries. Weights may be
#'   positive or negative. Internally, base weights are rescaled to unit RMS
#'   (`sqrt(mean(w^2)) = 1`) for numerical stability and better identifiability.
#'   When NULL, equal base weights are used.
#' @param estimate_marker_weights Logical; if TRUE, estimate marker weights via a
#'   signed additive perturbation around base weights. The effective weights are
#'   RMS-stabilized so composition can change while overall scale remains stable.
#'   Latent perturbations follow Normal or Laplace priors depending on `shrinkage`.
#' @param marker_weight_scale Non-negative prior scale for marker-weight shrinkage.
#' @param flag_resid_dim Integer flag to include residual dimension checks.
#' @param basehaz Baseline hazard basis type: "bs", "ns", or "formula".
#' @param n_knots Number of internal knots for spline baseline hazards.
#' @param basehaz_degree Degree of spline basis for baseline hazard.
#' @param basehaz_formula Formula for baseline hazard when `basehaz = "formula"`.
#' @param tau_spline Prior scale for spline coefficients (penalized spline).
#' @param vcov_diag_link Link for covariance regression diagonals: "softplus" or "exp".
#' @param tau_sde_fixed Optional fixed tau for skew-double-exponential (0 < tau < 1).
#' @param seed Optional random seed for deterministic components of standata.
#' @export
# File overview:
# - Parse formulas, build design matrices, and scale time.
# - Assemble longitudinal, survival, and association data blocks for Stan.
joinme_standata <- function(
  formulaLong,
  dataLong,
  formulaEvent,
  dataEvent,
  formulaVcov = ~1,
  formulaDist = NULL,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  y_var = "y",
  event_time_var = "time",
  event_var = "event",
  event_type_var = "event_type",
  stage_from_var = "state_from",
  stage_to_var = "state_to",
  eps_fd = 1e-2,
  assoc = c("cv_mean"),
  families = NULL,
  transforms = NULL,
  beta_prior = NULL,
  alpha_prior = NULL,
  lkj_prior = NULL,
  allow_marker_crosscorr = 1L,
  shrinkage = 2L,
  marker_weights = NULL,
  estimate_marker_weights = TRUE,
  marker_weight_scale = 1,
  flag_resid_dim = 0L,
  basehaz = c("bs", "ns", "formula"),
  n_knots = 5L,
  basehaz_degree = 3L,
  basehaz_formula = ~ 1 + time,
  tau_spline = 0.4,
  vcov_diag_link = c("softplus", "exp"),
  tau_sde_fixed = NULL,
  seed = NULL
) {
  assertthat::assert_that(is.data.frame(dataLong), msg = "dataLong must be a data.frame")
  assertthat::assert_that(is.data.frame(dataEvent), msg = "dataEvent must be a data.frame")
  basehaz <- match.arg(basehaz)
  vcov_diag_link <- match.arg(vcov_diag_link)

  # Workflow: parse formulas -> build matrices -> assemble survival basis -> pack list
  # Handle mixed families
  priors <- .build_priors(beta_prior = beta_prior, alpha_prior = alpha_prior, lkj_prior = lkj_prior)
  
  # Parse lme4 bars
  f_exp <- reformulas::expandDoubleVerts(formulaLong)
  bars <- reformulas::findbars(f_exp)
  f_fix <- reformulas::nobars(f_exp)
  if (length(bars) == 0) {
     cli::cli_abort(c(
      i = "formulaLong must include (...| {.arg id_var }) and a marker term ( ... | {.arg marker_var}.",
      x = "Missing individual random-effects terms."
    ))
  }

  grp <- vapply(bars, function(b) as.character(b[[3]]), character(1))
  id_idx <- which(grp == id_var)
  if (length(id_idx) == 0) {
    cli::cli_abort(c(
      i = "formulaLong must include (...| {.arg id_var }).",
      x = "Missing random-effects term for {.arg id_var}."
    ))
  }

  nested <- .extract_nested_marker_terms(formulaLong, marker_var = marker_var, id_var = id_var)
  mk_rhs_list <- nested$mk_rhs_list
  idm_rhs_list <- nested$idm_rhs_list
  if (length(nested$marker_terms) == 0) {
    cli::cli_abort(c(
      i = "formulaLong must include a marker block ( ... | {.arg marker_var} ).",
      x = "Missing random-effects term for {.arg marker_var}."
    ))
  }
  has_idm <- length(idm_rhs_list) > 0
  if (!has_idm && "vcov" %in% assoc) {
    cli::cli_abort(c(
      x = "Association {.arg vcov} requires marker-by-id random effects.",
      i = "Add an inner ( ... | {.arg id_var} ) term inside the marker block, or remove {.arg vcov} from {.arg assoc}."
    ))
  }

  # Fixed RHS only
  fixed_rhs <- stats::update(f_fix, . ~ .)
  fixed_rhs[[2]] <- NULL

  # Map ids by event data order
  ids <- sort(unique(dataEvent[[id_var]]))
  id_map <- setNames(seq_along(ids), ids)
  dataEvent <- dataEvent[match(ids, dataEvent[[id_var]]), , drop = FALSE]
  dataEvent$id_int <- as.integer(id_map[as.character(dataEvent[[id_var]])])
  dataLong$id_int <- as.integer(id_map[as.character(dataLong[[id_var]])])

  dataLong <- dataLong |>
    dplyr::filter(
      !is.na(.data$id_int),
      is.finite(.data[[time_var]]),
      is.finite(.data[[y_var]]),
      !is.na(.data[[marker_var]])
    ) |>
    dplyr::arrange(.data$id_int, .data[[marker_var]], .data[[time_var]])

  dataLong[[marker_var]] <- factor(dataLong[[marker_var]])
  marker_levels <- levels(dataLong[[marker_var]])
  D <- length(marker_levels)
  dataLong$marker_int <- as.integer(dataLong[[marker_var]])

  # Marker weights (for weighted means in CV/CS association components)
  if (is.null(marker_weights)) {
    marker_weights <- rep(1, D)
  } else {
    if (!is.numeric(marker_weights)) {
      cli::cli_abort(c(
        x = "{.arg marker_weights} must be numeric.",
        i = "Provide a numeric vector with length D or named by marker levels."
      ))
    }
    if (!is.null(names(marker_weights))) {
      marker_weights <- marker_weights[marker_levels]
    }
    if (length(marker_weights) != D) {
      cli::cli_abort(c(
        x = "{.arg marker_weights} must have length D={D}.",
        i = "Use one weight per marker level."
      ))
    }
    if (any(!is.finite(marker_weights))) {
      cli::cli_abort(c(
        x = "{.arg marker_weights} must be finite.",
        i = "Weights may be positive or negative, but cannot be NA/Inf."
      ))
    }
    if (all(abs(marker_weights) < 1e-12)) {
      cli::cli_abort(c(
        x = "{.arg marker_weights} cannot be all zeros.",
        i = "Provide at least one non-zero weight."
      ))
    }
  }

  # Normalize base weights to unit RMS so marker-side scale is identifiable.
  # This preserves signs and relative patterns while keeping the global magnitude stable.
  mw_rms <- sqrt(mean(marker_weights^2))
  if (!is.finite(mw_rms) || mw_rms <= 0) {
    cli::cli_abort(c(
      x = "Failed to scale {.arg marker_weights}.",
      i = "Provide finite weights with at least one non-zero value."
    ))
  }
  marker_weights <- marker_weights / mw_rms

  # Marker-weight estimation controls
  # - estimate_marker_weights toggles whether shrinkage is applied in Stan
  estimate_marker_weights <- isTRUE(estimate_marker_weights)
    # Independence flags from lme4-style double-bar syntax
    # - (... || id) enforces diagonal id RE covariance
    # - (... || marker) enforces diagonal marker RE covariance
    # - Nested (... || id) inside marker block enforces diagonal marker-by-id covariance
    indep_flags <- .resolve_re_independence(formulaLong, marker_var = marker_var, id_var = id_var)
  if (!is.numeric(marker_weight_scale) || length(marker_weight_scale) != 1 || !is.finite(marker_weight_scale) || marker_weight_scale < 0) {
    cli::cli_abort(c(
      x = "{.arg marker_weight_scale} must be a single non-negative numeric value.",
      i = "Example: marker_weight_scale = 1."
    ))
  }
  
  # Process mixed families (optional)
  # - family_codes drive distributional parameter availability
  if (!is.null(families)) {
    # If only length 1, duplicate to all markers
    if (length(families) == 1) {
      families <- rep(families, D)
    }
    family_codes <- .validate_family_list(families, D, dataLong, marker_var, y_var)
  } else {
    family_codes <- rep(1L, D)  # Default: all gaussian
  }
  family_names <- vapply(family_codes, .family_code_to_name, character(1))

  # Family-level distributional parameter indexing.
  #
  # Goal:
  # - for each distributional parameter (sigma, nu, phi, alpha, phi_beta,
  #   tau_sde), build one shared parameter per unique family that requires it,
  # - map each marker to its family-level parameter index,
  # - use index 0 for markers/families that do not use that parameter.
  .build_family_param_index <- function(param_name) {
    need_param <- vapply(family_codes, function(fc) {
      param_name %in% .family_distrib_params(fc)
    }, logical(1))
    fam_codes <- sort(unique(as.integer(family_codes[need_param])))
    map <- integer(D)
    if (length(fam_codes) > 0) {
      for (d in seq_len(D)) {
        if (need_param[d]) {
          map[d] <- match(as.integer(family_codes[d]), fam_codes)
        }
      }
    }
    list(
      n = as.integer(length(fam_codes)),
      marker_to = as.integer(map),
      family_codes = as.integer(fam_codes),
      family_names = if (length(fam_codes) > 0) {
        vapply(fam_codes, .family_code_to_name, character(1))
      } else {
        character(0)
      }
    )
  }

  fam_sigma <- .build_family_param_index("sigma")
  fam_nu <- .build_family_param_index("nu")
  fam_phi <- .build_family_param_index("phi")
  fam_alpha <- .build_family_param_index("alpha")
  fam_phi_beta <- .build_family_param_index("phi_beta")
  fam_tau_sde <- .build_family_param_index("tau_sde")

  vcov_diag_link_code <- if (vcov_diag_link == "exp") 1L else 0L

  use_tau_sde_fixed <- 0L
  tau_sde_fixed_value <- 0.5
  if (!is.null(tau_sde_fixed)) {
    if (!is.numeric(tau_sde_fixed) || length(tau_sde_fixed) != 1 || !is.finite(tau_sde_fixed)) {
      cli::cli_abort(c(
        x = "{.arg tau_sde_fixed} must be a single finite numeric value.",
        i = "Provide a number strictly between 0 and 1."
      ))
    }
    tau_sde_fixed_value <- as.numeric(tau_sde_fixed)
    if (tau_sde_fixed_value <= 0 || tau_sde_fixed_value >= 1) {
      cli::cli_abort(c(
        x = "{.arg tau_sde_fixed} must be between 0 and 1.",
        i = "Provide a number strictly between 0 and 1."
      ))
    }
    if (any(family_codes == 9L)) {
      use_tau_sde_fixed <- 1L
    } else {
      cli::cli_warn(c(
        x = "{.arg tau_sde_fixed} is ignored because no skew_double_exponential family is present.",
        i = "Remove {.arg tau_sde_fixed} or include the skew_double_exponential family."
      ))
    }
  }
  
  dist_formulas <- .normalize_formula_dist(formulaDist)
  family_names_present <- vapply(family_codes, .family_code_to_name, character(1))
  .validate_dist_formula_scopes(dist_formulas, family_names_present)
  allowed_dist <- unique(unlist(lapply(family_codes, .family_distrib_params)))
  if (length(dist_formulas) > 0) {
    bad <- setdiff(names(dist_formulas), allowed_dist)
    if (length(bad) > 0) {
      cli::cli_abort(c(
        x = "Distributional regression specified for unsupported parameters: {paste(bad, collapse = ', ')}.",
        i = "Allowed for the selected families: {paste(allowed_dist, collapse = ', ')}."
      ))
    }
  }
  tf_compositions <- list(
    tf_cv_mean_comp = list(codes = c(1L, rep(0L, 6)), n_codes = 1L),
    tf_cv_marker_comp = list(codes = c(1L, rep(0L, 6)), n_codes = 1L),
    tf_cs_mean_comp = list(codes = c(1L, rep(0L, 6)), n_codes = 1L),
    tf_cs_marker_comp = list(codes = c(1L, rep(0L, 6)), n_codes = 1L)
  )

  # Scale time to [0,1] by max event time
  # - tmax is also returned for coefficient rescaling in Stan
  tmax <- max(dataEvent[[event_time_var]])
  if (!is.finite(tmax) || tmax <= 0) stop("Invalid max event time.")
  dataLong$t_scaled <- dataLong[[time_var]] / tmax
  dataEvent$S_scaled <- dataEvent[[event_time_var]] / tmax

  # Build obs matrices on scaled time
  # - distributional regression matrices follow the same scaled timeline
  dl <- dataLong
  dl[[time_var]] <- dl$t_scaled
  family_by_row <- .family_by_row_from_marker(
    marker_values = dl[[marker_var]],
    marker_levels = marker_levels,
    family_codes = family_codes
  )
  attr(dl, "joinme_family_by_row") <- family_by_row

  dist_sigma <- .build_dist_matrix(dist_formulas$sigma, dl, family_by_row = family_by_row)
  dist_nu <- .build_dist_matrix(dist_formulas$nu, dl, family_by_row = family_by_row)
  dist_phi <- .build_dist_matrix(dist_formulas$phi, dl, family_by_row = family_by_row)
  dist_alpha <- .build_dist_matrix(dist_formulas$alpha, dl, family_by_row = family_by_row)
  dist_phi_beta <- .build_dist_matrix(dist_formulas$phi_beta, dl, family_by_row = family_by_row)
  dist_tau_sde <- .build_dist_matrix(dist_formulas$tau_sde, dl, family_by_row = family_by_row)

  re_sigma <- .pad_re_terms(.build_dist_re_terms(dist_formulas$sigma, dl), nrow(dl))
  re_nu <- .pad_re_terms(.build_dist_re_terms(dist_formulas$nu, dl), nrow(dl))
  re_phi <- .pad_re_terms(.build_dist_re_terms(dist_formulas$phi, dl), nrow(dl))
  re_alpha <- .pad_re_terms(.build_dist_re_terms(dist_formulas$alpha, dl), nrow(dl))
  re_phi_beta <- .pad_re_terms(.build_dist_re_terms(dist_formulas$phi_beta, dl), nrow(dl))
  re_tau_sde <- .pad_re_terms(.build_dist_re_terms(dist_formulas$tau_sde, dl), nrow(dl))

  X_obs <- .mm(fixed_rhs, dl)
  P <- ncol(X_obs)

  # id block
  # - always required (subject random effects)
  id_rhs_list <- .bar_terms_to_rhs_list(bars[id_idx])
  Z_id_obs <- do.call(cbind, lapply(id_rhs_list, function(rhs) .mm(rhs, dl)))
  R_id <- ncol(Z_id_obs)

  # marker-only optional
  # - empty marker block yields zero columns
  if (length(mk_rhs_list) == 0) {
    mk0 <- .zero_marker_block(N = nrow(dl), n_id = nrow(dataEvent))
    R_mk <- mk0$R_mk
    Z_mk_obs <- mk0$Z_mk_obs
  } else {
    Z_mk_obs <- do.call(cbind, lapply(mk_rhs_list, function(rhs) .mm(rhs, dl)))
    R_mk <- ncol(Z_mk_obs)
  }

  # marker-by-id basis
  # - enables vcov association features
  if (!has_idm) {
    Z_idm_obs <- matrix(0.0, nrow(dl), 0)
    Q_idm <- 0L
  } else {
    Z_idm_obs <- do.call(cbind, lapply(idm_rhs_list, function(rhs) .mm(rhs, dl)))
    Q_idm <- ncol(Z_idm_obs)
  }

  # Time presence checks + CS warning
  # - disable slope associations if no time terms exist
  id_has_time <- any(vapply(bars[id_idx], function(bt) .expr_has_time(bt[[2]], time_var), logical(1)))
  mk_has_time <- any(vapply(mk_rhs_list, function(rhs) .expr_has_time(rhs, time_var), logical(1)))
  idm_has_time <- any(vapply(idm_rhs_list, function(rhs) .expr_has_time(rhs, time_var), logical(1)))

  if (!id_has_time && !mk_has_time && !idm_has_time) {
    stop("Neither id block nor marker block contains '", time_var, "'. At least one must include time.")
  }
  if (!id_has_time) {
    assoc <- .disable_assoc(
      assoc, "cs_mean",
      sprintf("no '%s' in ( ... | %s ), so cs_mean is time-invariant", time_var, id_var)
    )
  }
  if (!mk_has_time && !idm_has_time) {
    assoc <- .disable_assoc(
      assoc, "cs_marker",
      sprintf("no '%s' in marker block, so cs_marker is time-invariant", time_var)
    )
  }
  if (!id_has_time && !mk_has_time && !idm_has_time) {
    assoc <- .disable_assoc(
      assoc, "cs_total",
      sprintf("no '%s' in id or marker terms, so cs_total is time-invariant", time_var)
    )
  }

  # Prepare outcomes (mixed families)
  # - y_real/y_int mirror Stan data requirements
  y_real <- as.numeric(dl[[y_var]])
  y_int <- as.integer(dl[[y_var]])
  trials <- rep.int(1L, nrow(dl))
  if ("trials" %in% colnames(dl)) {
    trials <- as.integer(dl$trials)
  }

  # Hazard covariates W
  # - baseline covariates only, handled in survival submodel
  fe_rhs <- stats::update(formulaEvent, . ~ .)
  fe_rhs[[2]] <- NULL
  W <- .mm(fe_rhs, dataEvent)
  p_w <- ncol(W)

  # Covariance covariates Xcov (intercept-only -> Xcov=0)
  # - used to build subject-specific L_i in Stan
  if (length(reformulas::findbars(formulaVcov)) > 0) {
    cli::cli_abort(c(
      x = "{.arg formulaVcov} does not support random-effects terms.",
      i = "Remove all ( ... | ... ) terms from {.arg formulaVcov}."
    ))
  }
  fv_rhs <- stats::update(formulaVcov, . ~ .)
  fv_rhs[[2]] <- NULL
  if (length(fv_rhs) >= 3 && .expr_has_time(fv_rhs[[3]], time_var)) {
    cli::cli_abort(c(
      x = "{.arg formulaVcov} cannot include the time variable {.arg {time_var}}.",
      i = "Remove time from {.arg formulaVcov} or move it to longitudinal formulas."
    ))
  }
  Xtmp <- .mm(fv_rhs, dataEvent)
  if (ncol(Xtmp) == 1 && colnames(Xtmp)[1] == "(Intercept)") {
    K_cov <- 1L
    Xcov <- matrix(0.0, nrow(dataEvent), 1)
  } else {
    if ("(Intercept)" %in% colnames(Xtmp)) Xtmp <- Xtmp[, colnames(Xtmp) != "(Intercept)", drop = FALSE]
    if (ncol(Xtmp) < 1) {
      K_cov <- 1L
      Xcov <- matrix(0.0, nrow(dataEvent), 1)
    } else {
      K_cov <- ncol(Xtmp)
      Xcov <- Xtmp
    }
  }

  # Survival outcomes (scaled)
  # - S_event is on [0,1] after scaling by tmax
  n_id <- nrow(dataEvent)
  S_event <- as.numeric(dataEvent$S_scaled)
  d_event <- as.integer(dataEvent[[event_var]])

  # Competing risks / multi-stage event type
  # - event_type collapses stages into a single factor
  if (all(c(stage_from_var, stage_to_var) %in% colnames(dataEvent))) {
    event_type_raw <- interaction(
      dataEvent[[stage_from_var]],
      dataEvent[[stage_to_var]],
      drop = TRUE
    )
  } else if (event_type_var %in% colnames(dataEvent)) {
    event_type_raw <- dataEvent[[event_type_var]]
  } else {
    event_type_raw <- rep.int(1L, n_id)
  }
  event_type_fac <- factor(event_type_raw)
  event_type <- as.integer(event_type_fac)
  K_event <- nlevels(event_type_fac)

  # GK times now/fwd on scaled domain
  # - u_now/u_fwd feed CV/CS feature computations
  gk <- .gk15_nodes()
  u_now <- matrix(NA_real_, n_id, 15)
  u_fwd <- matrix(NA_real_, n_id, 15)
  for (i in seq_len(n_id)) {
    u <- 0.5 * S_event[i] * (gk + 1)
    uf <- pmin(u + eps_fd, S_event[i])
    u_now[i, ] <- u
    u_fwd[i, ] <- uf
  }

  # Baseline hazard basis
  if (basehaz == "formula") {
    de_time <- dataEvent
    de_time[[time_var]] <- S_event
    Bs_event_raw <- .mm(basehaz_formula, de_time)
    Kbs <- ncol(Bs_event_raw)
    Bs_gk_raw <- .eval_rhs_list_on_times(list(basehaz_formula), dataEvent, time_var, u_now)
  } else {
    probs <- (seq_len(n_knots)) / (n_knots + 1)
    knots <- as.numeric(stats::quantile(S_event, probs = probs, names = FALSE, type = 7))
    knots <- pmin(pmax(knots, 1e-6), 1 - 1e-6)
    knots <- unique(knots)

    Bs_obj <- .make_basehaz_basis(
      x = S_event, basis = basehaz, knots = knots, degree = basehaz_degree, boundary = c(0, 1)
    )
    Bs_event_raw <- as.matrix(Bs_obj)
    Kbs <- ncol(Bs_event_raw)

    u_vec <- as.vector(t(u_now))
    B_now <- as.matrix(predict(Bs_obj, newx = u_vec))
    Bs_gk_raw <- array(B_now, dim = c(15, n_id, Kbs))
    Bs_gk_raw <- aperm(Bs_gk_raw, c(2, 1, 3))
  }

  centered <- .center_baseline(Bs_event_raw, Bs_gk_raw)
  Bs_event_c <- centered$Bs_event_c
  Bs_gk_c <- centered$Bs_gk_c

  # Association designs at GK and event times
  S_now <- S_event
  S_fwd <- pmin(S_event + eps_fd, 1)

  de_now <- dataEvent
  de_now[[time_var]] <- S_now
  de_fwd <- dataEvent
  de_fwd[[time_var]] <- S_fwd

  # Mean designs
  X_gk_now <- .eval_rhs_list_on_times(list(fixed_rhs), dataEvent, time_var, u_now)
  X_gk_fwd <- .eval_rhs_list_on_times(list(fixed_rhs), dataEvent, time_var, u_fwd)
  X_event_now <- .mm(fixed_rhs, de_now)
  X_event_fwd <- .mm(fixed_rhs, de_fwd)

  Z_id_gk_now <- .eval_rhs_list_on_times(id_rhs_list, dataEvent, time_var, u_now)
  Z_id_gk_fwd <- .eval_rhs_list_on_times(id_rhs_list, dataEvent, time_var, u_fwd)
  Z_id_event_now <- .eval_rhs_list_at_event(id_rhs_list, dataEvent, time_var, S_now)
  Z_id_event_fwd <- .eval_rhs_list_at_event(id_rhs_list, dataEvent, time_var, S_fwd)

  # Marker-only designs (optional)
  if (R_mk == 0) {
    mk0 <- .zero_marker_block(N = nrow(dl), n_id = n_id)
    Z_mk_gk_now <- mk0$Z_mk_gk_now
    Z_mk_gk_fwd <- mk0$Z_mk_gk_fwd
    Z_mk_event_now <- mk0$Z_mk_event_now
    Z_mk_event_fwd <- mk0$Z_mk_event_fwd
  } else {
    Z_mk_gk_now <- .eval_rhs_list_on_times(mk_rhs_list, dataEvent, time_var, u_now)
    Z_mk_gk_fwd <- .eval_rhs_list_on_times(mk_rhs_list, dataEvent, time_var, u_fwd)
    Z_mk_event_now <- .eval_rhs_list_at_event(mk_rhs_list, dataEvent, time_var, S_now)
    Z_mk_event_fwd <- .eval_rhs_list_at_event(mk_rhs_list, dataEvent, time_var, S_fwd)
  }

  # Marker-by-id designs
  Z_idm_gk_now <- .eval_rhs_list_on_times(idm_rhs_list, dataEvent, time_var, u_now)
  Z_idm_gk_fwd <- .eval_rhs_list_on_times(idm_rhs_list, dataEvent, time_var, u_fwd)
  Z_idm_event_now <- .eval_rhs_list_at_event(idm_rhs_list, dataEvent, time_var, S_now)
  Z_idm_event_fwd <- .eval_rhs_list_at_event(idm_rhs_list, dataEvent, time_var, S_fwd)

  # Association flags
  af <- .parse_assoc(assoc)

  # --------------------------
  # Time-index metadata
  # --------------------------
  x_cols <- colnames(X_obs)
  zid_cols <- colnames(Z_id_obs)
  zmk_cols <- if (R_mk > 0) colnames(Z_mk_obs) else character(0)
  zidm_cols <- colnames(Z_idm_obs)
  if (is.null(zidm_cols)) zidm_cols <- character(0)

  time_meta <- .make_time_index_metadata(
    x_cols = x_cols,
    zid_cols = zid_cols,
    zmk_cols = zmk_cols,
    zidm_cols = zidm_cols,
    time_var = time_var
  )

  # Prior scales (align to P)
  beta_scale <- as.numeric(priors$beta_scale)
  if (length(beta_scale) < P) {
    beta_scale <- c(beta_scale, rep(2.0, P - length(beta_scale)))
  }
  beta_scale <- beta_scale[seq_len(P)]

  # Family codes (vector by marker)
  family_long <- as.integer(family_codes)

  # Ordinal (cumulative logit) category count
  has_ordinal <- any(family_long == 11L)
  if (has_ordinal) {
    y_ord_vals <- y_int[family_long[dl$marker_int] == 11L]
    K_ord <- max(y_ord_vals, na.rm = TRUE)
    if (!is.finite(K_ord) || K_ord < 2) {
      stop("Cumulative logit requires at least 2 categories.")
    }
  } else {
    K_ord <- 2L
  }

  # Build arbitrary transformations (optional)
  arbitrary_tf_data <- build_standata_transforms(transforms)

  if (af$assoc_cv_total == 0 && af$assoc_cv_mean == 1 && af$assoc_cv_marker == 1 &&
      arbitrary_tf_data$tf_mode_cv_mean == 0 && arbitrary_tf_data$tf_mode_cv_marker == 0) {
    warning("cv_mean and cv_marker are both identity; consider using cv_total instead.", call. = FALSE)
  }
  if (af$assoc_cv_total == 1 && af$assoc_cv_mean == 1 && af$assoc_cv_marker == 1 &&
      arbitrary_tf_data$tf_mode_cv_tot == 0 &&
      arbitrary_tf_data$tf_mode_cv_mean == 0 && arbitrary_tf_data$tf_mode_cv_marker == 0) {
    warning("cv_total, cv_mean, and cv_marker are all identity; consider using only cv_total to avoid redundant association terms.", call. = FALSE)
  }
  if (af$assoc_cs_total == 0 && af$assoc_cs_mean == 1 && af$assoc_cs_marker == 1 &&
      arbitrary_tf_data$tf_mode_cs_mean == 0 && arbitrary_tf_data$tf_mode_cs_marker == 0) {
    warning("cs_mean and cs_marker are both identity; consider using cs_total instead.", call. = FALSE)
  }
  if (af$assoc_cs_total == 1 && af$assoc_cs_mean == 1 && af$assoc_cs_marker == 1 &&
      arbitrary_tf_data$tf_mode_cs_tot == 0 &&
      arbitrary_tf_data$tf_mode_cs_mean == 0 && arbitrary_tf_data$tf_mode_cs_marker == 0) {
    warning("cs_total, cs_mean, and cs_marker are all identity; consider using only cs_total to avoid redundant association terms.", call. = FALSE)
  }

  # Return Stan data
  standata_base <- list(
    n_id = as.integer(n_id),
    N = as.integer(nrow(dl)),
    id = as.integer(dl$id_int),
    marker = as.integer(dl$marker_int),
    D = as.integer(D),
    y_real = as.numeric(y_real),
    y_int = as.integer(y_int),
    trials = as.integer(trials),
    family_long = as.integer(family_long),
    n_family_sigma = fam_sigma$n,
    marker_to_sigma_family = fam_sigma$marker_to,
    n_family_nu = fam_nu$n,
    marker_to_nu_family = fam_nu$marker_to,
    n_family_phi = fam_phi$n,
    marker_to_phi_family = fam_phi$marker_to,
    n_family_alpha = fam_alpha$n,
    marker_to_alpha_family = fam_alpha$marker_to,
    n_family_phi_beta = fam_phi_beta$n,
    marker_to_phi_beta_family = fam_phi_beta$marker_to,
    n_family_tau_sde = fam_tau_sde$n,
    marker_to_tau_sde_family = fam_tau_sde$marker_to,
    P = as.integer(P),
    X_obs = X_obs,
    R_id = as.integer(R_id),
    Z_id_obs = Z_id_obs,
    R_mk = as.integer(R_mk),
    Z_mk_obs = Z_mk_obs,
    Q_idm = as.integer(Q_idm),
    Z_idm_obs = Z_idm_obs,
    flag_resid_dim = as.integer(flag_resid_dim),
    indep_id_re = as.integer(indep_flags$indep_id_re),
    indep_marker_re = as.integer(indep_flags$indep_marker_re),
    indep_marker_byid_latent_re = as.integer(1L), # keep as in your Page; extend if needed
    indep_idmarker_cov = as.integer(indep_flags$indep_idmarker_cov),
    allow_marker_crosscorr = as.integer(allow_marker_crosscorr),
    vcov_diag_link = as.integer(vcov_diag_link_code),
    use_tau_sde_fixed = as.integer(use_tau_sde_fixed),
    tau_sde_fixed = as.numeric(tau_sde_fixed_value),
    p_w = as.integer(p_w),
    W = W,
    K_cov = as.integer(K_cov),
    Xcov = Xcov,
    Kbs = as.integer(Kbs),
    Bs_event_c = Bs_event_c,
    Bs_gk_c = Bs_gk_c,
    tau_spline = tau_spline,
    S_event = as.numeric(S_event),
    d_event = as.integer(d_event),
    K_event = as.integer(K_event),
    event_type = as.integer(event_type),
    eps_fd = eps_fd,

    # Distributional regression
    P_sigma = as.integer(dist_sigma$P),
    X_sigma = dist_sigma$X,
    n_re_sigma = as.integer(re_sigma$n_re),
    K_sigma = as.array(as.integer(re_sigma$K)),
    G_sigma = as.array(as.integer(re_sigma$G)),
    K_sigma_max = as.integer(re_sigma$K_max),
    G_sigma_max = as.integer(re_sigma$G_max),
    Z_sigma = re_sigma$Z,
    J_sigma = re_sigma$J_mat,
    P_nu = as.integer(dist_nu$P),
    X_nu = dist_nu$X,
    n_re_nu = as.integer(re_nu$n_re),
    K_nu = as.array(as.integer(re_nu$K)),
    G_nu = as.array(as.integer(re_nu$G)),
    K_nu_max = as.integer(re_nu$K_max),
    G_nu_max = as.integer(re_nu$G_max),
    Z_nu = re_nu$Z,
    J_nu = re_nu$J_mat,
    P_phi = as.integer(dist_phi$P),
    X_phi = dist_phi$X,
    n_re_phi = as.integer(re_phi$n_re),
    K_phi = as.array(as.integer(re_phi$K)),
    G_phi = as.array(as.integer(re_phi$G)),
    K_phi_max = as.integer(re_phi$K_max),
    G_phi_max = as.integer(re_phi$G_max),
    Z_phi = re_phi$Z,
    J_phi = re_phi$J_mat,
    P_alpha = as.integer(dist_alpha$P),
    X_alpha = dist_alpha$X,
    n_re_alpha = as.integer(re_alpha$n_re),
    K_alpha = as.array(as.integer(re_alpha$K)),
    G_alpha = as.array(as.integer(re_alpha$G)),
    K_alpha_max = as.integer(re_alpha$K_max),
    G_alpha_max = as.integer(re_alpha$G_max),
    Z_alpha = re_alpha$Z,
    J_alpha = re_alpha$J_mat,
    P_phi_beta = as.integer(dist_phi_beta$P),
    X_phi_beta = dist_phi_beta$X,
    n_re_phi_beta = as.integer(re_phi_beta$n_re),
    K_phi_beta = as.array(as.integer(re_phi_beta$K)),
    G_phi_beta = as.array(as.integer(re_phi_beta$G)),
    K_phi_beta_max = as.integer(re_phi_beta$K_max),
    G_phi_beta_max = as.integer(re_phi_beta$G_max),
    Z_phi_beta = re_phi_beta$Z,
    J_phi_beta = re_phi_beta$J_mat,
    P_tau_sde = as.integer(dist_tau_sde$P),
    X_tau_sde = dist_tau_sde$X,
    n_re_tau_sde = as.integer(re_tau_sde$n_re),
    K_tau_sde = as.array(as.integer(re_tau_sde$K)),
    G_tau_sde = as.array(as.integer(re_tau_sde$G)),
    K_tau_sde_max = as.integer(re_tau_sde$K_max),
    G_tau_sde_max = as.integer(re_tau_sde$G_max),
    Z_tau_sde = re_tau_sde$Z,
    J_tau_sde = re_tau_sde$J_mat,
    K_ord = as.integer(K_ord),

    # Mean association designs
    X_gk_now = X_gk_now,
    X_gk_fwd = X_gk_fwd,
    Z_id_gk_now = Z_id_gk_now,
    Z_id_gk_fwd = Z_id_gk_fwd,
    X_event_now = X_event_now,
    X_event_fwd = X_event_fwd,
    Z_id_event_now = Z_id_event_now,
    Z_id_event_fwd = Z_id_event_fwd,

    # Marker association designs
    Z_mk_gk_now = Z_mk_gk_now,
    Z_mk_gk_fwd = Z_mk_gk_fwd,
    Z_idm_gk_now = Z_idm_gk_now,
    Z_idm_gk_fwd = Z_idm_gk_fwd,
    Z_mk_event_now = Z_mk_event_now,
    Z_mk_event_fwd = Z_mk_event_fwd,
    Z_idm_event_now = Z_idm_event_now,
    Z_idm_event_fwd = Z_idm_event_fwd,
    assoc_cv_total = af$assoc_cv_total,
    assoc_cv_mean = af$assoc_cv_mean,
    assoc_cv_marker = af$assoc_cv_marker,
    assoc_cs_total = af$assoc_cs_total,
    assoc_cs_mean = af$assoc_cs_mean,
    assoc_cs_marker = af$assoc_cs_marker,
    assoc_vcov = af$assoc_vcov,
    shrinkage = as.integer(shrinkage),

    # --------------------------
    # TRANSFORMATIONS
    # --------------------------
    # --------------------------
    # PRIOR SCALES 
    # --------------------------
    beta_scale = as.numeric(beta_scale),
    alpha_scale = as.numeric(priors$alpha_scale),
    lkj_eta = as.numeric(priors$lkj_eta),

    # --------------------------
    # Time metadata
    # --------------------------
    tmax = as.numeric(tmax),
    n_time_beta = as.integer(time_meta$n_time_beta),
    idx_time_beta = as.array(as.integer(time_meta$idx_time_beta)),
    n_time_uid = as.integer(time_meta$n_time_uid),
    idx_time_uid = as.array(as.integer(time_meta$idx_time_uid)),
    n_time_vmk = as.integer(time_meta$n_time_vmk),
    idx_time_vmk = as.array(as.integer(time_meta$idx_time_vmk)),
    n_time_widm = as.integer(time_meta$n_time_widm),
    idx_time_widm = as.array(as.integer(time_meta$idx_time_widm)),

    # --------------------------
    # Other metadata
    # --------------------------
    x_cols = x_cols,
    w_cols = colnames(W),
    zid_cols = zid_cols,
    zmk_cols = zmk_cols,
    zidm_cols = zidm_cols,
    time_var = time_var,
    marker_levels = marker_levels,
    marker_weights = as.numeric(marker_weights),
    estimate_marker_weights = as.integer(estimate_marker_weights),
    marker_weight_scale = as.numeric(marker_weight_scale),
    family_codes = family_codes,
    family_names = family_names,
    family_sigma_codes = fam_sigma$family_codes,
    family_sigma_names = fam_sigma$family_names,
    family_nu_codes = fam_nu$family_codes,
    family_nu_names = fam_nu$family_names,
    family_phi_codes = fam_phi$family_codes,
    family_phi_names = fam_phi$family_names,
    family_alpha_codes = fam_alpha$family_codes,
    family_alpha_names = fam_alpha$family_names,
    family_phi_beta_codes = fam_phi_beta$family_codes,
    family_phi_beta_names = fam_phi_beta$family_names,
    family_tau_sde_codes = fam_tau_sde$family_codes,
    family_tau_sde_names = fam_tau_sde$family_names,
    tf_compositions = NULL,
    basehaz = basehaz,
    n_knots = n_knots,
    basehaz_degree = basehaz_degree,
    Bs_obj = if (basehaz != "formula") Bs_obj else NULL,
    dist_cols = list(
      sigma = dist_sigma$cols,
      nu = dist_nu$cols,
      phi = dist_phi$cols,
      alpha = dist_alpha$cols,
      phi_beta = dist_phi_beta$cols,
      tau_sde = dist_tau_sde$cols
    ),
    dist_re_terms = list(
      sigma = re_sigma$terms,
      nu = re_nu$terms,
      phi = re_phi$terms,
      alpha = re_alpha$terms,
      phi_beta = re_phi_beta$terms,
      tau_sde = re_tau_sde$terms
    ),
    dist_formulas = dist_formulas
  )

  c(standata_base, arbitrary_tf_data)
}
