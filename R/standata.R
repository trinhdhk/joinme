#' Prepare standata for Joint Nested Mixed-effects (JoiNMe) model.
#'
#' @importFrom stats predict setNames
#'
#' @description
#' This function process the input from `fit` to prepare standata for Joint Nested Mixed-effects (JoiNMe) model.
#'
#' @details
#' This builder performs three key steps:
#' 1. Constructs fixed-effect and random-effect design matrices on a scaled time axis.
#' 2. Creates distributional regression matrices (including optional random effects).
#' 3. Assembles spline bases and Gauss-Kronrod nodes for survival integration.
#'    Only the node count is passed to Stan; nodes/weights are hardcoded in the
#'    Stan functions that evaluate the cumulative hazard.
#'
#' Time is scaled internally as `t_scaled = t/t_max` for numerical stability; indices
#' are recorded to rescale time-associated coefficients back to the original units
#' inside Stan.
#'
#' @param formulaLong Longitudinal formula defining fixed effects, id-level effects,
#'   and the marker block. The marker block may optionally include an inner
#'   `( ... | id )` term for marker-by-id random effects. When the inner term is
#'   omitted, marker-by-id random effects are disabled (Q_idm = 0). Use `||` to
#'   enforce independence. In nested marker terms, outer `( ... || marker )`
#'   keeps marker-only and marker-by-id blocks independent, while inner
#'   `( ... || id )` makes the marker-by-id covariance diagonal. Grouping
#'   terms may use `weighted(group, weights = <column>)`
#'   to define formula-scoped positive subject/group weights.
#' @param dataLong Long-format longitudinal data with columns for id, marker, time,
#'   outcome, and covariates referenced in `formulaLong`.
#' @param formulaEvent Survival formula for baseline covariates and event model.
#'   Supported LHS forms are:
#'   - `survival::Surv(time, status)`
#'   - `survival::Surv(start, stop, status)`
#'   - `survival::Surv(time, status, type = "left")`
#'   - `survival::Surv(time1, time2, type = "interval2")`
#'   The legacy `type = "interval"` form is intentionally rejected.
#' @param dataEvent Event-process data with either one row per id
#'   (`Surv(time, status)`) or multiple time-split rows per id
#'   (`Surv(start, stop, status)`).
#'   Covariates referenced in `formulaEvent` may vary by interval in the split form.
#' @param formulaVCov Covariance regression formula for id-specific marker-by-id effects.
#'   If the marker block omits the inner `( ... | id )`, then marker-by-id effects
#'   are absent and covariance-style associations (`corr`, `vcov`) are not allowed.
#'   Downstream, `corr` uses the off-diagonal entries of the subject-specific
#'   Cholesky-correlation factor `K`, whereas `vcov` uses those off-diagonal `K`
#'   entries together with the subject-specific standard deviations. If both
#'   `corr` and `vcov` are requested, `vcov` is kept and `corr` is ignored with a
#'   warning. The default `~ 1` is valid and yields an intercept-only covariance
#'   regression with no subject-level slope columns.
#' @param formulaDist Optional list of formulas for distributional regression.
#'   Supported parameters are \code{sigma}, \code{nu}, \code{phi},
#'   \code{alpha}, \code{kappa}, and \code{tau}. Here \code{nu} is reserved
#'   for Student-t degrees of freedom; \code{kappa} is the positive Beta
#'   sample-size parameter defining shapes \eqn{\mu\kappa} and
#'   \eqn{(1-\mu)\kappa}; and \code{tau} is the skew-double-exponential
#'   quantile/asymmetry parameter in \eqn{(0,1)}.
#'
#'   Three forms are supported:
#'   1. Named list with RHS-only formulas, e.g. `list(sigma = ~ 1 + time)`.
#'   2. Unnamed list with LHS parameter names, e.g. `list(sigma ~ 1 + time)`.
#'   Random-effects terms with `|` are supported; nested random-effects formulas are not.
#'   Grouping factors in distributional random-effects terms may also use
#'   `weighted(group, weights = <column>)`.
#' @param id_var Column name for subject id in both longitudinal and event data.
#' @param marker_var Column name for marker/biomarker id in longitudinal data.
#' @param time_var Column name for longitudinal time in `dataLong`.
#' @param eps_fd Positive finite-difference step for association derivatives.
#' @param assoc Character vector specifying association components
#'   (e.g., "cv_mean", "cs_total", "corr", "vcov"). The `corr` channel
#'   targets off-diagonal entries of the subject-specific Cholesky-correlation
#'   factor `K`; the `vcov` channel targets those off-diagonal `K` entries
#'   together with the subject-specific standard deviations.
#' @param families Optional marker-specific family specification. Can be
#'   character family names or `jm_family(...)` entries with per-marker links.
#'   If NULL, all markers use Gaussian responses with identity link.
#'   Supported named forward links are `identity`, `log`, `logit`, `probit`,
#'   and `exp`; `jm_family()` also accepts an invertible formula link or a
#'   directly specified inverse-link formula. In formula syntax,
#'   `inv_Phi`/`qnorm`/`probit` denote the standard normal quantile and
#'   `Phi`/`pnorm` denote the standard normal CDF. Thus a probit forward link
#'   is inverted to `Phi` before Stan data are constructed.
#' @param transforms Optional list specifying transformations for association terms.
#'   Each element (cv_total, cs_total, corr, vcov) is a list with a `type` and fields
#'   required by that type (see `build_standata_transforms()`). For covariance-style
#'   terms, transforms apply either to off-diagonal `K` features (`corr`) or to the
#'   combined off-diagonal `K` plus subject-specific SD features (`vcov`).
#'   Functional transforms may request fit-only affine-shift parameters through
#'   `intercept = TRUE` and/or `slope = TRUE` inside the transform formula.
#'   Ordered piecewise-linear fits use
#'   `list(type = "pwlin", knots = ..., direction = "increasing")` (or
#'   `"decreasing"`). Their knot ordinates are estimated from simplex
#'   increments; legacy `x` is accepted as an alias for `knots`, while
#'   legacy `y` no longer fixes the fitted curve.
#' @param beta_prior Prior specification for longitudinal fixed effects.
#' @param alpha_prior Prior specification for association parameters.
#' @param iota_prior Prior specification for fit-only affine-shift intercept and
#'   slope parameters used by functional association transforms.
#' @param lkj_prior Prior specification for correlation structures.
#' @param allow_marker_crosscorr Integer flag; 1 allows cross-marker correlation in marker RE.
#' @param shrinkage Integer flag controlling shrinkage behaviour for marker-by-id effects. 
#'   0 = student_t(6, 0, 1), 1 = double_exponential(0, 1), 2 = std_normal();
#' @param marker_weights Optional base weights for marker-specific association
#'   components. If `shared_marker_weights = TRUE`, provide one numeric vector of
#'   length D (or one named numeric vector keyed by marker level) that is shared
#'   across all active weighted association terms. If
#'   `shared_marker_weights = FALSE`, you may instead provide a named list with
#'   entries `cv_total`, `cs_total`, `cv_marker`, and `cs_marker`. Each supplied
#'   entry is aligned to marker order and used only for the matching active
#'   weighted association term.
#' @param fixed_marker_weights Logical; if TRUE, marker weights are kept fixed
#'   at `marker_weights` (no perturbation). If FALSE, marker weights are estimated
#'   via signed additive perturbations,
#'   `marker_weights + z_marker_weights`. The standardized latent family is
#'   selected by `shrinkage`: Student-t(6), Laplace, or Normal for 0, 1, or 2.
#' @param shared_marker_weights Logical; if TRUE, all active weighted
#'   marker-based association terms share one marker-weight structure. If FALSE,
#'   each active weighted marker-based association term gets its own marker-weight
#'   structure.
#' @param flag_resid_dim Integer flag to include residual dimension checks.
#' @param basehaz Baseline hazard basis type: "bs", "ns", or "formula".
#' @param basehaz_n_knots Number of internal knots for spline baseline hazards.
#' @param basehaz_knots Optional numeric vector of internal knots for spline baseline hazards.
#' @param basehaz_degree Degree of spline basis for baseline hazard.
#' @param basehaz_formula Formula for baseline hazard when `basehaz = "formula"`.
#' @param tau_spline Prior scale for spline coefficients (penalised spline).
#' @param quadrature_nodes Optional positive integer target for total quadrature
#'   points. Allowed values are exactly `7`, `15`, `31`, `41`, `51`, and `61`.
#'   Only the node count is passed to Stan; GK nodes/weights are fixed in Stan.
#' @param vcov_diag_link Link for the subject-specific standard deviation regression:
#'   "softplus" or "exp".
#' @param tau_fixed Optional fixed quantile/asymmetry parameter for the skew
#'   double exponential family. It must lie strictly between zero and one; the
#'   value \eqn{0.5} gives the symmetric double exponential distribution. It
#'   cannot be combined with a \code{tau} distributional regression.
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
  formulaVCov = ~1,
  formulaDist = NULL,
  id_var = "id",
  marker_var = "marker",
  time_var = "time",
  eps_fd = 1e-2,
  assoc = c("cv_mean"),
  families = NULL,
  transforms = NULL,
  beta_prior = NULL,
  alpha_prior = NULL,
  iota_prior = NULL,
  lkj_prior = NULL,
  allow_marker_crosscorr = 1L,
  shrinkage = 0L,
  marker_weights = NULL,
  fixed_marker_weights = FALSE,
  shared_marker_weights = TRUE,
  # flag_resid_dim = 0L,
  basehaz = c("bs", "ns", "formula"),
  basehaz_n_knots = 5L,
  basehaz_knots = NULL,
  basehaz_degree = 3L,
  basehaz_formula = ~ 1 + time,
  tau_spline = 0.4,
  quadrature_nodes = NULL,
  vcov_diag_link = c("softplus", "exp"),
  tau_fixed = NULL,
  seed = .Random.seed[[1]]
) {
  assertthat::assert_that(is.data.frame(dataLong), msg = "dataLong must be a data.frame")
  assertthat::assert_that(is.data.frame(dataEvent), msg = "dataEvent must be a data.frame")
  assoc <- .validate_assoc_channels(assoc, context = "joinme_standata()")
  y_var <- .resolve_response_var(
    formulaLong = formulaLong,
    dataLong = dataLong,
    context = "joinme_standata()"
  )
  event_vars <- .resolve_event_model_vars(
    formulaEvent = formulaEvent,
    dataEvent = dataEvent,
    context = "joinme_standata()"
  )
  event_start <- as.numeric(event_vars$event_start)
  event_stop <- as.numeric(event_vars$event_stop)
  event_status <- event_vars$event_status
  event_censor_type <- as.integer(event_vars$event_censor_type %||% rep(0L, length(event_stop)))
  surv_type <- as.character(event_vars$surv_type %||% "right")
  basehaz <- match.arg(basehaz)
  vcov_diag_link <- match.arg(vcov_diag_link)
  quad_req <- .resolve_gk_request(nodes = quadrature_nodes %||% 15L)

  if (!is.logical(fixed_marker_weights) || length(fixed_marker_weights) != 1L || is.na(fixed_marker_weights)) {
    cli::cli_abort(c(
      x = "{.arg fixed_marker_weights} must be TRUE/FALSE.",
      i = "Use TRUE to keep supplied marker weights fixed; FALSE to estimate marker-weight perturbations."
    ))
  }
  if (!is.logical(shared_marker_weights) || length(shared_marker_weights) != 1L || is.na(shared_marker_weights)) {
    cli::cli_abort(c(
      x = "{.arg shared_marker_weights} must be TRUE/FALSE.",
      i = "Use TRUE to share one marker-weight structure across weighted association terms, or FALSE to use one structure per active weighted association term."
    ))
  }

  # Workflow: parse formulas -> build matrices -> assemble survival basis -> pack list
  # Handle mixed families
  priors <- .build_priors(beta_prior = beta_prior, alpha_prior = alpha_prior, iota_prior = iota_prior, lkj_prior = lkj_prior)
  
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

  grp <- vapply(bars, function(b) .group_name_from_expr(b[[3]]), character(1))
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
  if (!has_idm && any(c("corr", "vcov") %in% assoc)) {
    cli::cli_abort(c(
      x = "Associations {.arg corr} and {.arg vcov} require marker-by-id random effects.",
      i = "Add an inner ( ... | {.arg id_var} ) term inside the marker block, or remove covariance-style association terms from {.arg assoc}."
    ))
  }

  # Fixed RHS only
  fixed_rhs <- stats::update(f_fix, . ~ .)
  fixed_rhs[[2]] <- NULL

  # Map ids by event process data order
  ids <- sort(unique(dataEvent[[id_var]]))
  id_map <- setNames(seq_along(ids), ids)
  dataEvent$id_int <- as.integer(id_map[as.character(dataEvent[[id_var]])])

  # Attach resolved event interval columns and order rows by (id, stop, start)
  dataEvent$event_start <- event_start
  dataEvent$event_stop <- event_stop
  dataEvent$event_status <- event_status
  dataEvent$event_censor_type <- event_censor_type
  ord_event <- order(dataEvent$id_int, dataEvent$event_stop, dataEvent$event_start)
  dataEvent <- dataEvent[ord_event, , drop = FALSE]
  rownames(dataEvent) <- NULL

  # Subject-level snapshot (first interval row per id) for id-level constructs
  first_idx_by_id <- match(seq_along(ids), dataEvent$id_int)
  dataEvent_id <- dataEvent[first_idx_by_id, , drop = FALSE]

  # Interval-level vectors after sorting
  event_start <- as.numeric(dataEvent$event_start)
  event_stop <- as.numeric(dataEvent$event_stop)
  event_status <- dataEvent$event_status
  event_censor_type <- as.integer(dataEvent$event_censor_type)

  # Validate interval structure within each id
  n_id <- length(ids)
  interval_index <- .build_event_interval_index(
    id_int = dataEvent$id_int,
    n_id = n_id,
    context = "joinme_standata()"
  )
  for (i in seq_len(n_id)) {
    ii <- interval_index$event_start_idx[i]:interval_index$event_end_idx[i]
    starts_i <- event_start[ii]
    stops_i <- event_stop[ii]
    if (is.unsorted(stops_i, strictly = FALSE)) {
      cli::cli_abort(c(
        x = "joinme_standata(): event stop times must be non-decreasing within each id.",
        i = "Check interval ordering for subject index {i}."
      ))
    }
    if (any(starts_i[-1] < stops_i[-length(stops_i)])) {
      cli::cli_abort(c(
        x = "joinme_standata(): overlapping event intervals are not supported.",
        i = "Each interval must start at or after the previous interval stop for subject index {i}."
      ))
    }
  }

  if (surv_type %in% c("right", "counting")) {
    event_outcomes_all <- .derive_event_outcomes(
      status_raw = event_status,
      context = "joinme_standata()"
    )
    d_event <- event_outcomes_all$d_event
    event_type <- event_outcomes_all$event_type
    K_event <- event_outcomes_all$K_event
  } else {
    d_event <- as.integer(event_censor_type == 1L)
    event_type <- rep.int(1L, length(event_censor_type))
    K_event <- 1L
  }

  for (i in seq_len(n_id)) {
    ii <- interval_index$event_start_idx[i]:interval_index$event_end_idx[i]
    event_rows <- which(d_event[ii] == 1L)
    if (length(event_rows) > 1L) {
      cli::cli_abort(c(
        x = "joinme_standata(): each id can have at most one event interval.",
        i = "Subject index {i} has multiple rows with event status != 0."
      ))
    }
    if (length(event_rows) == 1L && event_rows != length(ii)) {
      cli::cli_abort(c(
        x = "joinme_standata(): event row must be the final interval for each id.",
        i = "Move the event indicator to the last interval row for subject index {i}."
      ))
    }
  }

  dataLong$id_int <- as.integer(id_map[as.character(dataLong[[id_var]])])

  keep_long <- !is.na(dataLong$id_int) &
    is.finite(dataLong[[time_var]]) &
    is.finite(dataLong[[y_var]]) &
    !is.na(dataLong[[marker_var]])
  dataLong <- dataLong[keep_long, , drop = FALSE]
  if (nrow(dataLong) > 0L) {
    ord_long <- order(dataLong$id_int, dataLong[[marker_var]], dataLong[[time_var]])
    dataLong <- dataLong[ord_long, , drop = FALSE]
  }

  dataLong[[marker_var]] <- factor(dataLong[[marker_var]])
  marker_levels <- levels(dataLong[[marker_var]])
  D <- length(marker_levels)
  dataLong$marker_int <- as.integer(dataLong[[marker_var]])

  .resolve_group_weights <- function(group_exprs, data, group_index, n_groups, context) {
    if (length(group_exprs) == 0) {
      return(rep(1.0, n_groups))
    }

    parsed <- lapply(group_exprs, function(expr) {
      .parse_weighted_group_expr(expr, data = data, context = context)
    })
    weighted <- parsed[vapply(parsed, function(x) isTRUE(x$is_weighted), logical(1))]
    if (length(weighted) == 0) {
      return(rep(1.0, n_groups))
    }

    ref <- weighted[[1]]$weights
    if (length(weighted) > 1) {
      for (k in 2:length(weighted)) {
        if (max(abs(weighted[[k]]$weights - ref)) > 1e-8) {
          cli::cli_abort(c(
            x = "Inconsistent weighted-group declarations in {context}.",
            i = "Use the same {.code weights = ...} specification across random-effects terms sharing a grouping factor."
          ))
        }
      }
    }

    out <- as.numeric(tapply(ref, group_index, mean))
    if (length(out) != n_groups || any(!is.finite(out)) || any(out <= 0)) {
      cli::cli_abort(c(
        x = "Invalid aggregated group weights in {context}.",
        i = "Weights must map to one finite, strictly positive value per grouping level."
      ))
    }
    out
  }

  n_id <- nrow(dataEvent_id)
  id_group_exprs <- lapply(bars[id_idx], function(bt) bt[[3]])
  marker_idx <- which(grp == marker_var)
  marker_group_exprs <- lapply(bars[marker_idx], function(bt) bt[[3]])
  idm_group_exprs <- nested$idm_group_exprs

  re_weight_id <- .resolve_group_weights(
    group_exprs = id_group_exprs,
    data = dataLong,
    group_index = dataLong$id_int,
    n_groups = n_id,
    context = "formulaLong id-level random effects"
  )
  re_weight_marker <- .resolve_group_weights(
    group_exprs = marker_group_exprs,
    data = dataLong,
    group_index = dataLong$marker_int,
    n_groups = D,
    context = "formulaLong marker-level random effects"
  )
  re_weight_idm <- .resolve_group_weights(
    group_exprs = idm_group_exprs,
    data = dataLong,
    group_index = dataLong$id_int,
    n_groups = n_id,
    context = "formulaLong marker-by-id random effects"
  )

  has_id_weight <- any(abs(re_weight_id - 1.0) > 1e-8)
  has_idm_weight <- any(abs(re_weight_idm - 1.0) > 1e-8)

  if (has_id_weight && has_idm_weight && max(abs(re_weight_id - re_weight_idm)) > 1e-8) {
    cli::cli_abort(c(
      x = "Conflicting id weights between top-level id and nested marker-by-id terms.",
      i = "When both are weighted, use the same subject weights for {.code (.. | id)} and nested {.code (.. | id)} inside marker blocks."
    ))
  }

  subject_weights <- if (has_id_weight) {
    re_weight_id
  } else if (has_idm_weight) {
    re_weight_idm
  } else {
    rep(1.0, n_id)
  }

  # Covariance-regression latent effects are id-specific and should follow
  # nested marker-by-id weighting whenever it is supplied.
  re_weight_L <- if (has_idm_weight) re_weight_idm else subject_weights

  # Marker-weighted association channels use marker weights only for:
  # cv_total, cv_marker, cs_total, cs_marker.
  active_weight_terms <- .active_weighted_assoc_terms(assoc)
  marker_weight_assoc_active <- as.integer(length(active_weight_terms) > 0L)
  estimate_marker_weights_active <- as.integer(!isTRUE(fixed_marker_weights) && marker_weight_assoc_active == 1L)
  marker_weight_spec <- .resolve_marker_weight_structure(
    marker_weights = marker_weights,
    marker_levels = marker_levels,
    active_terms = active_weight_terms,
    shared_marker_weights = shared_marker_weights,
    estimate_marker_weights = as.logical(estimate_marker_weights_active),
    context = "joinme_standata()"
  )

  # Marker-weight estimation controls
  # - estimate_marker_weights_active toggles whether shrinkage is applied in Stan
    # Independence flags from lme4-style double-bar syntax
    # - (... || id) enforces diagonal id RE covariance
    # - (... || marker) enforces diagonal marker RE covariance
    # - Nested (... || id) inside marker block enforces diagonal marker-by-id covariance
    indep_flags <- .resolve_re_independence(formulaLong, marker_var = marker_var, id_var = id_var)
  # Process mixed families (optional)
  # - family_codes drive distributional parameter availability
  # - link_codes drive family-specific inverse-link mapping in Stan
  if (!is.null(families)) {
    # Treat a single family spec object/list as scalar and recycle to D markers.
    if (inherits(families, "JoiNMe_family_spec") || (is.list(families) && !is.null(families$family))) {
      families <- rep(list(families), D)
    }
    # If only length 1, duplicate to all markers
    if (length(families) == 1) {
      families <- rep(families, D)
    }
    family_spec <- .validate_family_list(families, D, dataLong, marker_var, y_var)
    family_codes <- as.integer(family_spec$family_codes)
    link_codes <- as.integer(family_spec$link_codes)
    link_names <- as.character(family_spec$link_names)
    inv_link_n_ops <- as.integer(family_spec$inv_link_n_ops)
    inv_link_ops <- family_spec$inv_link_ops
    inv_link_n_const <- as.integer(family_spec$inv_link_n_const)
    inv_link_const <- family_spec$inv_link_const
    max_inv_link_ops <- ncol(inv_link_ops)
    max_inv_link_const <- ncol(inv_link_const)
  } else {
    family_codes <- rep(1L, D)  # Default: all gaussian
    link_names <- vapply(family_codes, .default_link_for_family, character(1))
    link_codes <- as.integer(vapply(link_names, .link_code_from_name, integer(1)))
    inv_link_specs <- lapply(link_names, .inv_link_bc_from_name)
    inv_link_n_ops <- as.integer(vapply(inv_link_specs, function(x) length(x$bytecode), integer(1)))
    inv_link_n_const <- as.integer(vapply(inv_link_specs, function(x) length(x$const_data), integer(1)))
    max_inv_link_ops <- max(inv_link_n_ops, 1L)
    max_inv_link_const <- max(inv_link_n_const, 1L)
    inv_link_ops <- matrix(0L, nrow = D, ncol = max_inv_link_ops)
    inv_link_const <- matrix(0.0, nrow = D, ncol = max_inv_link_const)
    for (d in seq_len(D)) {
      if (inv_link_n_ops[d] > 0L) {
        inv_link_ops[d, seq_len(inv_link_n_ops[d])] <- as.integer(inv_link_specs[[d]]$bytecode)
      }
      if (inv_link_n_const[d] > 0L) {
        inv_link_const[d, seq_len(inv_link_n_const[d])] <- as.numeric(inv_link_specs[[d]]$const_data)
      }
    }
  }

  family_names <- vapply(family_codes, .family_code_to_name, character(1))

  # Family-level distributional parameter indexing.
  #
  # Goal:
  # - for each distributional parameter (sigma, nu, phi, alpha, kappa,
  #   tau), build one shared parameter per unique family that requires it,
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
  fam_kappa <- .build_family_param_index("kappa")
  fam_tau <- .build_family_param_index("tau")

  vcov_diag_link_code <- if (vcov_diag_link == "exp") 1L else 0L

  use_tau_fixed <- 0L
  tau_fixed_value <- 0.5
  if (!is.null(tau_fixed)) {
    if (!is.numeric(tau_fixed) || length(tau_fixed) != 1 || !is.finite(tau_fixed)) {
      cli::cli_abort(c(
        x = "{.arg tau_fixed} must be a single finite numeric value.",
        i = "Provide a number strictly between 0 and 1."
      ))
    }
    tau_fixed_value <- as.numeric(tau_fixed)
    if (tau_fixed_value <= 0 || tau_fixed_value >= 1) {
      cli::cli_abort(c(
        x = "{.arg tau_fixed} must be between 0 and 1.",
        i = "Provide a number strictly between 0 and 1."
      ))
    }
    if (any(family_codes == 9L)) {
      use_tau_fixed <- 1L
    } else {
      cli::cli_warn(c(
        x = "{.arg tau_fixed} is ignored because no skew_double_exponential family is present.",
        i = "Remove {.arg tau_fixed} or include the skew_double_exponential family."
      ))
    }
  }
  
  dist_formulas <- .normalize_formula_dist(formulaDist)

  # A fixed tau and a tau regression describe competing quantile/asymmetry
  # models. Require exactly one specification so the posterior never contains
  # regression coefficients that are disconnected from the likelihood.
  if (use_tau_fixed == 1L && !is.null(dist_formulas$tau)) {
    cli::cli_abort(c(
      x = "Specify either {.arg tau_fixed} or a {.code tau ~ ...} distributional regression, not both.",
      i = "Remove {.arg tau_fixed} to estimate tau, or remove the tau formula to keep it fixed."
    ))
  }

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

  # Scale time to [0,1] by max event stop time
  # - tmax is also returned for coefficient rescaling in Stan
  tmax <- max(event_stop, na.rm = TRUE)
  if (!is.finite(tmax) || tmax <= 0) stop("Invalid max event time.")
  dataLong$t_scaled <- dataLong[[time_var]] / tmax
  dataEvent$S_entry_scaled <- event_start / tmax
  dataEvent$S_scaled <- event_stop / tmax
  dataEvent_id$S_entry_scaled <- dataEvent_id$event_start / tmax
  dataEvent_id$S_scaled <- dataEvent_id$event_stop / tmax

  # Build obs matrices on scaled time
  # - distributional regression matrices follow the same scaled timeline
  dl <- dataLong
  dl[[time_var]] <- dl$t_scaled
  family_by_row <- .family_by_row_from_marker(
    marker_values = dl[[marker_var]],
    marker_levels = marker_levels,
    family_codes = family_codes
  )
  attr(dl, "JoiNMe_family_by_row") <- family_by_row

  dist_sigma <- .build_dist_matrix(dist_formulas$sigma, dl, family_by_row = family_by_row)
  dist_nu <- .build_dist_matrix(dist_formulas$nu, dl, family_by_row = family_by_row)
  dist_phi <- .build_dist_matrix(dist_formulas$phi, dl, family_by_row = family_by_row)
  dist_alpha <- .build_dist_matrix(dist_formulas$alpha, dl, family_by_row = family_by_row)
  dist_kappa <- .build_dist_matrix(dist_formulas$kappa, dl, family_by_row = family_by_row)
  dist_tau <- .build_dist_matrix(dist_formulas$tau, dl, family_by_row = family_by_row)

  re_sigma <- .pad_re_terms(.build_dist_re_terms(dist_formulas$sigma, dl), nrow(dl))
  re_nu <- .pad_re_terms(.build_dist_re_terms(dist_formulas$nu, dl), nrow(dl))
  re_phi <- .pad_re_terms(.build_dist_re_terms(dist_formulas$phi, dl), nrow(dl))
  re_alpha <- .pad_re_terms(.build_dist_re_terms(dist_formulas$alpha, dl), nrow(dl))
  re_kappa <- .pad_re_terms(.build_dist_re_terms(dist_formulas$kappa, dl), nrow(dl))
  re_tau <- .pad_re_terms(.build_dist_re_terms(dist_formulas$tau, dl), nrow(dl))

  fixed_template <- .make_model_matrix_template(
    fixed_rhs,
    dl,
    boundary_var = time_var,
    boundary_values = c(0, 1)
  )
  id_rhs_list <- .bar_terms_to_rhs_list(bars[id_idx])
  id_templates <- lapply(id_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, dl, boundary_var = time_var, boundary_values = c(0, 1))
  })
  mk_templates <- lapply(mk_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, dl, boundary_var = time_var, boundary_values = c(0, 1))
  })
  idm_templates <- lapply(idm_rhs_list, function(rhs) {
    .make_model_matrix_template(rhs, dl, boundary_var = time_var, boundary_values = c(0, 1))
  })

  X_obs <- .mm(fixed_template, dl)
  P <- ncol(X_obs)

  # id block
  # - always required (subject random effects)
  Z_id_obs <- do.call(cbind, lapply(id_templates, function(rhs) .mm(rhs, dl)))
  R_id <- ncol(Z_id_obs)

  # marker-only optional
  # - empty marker block yields zero columns
  if (length(mk_rhs_list) == 0) {
    mk0 <- .zero_marker_block(N = nrow(dl), n_id = n_id)
    R_mk <- mk0$R_mk
    Z_mk_obs <- mk0$Z_mk_obs
  } else {
    Z_mk_obs <- do.call(cbind, lapply(mk_templates, function(rhs) .mm(rhs, dl)))
    R_mk <- ncol(Z_mk_obs)
  }

  # marker-by-id basis
  # - enables covariance-style association features (corr, vcov)
  if (!has_idm) {
    Z_idm_obs <- matrix(0.0, nrow(dl), 0)
    Q_idm <- 0L
  } else {
    Z_idm_obs <- do.call(cbind, lapply(idm_templates, function(rhs) .mm(rhs, dl)))
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

  .validate_assoc_cov_structure(
    assoc = assoc,
    q_idm = Q_idm,
    diagonal_only = indep_flags$indep_idmarker_cov,
    context = "joinme_standata()"
  )

  # Prepare outcomes (mixed families)
  # - y_real/y_int mirror Stan data requirements
  y_real <- as.numeric(dl[[y_var]])
  y_int <- as.integer(dl[[y_var]])
  trials <- rep.int(1L, nrow(dl))
  if ("trials" %in% colnames(dl)) {
    trials <- as.integer(dl$trials)
  }

  # Hazard covariates W
  # - interval covariates for piecewise survival contributions
  W <- .mm_event(formulaEvent, dataEvent)
  p_w <- ncol(W)

  # Covariance covariates Xcov
  # - used to build subject-specific lower-triangular factors L_i in Stan
  formulaVCov <- .resolve_vcov_formula(
    formulaVCov = formulaVCov,
    default = ~ 1,
    context = "joinme_standata()"
  )
  vcov_design <- .build_vcov_design(
    formulaVCov = formulaVCov,
    dataEvent = dataEvent_id,
    time_var = time_var,
    context = "joinme_standata()"
  )
  K_cov <- vcov_design$K_cov
  Xcov <- vcov_design$Xcov

  # Survival outcomes (scaled)
  # - S_entry/S_event are interval bounds on [0,1] after scaling by tmax
  S_entry <- as.numeric(dataEvent$S_entry_scaled)
  S_event <- as.numeric(dataEvent$S_scaled)
  n_event <- nrow(dataEvent)

  # GK times now/fwd on scaled domain
  # - u_now/u_fwd feed CV/CS feature computations
  # - only n_gk is passed to Stan; nodes/weights are hardcoded in Stan
  quad <- .gk_single_panel(rule = quad_req$rule)
  n_gk <- quad$n_gk
  gk_nodes <- quad$nodes
  u_now <- matrix(NA_real_, n_event, n_gk)
  u_fwd <- matrix(NA_real_, n_event, n_gk)
  for (e in seq_len(n_event)) {
    delta_e <- S_event[e] - S_entry[e]
    u <- S_entry[e] + delta_e * gk_nodes
    uf <- pmin(u + eps_fd, 1)
    u_now[e, ] <- u
    u_fwd[e, ] <- uf
  }

  # Baseline hazard basis
  if (basehaz == "formula") {
    de_time <- dataEvent
    de_time[[time_var]] <- S_event
    Bs_event_raw <- .mm(basehaz_formula, de_time)
    Kbs <- ncol(Bs_event_raw)
    Bs_gk_raw <- .eval_rhs_list_on_times(list(basehaz_formula), dataEvent, time_var, u_now)
  } else {
    if (is.null(basehaz_knots)) {
      probs <- (seq_len(basehaz_n_knots)) / (basehaz_n_knots + 1)
      knots <- as.numeric(stats::quantile(S_event, probs = probs, names = FALSE, type = 7))
      knots <- pmin(pmax(knots, 1e-6), 1 - 1e-6)
      knots  <- unique(knots)
      basehaz_knots <- knots * tmax
    } else {
      knots <- unique(as.numeric(basehaz_knots))/tmax
    }
    

    Bs_obj <- .make_basehaz_basis(
      x = S_event, basis = basehaz, knots = knots, degree = basehaz_degree, boundary = c(0, 1)
    )
    Bs_event_raw <- as.matrix(Bs_obj)
    Kbs <- ncol(Bs_event_raw)

    u_vec <- as.vector(t(u_now))
    B_now <- as.matrix(predict(Bs_obj, newx = u_vec))
    Bs_gk_raw <- array(B_now, dim = c(n_gk, n_event, Kbs))
    Bs_gk_raw <- aperm(Bs_gk_raw, c(2, 1, 3))
  }

  centered <- .center_baseline(Bs_event_raw, Bs_gk_raw)
  Bs_event_c <- centered$Bs_event_c
  Bs_gk_c <- centered$Bs_gk_c
  # Retain the exact offsets subtracted from non-constant basis columns.  The
  # dynamic-prediction program must apply these same offsets to its quadrature
  # basis; recomputing them from a different spline type changes both the
  # baseline-hazard dimension and its interpretation.
  basehaz_col_means <- colMeans(Bs_event_raw - Bs_event_c)
  basehaz_cols <- colnames(Bs_event_raw) %||% paste0("basis_", seq_len(Kbs))
  # Association designs at GK and event times
  S_now <- S_event
  S_fwd <- pmin(S_event + eps_fd, 1)

  de_now <- dataEvent
  de_now[[time_var]] <- S_now
  de_fwd <- dataEvent
  de_fwd[[time_var]] <- S_fwd

  # Mean designs
  X_gk_now <- .eval_rhs_list_on_times(list(fixed_template), dataEvent, time_var, u_now)
  X_gk_fwd <- .eval_rhs_list_on_times(list(fixed_template), dataEvent, time_var, u_fwd)
  X_event_now <- .mm(fixed_template, de_now)
  X_event_fwd <- .mm(fixed_template, de_fwd)

  Z_id_gk_now <- .eval_rhs_list_on_times(id_templates, dataEvent, time_var, u_now)
  Z_id_gk_fwd <- .eval_rhs_list_on_times(id_templates, dataEvent, time_var, u_fwd)
  Z_id_event_now <- .eval_rhs_list_at_event(id_templates, dataEvent, time_var, S_now)
  Z_id_event_fwd <- .eval_rhs_list_at_event(id_templates, dataEvent, time_var, S_fwd)

  # Marker-only designs (optional)
  if (R_mk == 0) {
    mk0 <- .zero_marker_block(N = nrow(dl), n_id = n_event, n_gk = n_gk)
    Z_mk_gk_now <- mk0$Z_mk_gk_now
    Z_mk_gk_fwd <- mk0$Z_mk_gk_fwd
    Z_mk_event_now <- mk0$Z_mk_event_now
    Z_mk_event_fwd <- mk0$Z_mk_event_fwd
  } else {
    Z_mk_gk_now <- .eval_rhs_list_on_times(mk_templates, dataEvent, time_var, u_now)
    Z_mk_gk_fwd <- .eval_rhs_list_on_times(mk_templates, dataEvent, time_var, u_fwd)
    Z_mk_event_now <- .eval_rhs_list_at_event(mk_templates, dataEvent, time_var, S_now)
    Z_mk_event_fwd <- .eval_rhs_list_at_event(mk_templates, dataEvent, time_var, S_fwd)
  }

  # Marker-by-id designs
  Z_idm_gk_now <- .eval_rhs_list_on_times(idm_templates, dataEvent, time_var, u_now)
  Z_idm_gk_fwd <- .eval_rhs_list_on_times(idm_templates, dataEvent, time_var, u_fwd)
  Z_idm_event_now <- .eval_rhs_list_at_event(idm_templates, dataEvent, time_var, S_now)
  Z_idm_event_fwd <- .eval_rhs_list_at_event(idm_templates, dataEvent, time_var, S_fwd)

  # Association flags
  af <- .parse_assoc(assoc)
  M_corr_tf <- .assoc_transform_component_count("corr", Q_idm, diagonal_only = FALSE)
  M_vcov_tf <- .assoc_transform_component_count(
    "vcov",
    Q_idm,
    diagonal_only = as.integer(indep_flags$indep_idmarker_cov %||% 0L) == 1L
  )

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
    time_var = time_var,
    fixed_design = fixed_template,
    id_design = id_templates,
    marker_design = mk_templates,
    idm_design = idm_templates
  )

  # Prior scales (align to P)
  beta_scale <- as.numeric(priors$beta_scale)
  if (length(beta_scale) == 0) beta_scale <- 2.5
  if (length(beta_scale) < P) {
    # beta_scale <- c(beta_scale, rep(2.0, P - length(beta_scale)))
    beta_scale <- rep(beta_scale, length.out = P)
  }
  beta_scale <- beta_scale[seq_len(P)]

  allow_marker_crosscorr <- as.integer(allow_marker_crosscorr)
  if (!allow_marker_crosscorr %in% c(0L, 1L)) {
    cli::cli_abort(c(
      x = "{.arg allow_marker_crosscorr} must be binary.",
      i = "Set it to 1 to allow marker-to-marker-by-id cross-correlation when the formula structure permits it."
    ))
  }
  # if (as.integer(indep_flags$indep_marker_id_crosscorr %||% 0L) == 1L && allow_marker_crosscorr == 1L) {
  #   cli::cli_warn(c(
  #     i = "Marker-to-marker-by-id cross-correlation is disabled by nested {.code || marker} syntax.",
  #     v = "The outer marker double-bar makes the marker-only and marker-by-id blocks independent, so {.arg allow_marker_crosscorr} is being set to 0."
  #   ))
  #   allow_marker_crosscorr <- 0L
  # }
  if (as.integer(indep_flags$indep_marker_id_crosscorr %||% 0L) == 1L) {
    allow_marker_crosscorr <- 0L
  }

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
  functional_tf_data <- build_standata_transforms(
    transforms,
    n_corr_components = M_corr_tf,
    n_vcov_components = M_vcov_tf
  )

  # if (af$assoc_cv_total == 0 && af$assoc_cv_mean == 1 && af$assoc_cv_marker == 1 &&
  #     functional_tf_data$tf_mode_cv_mean == 0 && functional_tf_data$tf_mode_cv_marker == 0) {
  #   warning("cv_mean and cv_marker are both identity; consider using cv_total instead.", call. = FALSE)
  # }
  # if (af$assoc_cv_total == 1 && af$assoc_cv_mean == 1 && af$assoc_cv_marker == 1 &&
  #     functional_tf_data$tf_mode_cv_tot == 0 &&
  #     functional_tf_data$tf_mode_cv_mean == 0 && functional_tf_data$tf_mode_cv_marker == 0) {
  #   warning("cv_total, cv_mean, and cv_marker are all identity; consider using only cv_total to avoid redundant association terms.", call. = FALSE)
  # }
  # if (af$assoc_cs_total == 0 && af$assoc_cs_mean == 1 && af$assoc_cs_marker == 1 &&
  #     functional_tf_data$tf_mode_cs_mean == 0 && functional_tf_data$tf_mode_cs_marker == 0) {
  #   warning("cs_mean and cs_marker are both identity; consider using cs_total instead.", call. = FALSE)
  # }
  # if (af$assoc_cs_total == 1 && af$assoc_cs_mean == 1 && af$assoc_cs_marker == 1 &&
  #     functional_tf_data$tf_mode_cs_tot == 0 &&
  #     functional_tf_data$tf_mode_cs_mean == 0 && functional_tf_data$tf_mode_cs_marker == 0) {
  #   warning("cs_total, cs_mean, and cs_marker are all identity; consider using only cs_total to avoid redundant association terms.", call. = FALSE)
  # }

  # Return Stan data
  standata_base <- list(
    n_id = as.integer(n_id),
    N_event = as.integer(n_event),
    N = as.integer(nrow(dl)),
    id = as.integer(dl$id_int),
    event_id = as.integer(dataEvent$id_int),
    event_start_idx = as.integer(interval_index$event_start_idx),
    event_end_idx = as.integer(interval_index$event_end_idx),
    marker = as.integer(dl$marker_int),
    D = as.integer(D),
    subject_weights = as.numeric(subject_weights),
    y_real = as.numeric(y_real),
    y_int = as.integer(y_int),
    trials = as.integer(trials),
    family_long = as.integer(family_long),
    link_long = as.integer(link_codes),
    max_inv_link_ops = as.integer(max_inv_link_ops),
    inv_link_n_ops = as.integer(inv_link_n_ops),
    inv_link_ops = matrix(as.integer(inv_link_ops), nrow = nrow(inv_link_ops), ncol = ncol(inv_link_ops)),
    max_inv_link_const = as.integer(max_inv_link_const),
    inv_link_n_const = as.integer(inv_link_n_const),
    inv_link_const = matrix(as.numeric(inv_link_const), nrow = nrow(inv_link_const), ncol = ncol(inv_link_const)),
    n_family_sigma = fam_sigma$n,
    marker_to_sigma_family = fam_sigma$marker_to,
    n_family_nu = fam_nu$n,
    marker_to_nu_family = fam_nu$marker_to,
    n_family_phi = fam_phi$n,
    marker_to_phi_family = fam_phi$marker_to,
    n_family_alpha = fam_alpha$n,
    marker_to_alpha_family = fam_alpha$marker_to,
    n_family_kappa = fam_kappa$n,
    marker_to_kappa_family = fam_kappa$marker_to,
    n_family_tau = fam_tau$n,
    marker_to_tau_family = fam_tau$marker_to,
    P = as.integer(P),
    X_obs = X_obs,
    R_id = as.integer(R_id),
    Z_id_obs = Z_id_obs,
    re_weight_id = as.numeric(re_weight_id),
    R_mk = as.integer(R_mk),
    Z_mk_obs = Z_mk_obs,
    re_weight_marker = as.numeric(re_weight_marker),
    Q_idm = as.integer(Q_idm),
    Z_idm_obs = Z_idm_obs,
    re_weight_idm = as.numeric(re_weight_idm),
    re_weight_L = as.numeric(re_weight_L),
    # flag_resid_dim = as.integer(flag_resid_dim),
    indep_id_re = as.integer(indep_flags$indep_id_re),
    indep_marker_re = as.integer(indep_flags$indep_marker_re),
    indep_marker_byid_latent_re = as.integer(1L), # keep as in your Page; extend if needed
    indep_idmarker_cov = as.integer(indep_flags$indep_idmarker_cov),
    M_corr_tf = as.integer(M_corr_tf),
    M_vcov_tf = as.integer(M_vcov_tf),
    allow_marker_crosscorr = as.integer(allow_marker_crosscorr),
    vcov_diag_link = as.integer(vcov_diag_link_code),
    use_tau_fixed = as.integer(use_tau_fixed),
    tau_fixed = as.numeric(tau_fixed_value),
    p_w = as.integer(p_w),
    W = W,
    K_cov = as.integer(K_cov),
    Xcov = Xcov,
    Kbs = as.integer(Kbs),
    Bs_event_c = Bs_event_c,
    Bs_gk_c = Bs_gk_c,
    tau_spline = tau_spline,
    S_entry = as.numeric(S_entry),
    S_event = as.numeric(S_event),
    d_event = as.integer(d_event),
    event_censor_type = as.integer(event_censor_type),
    surv_type = as.character(surv_type),
    event_censor_types_present = as.integer(sort(unique(event_censor_type))),
    K_event = as.integer(K_event),
    event_type = as.integer(event_type),
    eps_fd = eps_fd,
    n_gk = as.integer(n_gk),

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
    re_weight_sigma = re_sigma$W,
    P_nu = as.integer(dist_nu$P),
    X_nu = dist_nu$X,
    n_re_nu = as.integer(re_nu$n_re),
    K_nu = as.array(as.integer(re_nu$K)),
    G_nu = as.array(as.integer(re_nu$G)),
    K_nu_max = as.integer(re_nu$K_max),
    G_nu_max = as.integer(re_nu$G_max),
    Z_nu = re_nu$Z,
    J_nu = re_nu$J_mat,
    re_weight_nu = re_nu$W,
    P_phi = as.integer(dist_phi$P),
    X_phi = dist_phi$X,
    n_re_phi = as.integer(re_phi$n_re),
    K_phi = as.array(as.integer(re_phi$K)),
    G_phi = as.array(as.integer(re_phi$G)),
    K_phi_max = as.integer(re_phi$K_max),
    G_phi_max = as.integer(re_phi$G_max),
    Z_phi = re_phi$Z,
    J_phi = re_phi$J_mat,
    re_weight_phi = re_phi$W,
    P_alpha = as.integer(dist_alpha$P),
    X_alpha = dist_alpha$X,
    n_re_alpha = as.integer(re_alpha$n_re),
    K_alpha = as.array(as.integer(re_alpha$K)),
    G_alpha = as.array(as.integer(re_alpha$G)),
    K_alpha_max = as.integer(re_alpha$K_max),
    G_alpha_max = as.integer(re_alpha$G_max),
    Z_alpha = re_alpha$Z,
    J_alpha = re_alpha$J_mat,
    re_weight_alpha = re_alpha$W,
    P_kappa = as.integer(dist_kappa$P),
    X_kappa = dist_kappa$X,
    n_re_kappa = as.integer(re_kappa$n_re),
    K_kappa = as.array(as.integer(re_kappa$K)),
    G_kappa = as.array(as.integer(re_kappa$G)),
    K_kappa_max = as.integer(re_kappa$K_max),
    G_kappa_max = as.integer(re_kappa$G_max),
    Z_kappa = re_kappa$Z,
    J_kappa = re_kappa$J_mat,
    re_weight_kappa = re_kappa$W,
    P_tau = as.integer(dist_tau$P),
    X_tau = dist_tau$X,
    n_re_tau = as.integer(re_tau$n_re),
    K_tau = as.array(as.integer(re_tau$K)),
    G_tau = as.array(as.integer(re_tau$G)),
    K_tau_max = as.integer(re_tau$K_max),
    G_tau_max = as.integer(re_tau$G_max),
    Z_tau = re_tau$Z,
    J_tau = re_tau$J_mat,
    re_weight_tau = re_tau$W,
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
    assoc_corr = af$assoc_corr,
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
    iota_scale = as.numeric(priors$iota_scale),
    lkj_eta = as.numeric(priors$lkj_eta),

    # --------------------------
    # Time metadata
    # --------------------------
    tmax = as.numeric(tmax),
    # tmax_internal = as.numeric(tmax),
    quadrature_nodes = as.integer(n_gk),
    n_time_beta = as.integer(time_meta$n_time_beta),
    idx_time_beta = as.array(as.integer(time_meta$idx_time_beta)),
    n_time_uid = as.integer(time_meta$n_time_uid),
    idx_time_uid = as.array(as.integer(time_meta$idx_time_uid)),
    n_time_vmk = as.integer(time_meta$n_time_vmk),
    idx_time_vmk = as.array(as.integer(time_meta$idx_time_vmk)),
    n_time_idm = as.integer(time_meta$n_time_idm),
    idx_time_idm = as.array(as.integer(time_meta$idx_time_idm)),

    # --------------------------
    # Other metadata
    # --------------------------
    x_cols = x_cols,
    w_cols = colnames(W),
    zid_cols = zid_cols,
    zmk_cols = zmk_cols,
    zidm_cols = zidm_cols,
    design_templates = list(
      fixed = fixed_template,
      id = id_templates,
      marker = mk_templates,
      idm = idm_templates
    ),
    time_var = time_var,
    marker_levels = marker_levels,
    marker_weights = if (isTRUE(shared_marker_weights)) {
      as.numeric(marker_weight_spec$base_matrix[1, ])
    } else {
      NULL
    },
    marker_weights_by_term = marker_weight_spec$base_by_term,
    marker_weights_cv_total = as.numeric(marker_weight_spec$base_by_term$cv_total),
    marker_weights_cs_total = as.numeric(marker_weight_spec$base_by_term$cs_total),
    marker_weights_cv_marker = as.numeric(marker_weight_spec$base_by_term$cv_marker),
    marker_weights_cs_marker = as.numeric(marker_weight_spec$base_by_term$cs_marker),
    shared_marker_weights = as.integer(isTRUE(marker_weight_spec$shared_marker_weights)),
    n_marker_weight_sets = as.integer(marker_weight_spec$n_sets),
    marker_weight_set_cv_total = as.integer(marker_weight_spec$set_index[["cv_total"]]),
    marker_weight_set_cs_total = as.integer(marker_weight_spec$set_index[["cs_total"]]),
    marker_weight_set_cv_marker = as.integer(marker_weight_spec$set_index[["cv_marker"]]),
    marker_weight_set_cs_marker = as.integer(marker_weight_spec$set_index[["cs_marker"]]),
    estimate_marker_weights = as.integer(estimate_marker_weights_active),
    fixed_marker_weights = as.integer(!as.logical(estimate_marker_weights_active)),
    use_marker_weight_assoc = as.integer(marker_weight_assoc_active),
    family_codes = family_codes,
    family_names = family_names,
    link_names = link_names,
    family_sigma_codes = fam_sigma$family_codes,
    family_sigma_names = fam_sigma$family_names,
    family_nu_codes = fam_nu$family_codes,
    family_nu_names = fam_nu$family_names,
    family_phi_codes = fam_phi$family_codes,
    family_phi_names = fam_phi$family_names,
    family_alpha_codes = fam_alpha$family_codes,
    family_alpha_names = fam_alpha$family_names,
    family_kappa_codes = fam_kappa$family_codes,
    family_kappa_names = fam_kappa$family_names,
    family_tau_codes = fam_tau$family_codes,
    family_tau_names = fam_tau$family_names,
    tf_compositions = NULL,
    basehaz = basehaz,
    basehaz_n_knots = basehaz_n_knots,
    basehaz_knots = basehaz_knots,
    basehaz_degree = basehaz_degree,
    basehaz_formula = if (basehaz == "formula") basehaz_formula else NULL,
    basehaz_col_means = as.numeric(basehaz_col_means),
    basehaz_cols = as.character(basehaz_cols),
    Bs_obj = if (basehaz != "formula") Bs_obj else NULL,
    dist_cols = list(
      sigma = dist_sigma$cols,
      nu = dist_nu$cols,
      phi = dist_phi$cols,
      alpha = dist_alpha$cols,
      kappa = dist_kappa$cols,
      tau = dist_tau$cols
    ),
    dist_re_terms = list(
      sigma = re_sigma$terms,
      nu = re_nu$terms,
      phi = re_phi$terms,
      alpha = re_alpha$terms,
      kappa = re_kappa$terms,
      tau = re_tau$terms
    ),
    dist_formulas = dist_formulas
  )

  c(standata_base, functional_tf_data)
}
