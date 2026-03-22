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

# File overview:
# - Common helpers for matrix construction, formula parsing, and validation.
# - Shared utilities used by standata, simulation, and diagnostics.

suppressPackageStartupMessages({
  library(lme4)
  library(reformulas)
  library(splines)
  library(posterior)
})

#' Resolve longitudinal response column name from formula LHS
#'
#' @description
#' Determine the response column used by the longitudinal model by reading the
#' left-hand side of `formulaLong` and validating it against `dataLong`.
#'
#' @param formulaLong Longitudinal formula.
#' @param dataLong Longitudinal dataset.
#' @param context Character label for error messages.
#'
#' @return Character scalar naming the resolved response column.
#' @keywords internal
.resolve_response_var <- function(formulaLong, dataLong, context = "joinme") {
  if (!inherits(formulaLong, "formula") || length(formulaLong) < 3) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaLong} must be a two-sided formula with a response on the left-hand side.",
      i = "Example: {.code y ~ time + x1 + (1 + time | id)}"
    ))
  }

  lhs_vars <- all.vars(formulaLong[[2]])
  if (length(lhs_vars) != 1L || !nzchar(lhs_vars[[1]])) {
    cli::cli_abort(c(
      x = "{context}: could not uniquely determine the longitudinal response variable from formula LHS.",
      i = "Use a single response variable on the left-hand side, e.g. {.code vol ~ ...}."
    ))
  }

  response_var <- lhs_vars[[1]]
  if (!(response_var %in% names(dataLong))) {
    cli::cli_abort(c(
      x = "{context}: formulaLong response column {.val {response_var}} was not found in {.arg dataLong}.",
      i = "Available columns include: {paste(utils::head(names(dataLong), 20), collapse = ', ')}{if (ncol(dataLong) > 20) ', ...' else ''}."
    ))
  }

  response_var
}

#' Resolve event outcomes from formulaEvent via Surv model response
#'
#' @param formulaEvent Event formula with `Surv(...)` on the left-hand side.
#' @param dataEvent Event dataset.
#' @param context Character label for error messages.
#'
#' @return Named list with extracted `event_time`, `event_status`, and `surv_type`.
#' @keywords internal
.resolve_event_model_vars <- function(formulaEvent, dataEvent, context = "joinme") {
  if (!inherits(formulaEvent, "formula") || length(formulaEvent) < 3) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaEvent} must be a two-sided survival formula.",
      i = "Example: {.code survival::Surv(event_time, event_status) ~ x1 + x2}."
    ))
  }

  mf <- tryCatch(
    stats::model.frame(formulaEvent, data = dataEvent, na.action = stats::na.pass),
    error = function(e) {
      cli::cli_abort(c(
        x = "{context}: failed to evaluate {.arg formulaEvent} on {.arg dataEvent}.",
        i = "Error: {e$message}"
      ))
    }
  )
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv")) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaEvent} response must evaluate to a {.code Surv} object.",
      i = "Use {.code survival::Surv(time, status)} or {.code survival::Surv(start, stop, status)} on the left-hand side."
    ))
  }

  surv_mat <- unclass(y)
  if (!is.matrix(surv_mat) || ncol(surv_mat) < 2L) {
    cli::cli_abort(c(
      x = "{context}: unsupported {.code Surv} response format.",
      i = "Expected at least time and status columns from survival::Surv()."
    ))
  }

  if (ncol(surv_mat) == 2L) {
    event_time <- surv_mat[, 1]
    event_status <- surv_mat[, 2]
  } else {
    event_time <- surv_mat[, 2]
    event_status <- surv_mat[, 3]
  }

  list(
    event_time = as.numeric(event_time),
    event_status = event_status,
    surv_type = as.character(attr(y, "type") %||% "right")
  )
}

#' Resolve the covariance-regression formula
#'
#' @param formulaVCov Covariance-regression formula.
#' @param default Default formula when the input is not supplied.
#' @param context Character label for error messages.
#'
#' @return A formula for the subject-level covariance regression.
#' @keywords internal
.resolve_vcov_formula <- function(formulaVCov = NULL,
                                  default = ~ 1,
                                  context = "joinme") {
  chosen <- formulaVCov %||% default
  if (!inherits(chosen, "formula")) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaVCov} must be a formula.",
      i = "Example: {.code ~ x1 + x2}."
    ))
  }
  chosen
}

#' Build the subject-level covariance-regression design matrix
#'
#' @param formulaVCov Covariance-regression formula.
#' @param dataEvent Event-level data with one row per subject.
#' @param time_var Longitudinal time variable name, forbidden in `formulaVCov`.
#' @param context Character label for error messages.
#'
#' @return Named list with `formulaVCov`, `K_cov`, and `Xcov`.
#' @keywords internal
.build_vcov_design <- function(formulaVCov,
                               dataEvent,
                               time_var,
                               context = "joinme") {
  if (length(reformulas::findbars(formulaVCov)) > 0) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaVCov} does not support random-effects terms.",
      i = "Remove all ( ... | ... ) terms from {.arg formulaVCov}."
    ))
  }

  fv_rhs <- stats::update(formulaVCov, . ~ .)
  fv_rhs[[2]] <- NULL
  if (length(fv_rhs) >= 3 && .expr_has_time(fv_rhs[[3]], time_var)) {
    cli::cli_abort(c(
      x = "{context}: {.arg formulaVCov} cannot include the time variable {.arg {time_var}}.",
      i = "Remove time from {.arg formulaVCov} or move it to longitudinal formulas."
    ))
  }

  Xtmp <- .mm(fv_rhs, dataEvent)
  if ("(Intercept)" %in% colnames(Xtmp)) {
    Xtmp <- Xtmp[, colnames(Xtmp) != "(Intercept)", drop = FALSE]
  }

  if (ncol(Xtmp) < 1L) {
    return(list(
      formulaVCov = formulaVCov,
      K_cov = 0L,
      Xcov = matrix(0.0, nrow(dataEvent), 0)
    ))
  }

  list(
    formulaVCov = formulaVCov,
    K_cov = as.integer(ncol(Xtmp)),
    Xcov = Xtmp
  )
}

#' Decode event status into event indicator and cause/type index
#'
#' @param status_raw Raw event-status vector.
#' @param context Character label for error messages.
#'
#' @return Named list with `d_event`, `event_type`, and `K_event`.
#' @keywords internal
.derive_event_outcomes <- function(status_raw, context = "joinme") {
  n <- length(status_raw)

  if (is.numeric(status_raw) || is.integer(status_raw) || is.logical(status_raw)) {
    status_num <- suppressWarnings(as.numeric(status_raw))
    status_chr <- as.character(status_raw)
    numeric_like <- TRUE
  } else {
    status_chr <- trimws(as.character(status_raw))
    status_num <- suppressWarnings(as.numeric(status_chr))
    numeric_like <- all(is.na(status_chr) == is.na(status_num))
  }

  if (numeric_like) {
    d_event <- as.integer(!is.na(status_num) & status_num != 0)
    event_type_raw <- status_num
  } else {
    status_chr_norm <- trimws(tolower(status_chr))
    is_censored <- status_chr_norm %in% c("", "0", "false", "no", "censor", "censored", "none", "na")
    d_event <- as.integer(!is.na(status_raw) & !is_censored)
    event_type_raw <- status_chr
  }

  event_type <- rep.int(1L, n)
  observed_types <- unique(event_type_raw[d_event == 1L & !is.na(event_type_raw)])
  observed_types <- observed_types[nzchar(as.character(observed_types))]

  if (length(observed_types) > 1L) {
    observed_types_chr <- as.character(observed_types)
    event_type[d_event == 1L] <- match(as.character(event_type_raw[d_event == 1L]), observed_types_chr)
    K_event <- length(observed_types_chr)
  } else {
    K_event <- 1L
  }

  if (any(!is.finite(d_event))) {
    cli::cli_abort(c(
      x = "{context}: failed to derive finite event indicators from status values.",
      i = "Ensure event status encodes censoring as 0/FALSE and events as non-zero values."
    ))
  }

  list(
    d_event = as.integer(d_event),
    event_type = as.integer(event_type),
    K_event = as.integer(K_event)
  )
}

#' Gauss-Kronrod quadrature grid helper
#'
#' Build a quadrature grid on `[0, 1]` for survival integration.
#'
#' @param nodes Positive integer node count. Allowed values are exactly
#'   `7`, `15`, `31`, `41`, `51`, and `61`. Defaults to `15`.
#'
#' @return A named list with components:
#' \\describe{
#'   \\item{n_gk}{Integer total node count.}
#'   \\item{nodes}{Numeric vector of nodes on `[0, 1]`.}
#'   \\item{weights}{Numeric vector of normalised weights summing to `1`.}
#'   \\item{panels}{Always `1L` (single fixed rule).}
#'   \\item{rule}{Character rule label (`"gk7"`, `"gk15"`, `"gk31"`, `"gk41"`, `"gk51"`, or `"gk61"`).}
#' }
#'
#' @details
#' Only fixed single-panel rules are supported.
#'
#' @export
gk_quadrature <- function(nodes = 15L) {
  spec <- .resolve_gk_request(nodes = nodes)
  .gk_single_panel(rule = spec$rule)
}

#' @keywords internal
.resolve_gk_request <- function(nodes = 15L) {
  supported_single <- c(7L, 15L, 31L, 41L, 51L, 61L)

  nodes <- as.integer(nodes)
  if (!is.finite(nodes) || length(nodes) != 1L || nodes < 1L) {
    cli::cli_abort(c(
      x = "{.arg nodes} must be a positive integer.",
      i = "Allowed values are {.code 7}, {.code 15}, {.code 31}, {.code 41}, {.code 51}, and {.code 61}."
    ))
  }

  if (!nodes %in% supported_single) {
    cli::cli_abort(c(
      x = "Unsupported quadrature node count: {.val {nodes}}.",
      i = "Allowed values are {.code 7}, {.code 15}, {.code 31}, {.code 41}, {.code 51}, and {.code 61}."
    ))
  }

  list(kind = "single", rule = nodes, nodes = nodes, panels = 1L, rounded = FALSE)
}

#' @keywords internal
.gk_nodes <- function(rule) {
  rule <- as.integer(rule)
  if (!rule %in% c(7L, 15L, 31L, 41L, 51L, 61L)) {
    cli::cli_abort(c(
      x = "Unsupported Gauss-Kronrod rule: {.val {rule}}.",
      i = "Supported single-panel rules are 7, 15, 31, 41, 51, and 61."
    ))
  }

  if (rule == 7L) {
    return(list(
      xi = c(
        0.0000000000000000,
        0.4342437493468026,
        0.7745966692414834,
        0.9604912687080203
      ),
      wi = c(
        0.4509165386584741,
        0.4013974147759622,
        0.2684880898683334,
        0.1046562260264673
      )
    ))
  }

  if (rule == 15L) {
    return(list(
      xi = c(
        0.0000000000000000,
        0.2077849550078985,
        0.4058451513773972,
        0.5860872354676911,
        0.7415311855993945,
        0.8648644233597691,
        0.9491079123427585,
        0.9914553711208126
      ),
      wi = c(
        0.2094821410847278,
        0.2044329400752989,
        0.1903505780647854,
        0.1690047266392679,
        0.1406532597155259,
        0.1047900103222502,
        0.06309209262997855,
        0.02293532201052922
      )
    ))
  }

  if (rule == 31L) {
    return(list(
      xi = c(
        0.0000000000000000,
        0.1011420669187175,
        0.2011940939974345,
        0.2991800071531688,
        0.3941513470775634,
        0.4850818636402397,
        0.5709721726085388,
        0.6509967412974170,
        0.7244177313601700,
        0.7904185014424659,
        0.8482065834104272,
        0.8972645323440819,
        0.9372733924007059,
        0.9677390756791391,
        0.9879925180204854,
        0.9980022986933971
      ),
      wi = c(
        0.1013300070147915,
        0.1007698455238756,
        0.0991735987217920,
        0.09664272698362368,
        0.09312659817082532,
        0.08856444305621177,
        0.08308050282313302,
        0.07684968075772038,
        0.06985412131872826,
        0.06200956780067064,
        0.05348152469092809,
        0.04458975132476488,
        0.03534636079137585,
        0.02546084732671532,
        0.01500794732931612,
        0.005377479872923349
      )
    ))
  }

  if (rule == 41L) {
    return(list(
      xi = c(
        0.0000000000000000,
        0.07652652113349733,
        0.1526054652409227,
        0.2277858511416451,
        0.3016278681149130,
        0.3737060887154196,
        0.4435931752387251,
        0.5108670019508271,
        0.5751404468197103,
        0.6360536807265150,
        0.6932376563347514,
        0.7463319064601508,
        0.7950414288375512,
        0.8391169718222188,
        0.8782768112522820,
        0.9122344282513259,
        0.9408226338317548,
        0.9639719272779138,
        0.9815078774502503,
        0.9931285991850949,
        0.9988590315882777
      ),
      wi = c(
        0.07660071191799966,
        0.07637786767208074,
        0.07570449768455667,
        0.07458287540049919,
        0.07303069033278667,
        0.07105442355344407,
        0.06864867292852162,
        0.06583459713361842,
        0.06265323755478117,
        0.05911140088063957,
        0.05519510534828600,
        0.05094457392372869,
        0.04643482186749767,
        0.04166887332797369,
        0.03660016975820080,
        0.03128730677703280,
        0.02588213360495116,
        0.02038837346126652,
        0.01462616925697125,
        0.008600269855642943,
        0.003073583718520531
      )
    ))
  }

  if (rule == 51L) {
    return(list(
      xi = c(
        0.0000000000000000,
        0.06154448300568508,
        0.1228646926107104,
        0.1837189394210489,
        0.2438668837209884,
        0.3030895389311078,
        0.3611723058093878,
        0.4178853821930377,
        0.4730027314457150,
        0.5263252843347192,
        0.5776629302412230,
        0.6268100990103174,
        0.6735663684734684,
        0.7177664068130844,
        0.7592592630373576,
        0.7978737979985001,
        0.8334426287608340,
        0.8658470652932756,
        0.8949919978782754,
        0.9207471152817016,
        0.9429745712289743,
        0.9616149864258425,
        0.9766639214595175,
        0.9880357945340773,
        0.9955569697904981,
        0.9992621049926098
      ),
      wi = c(
        0.06158081806783294,
        0.06147118987142532,
        0.06112850971705305,
        0.06053945537604586,
        0.05972034032417406,
        0.05868968002239421,
        0.05743711636156783,
        0.05595081122041232,
        0.05425112988854549,
        0.05236288580640748,
        0.05027767908071567,
        0.04798253713883671,
        0.04550291304992179,
        0.04287284502017005,
        0.04008382550403238,
        0.03711627148341554,
        0.03400213027432934,
        0.03079230016738749,
        0.02747531758785174,
        0.02400994560695322,
        0.02043537114588284,
        0.01684781770912830,
        0.01323622919557167,
        0.009473973386174152,
        0.005561932135356714,
        0.001987383892330316
      )
    ))
  }

  list(
    xi = c(
      0.0000000000000000,
      0.05147184255531770,
      0.1028069379667370,
      0.1538699136085835,
      0.2045251166823099,
      0.2546369261678898,
      0.3040732022736251,
      0.3527047255308781,
      0.4004012548303944,
      0.4470337695380892,
      0.4924804678617786,
      0.5366241481420199,
      0.5793452358263617,
      0.6205261829892429,
      0.6600610641266270,
      0.6978504947933158,
      0.7337900624532268,
      0.7677774321048262,
      0.7997278358218391,
      0.8295657623827684,
      0.8572052335460611,
      0.8825605357920527,
      0.9055733076999078,
      0.9262000474292743,
      0.9443744447485600,
      0.9600218649683075,
      0.9731163225011263,
      0.9836681232797472,
      0.9916309968704046,
      0.9968934840746495,
      0.9994844100504906
    ),
    wi = c(
      0.05149472942945157,
      0.05142612853745903,
      0.05122154784925877,
      0.05088179589874961,
      0.05040592140278235,
      0.04979568342707421,
      0.04905543455502978,
      0.04818586175708713,
      0.04718554656929915,
      0.04605923827100700,
      0.04481480013316266,
      0.04345253970135607,
      0.04196981021516425,
      0.04037453895153596,
      0.03867894562472759,
      0.03688236465182123,
      0.03497933802806002,
      0.03298144705748373,
      0.03090725756238776,
      0.02875404876504129,
      0.02650995488233310,
      0.02419116207808060,
      0.02182803582160919,
      0.01941414119394238,
      0.01692088918905327,
      0.01436972950704580,
      0.01182301525349634,
      0.009273279659517763,
      0.006630703915931292,
      0.003890461127099884,
      0.001389013698677008
    )
  )
}

#' @keywords internal
.gk_expand_rule <- function(xi_half, wi_half) {
  c(
    -rev(xi_half[-1L]),
    xi_half
  ) -> xi
  c(
    rev(wi_half[-1L]),
    wi_half
  ) -> wi
  list(xi = xi, wi = wi)
}

#' @keywords internal
.gk_single_panel <- function(rule = 15L) {
  tab <- .gk_nodes(rule)
  expanded <- .gk_expand_rule(tab$xi, tab$wi)

  nodes <- 0.5 * (expanded$xi + 1.0)
  weights <- 0.5 * expanded$wi

  list(
    n_gk = as.integer(length(nodes)),
    nodes = as.numeric(nodes),
    weights = as.numeric(weights / sum(weights)),
    panels = 1L,
    rule = paste0("gk", as.integer(rule))
  )
}

#' Safe model.matrix wrapper
#' @keywords internal
.mm <- function(formula, data) {
  # Always return double matrices for Stan compatibility
  X <- stats::model.matrix(formula, data = data)
  storage.mode(X) <- "double"
  X
}

#' Survival/event model matrix
#' @keywords internal
.mm_event <- function(formulaEvent, data) {
  rhs <- stats::delete.response(stats::terms(formulaEvent))
  X <- stats::model.matrix(rhs, data = data)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  storage.mode(X) <- "double"
  X
}

#' @keywords internal
.new_dist_scope <- function() {
  structure(list(default = NULL, by_family = list()), class = "joinme_dist_scope")
}

#' @keywords internal
.is_dist_scope <- function(x) {
  is.list(x) && all(c("default", "by_family") %in% names(x))
}

#' @keywords internal
.canonical_dist_param <- function(param_raw) {
  param <- tolower(trimws(as.character(param_raw)))
  if (param == "skew") param <- "alpha"
  if (param %in% c("alpha_skew", "skew_alpha")) param <- "alpha"
  if (param %in% c("phi_beta", "beta_phi", "precision")) param <- "phi_beta"
  if (param %in% c("tau", "tau_sde", "skew_sde")) param <- "tau_sde"
  param
}

#' @keywords internal
.canonical_family_name <- function(family_raw) {
  fam <- tolower(trimws(as.character(family_raw)))
  fam <- gsub("['\"]", "", fam)
  fam <- switch(fam,
    normal = "gaussian",
    gaussian = "gaussian",
    student = "student_t",
    "student-t" = "student_t",
    fam
  )
  fam_code <- tryCatch(.parse_family(fam), error = function(e) NA_integer_)
  if (!is.finite(fam_code)) {
    cli::cli_abort(c(
      x = "Unknown family in distributional formula scope: {.val {family_raw}}.",
      i = "Use a supported family name, e.g. gaussian, student_t, negbin2, skew_normal, skew_double_exponential."
    ))
  }
  .family_code_to_name(fam_code)
}

#' @keywords internal
.parse_dist_lhs <- function(lhs_expr) {
  lhs_txt <- paste(deparse(lhs_expr, width.cutoff = 500L), collapse = " ")
  lhs_compact <- gsub("\\s+", "", lhs_txt)

  m <- regexec("^([A-Za-z_][A-Za-z0-9_]*)(?:\\[(.*)\\])?$", lhs_compact)
  cap <- regmatches(lhs_compact, m)[[1]]
  if (length(cap) == 0) {
    cli::cli_abort(c(
      x = "Invalid distributional regression LHS: {.val {lhs_txt}}.",
      i = "Use {.code param ~ ...} or {.code param[family=name] ~ ...}."
    ))
  }

  param <- .canonical_dist_param(cap[2])
  if (!param %in% c("sigma", "nu", "phi", "alpha", "phi_beta", "tau_sde")) {
    cli::cli_abort(c(
      x = "Unknown distributional parameter: {.val {cap[2]}}.",
      i = "Allowed parameters: sigma, nu, phi, alpha (or alpha_skew/skew), phi_beta, tau_sde."
    ))
  }

  family_name <- NULL
  scope_txt <- cap[3]
  if (!is.na(scope_txt) && nzchar(scope_txt)) {
    scope_m <- regexec("^(family|fam)=([A-Za-z0-9_.-]+)$", scope_txt)
    scope_cap <- regmatches(scope_txt, scope_m)[[1]]
    if (length(scope_cap) == 0) {
      cli::cli_abort(c(
        x = "Invalid distributional formula scope: {.val [{scope_txt}]}",
        i = "Use {.code [family=<name>]} (example: {.code sigma[family=student_t] ~ 1 + time})."
      ))
    }
    family_name <- .canonical_family_name(scope_cap[3])
  }

  list(param = param, family = family_name)
}

#' Normalise distributional formula input
#' @keywords internal
.normalize_formula_dist <- function(formulaDist) {
  # Normalise list input to named distributional formulas
  if (is.null(formulaDist)) return(list())
  if (is.list(formulaDist) && length(formulaDist) > 0 &&
      all(vapply(formulaDist, .is_dist_scope, logical(1)))) {
    return(formulaDist)
  }
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

    lhs_info <- .parse_dist_lhs(f[[2]])
    param <- lhs_info$param
    family_name <- lhs_info$family

    if (has_name) {
      name_lhs <- stats::as.formula(paste0(nm, " ~ 1"))[[2]]
      name_info <- .parse_dist_lhs(name_lhs)
      if (!identical(name_info$param, param) || !identical(name_info$family, family_name)) {
        cli::cli_abort(c(
          x = "Distributional formula name '{nm}' does not match LHS '{paste(deparse(f[[2]]), collapse = ' ')}'.",
          i = "Use matching list names (including optional family scope), or leave entries unnamed."
        ))
      }
    }

    if (is.null(out[[param]])) out[[param]] <- .new_dist_scope()
    if (is.null(family_name)) {
      out[[param]]$default <- f
    } else {
      out[[param]]$by_family[[family_name]] <- f
    }
  }

  out
}

#' @keywords internal
.validate_dist_formula_scopes <- function(dist_formulas, family_names_present) {
  if (length(dist_formulas) == 0) return(invisible(TRUE))
  family_names_present <- unique(as.character(family_names_present %||% character(0)))

  for (param_name in names(dist_formulas)) {
    spec <- dist_formulas[[param_name]]
    if (!.is_dist_scope(spec)) next
    scoped_families <- names(spec$by_family %||% list())
    if (length(scoped_families) == 0) next

    missing_families <- setdiff(scoped_families, family_names_present)
    if (length(missing_families) > 0) {
      cli::cli_abort(c(
        x = "Distributional formula scope for {.val {param_name}} references family/families not present in data: {.val {paste(missing_families, collapse = ', ')}}.",
        i = "Available families: {.val {paste(family_names_present, collapse = ', ')}}."
      ))
    }

    bad_param_families <- scoped_families[!vapply(scoped_families, function(fm) {
      param_name %in% .family_distrib_params(.parse_family(fm))
    }, logical(1))]

    if (length(bad_param_families) > 0) {
      cli::cli_abort(c(
        x = "Distributional parameter {.val {param_name}} is not used by family/families: {.val {paste(bad_param_families, collapse = ', ')}}.",
        i = "Use supported combinations only (e.g., nu for student_t, phi for negbin2, phi_beta for beta)."
      ))
    }
  }

  invisible(TRUE)
}

#' @keywords internal
.family_by_row_from_marker <- function(marker_values, marker_levels, family_codes) {
  marker_values <- as.character(marker_values)
  marker_levels <- as.character(marker_levels)
  family_codes <- as.integer(family_codes)
  if (length(marker_levels) != length(family_codes)) {
    cli::cli_abort(c(
      x = "Marker levels and family code lengths do not match.",
      i = "Cannot derive row-level family labels for distributional formulas."
    ))
  }
  family_names <- vapply(family_codes, .family_code_to_name, character(1))
  fam_map <- setNames(family_names, marker_levels)
  fam <- unname(fam_map[marker_values])
  if (any(is.na(fam))) {
    cli::cli_abort(c(
      x = "Cannot map marker value(s) to family labels: {.val {paste(sort(unique(marker_values[is.na(fam)])), collapse = ', ')}}.",
      i = "Ensure marker levels align with fitted marker-family mapping."
    ))
  }
  fam
}

#' Build distributional design matrix
#' @keywords internal
.build_dist_matrix <- function(formula, data, family_by_row = NULL) {
  # Construct fixed-effect matrix for distributional regression
  if (is.null(formula)) {
    return(list(P = 0L, X = matrix(0.0, nrow(data), 0), cols = character(0)))
  }

  if (.is_dist_scope(formula)) {
    x_parts <- list()
    col_parts <- character(0)

    if (!is.null(formula$default)) {
      base_default <- .build_dist_matrix(formula$default, data)
      if (base_default$P > 0) {
        Xd <- base_default$X
        colnames(Xd) <- paste0("all::", base_default$cols)
        x_parts <- c(x_parts, list(Xd))
        col_parts <- c(col_parts, colnames(Xd))
      }
    }

    fam_specs <- formula$by_family %||% list()
    if (length(fam_specs) > 0) {
      if (is.null(family_by_row) || length(family_by_row) != nrow(data)) {
        cli::cli_abort(c(
          x = "Family-scoped distributional formulas require {.arg family_by_row} aligned with data rows.",
          i = "Provide one family label per row when building scoped distributional matrices."
        ))
      }
      family_by_row <- as.character(family_by_row)

      for (family_name in names(fam_specs)) {
        base_fam <- .build_dist_matrix(fam_specs[[family_name]], data)
        if (base_fam$P == 0) next
        gate <- as.numeric(family_by_row == family_name)
        Xf <- base_fam$X * gate
        colnames(Xf) <- paste0("family=", family_name, "::", base_fam$cols)
        x_parts <- c(x_parts, list(Xf))
        col_parts <- c(col_parts, colnames(Xf))
      }
    }

    if (length(x_parts) == 0) {
      return(list(P = 0L, X = matrix(0.0, nrow(data), 0), cols = character(0)))
    }
    X <- do.call(cbind, x_parts)
    storage.mode(X) <- "double"
    return(list(P = ncol(X), X = X, cols = col_parts))
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
  # Parse random-effect terms for distributional regression
  if (is.null(formula)) {
    return(list(n_re = 0L, K = integer(0), G = integer(0), Z = list(), J = list(), W = list(), terms = character(0)))
  }

  if (.is_dist_scope(formula)) {
    merge_re <- function(base, add) {
      if (is.null(add) || add$n_re == 0L) return(base)
      list(
        n_re = base$n_re + add$n_re,
        K = c(base$K, add$K),
        G = c(base$G, add$G),
        Z = c(base$Z, add$Z),
        J = c(base$J, add$J),
        W = c(base$W, add$W),
        terms = c(base$terms, add$terms)
      )
    }

    out <- list(n_re = 0L, K = integer(0), G = integer(0), Z = list(), J = list(), W = list(), terms = character(0))

    if (!is.null(formula$default)) {
      out <- merge_re(out, .build_dist_re_terms(formula$default, data))
    }

    fam_specs <- formula$by_family %||% list()
    if (length(fam_specs) > 0) {
      family_by_row <- attr(data, "joinme_family_by_row", exact = TRUE)
      if (is.null(family_by_row) || length(family_by_row) != nrow(data)) {
        cli::cli_abort(c(
          x = "Family-scoped distributional random effects require row-level family labels.",
          i = "Attach {.code attr(data, 'joinme_family_by_row')} before parsing random effects."
        ))
      }
      family_by_row <- as.character(family_by_row)

      for (family_name in names(fam_specs)) {
        cur <- .build_dist_re_terms(fam_specs[[family_name]], data)
        if (cur$n_re == 0L) next
        gate <- as.numeric(family_by_row == family_name)
        cur$Z <- lapply(cur$Z, function(Zm) {
          Zm * gate
        })
        cur$terms <- paste0("family=", family_name, "::", cur$terms)
        out <- merge_re(out, cur)
      }
    }
    return(out)
  }

  if (!inherits(formula, "formula")) formula <- stats::as.formula(formula)

  f_exp <- reformulas::expandDoubleVerts(formula)
  bars <- reformulas::findbars(f_exp)
  if (length(bars) == 0) {
    return(list(n_re = 0L, K = integer(0), G = integer(0), Z = list(), J = list(), W = list(), terms = character(0)))
  }

  Z_list <- list()
  J_list <- list()
  W_list <- list()
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

    grp_info <- .parse_weighted_group_expr(
      grp_expr = grp_expr,
      data = data,
      context = "formulaDist random-effects term"
    )
    grp_eval_expr <- grp_info$group_expr
    obs_weights <- grp_info$weights

    grp_df <- stats::model.frame(
      stats::as.formula(paste0("~ ", paste(deparse(grp_eval_expr, width.cutoff = 500L), collapse = " "))),
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
    W <- as.numeric(tapply(obs_weights, grp, mean))
    if (length(W) != nlevels(grp) || any(!is.finite(W)) || any(W <= 0)) {
      cli::cli_abort(c(
        x = "Invalid group weights in distributional random-effects term.",
        i = "Weights declared via {.code weighted(..., weights = ...)} must be finite and strictly positive within each grouping level."
      ))
    }

    Z_list[[length(Z_list) + 1L]] <- Z
    J_list[[length(J_list) + 1L]] <- J
    W_list[[length(W_list) + 1L]] <- W
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
    W = W_list,
    terms = term_names
  )
}

#' Pad random-effects term matrices to max dimensions
#' @keywords internal
.pad_re_terms <- function(re_terms, n_rows) {
  # Pad variable-size RE terms into aligned arrays for Stan
  if (re_terms$n_re == 0) {
    return(list(
      n_re = 0L,
      K = integer(0),
      G = integer(0),
      K_max = 0L,
      G_max = 0L,
      Z = list(),
      J = list(),
      W = list(),
      J_mat = matrix(0L, nrow = 0, ncol = n_rows),
      W_mat = matrix(1.0, nrow = 0, ncol = 0),
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
  W_padded <- lapply(seq_len(re_terms$n_re), function(i) {
    w <- as.numeric(re_terms$W[[i]])
    if (length(w) != re_terms$G[[i]] || any(!is.finite(w)) || any(w <= 0)) {
      cli::cli_abort(c(
        x = "Invalid group-weight vector in random-effects terms.",
        i = "Each random-effects term must have one finite, strictly positive weight per group level."
      ))
    }
    if (length(w) < G_max) {
      w <- c(w, rep(1.0, G_max - length(w)))
    }
    w
  })
  J_mat <- do.call(rbind, lapply(re_terms$J, function(j) as.integer(j)))
  W_mat <- do.call(rbind, lapply(W_padded, as.numeric))
  list(
    n_re = re_terms$n_re,
    K = re_terms$K,
    G = re_terms$G,
    K_max = K_max,
    G_max = G_max,
    Z = Z_padded,
    J = re_terms$J,
    W = W_padded,
    J_mat = J_mat,
    W_mat = W_mat,
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
    assoc_corr      = as.integer("corr" %in% assoc),
    assoc_vcov      = as.integer("vcov" %in% assoc)
  )
}

#' Validate association channel combinations
#' @keywords internal
.validate_assoc_channels <- function(assoc, context = "association specification") {
  assoc <- unique(as.character(assoc %||% character(0)))
  if (all(c("corr", "vcov") %in% assoc)) {
    cli::cli_abort(c(
      x = "{.arg corr} and {.arg vcov} cannot be used together in {.field {context}}.",
      i = "Choose {.arg corr} for off-diagonal correlation features or {.arg vcov} for lower-triangular Cholesky-factor features from {.arg L}."
    ))
  }
  assoc
}

#' Covariance-association feature count
#' @keywords internal
.assoc_cov_feature_count <- function(q_idm, include_diag = FALSE) {
  q_idm <- as.integer(q_idm %||% 0L)
  if (q_idm <= 0L) {
    return(0L)
  }
  if (isTRUE(include_diag)) {
    return(as.integer((q_idm * (q_idm + 1L)) %/% 2L))
  }
  if (q_idm < 2L) {
    return(0L)
  }
  as.integer((q_idm * (q_idm - 1L)) %/% 2L)
}

#' Association transform component count
#' @keywords internal
.assoc_transform_component_count <- function(term_key, q_idm, diagonal_only = FALSE) {
  term_key <- as.character(term_key %||% "")[1]
  q_idm <- as.integer(q_idm %||% 0L)
  diagonal_only <- isTRUE(diagonal_only)

  if (identical(term_key, "corr")) {
    return(.assoc_cov_feature_count(q_idm, include_diag = FALSE))
  }
  if (identical(term_key, "vcov")) {
    if (q_idm <= 0L) {
      return(0L)
    }
    if (diagonal_only) {
      return(q_idm)
    }
    return(.assoc_cov_feature_count(q_idm, include_diag = TRUE))
  }
  1L
}

#' Association transform component labels
#' @keywords internal
.assoc_transform_component_labels <- function(term_key, n_components) {
  term_key <- as.character(term_key %||% "")[1]
  n_components <- as.integer(n_components %||% 0L)
  if (!(term_key %in% c("corr", "vcov")) || n_components <= 0L) {
    return(term_key)
  }
  paste0(term_key, "[", seq_len(n_components), "]")
}

#' Covariance-association feature map
#' @keywords internal
.assoc_cov_feature_map <- function(q_idm, diagonal_only = FALSE, include_diag = TRUE) {
  q_idm <- as.integer(q_idm %||% 0L)
  if (q_idm <= 0L) {
    return(matrix(integer(0), ncol = 2L))
  }
  if (isTRUE(diagonal_only)) {
    idx <- seq_len(q_idm)
    return(cbind(row = idx, col = idx))
  }

  rows <- list()
  pos <- 1L
  for (r in seq_len(q_idm)) {
    c_start <- if (isTRUE(include_diag)) 1L else 1L
    c_end <- if (isTRUE(include_diag)) r else r - 1L
    if (c_end < c_start) next
    for (c in seq.int(c_start, c_end)) {
      rows[[pos]] <- c(r, c)
      pos <- pos + 1L
    }
  }
  do.call(rbind, rows)
}

#' Extract raw vcov-association features from a Cholesky factor
#' @keywords internal
.assoc_vcov_features_from_chol <- function(L_i, diagonal_only = FALSE) {
  L_i <- as.matrix(L_i)
  q_idm <- nrow(L_i)
  if (!q_idm || ncol(L_i) != q_idm) {
    return(numeric(0))
  }
  if (isTRUE(diagonal_only)) {
    return(diag(L_i))
  }

  out <- numeric(.assoc_cov_feature_count(q_idm, include_diag = TRUE))
  pos <- 1L
  for (r in seq_len(q_idm)) {
    for (c in seq_len(r)) {
      out[pos] <- L_i[r, c]
      pos <- pos + 1L
    }
  }
  out
}

#' Validate transform flags
#'
#' Transform codes: 0 id, 1 exp, 2 log, 3 inv_logit, 4 probit, 5 sqrt, 6 cbrt.
#' Restrictions: log and sqrt are restricted to cv_mean (legacy constraint).
#'
#' @keywords internal
.validate_tf <- function(tf) {
  nm <- c("cv_mean", "cv_marker", "cs_mean", "cs_marker", "corr", "vcov")
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
  bad <- c(tf$cv_marker, tf$cs_mean, tf$cs_marker, tf$corr, tf$vcov)
  if (any(bad %in% c(2L, 5L))) {
    cli::cli_abort(c(
      x = "log/sqrt transforms are restricted to cv_mean only.",
      i = "Use identity or other transforms for cs_mean/cs_marker/corr/vcov."
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

#' Centre baseline hazard basis columns
#' @keywords internal
.center_baseline <- function(Bs_event_raw, Bs_gk_raw) {
  n_id <- nrow(Bs_event_raw)
  Kbs <- ncol(Bs_event_raw)
  n_gk <- dim(Bs_gk_raw)[2]
  big <- rbind(
    Bs_event_raw,
    matrix(aperm(Bs_gk_raw, c(1, 3, 2)), nrow = n_id * n_gk, ncol = Kbs)
  )
  colm <- colMeans(big)
  # Keep constant (intercept) columns uncentred so the baseline level is preserved.
  col_sd <- apply(big, 2, stats::sd)
  center_mask <- !(col_sd < 1e-8 | !is.finite(col_sd))
  colm_use <- ifelse(center_mask, colm, 0)
  Bs_event_c <- sweep(Bs_event_raw, 2, colm_use, "-")
  Bs_gk_c <- sweep(Bs_gk_raw, MARGIN = 3, STATS = colm_use, FUN = "-")
  list(Bs_event_c = Bs_event_c, Bs_gk_c = Bs_gk_c, col_means = colm)
}

#' Construct a zero-dimension marker-only block
#' @keywords internal
.zero_marker_block <- function(N, n_id, n_gk = 15L) {
  n_gk <- as.integer(n_gk)
  list(
    R_mk = 0L,
    Z_mk_obs = matrix(0.0, N, 0),
    Z_mk_gk_now = array(0.0, dim = c(n_id, n_gk, 0)),
    Z_mk_gk_fwd = array(0.0, dim = c(n_id, n_gk, 0)),
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
    suppressWarnings(
      splines::bs(x, knots = knots, Boundary.knots = boundary, degree = degree, intercept = TRUE)
    )
  } else {
    suppressWarnings(
      splines::ns(x, knots = knots, Boundary.knots = boundary, intercept = TRUE)
    )
  }
}

#' Resolve grouping factor names from random-effect terms
#'
#' @param grp_expr Grouping expression from a random-effects term.
#' @return Character scalar with the base grouping name.
#'
#' @keywords internal
.group_name_from_expr <- function(grp_expr) {
  if (is.call(grp_expr) && identical(as.character(grp_expr[[1]]), "weighted")) {
    if (length(grp_expr) < 2) {
      cli::cli_abort(c(
        x = "Invalid weighted grouping expression.",
        i = "Use {.code weighted(group, weights = weight_column)}."
      ))
    }
    return(.group_name_from_expr(grp_expr[[2]]))
  }
  if (is.call(grp_expr) && as.character(grp_expr)[[1]] %in% c(":", "*")) {
    # Handle marker:id or marker*id patterns - extract first part
    return(deparse(grp_expr[[2]], width.cutoff = 500L)[1])
  }
  # Simple symbol or name
  as.character(grp_expr)[[1]]
}

#' Parse weighted grouping expressions used in random-effect terms
#' @keywords internal
.parse_weighted_group_expr <- function(grp_expr, data, context = "random-effects term") {
  if (!(is.call(grp_expr) && identical(as.character(grp_expr[[1]]), "weighted"))) {
    return(list(
      group_expr = grp_expr,
      is_weighted = FALSE,
      weights = rep(1.0, nrow(data)),
      weight_label = NULL
    ))
  }

  if (length(grp_expr) < 2) {
    cli::cli_abort(c(
      x = "Invalid weighted grouping expression in {context}.",
      i = "Use {.code weighted(group, weights = weight_column)}."
    ))
  }

  args <- as.list(grp_expr)
  arg_names <- names(args)
  weight_pos <- which(arg_names == "weights")
  if (length(weight_pos) != 1) {
    cli::cli_abort(c(
      x = "Invalid weighted grouping expression in {context}.",
      i = "Provide exactly one {.code weights = <column>} argument."
    ))
  }

  base_group_expr <- args[[2]]
  weight_expr <- args[[weight_pos]]
  weight_formula <- stats::as.formula(paste0("~ ", paste(deparse(weight_expr, width.cutoff = 500L), collapse = " ")))
  weight_df <- stats::model.frame(weight_formula, data = data, na.action = stats::na.fail)
  if (ncol(weight_df) != 1) {
    cli::cli_abort(c(
      x = "Invalid weight specification in weighted grouping expression.",
      i = "Weight must evaluate to a single numeric column."
    ))
  }

  weight_vec <- as.numeric(weight_df[[1]])
  if (length(weight_vec) != nrow(data) || any(!is.finite(weight_vec)) || any(weight_vec <= 0)) {
    cli::cli_abort(c(
      x = "Invalid weights in weighted grouping expression.",
      i = "Weights must be finite, strictly positive, and aligned with data rows."
    ))
  }

  list(
    group_expr = base_group_expr,
    is_weighted = TRUE,
    weights = weight_vec,
    weight_label = paste(deparse(weight_expr, width.cutoff = 500L), collapse = " ")
  )
}

#' Resolve independence flags from double-bar syntax
#'
#' @description
#' Uses lme4-style `||` parsing to determine which random-effects blocks should
#' be modelled with diagonal covariance structures. The resolution rules are:
#' - Top-level `(... || id)` sets `indep_id_re = 1`.
#' - `(... || marker)` sets `indep_marker_re = 1`.
#' - Nested `(... || id)` inside a marker block sets `indep_idmarker_cov = 1`.
#'
#' @param formulaLong Longitudinal formula with random-effects terms.
#' @param marker_var Marker grouping variable name.
#' @param id_var Subject grouping variable name.
#' @return Named list with `indep_id_re`, `indep_marker_re`, and `indep_idmarker_cov`.
#'
#' @keywords internal
.resolve_re_independence <- function(formulaLong, marker_var, id_var) {
  # Collect random-effects terms with their operator and context so that
  # nested marker-by-id specifications can be distinguished from top-level id blocks.
  .collect_re_terms <- function(expr, context = NULL) {
    terms <- list()
    if (!is.call(expr)) {
      return(terms)
    }

    op <- as.character(expr[[1]])
    if (op %in% c("|", "||")) {
      grp <- .group_name_from_expr(expr[[3]])
      terms <- c(terms, list(list(op = op, group = grp, context = context)))
      # Traverse the left-hand side within the current group context.
      return(c(terms, .collect_re_terms(expr[[2]], context = grp)))
    }

    for (i in seq_along(expr)[-1]) {
      terms <- c(terms, .collect_re_terms(expr[[i]], context = context))
    }
    terms
  }

  rhs <- formulaLong[[3]]
  terms <- .collect_re_terms(rhs, context = NULL)
  if (length(terms) == 0) {
    return(list(indep_id_re = 0L, indep_marker_re = 0L, indep_idmarker_cov = 0L))
  }

  indep_id_re <- any(vapply(
    terms,
    function(t) t$op == "||" && t$group == id_var && (is.null(t$context) || t$context != marker_var),
    logical(1)
  ))
  indep_marker_re <- any(vapply(
    terms,
    function(t) t$op == "||" && t$group == marker_var,
    logical(1)
  ))
  indep_idmarker_cov <- any(vapply(
    terms,
    function(t) t$op == "||" && t$group == id_var && !is.null(t$context) && t$context == marker_var,
    logical(1)
  ))

  list(
    indep_id_re = as.integer(indep_id_re),
    indep_marker_re = as.integer(indep_marker_re),
    indep_idmarker_cov = as.integer(indep_idmarker_cov)
  )
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
    return(list(mk_rhs_list = list(), idm_rhs_list = list(), idm_group_exprs = list(), marker_terms = list()))
  }

  expr_has_literal <- function(expr, value) {
    if (is.null(expr)) return(FALSE)
    if (is.numeric(expr) && length(expr) == 1) return(identical(as.numeric(expr), value))
    if (!is.call(expr)) return(FALSE)
    any(vapply(as.list(expr)[-1], expr_has_literal, logical(1), value = value))
  }

  expr_is_zero <- function(expr) {
    is.numeric(expr) && length(expr) == 1 && identical(as.numeric(expr), 0)
  }
  
  # Extract grouping variable name - handle both symbols and compound expressions
  grp <- vapply(bars, function(b) .group_name_from_expr(b[[3]]), character(1))
  
  mk_idx <- which(grp == marker_var)
  if (length(mk_idx) == 0) {
    return(list(mk_rhs_list = list(), idm_rhs_list = list(), idm_group_exprs = list(), marker_terms = list()))
  }

  mk_rhs_list <- list()
  idm_rhs_list <- list()
  idm_group_exprs <- list()
  marker_terms <- list()

  for (j in mk_idx) {
    bt <- bars[[j]]
    outer_expr <- bt[[2]]
    marker_terms <- c(marker_terms, list(outer_expr))

    f_inner <- stats::as.formula(paste0("~ ", paste(deparse(outer_expr, width.cutoff = 500L), collapse = " ")))
    f_inner <- reformulas::expandDoubleVerts(f_inner)

    inner_bars <- reformulas::findbars(f_inner)
    if (length(inner_bars) > 0) {
      inner_grp <- vapply(inner_bars, function(b) .group_name_from_expr(b[[3]]), character(1))
      id_inner_idx <- which(inner_grp == id_var)
      if (length(id_inner_idx) > 0) {
        idm_rhs_list <- c(idm_rhs_list, .bar_terms_to_rhs_list(inner_bars[id_inner_idx]))
        idm_group_exprs <- c(idm_group_exprs, lapply(inner_bars[id_inner_idx], function(b) b[[3]]))
      }
    }

    f_mk <- reformulas::nobars(f_inner)
    rhs_terms <- terms(f_mk)
    term_labels <- attr(rhs_terms, "term.labels")
    has_only_intercept <- (length(term_labels) == 0) && (attr(rhs_terms, "intercept") == 1)

    user_explicit_intercept <- expr_has_literal(outer_expr, 1)

    if (has_only_intercept && !user_explicit_intercept) next

    rhs_expr <- f_mk[[2]]
    if (expr_is_zero(rhs_expr)) next

    mk_rhs_list <- c(mk_rhs_list, list(f_mk))
  }

  list(mk_rhs_list = mk_rhs_list, idm_rhs_list = idm_rhs_list, idm_group_exprs = idm_group_exprs, marker_terms = marker_terms)
}

# ---- time-index metadata for internal scaling in Stan -----------------

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

#' Materialize CmdStanR fit data in memory
#' @keywords internal
.materialize_cmdstanr_fit <- function(fit) {
  if (!.is_cmdstanr_fit(fit)) {
    return(fit)
  }

  # Mirror cmdstanr::save_object() so a later saveRDS() does not depend on the
  # original CSV files still being present.
  fit$draws()
  tryCatch(fit$sampler_diagnostics(), error = function(e) NULL)
  tryCatch(fit$init(), error = function(e) NULL)
  tryCatch(fit$profiles(), error = function(e) NULL)
  fit
}

#' Unified draws object helper
#' @keywords internal
.get_draws_obj <- function(fit, variables = NULL, draws = NULL, seed = 1, keep_chains = FALSE) {
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
  if (!is.null(draws) && is.finite(draws) && !isTRUE(keep_chains)) {
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

#' Unified draws array helper (keeps chains separate)
#' @keywords internal
.get_draws_array <- function(fit, variables = NULL, draws = NULL, seed = 1) {
  d <- .get_draws_obj(fit, variables = variables, draws = NULL, seed = seed, keep_chains = TRUE)
  arr <- posterior::as_draws_array(d)
  if (!is.null(draws) && is.finite(draws)) {
    n_iter <- dim(arr)[1]
    if (draws < n_iter) {
      set.seed(seed)
      idx <- sample.int(n_iter, size = draws)
      arr <- arr[idx, , , drop = FALSE]
    }
  }
  arr
}

# ---- Posterior draw helpers for methods -----------------------------------

#' Subset draws safely before conversion to draws_df
#' @keywords internal
.get_draws_df <- function(fit, variables, draws = NULL, seed = 1) {
  d <- .get_draws_obj(fit, variables = variables, draws = draws, seed = seed)
  posterior::as_draws_df(d)
}

#' Summarise a numeric vector of draws
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
    if (grepl("^[[:space:]]*#include[[:space:]]+", line)) {
      inc <- sub("^[[:space:]]*#include[[:space:]]+\"?([^\"[:space:]]+)\"?.*$", "\\1", line)
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
  sd$re_weight_sigma <- make_num_array(sd$re_weight_sigma, c(sd$n_re_sigma %||% 0L, sd$G_sigma_max %||% 0L))
  sd$Z_nu <- make_num_array(sd$Z_nu, c(sd$n_re_nu %||% 0L, N, sd$K_nu_max %||% 0L))
  sd$J_nu <- make_int_array(sd$J_nu, c(sd$n_re_nu %||% 0L, N))
  sd$re_weight_nu <- make_num_array(sd$re_weight_nu, c(sd$n_re_nu %||% 0L, sd$G_nu_max %||% 0L))
  sd$Z_phi <- make_num_array(sd$Z_phi, c(sd$n_re_phi %||% 0L, N, sd$K_phi_max %||% 0L))
  sd$J_phi <- make_int_array(sd$J_phi, c(sd$n_re_phi %||% 0L, N))
  sd$re_weight_phi <- make_num_array(sd$re_weight_phi, c(sd$n_re_phi %||% 0L, sd$G_phi_max %||% 0L))
  sd$Z_alpha <- make_num_array(sd$Z_alpha, c(sd$n_re_alpha %||% 0L, N, sd$K_alpha_max %||% 0L))
  sd$J_alpha <- make_int_array(sd$J_alpha, c(sd$n_re_alpha %||% 0L, N))
  sd$re_weight_alpha <- make_num_array(sd$re_weight_alpha, c(sd$n_re_alpha %||% 0L, sd$G_alpha_max %||% 0L))
  sd$Z_phi_beta <- make_num_array(sd$Z_phi_beta, c(sd$n_re_phi_beta %||% 0L, N, sd$K_phi_beta_max %||% 0L))
  sd$J_phi_beta <- make_int_array(sd$J_phi_beta, c(sd$n_re_phi_beta %||% 0L, N))
  sd$re_weight_phi_beta <- make_num_array(sd$re_weight_phi_beta, c(sd$n_re_phi_beta %||% 0L, sd$G_phi_beta_max %||% 0L))
  sd$Z_tau_sde <- make_num_array(sd$Z_tau_sde, c(sd$n_re_tau_sde %||% 0L, N, sd$K_tau_sde_max %||% 0L))
  sd$J_tau_sde <- make_int_array(sd$J_tau_sde, c(sd$n_re_tau_sde %||% 0L, N))
  sd$re_weight_tau_sde <- make_num_array(sd$re_weight_tau_sde, c(sd$n_re_tau_sde %||% 0L, sd$G_tau_sde_max %||% 0L))
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

#' Safe with progress  helper
#' @keywords internal
.safe_progress <- function(show_progress = FALSE, expr) {
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  if (has_progressr && show_progress) {
    with_progress_f <- get("with_progress", envir = asNamespace("progressr"))
    args <- list(expr = substitute(expr))
    do.call(with_progress_f, args, envir = parent.frame())
  } else {
    eval(expr, envir = parent.frame())
  }
}
