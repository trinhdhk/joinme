#' Family and Transformation Utilities for Mixed-Family Joint Models
#'
#' @description
#' Helper functions for specifying mixed families per marker and transformation compositions.
#'
#' @name family_utils
#' @noRd
NULL

#' JoiNMe longitudinal family specification
#'
#' @description
#' Creates a marker-specific family specification for
#' `joinme(..., families = ...)`. A specification may define a custom
#' longitudinal link (or inverse link) and, for the skew-Laplace family, a
#' fixed quantile parameter.
#'
#' Character values supplied through `link` name the *forward* link
#' \eqn{g(\mu)}.  Formula values supplied through `link` are inverted
#' symbolically and then compiled by the same functional bytecode evaluator
#' used for custom association transformations.  In contrast, `inv_link`
#' directly specifies \eqn{g^{-1}(\eta)}.
#'
#' The named forward links and their inverse links are:
#' \describe{
#'   \item{`"identity"`}{\eqn{g(\mu)=\mu} and
#'     \eqn{g^{-1}(\eta)=\eta}.}
#'   \item{`"log"`}{\eqn{g(\mu)=\log(\mu)} and
#'     \eqn{g^{-1}(\eta)=\exp(\eta)}.}
#'   \item{`"logit"`}{\eqn{g(\mu)=\operatorname{logit}(\mu)} and
#'     \eqn{g^{-1}(\eta)=\operatorname{logit}^{-1}(\eta)}.}
#'   \item{`"probit"`}{\eqn{g(\mu)=\Phi^{-1}(\mu)} and
#'     \eqn{g^{-1}(\eta)=\Phi(\eta)}.}
#'   \item{`"exp"`}{\eqn{g(\mu)=\exp(\mu)} and
#'     \eqn{g^{-1}(\eta)=\log(\eta)}.}
#' }
#'
#' @param name Character scalar family name, e.g. `"student_t"`, `"bernoulli"`.
#' @param link Optional character scalar naming a forward link, or a one-sided
#'   formula defining a forward link in `x`, for example `~ log(x + 1)`.
#'   Formula links must be one-to-one compositions that can be inverted
#'   algebraically by [invert_transform_expr()].
#' @param inv_link Optional one-sided formula `~ ...` or expression string
#'   using `x` (e.g. `~ exp(x)`, `~ inv_logit(x)`, or `~ Phi(x)`).
#' @param tau Optional fixed quantile/asymmetry parameter for
#'   `"skew_laplace"` or `"skew_double_exponential"`. It must be a finite
#'   scalar strictly between zero and one. Under Stan's quantile
#'   parameterisation, `tau = 0.5` is the symmetric Laplace distribution.
#'   Values zero and one are excluded because they yield a degenerate,
#'   non-normalisable limiting distribution. When omitted, `tau` is estimated
#'   from a distributional regression or as a family-level parameter.
#'
#' @return Object of class `"JoiNMe_family_spec"`.
#' @examples
#' jm_family("poisson", link = "log")
#' jm_family("poisson", link = ~ log(x))
#' jm_family("bernoulli", inv_link = ~ inv_logit(x))
#' jm_family("bernoulli", link = ~ inv_Phi(x))
#' jm_family("bernoulli", inv_link = ~ Phi(x))
#' jm_family("skew_laplace", tau = 0.8)
#' @aliases jm_family
#' @export
joinme_family <- function(name, link = NULL, inv_link = NULL, tau = NULL) {
  fam_code <- .parse_family(name)
  fam_name <- .family_code_to_name(fam_code)
  tau_fixed <- .validate_family_tau(
    tau = tau,
    family_code = fam_code,
    context = "jm_family()"
  )

  if (!is.null(link) && !is.null(inv_link)) {
    cli::cli_abort(c(
      x = "Specify only one of {.arg link} or {.arg inv_link}.",
      i = "Use {.code link = 'logit'} or {.code inv_link = ~ inv_logit(x)}."
    ))
  }

  inv_link_bc <- if (!is.null(inv_link)) {
    parse_transform_expr(inv_link)
  } else if (!is.null(link)) {
    .inv_link_bc_from_link(link)
  } else {
    .inv_link_bc_from_name(.default_link_for_family(fam_code))
  }

  if (!is.null(inv_link)) {
    .validate_noninvertible(
      inv_link_bc,
      context = paste0("jm_family(family = '", fam_name, "')")
    )
  }

  link_name <- .canonical_link_from_inv_link_bc(inv_link_bc)

  structure(
    list(
      family = fam_name,
      link = link_name,
      inv_link = inv_link_bc,
      tau = tau_fixed
    ),
    class = "JoiNMe_family_spec"
  )
}

#' @rdname joinme_family
#' @export
jm_family <- joinme_family


#' Method for printing a `JoiNMe_family_spec` object
#' 
#' @description Print out the spec details of a `JoiNMe_family_spec` object.
#' @param x An object of class `JoiNMe_family_spec`.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the input object `x`.
#' @method print JoiNMe_family_spec
#' @export
print.JoiNMe_family_spec <- function(x, ...){
  cat('Joint Nested Mixed-effects Longitudinal Family Specification:\n')
  cat('  Family: ', x$family, '\n', sep = '')
  cat('  Link: ', x$link, '\n', sep = '')
  cat('  Inverse Link: ', x$inv_link, '\n', sep = '')
  invisible(x)
}

#' Validate a fixed skew-Laplace quantile parameter
#'
#' @description
#' Validates the marker-specific constant used as the quantile/asymmetry
#' parameter of Stan's skew double exponential distribution. The validation is
#' centralised so formal `JoiNMe_family_spec` objects and compatible list-style
#' declarations obey identical statistical constraints.
#'
#' @param tau Optional numeric value supplied in a family declaration.
#' @param family_code Integer longitudinal-family code.
#' @param context Character description of the calling family declaration,
#'   used to produce a precise diagnostic.
#'
#' @return `NA_real_` when `tau` is absent; otherwise the validated numeric
#'   scalar.
#' @keywords internal
#' @noRd
.validate_family_tau <- function(tau, family_code, context = "family specification") {
  if (
    is.null(tau) ||
      (is.numeric(tau) && length(tau) == 1L && is.na(tau))
  ) {
    return(NA_real_)
  }

  if (!identical(as.integer(family_code), 9L)) {
    cli::cli_abort(c(
      x = "{.arg tau} is only defined for the skew-Laplace family in {context}.",
      i = "Use {.code jm_family('skew_laplace', tau = 0.5)}, or remove {.arg tau}."
    ))
  }

  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau)) {
    cli::cli_abort(c(
      x = "{.arg tau} in {context} must be one finite numeric value.",
      i = "Provide a quantile strictly between zero and one."
    ))
  }

  tau <- as.numeric(tau)
  if (tau <= 0 || tau >= 1) {
    cli::cli_abort(c(
      x = "{.arg tau} in {context} must lie strictly between zero and one.",
      i = "{.code tau = 1} is not a valid skew-Laplace distribution; {.code tau = 0.5} is symmetric."
    ))
  }

  tau
}

#' Normalise fixed skew-Laplace metadata by marker
#'
#' @description
#' Produces the marker-aligned fixed-`tau` vectors required by fitting and
#' dynamic prediction. Scalar values from older fitted objects are recycled so
#' that saved models remain usable after the family API became marker-specific.
#' Flags are cleared for non-skew-Laplace markers because `tau` has no
#' likelihood interpretation for those response families.
#'
#' @param use_tau_fixed Integer or logical flag, either scalar or one value per
#'   marker.
#' @param tau_fixed Numeric fixed value, either scalar or one value per marker.
#' @param family_codes Integer family code for each marker.
#'
#' @return A list with integer `use_tau_fixed` and numeric `tau_fixed`, each
#'   having one element per marker.
#' @keywords internal
#' @noRd
.normalise_fixed_tau_by_marker <- function(
  use_tau_fixed,
  tau_fixed,
  family_codes
) {
  family_codes <- as.integer(family_codes)
  n_markers <- length(family_codes)
  if (n_markers == 0L) {
    return(list(use_tau_fixed = integer(0), tau_fixed = numeric(0)))
  }

  use_tau_fixed <- as.integer(use_tau_fixed %||% 0L)
  tau_fixed <- as.numeric(tau_fixed %||% 0.5)
  if (length(use_tau_fixed) == 1L) {
    use_tau_fixed <- rep.int(use_tau_fixed, n_markers)
  }
  if (length(tau_fixed) == 1L) {
    tau_fixed <- rep.int(tau_fixed, n_markers)
  }
  if (length(use_tau_fixed) != n_markers || length(tau_fixed) != n_markers) {
    cli::cli_abort(c(
      x = "Fixed skew-Laplace metadata is not aligned with the fitted markers.",
      i = "Refit the model or supply one fixed-tau flag and value per marker."
    ))
  }
  if (anyNA(use_tau_fixed) || any(!use_tau_fixed %in% c(0L, 1L))) {
    cli::cli_abort("Fixed skew-Laplace flags must be zero or one.")
  }

  not_skew_laplace <- is.na(family_codes) | family_codes != 9L
  use_tau_fixed[not_skew_laplace] <- 0L
  fixed_index <- which(use_tau_fixed == 1L)
  if (
    length(fixed_index) > 0L &&
      (
        any(!is.finite(tau_fixed[fixed_index])) ||
          any(tau_fixed[fixed_index] <= 0) ||
          any(tau_fixed[fixed_index] >= 1)
      )
  ) {
    cli::cli_abort(
      "Every fixed skew-Laplace tau must lie strictly between zero and one."
    )
  }

  # Stan requires a valid bounded value even where the flag is zero. The
  # symmetric value is a neutral placeholder and is never read by the
  # likelihood for an estimated marker.
  tau_fixed[use_tau_fixed == 0L] <- 0.5
  list(
    use_tau_fixed = as.integer(use_tau_fixed),
    tau_fixed = as.numeric(tau_fixed)
  )
}

#' Identify a scalar family specification
#'
#' @description
#' Distinguishes one structured family declaration from a list containing
#' several marker-specific declarations. This distinction is required because
#' a `JoiNMe_family_spec` is itself represented as a list.
#'
#' @param x An object supplied through a `families` argument.
#'
#' @return `TRUE` when `x` is one formal or compatible list-style family
#'   specification; otherwise `FALSE`.
#' @keywords internal
#' @noRd
.is_single_family_spec <- function(x) {
  inherits(x, "JoiNMe_family_spec") ||
    (is.list(x) && !is.null(x$family))
}

#' Normalise a named forward link
#'
#' @param link Character scalar containing a supported forward-link name or
#'   recognised alias.
#'
#' @return One canonical name among `"identity"`, `"log"`, `"logit"`,
#'   `"probit"`, and `"exp"`.
#' @keywords internal
#' @noRd
.normalize_link_name <- function(link) {
  if (!is.character(link) || length(link) != 1L || !nzchar(link)) {
    cli::cli_abort(c(
      x = "{.arg link} must be a non-empty character scalar.",
      i = "Allowed values: identity, log, logit, probit, exp."
    ))
  }

  key <- tolower(trimws(link))
  key <- switch(
    key,
    id = "identity",
    inverse = "identity",
    inverse_logit = "logit",
    inv_logit = "logit",
    sigmoid = "logit",
    expit = "logit",
    key
  )

  allowed <- c("identity", "log", "logit", "probit", "exp")
  if (!key %in% allowed) {
    cli::cli_abort(c(
      x = "Unsupported link: {.val {link}}.",
      i = "Allowed values: identity, log, logit, probit, exp."
    ))
  }
  key
}

#' Encode a canonical forward link for Stan data
#'
#' @param link_name Character forward-link name accepted by
#'   `.normalize_link_name()`.
#'
#' @return Integer code in one through five.
#' @keywords internal
#' @noRd
.link_code_from_name <- function(link_name) {
  nm <- .normalize_link_name(link_name)
  switch(
    nm,
    identity = 1L,
    log = 2L,
    logit = 3L,
    probit = 4L,
    exp = 5L
  )
}

#' Compile the inverse of a named forward link
#'
#' @description
#' Maps a canonical forward link \eqn{g} to bytecode evaluating
#' \eqn{g^{-1}}.  In particular, a log forward link maps to exponential
#' bytecode and an exponential forward link maps to logarithm bytecode.
#'
#' @param link_name Character forward-link name accepted by
#'   `.normalize_link_name()`.
#'
#' @return A functional-bytecode list with operation and constant counts.
#' @keywords internal
#' @noRd
.inv_link_bc_from_name <- function(link_name) {
  nm <- .normalize_link_name(link_name)
  normal_ops <- .bytecode_normal_ops()
  bc <- switch(
    nm,
    identity = list(bytecode = c(0L), const_data = numeric(0)),
    log = list(bytecode = c(0L, 7L), const_data = numeric(0)),
    logit = list(bytecode = c(0L, 9L), const_data = numeric(0)),
    probit = list(
      bytecode = c(0L, unname(normal_ops[["PHI"]])),
      const_data = numeric(0)
    ),
    exp = list(bytecode = c(0L, 6L), const_data = numeric(0))
  )
  bc$n_ops <- length(bc$bytecode)
  bc$n_const <- length(bc$const_data)
  bc
}

#' Compile the inverse of a named or formula link
#'
#' @description
#' Character links are translated through the canonical link table.  A formula
#' is treated as a forward link \eqn{g(x)}, inverted symbolically, and compiled
#' to functional bytecode.
#'
#' @param link A supported character link name or an invertible one-sided
#'   formula in `x`.
#'
#' @return A validated functional-bytecode specification for the inverse link.
#' @keywords internal
#' @noRd
.inv_link_bc_from_link <- function(link) {
  if (is.character(link)) {
    return(.inv_link_bc_from_name(link))
  }
  inverse_expr <- invert_transform_expr(link)
  parse_transform_expr(inverse_expr)
}

#' Recognise the forward link represented by inverse-link bytecode
#'
#' @param inv_link_bc A functional-bytecode specification that evaluates an
#'   inverse link.
#'
#' @return A canonical *forward* link name, or `NA_character_` for a custom
#'   inverse link.
#' @keywords internal
#' @noRd
.canonical_link_from_inv_link_bc <- function(inv_link_bc) {
  ops <- as.integer(inv_link_bc$bytecode %||% integer(0))
  const_data <- as.numeric(inv_link_bc$const_data %||% numeric(0))
  normal_ops <- .bytecode_normal_ops()
  if (length(const_data) > 0L) {
    return(NA_character_)
  }
  if (identical(ops, c(0L))) {
    return("identity")
  }
  if (identical(ops, c(0L, 7L))) {
    return("log")
  }
  if (identical(ops, c(0L, 9L))) {
    return("logit")
  }
  if (identical(ops, c(0L, unname(normal_ops[["PHI"]])))) {
    return("probit")
  }
  if (identical(ops, c(0L, 6L))) {
    return("exp")
  }
  NA_character_
}

#' Detect operations that are commonly non-injective
#'
#' @param inv_link_bc A functional-bytecode inverse-link specification.
#'
#' @return `TRUE` when the program contains a periodic, absolute-value, square,
#'   or hyperbolic-cosine operation; otherwise `FALSE`.
#' @keywords internal
#' @noRd
.bytecode_has_potentially_noninvertible_ops <- function(inv_link_bc) {
  # Heuristic operators that are non-injective on their natural/global domains.
  risky_ops <- c(
    13L, # sin
    14L, # cos
    15L, # tan
    16L, # abs
    17L, # square
    # 18L, # sinh (strictly monotone but can induce extreme tails)
    19L # cosh
    # 20L  # tanh (bounded; may flatten heavily)
  )
  ops <- as.integer(inv_link_bc$bytecode %||% integer(0))
  any(ops %in% risky_ops)
}

#' Check numerical monotonicity of inverse-link bytecode
#'
#' @param inv_link_bc A functional-bytecode inverse-link specification.
#' @param grid Numeric reference grid on the linear-predictor scale.
#' @param tol Non-negative numerical tolerance for successive differences.
#'
#' @return `TRUE` when all finite evaluated differences have one direction;
#'   otherwise `FALSE`.
#' @keywords internal
#' @noRd
.validate_monotonic_bytecode <- function(
  inv_link_bc,
  grid = seq(-8, 8, length.out = 257L),
  tol = 1e-10
) {
  vals <- tryCatch(
    eval_bytecode_vector(
      grid,
      bytecode = inv_link_bc$bytecode,
      const_data = inv_link_bc$const_data %||% numeric(0)
    ),
    error = function(e) NULL
  )

  if (is.null(vals) || length(vals) != length(grid) || any(!is.finite(vals))) {
    return(FALSE)
  }

  diffs <- diff(as.numeric(vals))
  all(diffs >= -tol) || all(diffs <= tol)
}

#' Warn about a potentially non-injective transformation
#'
#' @param bc A functional-bytecode specification.
#' @param context Character description of the family declaration used in the
#'   warning.
#' @param mode Character scalar `"warning"` or `"error"`.
#'
#' @return `NULL`, invisibly.  A warning is issued when canonical recognition
#'   and numerical monotonicity checks do not support a one-to-one map.
#' @keywords internal
#' @noRd
.validate_noninvertible <- function(bc, context, mode = c('warning', 'error')) {
  # Issue an explicit identifiability warning when inverse-link program is not canonical
  # or uses operations that are frequently non-invertible in practical model domains.
  canonical_name <- .canonical_link_from_inv_link_bc(bc)
  has_risky_ops <- .bytecode_has_potentially_noninvertible_ops(bc)
  is_monotone <- .validate_monotonic_bytecode(bc)
  if ((!is.na(canonical_name) || is_monotone) && !has_risky_ops) {
    return(invisible(NULL))
  }
  mode <- match.arg(mode)
  mode_fn <- switch(
    mode,
    warning = cli::cli_warn,
    error = cli::cli_abort
  )

  mode_fn(c(
    x = "Bytecode in {context} may be non-invertible.",
    i = "This can make the joint model weakly identified or unidentifiable (e.g., periodic transforms like sin/cos).",
    i = "Prefer monotone one-to-one inverse links when possible (identity, exp, inv_logit, Phi, etc.)."
  ))
  invisible(NULL)
}

#' Decode a forward-link integer code
#'
#' @param link_code Integer code stored in Stan data.
#'
#' @return Canonical forward-link name, `"custom"` for zero, or `"identity"`
#'   as a defensive fallback.
#' @keywords internal
#' @noRd
.link_name_from_code <- function(link_code) {
  switch(
    as.character(as.integer(link_code)),
    "0" = "custom",
    "1" = "identity",
    "2" = "log",
    "3" = "logit",
    "4" = "probit",
    "5" = "exp",
    "identity"
  )
}

#' Apply one fitted inverse link to a matrix of linear predictors
#'
#' @description
#' Transforms posterior linear predictors to the response-parameter scale
#' using the same canonical link codes and custom bytecode used by Stan.
#' Centralising this calculation ensures that plotting and conditional-effects
#' summaries cannot assign different meanings to the same family declaration.
#'
#' Canonical forward-link codes are inverted as follows: identity returns
#' \eqn{\eta}; log returns \eqn{\exp(\eta)}; logit returns the logistic CDF;
#' probit returns the standard normal CDF \eqn{\Phi(\eta)}; and exponential
#' returns \eqn{\log(\eta)}. Code zero denotes a custom inverse-link program.
#'
#' @param eta Numeric matrix with posterior draws in rows and evaluation
#'   points in columns.
#' @param link_code Integer canonical forward-link code. Zero requests custom
#'   bytecode; codes one through five denote identity, log, logit, probit, and
#'   exponential forward links.
#' @param bytecode Integer custom inverse-link instruction sequence, used only
#'   when `link_code` is zero or otherwise non-canonical.
#' @param const_data Numeric constants consumed by custom bytecode.
#'
#' @return Numeric matrix with the same dimensions and dimension names as
#'   `eta`, evaluated on the response-parameter scale.
#' @keywords internal
#' @noRd
.apply_inverse_link_matrix <- function(
  eta,
  link_code,
  bytecode = integer(0),
  const_data = numeric(0)
) {
  eta <- as.matrix(eta)
  link_code <- as.integer(link_code)[1L]
  if (!is.finite(link_code)) {
    cli::cli_abort("The fitted inverse-link code must be a finite integer.")
  }

  transformed <- switch(
    as.character(link_code),
    "1" = eta,
    "2" = exp(eta),
    "3" = stats::plogis(eta),
    "4" = stats::pnorm(eta),
    "5" = log(eta),
    NULL
  )
  if (!is.null(transformed)) {
    dim(transformed) <- dim(eta)
    dimnames(transformed) <- dimnames(eta)
    return(transformed)
  }

  transformed <- matrix(
    NA_real_,
    nrow = nrow(eta),
    ncol = ncol(eta),
    dimnames = dimnames(eta)
  )
  for (column in seq_len(ncol(eta))) {
    transformed[, column] <- eval_bytecode_vector(
      eta[, column],
      bytecode = bytecode,
      const_data = const_data
    )
  }
  transformed
}

#' Select the default forward link for a response family
#'
#' @param family_code Integer longitudinal-family code.
#'
#' @return Canonical forward-link name.  Binary, binomial, beta, and ordinal
#'   families use logit; count families use log; remaining families use
#'   identity.
#' @keywords internal
#' @noRd
.default_link_for_family <- function(family_code) {
  switch(
    as.character(as.integer(family_code)),
    "3" = "logit", # bernoulli
    "4" = "logit", # binomial
    "5" = "log", # poisson
    "6" = "log", # negbin2
    "10" = "logit", # beta mean
    "11" = "logit", # cumulative logit latent scale
    "identity"
  )
}

#' Recognise a canonical link from an inverse-link expression
#'
#' @param inv_link Formula or expression directly defining an inverse link.
#'
#' @return Canonical forward-link name or `NA_character_`.
#' @keywords internal
#' @noRd
.link_name_from_inv_link_expr <- function(inv_link) {
  bc <- parse_transform_expr(inv_link)
  .canonical_link_from_inv_link_bc(bc)
}

#' Resolve one family and link declaration
#'
#' @description
#' Accepts a `JoiNMe_family_spec`, a list-style declaration, or a bare family
#' name and returns the common representation used by standata construction.
#'
#' @param x One marker-specific family declaration.
#'
#' @return A list containing integer `family_code`, canonical or missing
#'   `link_name`, compiled `inv_link_bc`, and `tau_fixed`, where the latter is
#'   `NA_real_` unless a fixed skew-Laplace quantile was requested.
#' @keywords internal
#' @noRd
.extract_family_and_link <- function(x) {
  if (inherits(x, "JoiNMe_family_spec")) {
    fam_code <- .parse_family(x$family)
    inv_link_bc <- x$inv_link
    if (is.null(inv_link_bc)) {
      inv_link_bc <- .inv_link_bc_from_name(
        x$link %||% .default_link_for_family(fam_code)
      )
    }
    link_name <- .canonical_link_from_inv_link_bc(inv_link_bc)
    tau_fixed <- .validate_family_tau(
      tau = x$tau,
      family_code = fam_code,
      context = "JoiNMe family specification"
    )
    return(list(
      family_code = as.integer(fam_code),
      link_name = link_name,
      inv_link_bc = inv_link_bc,
      tau_fixed = tau_fixed
    ))
  }

  if (is.list(x) && !is.null(x$family)) {
    fam_code <- .parse_family(x$family)
    if (!is.null(x$link) && !is.null(x$inv_link)) {
      cli::cli_abort(c(
        x = "Family specification cannot contain both {.arg link} and {.arg inv_link}.",
        i = "Keep only one of them."
      ))
    }
    inv_link_bc <- if (!is.null(x$inv_link)) {
      parse_transform_expr(x$inv_link)
    } else if (!is.null(x$link)) {
      .inv_link_bc_from_link(x$link)
    } else {
      .inv_link_bc_from_name(.default_link_for_family(fam_code))
    }
    if (!is.null(x$inv_link)) {
      .validate_noninvertible(
        inv_link_bc,
        context = paste0("families[[", x$family, "]]")
      )
    }
    link_name <- .canonical_link_from_inv_link_bc(inv_link_bc)
    tau_fixed <- .validate_family_tau(
      tau = x$tau,
      family_code = fam_code,
      context = "list-style family specification"
    )
    return(list(
      family_code = as.integer(fam_code),
      link_name = link_name,
      inv_link_bc = inv_link_bc,
      tau_fixed = tau_fixed
    ))
  }

  fam_code <- .parse_family(x)
  inv_link_bc <- .inv_link_bc_from_name(.default_link_for_family(fam_code))
  list(
    family_code = as.integer(fam_code),
    link_name = .canonical_link_from_inv_link_bc(inv_link_bc),
    inv_link_bc = inv_link_bc,
    tau_fixed = NA_real_
  )
}

#' Parse family string to numeric code
#' @param family String family name (e.g., "gaussian", "student_t", "binomial") or vector of names
#' @return Integer family code (1-11) or vector of codes
#' @keywords internal
#' @noRd
.parse_family <- function(family) {
  # Map human-readable family names to Stan integer codes
  family <- tolower(as.character(family))
  # Handle vector input
  if (length(family) > 1) {
    return(sapply(family, .parse_family))
  }

  switch(
    family,
    "gaussian" = 1L,
    "normal" = 1L,
    "student_t" = 2L,
    "student" = 2L,
    "student-t" = 2L,
    "bernoulli" = 3L,
    "binomial" = 4L,
    "poisson" = 5L,
    "negbin2" = 6L,
    "negative_binomial" = 6L,
    "skew_normal" = 7L,
    "skew-normal" = 7L,
    "skewnormal" = 7L,
    "double_exponential" = 8L,
    "double-exponential" = 8L,
    "laplace" = 8L,
    "skew_double_exponential" = 9L,
    "skew-double-exponential" = 9L,
    "skew_laplace" = 9L,
    "skew-laplace" = 9L,
    "beta" = 10L,
    "cumulative_logit" = 11L,
    "ordered_logistic" = 11L,
    "ordinal" = 11L,
    stop(
      "Unknown family: ",
      family,
      ". Must be one of: gaussian, student_t, bernoulli, binomial, poisson, negbin2, ",
      "skew_normal, double_exponential, skew_double_exponential, beta, cumulative_logit, ordinal"
    )
  )
}

#' Convert family code to string
#' @param fam Integer family code (1-11), can be vector
#' @return Character family name(s)
#' @keywords internal
#' @noRd
.family_code_to_name <- function(fam) {
  # Map Stan integer codes to canonical family names
  fam <- as.integer(fam)
  # Handle vector input
  if (length(fam) > 1) {
    return(sapply(fam, .family_code_to_name))
  }

  switch(
    as.character(fam),
    "1" = "gaussian",
    "2" = "student_t",
    "3" = "bernoulli",
    "4" = "binomial",
    "5" = "poisson",
    "6" = "negbin2",
    "7" = "skew_normal",
    "8" = "double_exponential",
    "9" = "skew_double_exponential",
    "10" = "beta",
    "11" = "cumulative_logit",
    stop("Unknown family code: ", fam)
  )
}

#' Check family is appropriate for outcome type
#' @param family Integer family code
#' @param y_range Numeric range of observed values
#' @keywords internal
#' @noRd
.validate_family_outcome <- function(family, y_range, name = "y") {
  # Validate observed outcome ranges by family requirements
  fam_name <- .family_code_to_name(family)
  if (family %in% c(1, 2, 7, 8, 9)) {
    # Gaussian, student_t, skew normal, double exponential, skew double exponential: any real
    return(TRUE)
  } else if (family == 10) {
    # Beta: must be strictly between 0 and 1
    if (any(y_range <= 0) || any(y_range >= 1)) {
      stop(
        "Beta family requires ",
        name,
        " in (0,1), but got values in [",
        min(y_range),
        ", ",
        max(y_range),
        "]"
      )
    }
  } else if (family == 11) {
    # Cumulative logit (ordinal): must be integer categories >= 1
    if (any(y_range < 1) || any(y_range != as.integer(y_range))) {
      stop(
        "Cumulative logit family requires ",
        name,
        " to be integer categories >= 1, but got values in [",
        min(y_range),
        ", ",
        max(y_range),
        "]"
      )
    }
  } else if (family == 3) {
    # Bernoulli: must be 0 or 1
    if (!all(y_range %in% c(0, 1))) {
      stop(
        "Bernoulli family requires ",
        name,
        " in {0,1}, but got values in [",
        min(y_range),
        ", ",
        max(y_range),
        "]"
      )
    }
  } else if (family == 4) {
    # Binomial: integers in [0,n]
    if (any(y_range < 0) || any(y_range != as.integer(y_range))) {
      stop(
        "Binomial family requires ",
        name,
        " to be non-negative integers, but got values in [",
        min(y_range),
        ", ",
        max(y_range),
        "]"
      )
    }
  } else if (family %in% c(5, 6)) {
    # Poisson, negbin2: non-negative integers
    if (any(y_range < 0) || any(y_range != as.integer(y_range))) {
      stop(
        fam_name,
        " family requires ",
        name,
        " to be non-negative integers, but got values in [",
        min(y_range),
        ", ",
        max(y_range),
        "]"
      )
    }
  }

  return(TRUE)
}

#' Validate family list for markers
#' @description
#' Check that family list has right structure and all families are valid.
#'
#' @param families List or character vector specifying family per marker
#' @param D Number of markers
#' @param dataLong Longitudinal data frame
#' @param marker_var Name of marker column
#' @param y_var Name of outcome column
#'
#' @return A list containing marker-aligned family and inverse-link metadata,
#'   together with `use_tau_fixed` and `tau_fixed` vectors. The fixed-`tau`
#'   flag is one only for skew-Laplace markers whose family declaration
#'   supplies a constant quantile.
#' @keywords internal
#' @noRd
.validate_family_list <- function(families, D, dataLong, marker_var, y_var) {
  # Validate per-marker families against observed data
  # Convert a scalar character or structured declaration to a marker list.
  if (.is_single_family_spec(families)) {
    families <- rep(list(families), D)
  } else if (length(families) == 1 && !is.list(families)) {
    families <- rep(list(families), D)
  } else if (!is.list(families)) {
    families <- as.list(families)
  }

  if (length(families) != D) {
    stop(
      "families must have length ",
      D,
      " (number of markers), got ",
      length(families)
    )
  }

  # Parse each family/link pair with strict, deterministic defaults.
  family_codes <- integer(D)
  link_names <- rep(NA_character_, D)
  link_codes <- integer(D)
  inv_link_specs <- vector("list", D)
  tau_fixed <- rep.int(0.5, D)
  use_tau_fixed <- integer(D)
  for (d in seq_along(families)) {
    spec_d <- .extract_family_and_link(families[[d]])
    family_codes[d] <- as.integer(spec_d$family_code)
    inv_link_specs[[d]] <- spec_d$inv_link_bc
    link_names[d] <- spec_d$link_name
    link_codes[d] <- if (is.na(spec_d$link_name)) {
      0L
    } else {
      .link_code_from_name(spec_d$link_name)
    }
    if (!is.na(spec_d$tau_fixed)) {
      use_tau_fixed[d] <- 1L
      tau_fixed[d] <- as.numeric(spec_d$tau_fixed)
    }
  }

  inv_link_n_ops <- as.integer(vapply(
    inv_link_specs,
    function(x) length(x$bytecode %||% integer(0)),
    integer(1)
  ))
  inv_link_n_const <- as.integer(vapply(
    inv_link_specs,
    function(x) length(x$const_data %||% numeric(0)),
    integer(1)
  ))
  max_inv_link_ops <- as.integer(max(inv_link_n_ops, 1L))
  max_inv_link_const <- as.integer(max(inv_link_n_const, 1L))
  inv_link_ops <- matrix(0L, nrow = D, ncol = max_inv_link_ops)
  inv_link_const <- matrix(0.0, nrow = D, ncol = max_inv_link_const)
  for (d in seq_len(D)) {
    ops_d <- as.integer(inv_link_specs[[d]]$bytecode %||% integer(0))
    const_d <- as.numeric(inv_link_specs[[d]]$const_data %||% numeric(0))
    if (length(ops_d) > 0L) {
      inv_link_ops[d, seq_along(ops_d)] <- ops_d
    }
    if (length(const_d) > 0L) inv_link_const[d, seq_along(const_d)] <- const_d
  }

  # Validate against actual data
  # Make sure marker_var is a factor
  if (!is.factor(dataLong[[marker_var]])) {
    dataLong[[marker_var]] <- factor(dataLong[[marker_var]])
  }

  marker_levels <- levels(dataLong[[marker_var]])
  if (length(marker_levels) != D) {
    stop(
      "Number of unique markers (",
      length(marker_levels),
      ") does not match D (",
      D,
      ")"
    )
  }

  for (d in seq_len(D)) {
    marker_name <- marker_levels[d]
    y_vals <- dataLong[[y_var]][dataLong[[marker_var]] == marker_name]
    if (length(y_vals) > 0) {
      .validate_family_outcome(
        family_codes[d],
        y_vals,
        name = paste0(y_var, "[marker=", marker_name, "]")
      )
    }
  }

  list(
    family_codes = as.integer(family_codes),
    link_codes = as.integer(link_codes),
    link_names = as.character(link_names),
    inv_link_n_ops = inv_link_n_ops,
    inv_link_ops = inv_link_ops,
    inv_link_n_const = inv_link_n_const,
    inv_link_const = inv_link_const,
    use_tau_fixed = as.integer(use_tau_fixed),
    tau_fixed = as.numeric(tau_fixed)
  )
}

#' Get required distributional parameters for a family
#' @param family Integer family code
#' @return Character vector of parameter names
#' @keywords internal
#' @noRd
.family_distrib_params <- function(family) {
  # Enumerate distributional regression parameters per family
  switch(
    as.integer(family),
    "1" = c("sigma"), # gaussian
    "2" = c("sigma", "nu"), # student_t
    "3" = c(), # bernoulli
    "4" = c(), # binomial (uses trials)
    "5" = c(), # poisson
    "6" = c("phi"), # negbin2
    "7" = c("sigma", "alpha"), # skew_normal (scale + skew)
    "8" = c("sigma"), # double_exponential (scale)
    "9" = c("sigma", "tau"), # skew_double_exponential (scale + skew)
    "10" = c("kappa"), # beta (precision)
    "11" = c(), # cumulative_logit (cutpoints)
    stop("Unknown family code: ", family)
  )
}

#' Build marker-to-family indexing for a distributional parameter
#'
#' @description
#' Identifies the distinct response families that require one distributional
#' parameter and maps each eligible marker to the corresponding family-level
#' parameter. A marker receives index zero when its family does not use the
#' parameter or when a marker-specific fixed value replaces estimation.
#'
#' @param family_codes Integer family code for every longitudinal marker.
#' @param parameter Character scalar naming a supported distributional
#'   parameter.
#' @param eligible Optional logical vector aligned with `family_codes`. `FALSE`
#'   excludes a marker from family-level estimation even when its family
#'   ordinarily requires the parameter.
#'
#' @return A list containing the number of estimated family-level parameters,
#'   the marker-to-parameter index, and the corresponding family codes and
#'   names.
#' @keywords internal
#' @noRd
.build_family_parameter_index <- function(family_codes, parameter, eligible = NULL) {
  family_codes <- as.integer(family_codes)
  n_markers <- length(family_codes)
  if (is.null(eligible)) {
    eligible <- rep.int(TRUE, n_markers)
  }
  if (!is.logical(eligible) || length(eligible) != n_markers || anyNA(eligible)) {
    cli::cli_abort(c(
      x = "{.arg eligible} must contain one non-missing logical value per marker.",
      i = "Align eligibility with the marker-specific family vector."
    ))
  }

  requires_parameter <- vapply(
    family_codes,
    function(family_code) parameter %in% .family_distrib_params(family_code),
    logical(1)
  )
  estimate_parameter <- requires_parameter & eligible
  parameter_families <- sort(unique(family_codes[estimate_parameter]))

  marker_to_family <- integer(n_markers)
  if (length(parameter_families) > 0L) {
    marker_to_family[estimate_parameter] <- match(
      family_codes[estimate_parameter],
      parameter_families
    )
  }

  list(
    n = as.integer(length(parameter_families)),
    marker_to = as.integer(marker_to_family),
    family_codes = as.integer(parameter_families),
    family_names = if (length(parameter_families) > 0L) {
      vapply(parameter_families, .family_code_to_name, character(1))
    } else {
      character(0)
    }
  )
}

#' Generate prior predictive for outcome given family
#' @description
#' Simple prior predictive sampling for validation/simulation
#'
#' @param n Sample size
#' @param mu Mean/probability on the response scale (i.e., inverse-link applied).
#' @param family Integer family code
#' @param sigma Standard deviation (if needed)
#' @param nu Degrees of freedom (if needed, student_t)
#' @param phi Dispersion (if needed, negbin2)
#' @param kappa Positive sample-size parameter for the Beta family. Under the
#'   mean/sample-size parameterisation, the two Beta shape parameters are
#'   \eqn{a = \mu\kappa} and \eqn{b = (1 - \mu)\kappa}.
#' @param tau Quantile/asymmetry parameter in \eqn{(0,1)} for the skew double
#'   exponential family. The value \eqn{\tau = 0.5} gives the symmetric case.
#' @param trials Number of trials (if needed, binomial)
#'
#' @return Numeric vector of samples
#' @keywords internal
#' @noRd
.sample_from_family <- function(
  n,
  mu,
  family,
  sigma = 1,
  nu = 4,
  phi = 1,
  kappa = 10,
  tau = 0.5,
  trials = 1,
  skew = 0
) {
  # Prior predictive helper for a single family
  family <- as.integer(family)
  if (is.null(mu)) {
    cli::cli_abort(
      "{.arg mu} must be provided as the inverse-link mean/probability."
    )
  }

  if (family == 1) {
    # Gaussian
    rnorm(n, mu, sigma)
  } else if (family == 2) {
    # Student-t
    rt(n, df = nu) * sigma + mu
  } else if (family == 3) {
    # Bernoulli
    rbinom(n, size = 1, prob = mu)
  } else if (family == 4) {
    # Binomial
    rbinom(n, size = trials, prob = mu)
  } else if (family == 5) {
    # Poisson
    rpois(n, lambda = mu)
  } else if (family == 6) {
    # Negative binomial 2
    rnbinom(n, size = phi, mu = pmax(mu, 1e-10))
  } else if (family == 7) {
    # Skew-normal
    brms::rskew_normal(n, mu = mu, sigma = sigma, alpha = skew)
  } else if (family == 8) {
    # Double exponential (Laplace)
    # Sample by drawing an exponential magnitude and random sign.
    sign <- sample(c(-1, 1), size = n, replace = TRUE)
    mu + sign * rexp(n, rate = 1 / sigma)
  } else if (family == 9) {
    # Skew double exponential (asymmetric Laplace)
    # Piecewise exponential: probability mass split by tau.
    u <- runif(n)
    e <- rexp(n, rate = 1)
    mu + ifelse(u < tau, e * sigma / tau, -e * sigma / (1 - tau))
  } else if (family == 10) {
    # Step 1: express both Beta shapes through the response-scale mean and the
    # positive sample-size parameter kappa.
    shape1 <- mu * kappa
    shape2 <- (1 - mu) * kappa

    # Step 2: sample from the resulting mean/sample-size Beta distribution.
    rbeta(n, shape1 = shape1, shape2 = shape2)
  } else if (family == 11) {
    # Cumulative logit: return category labels (1..K)
    # This is a placeholder sampler; use with a valid cutpoint structure externally.
    stop(
      "Sampling for cumulative logit requires cutpoints; use model-based sampling instead."
    )
  } else {
    stop("Unknown family: ", family)
  }
}

#' Encode one validated prior declaration for Stan data
#'
#' @param prior A normalised `joinme_prior_spec` object.
#'
#' @return A list containing the integer family code and fixed
#'   hyperparameters used by the reusable Stan prior module.
#' @keywords internal
#' @noRd
.encode_joinme_prior <- function(prior) {
  family_codes <- c(
    student_t = 1L,
    normal = 2L,
    laplace = 3L,
    horseshoe = 4L
  ) # stable R-to-Stan family-code dictionary
  list(
    family = unname(family_codes[[prior$family]]),
    mu = as.numeric(prior$mu),
    scale = as.numeric(prior$scale),
    df = if (is.finite(prior$df)) prior$df else 1,
    global_df = if (is.finite(prior$global_df)) prior$global_df else 1,
    global_scale = if (is.finite(prior$global_scale)) prior$global_scale else 1,
    slab_df = if (is.finite(prior$slab_df)) prior$slab_df else 4,
    slab_scale = if (is.finite(prior$slab_scale)) prior$slab_scale else 2
  )
}

#' Coefficient-prior declarations at fitted block dimensions
#'
#' @description
#' Scalar locations and scales are recycled only after formula parsing reveals
#' the exact number and order of coefficients. Non-scalar declarations must
#' match that number exactly, preventing silent partial recycling across
#' scientifically different parameters.
#'
#' @param priors Named list of encoded prior declarations, such as the
#'   family-only marker blocks or the class-membership regression block.
#' @param dimensions Named integer vector giving the fitted dimension of every
#'   prior block to assemble.
#'
#' @return Named Stan-data fields for every requested prior block.
#' @keywords internal
#' @noRd
.assemble_joinme_prior_data <- function(priors, dimensions) {
  expand_values <- function(values, number_parameters, component, field) {
    values <- as.numeric(values)
    if (number_parameters == 0L) return(numeric(0))
    if (length(values) == 1L) return(rep.int(values, number_parameters))
    if (length(values) != number_parameters) {
      cli::cli_abort(c(
        x = "Prior {.arg {component}} provides {length(values)} {field} values for {number_parameters} parameters.",
        i = "Supply one value or exactly one value per parameter in the fitted block."
      ))
    }
    values
  }

  output <- list()
  for (component in names(dimensions)) {
    specification <- priors[[component]]
    number_parameters <- as.integer(dimensions[[component]])
    prefix <- paste0("prior_", component, "_")
    output[[paste0(prefix, "family")]] <- as.integer(specification$family)
    output[[paste0(prefix, "mu")]] <- expand_values(
      specification$mu,
      number_parameters,
      component,
      "location"
    )
    output[[paste0(prefix, "scale")]] <- expand_values(
      specification$scale,
      number_parameters,
      component,
      "scale"
    )
    output[[paste0(prefix, "df")]] <- as.numeric(specification$df)
    output[[paste0(prefix, "global_df")]] <- as.numeric(specification$global_df)
    output[[paste0(prefix, "global_scale")]] <- as.numeric(specification$global_scale)
    output[[paste0(prefix, "slab_df")]] <- as.numeric(specification$slab_df)
    output[[paste0(prefix, "slab_scale")]] <- as.numeric(specification$slab_scale)
  }
  output
}

#' Resolve scoped priors for one distributional regression
#'
#' @description
#' Converts parameter-wide, family-scoped and marker-scoped declarations into
#' disjoint coefficient assignments. Family scope is read from the prefixes
#' created by `.build_dist_matrix()`. A marker selector is permitted only when
#' that marker is the sole marker using its family block; otherwise the fitted
#' coefficient is shared and a marker-specific prior would misstate the model.
#'
#' @param parameter Distributional parameter name.
#' @param columns Distributional model-matrix column labels.
#' @param roles Intercept or slope role for every column.
#' @param specification Normalised prior specification for `parameter`.
#' @param marker_levels Longitudinal marker labels in fitted order.
#' @param family_codes Response-family code for every marker.
#'
#' @return A block accepted by `.pack_regression_priors()`.
#' @keywords internal
#' @noRd
.distributional_prior_block <- function(
  parameter,
  columns,
  roles,
  specification,
  marker_levels,
  family_codes
) {
  columns <- as.character(columns) # fitted distributional coefficient labels
  roles <- as.character(roles) # intercept or slope role aligned with columns
  marker_levels <- as.character(marker_levels) # exact response-marker labels
  family_codes <- as.integer(family_codes) # marker-aligned response-family codes
  if (length(columns) != length(roles)) {
    cli::cli_abort("Internal distributional prior columns and roles have different lengths for {.field {parameter}}.")
  }
  if (length(columns) == 0L) {
    if (length(specification$by_family) > 0L || length(specification$by_marker) > 0L) {
      cli::cli_abort(c(
        x = "Scoped prior declarations were supplied for {.field {parameter}}, but its {.arg formulaDist} regression has no coefficients.",
        i = "Add the corresponding family-scoped distributional formula or remove the bracket selector."
      ))
    }
    return(list(roles = roles, priors = specification))
  }

  coefficient_source_type <- rep.int("default", length(columns)) # default, family, or marker selector governing each coefficient
  coefficient_source_value <- rep.int(NA_character_, length(columns)) # canonical family or exact marker label for scoped coefficients
  family_prefix <- ifelse(
    grepl("^family=[^:]+::", columns),
    sub("^family=([^:]+)::.*$", "\\1", columns),
    NA_character_
  ) # canonical family scope encoded in each coefficient label

  for (family_name in names(specification$by_family)) {
    positions <- which(family_prefix == family_name) # coefficients fitted only for this response family
    if (length(positions) == 0L) {
      cli::cli_abort(c(
        x = "Prior selector {.code {parameter}[family='{family_name}']} does not match a fitted distributional coefficient block.",
        i = "Use the same family scope on the left-hand side of {.arg formulaDist}."
      ))
    }
    coefficient_source_type[positions] <- "family"
    coefficient_source_value[positions] <- family_name
  }

  family_names_by_marker <- vapply(family_codes, .family_code_to_name, character(1)) # canonical family for every marker
  for (marker_name in names(specification$by_marker)) {
    marker_position <- match(marker_name, marker_levels) # selected response marker in fitted marker order
    if (is.na(marker_position)) {
      cli::cli_abort(c(
        x = "Unknown marker in prior selector {.code {parameter}[marker='{marker_name}']}.",
        i = "Available markers are: {.val {paste(marker_levels, collapse = ', ')}}."
      ))
    }
    family_name <- family_names_by_marker[[marker_position]] # family block used by the selected marker
    markers_sharing_family <- marker_levels[family_names_by_marker == family_name] # markers governed by the same coefficients
    if (length(markers_sharing_family) != 1L) {
      cli::cli_abort(c(
        x = "Prior selector {.code {parameter}[marker='{marker_name}']} cannot isolate one fitted coefficient block.",
        i = "Markers {.val {paste(markers_sharing_family, collapse = ', ')}} share the {.val {family_name}} distributional coefficients.",
        i = "Use {.code {parameter}[family='{family_name}']} or define a distributional design with distinct coefficients."
      ))
    }
    positions <- which(family_prefix == family_name) # uniquely marker-associated family-scoped coefficients
    if (length(positions) == 0L) {
      cli::cli_abort(c(
        x = "Prior selector {.code {parameter}[marker='{marker_name}']} does not match a fitted distributional coefficient block.",
        i = "The marker prior can be resolved only when {.arg formulaDist} contains {.code {parameter}[family={family_name}] ~ ...}."
      ))
    }
    coefficient_source_type[positions] <- "marker"
    coefficient_source_value[positions] <- marker_name
  }

  assignment_table <- unique(data.frame(
    source_type = coefficient_source_type,
    source_value = ifelse(is.na(coefficient_source_value), "", coefficient_source_value),
    role = roles,
    stringsAsFactors = FALSE
  )) # distinct selector-role groups without encoding user labels into a delimiter-separated string
  assignments <- lapply(seq_len(nrow(assignment_table)), function(assignment_index) {
    assignment_row <- assignment_table[assignment_index, , drop = FALSE] # one selector-role group
    positions <- which(
      (coefficient_source_type == assignment_row$source_type) &
        (if (nzchar(assignment_row$source_value)) {
          coefficient_source_value == assignment_row$source_value
        } else {
          is.na(coefficient_source_value)
        }) &
        (roles == assignment_row$role)
    ) # coefficient positions governed by this selector and role
    source_type <- assignment_row$source_type # default, family, or marker
    source_value <- assignment_row$source_value # canonical family or exact marker label
    role <- assignment_row$role # intercept or slope shared by this assignment
    source_specification <- switch(
      source_type,
      default = specification,
      family = specification$by_family[[source_value]],
      marker = specification$by_marker[[source_value]]
    ) # normalised role declarations for the selected scope
    list(
      label = if (identical(source_type, "default")) {
        paste0(parameter, "$", role)
      } else {
        paste0(parameter, "[", source_type, "=", source_value, "]$", role)
      },
      positions = positions,
      prior = source_specification[[role]],
      role = role
    )
  }) # complete non-overlapping prior assignments in fitted coefficient order

  list(roles = roles, priors = specification, assignments = assignments)
}

#' Pack role-specific regression priors into one reusable coefficient layout
#'
#' @description
#' The fitting programme uses one common representation for every ordinary
#' regression coefficient. Each coefficient retains its scientific component
#' and intercept/slope role through the declaration chosen in R. Packing the
#' coefficients once avoids repeating a separate Stan prior implementation for
#' longitudinal, event, covariance, association, functional and
#' distributional regressions.
#'
#' @param blocks Named list. Each element contains `roles`, one role label per
#'   coefficient, and `priors`, a named list of normalised prior declarations
#'   for those roles.
#'
#' @return Stan data for coefficient-wise prior transformation and a named
#'   vector giving the first packed position of every block.
#' @keywords internal
#' @noRd
.pack_regression_priors <- function(blocks) {
  family_code <- c(student_t = 1L, normal = 2L, laplace = 3L, horseshoe = 4L) # common distribution-family codes used by the Stan prior interpreter
  coefficient_family <- integer(0) # family code for every packed regression coefficient
  coefficient_mu <- numeric(0) # prior location for every packed regression coefficient
  coefficient_scale <- numeric(0) # ordinary scale for every packed regression coefficient
  coefficient_df <- numeric(0) # Student-t degrees of freedom, or horseshoe local-scale degrees of freedom
  horseshoe_local_index <- integer(0) # positive local-scale index for horseshoe coefficients and zero otherwise
  horseshoe_group_index <- integer(0) # positive shared global-scale group for horseshoe coefficients and zero otherwise
  horseshoe_local_df <- numeric(0) # local half-Student-t degrees of freedom in packed horseshoe order
  horseshoe_global_df <- numeric(0) # global half-Student-t degrees of freedom by scientific prior group
  horseshoe_global_scale <- numeric(0) # fixed global scale by scientific prior group
  horseshoe_slab_df <- numeric(0) # finite-slab degrees of freedom by scientific prior group
  horseshoe_slab_scale <- numeric(0) # finite-slab scale by scientific prior group
  block_start <- integer(length(blocks)) # one-based first coefficient position for every block, or zero for an empty block
  names(block_start) <- names(blocks)
  next_position <- 1L # next unused coefficient position in the common packed layout

  expand_role_values <- function(values, number_coefficients, component, role, field) {
    values <- as.numeric(values) # scalar or coefficient-specific values declared for this role
    if (number_coefficients == 0L) return(numeric(0))
    if (length(values) == 1L) return(rep.int(values, number_coefficients))
    if (length(values) != number_coefficients) {
      cli::cli_abort(c(
        x = "Prior {.arg {component}${role}} provides {length(values)} {field} values for {number_coefficients} coefficients.",
        i = "Supply one value or exactly one value for each {role} coefficient in this component."
      ))
    }
    values
  }

  for (component in names(blocks)) {
    block <- blocks[[component]] # coefficient roles and role-specific prior declarations for one scientific component
    roles <- as.character(block$roles %||% character(0)) # intercept/slope role in the actual fitted coefficient order
    number_coefficients <- length(roles) # dimension of this component in the fitted model
    block_start[[component]] <- if (number_coefficients > 0L) next_position else 0L
    if (number_coefficients == 0L) next
    if (any(!roles %in% c("intercept", "slope"))) {
      cli::cli_abort("Internal prior layout for {.field {component}} contains an unknown coefficient role.")
    }

    block_family <- integer(number_coefficients) # family codes restored to the component's fitted coefficient order
    block_mu <- numeric(number_coefficients) # locations restored to fitted coefficient order
    block_scale <- numeric(number_coefficients) # scales restored to fitted coefficient order
    block_df <- numeric(number_coefficients) # fixed degrees of freedom restored to fitted coefficient order
    block_local_index <- integer(number_coefficients) # horseshoe local-scale map in fitted coefficient order
    block_group_index <- integer(number_coefficients) # horseshoe global-group map in fitted coefficient order

    assignments <- block$assignments %||% lapply(unique(roles), function(role) {
      list(
        label = paste0(component, "$", role),
        positions = which(roles == role),
        prior = block$priors[[role]],
        role = role
      )
    }) # disjoint selector-role assignments, or ordinary intercept/slope assignments
    assigned_positions <- unlist(lapply(assignments, `[[`, "positions"), use.names = FALSE) # all coefficient positions covered by declarations
    if (!identical(sort(as.integer(assigned_positions)), seq_len(number_coefficients))) {
      cli::cli_abort("Internal prior assignments for {.field {component}} must cover every coefficient exactly once.")
    }

    complete_prior_positions <- function(declaration) {
      complete_prior_id <- attr(declaration, "complete_prior_id", exact = TRUE) # identifier of a bare prior spanning several coefficient roles
      if (is.null(complete_prior_id)) return(NULL)
      sort(unique(unlist(lapply(assignments, function(candidate) {
        candidate_id <- attr(candidate$prior, "complete_prior_id", exact = TRUE) # group identifier attached during prior normalisation
        if (identical(candidate_id, complete_prior_id)) candidate$positions else integer(0)
      }), use.names = FALSE)))
    }
    complete_horseshoe_groups <- list() # shared global-scale index for every bare horseshoe declaration spanning several roles

    for (assignment in assignments) {
      role <- assignment$role # scientific intercept or slope role
      role_positions <- as.integer(assignment$positions) # component positions governed by this scoped declaration
      declaration <- assignment$prior # checked Normal, Student-t, Laplace or horseshoe prior for this selector-role group
      if (is.null(declaration)) {
        cli::cli_abort("Internal prior layout for {.field {component}} is missing declaration {.field {assignment$label}}.")
      }
      block_family[role_positions] <- family_code[[declaration$family]]
      complete_positions <- complete_prior_positions(declaration) # full fitted positions governed by one bare component declaration
      expand_assignment_values <- function(values, field) {
        values <- as.numeric(values) # scalar, role-specific, or complete-component hyperparameters
        if (is.null(complete_positions) || length(values) == 1L) {
          return(expand_role_values(values, length(role_positions), component, assignment$label, field))
        }
        if (length(values) != length(complete_positions)) {
          cli::cli_abort(c(
            x = "Prior {.arg {component}} provides {length(values)} {field} values for {length(complete_positions)} coefficients.",
            i = "A bare component prior supplies one value or exactly one value per fitted coefficient; use intercept/slope lists for role-specific vectors."
          ))
        }
        values[match(role_positions, complete_positions)]
      }
      block_mu[role_positions] <- expand_assignment_values(declaration$mu, "location")
      block_scale[role_positions] <- expand_assignment_values(declaration$scale, "scale")
      block_df[role_positions] <- if (is.finite(declaration$df)) declaration$df else 1

      if (identical(declaration$family, "horseshoe")) {
        complete_prior_id <- attr(declaration, "complete_prior_id", exact = TRUE) # bare component horseshoe shares one global scale across its roles
        stored_group_index <- if (is.null(complete_prior_id)) NULL else complete_horseshoe_groups[[complete_prior_id]]
        group_index <- stored_group_index %||% (length(horseshoe_global_df) + 1L) # one shared global shrinkage scale for this declared prior group
        local_indices <- seq.int(length(horseshoe_local_df) + 1L, length.out = length(role_positions)) # coefficient-specific local scales
        block_local_index[role_positions] <- local_indices
        block_group_index[role_positions] <- group_index
        horseshoe_local_df <- c(horseshoe_local_df, rep(declaration$df, length(role_positions)))
        if (is.null(stored_group_index)) {
          horseshoe_global_df <- c(horseshoe_global_df, declaration$global_df)
          horseshoe_global_scale <- c(horseshoe_global_scale, declaration$global_scale)
          horseshoe_slab_df <- c(horseshoe_slab_df, declaration$slab_df)
          horseshoe_slab_scale <- c(horseshoe_slab_scale, declaration$slab_scale)
          if (!is.null(complete_prior_id)) complete_horseshoe_groups[[complete_prior_id]] <- group_index
        }
      }
    }

    coefficient_family <- c(coefficient_family, block_family)
    coefficient_mu <- c(coefficient_mu, block_mu)
    coefficient_scale <- c(coefficient_scale, block_scale)
    coefficient_df <- c(coefficient_df, block_df)
    horseshoe_local_index <- c(horseshoe_local_index, block_local_index)
    horseshoe_group_index <- c(horseshoe_group_index, block_group_index)
    next_position <- next_position + number_coefficients
  }

  list(
    n_regression_prior = as.integer(length(coefficient_family)),
    prior_regression_family = as.array(as.integer(coefficient_family)),
    prior_regression_mu = coefficient_mu,
    prior_regression_scale = coefficient_scale,
    prior_regression_df = coefficient_df,
    prior_regression_horseshoe_local_index = as.array(as.integer(horseshoe_local_index)),
    prior_regression_horseshoe_group_index = as.array(as.integer(horseshoe_group_index)),
    n_regression_horseshoe_local = as.integer(length(horseshoe_local_df)),
    prior_regression_horseshoe_local_df = horseshoe_local_df,
    n_regression_horseshoe_group = as.integer(length(horseshoe_global_df)),
    prior_regression_horseshoe_global_df = horseshoe_global_df,
    prior_regression_horseshoe_global_scale = horseshoe_global_scale,
    prior_regression_horseshoe_slab_df = horseshoe_slab_df,
    prior_regression_horseshoe_slab_scale = horseshoe_slab_scale,
    prior_regression_block_start = block_start
  )
}
