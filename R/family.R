#' Family and Transformation Utilities for Mixed-Family Joint Models
#'
#' @description
#' Helper functions for specifying mixed families per marker and transformation compositions.
#'
#' @name family_utils
NULL

#' JoiNMe longitudinal family specification
#'
#' @description
#' Creates a family specification object for `joinme(..., families = ...)` with an
#' optional custom longitudinal link (or inverse-link) per marker family.
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
#'   `Phi` is the standard normal CDF and is the inverse link for a probit
#'   model. In contrast, `inv_Phi`, `qnorm`, and `probit` denote the standard
#'   normal quantile function.
#'
#' @return Object of class `"JoiNMe_family_spec"`.
#' @examples
#' jm_family("poisson", link = "log")
#' jm_family("poisson", link = ~ log(x))
#' jm_family("bernoulli", inv_link = ~ inv_logit(x))
#' jm_family("bernoulli", link = ~ inv_Phi(x))
#' jm_family("bernoulli", inv_link = ~ Phi(x))
#' @aliases jm_family
#' @export
joinme_family <- function(name, link = NULL, inv_link = NULL) {
  fam_code <- .parse_family(name)
  fam_name <- .family_code_to_name(fam_code)

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
      inv_link = inv_link_bc
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
#' @param mode Character scalar `"warning"` or `"error"` 
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
#'   `link_name`, and compiled `inv_link_bc`.
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
    return(list(
      family_code = as.integer(fam_code),
      link_name = link_name,
      inv_link_bc = inv_link_bc
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
    return(list(
      family_code = as.integer(fam_code),
      link_name = link_name,
      inv_link_bc = inv_link_bc
    ))
  }

  fam_code <- .parse_family(x)
  inv_link_bc <- .inv_link_bc_from_name(.default_link_for_family(fam_code))
  list(
    family_code = as.integer(fam_code),
    link_name = .canonical_link_from_inv_link_bc(inv_link_bc),
    inv_link_bc = inv_link_bc
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
#' @return Integer vector of family codes (length D)
#' @keywords internal
#' @noRd
.validate_family_list <- function(families, D, dataLong, marker_var, y_var) {
  # Validate per-marker families against observed data
  # Convert to list if single value
  if (length(families) == 1 && !is.list(families)) {
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
    inv_link_const = inv_link_const
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

#' Build prior specification for JoiNMe model
#'
#' @description
#' Create prior specifications for fixed effects (beta), baseline hazard coefficients (alpha),
#' and correlation structure (LKJ prior).
#'
#' @param beta_prior Named list or vector. If vector: uniform scale for all beta.
#'   If list with 'scale' or 'sd': uses normal(0, sd) for all.
#'   Example: list(scale = 2) or list(scale = c(2, 1, 1))
#' @param alpha_prior Scale for baseline hazard coefficients. Default: normal(0, 2)
#' @param iota_prior Scale for fit-only affine-shift intercept and slope
#'   parameters in functional association transforms. Default: normal(0, 1)
#' @param lkj_prior Concentration parameter for LKJ correlation prior. Default: 1 (uniform)
#'
#' @return List with prior specifications compatible with Stan
#' @keywords internal
#' @noRd
.build_priors <- function(
  beta_prior = NULL,
  alpha_prior = NULL,
  iota_prior = NULL,
  lkj_prior = NULL
) {
  # Create full prior spec list consumed by joinme_standata()
  priors <- list()

  # Beta priors (fixed effects)
  if (!is.null(beta_prior)) {
    if (is.numeric(beta_prior)) {
      if (length(beta_prior) == 1) {
        priors$beta_scale <- rep(beta_prior, 100) # Will be truncated to P
      } else {
        priors$beta_scale <- as.numeric(beta_prior)
      }
    } else if (is.list(beta_prior)) {
      priors$beta_scale <- beta_prior$scale %||% beta_prior$sd %||% 2
      if (length(priors$beta_scale) == 1) {
        priors$beta_scale <- rep(priors$beta_scale, 100)
      }
    }
  } else {
    priors$beta_scale <- rep(2.0, 100) # Default: normal(0, 2)
  }

  # Alpha priors (baseline hazard)
  if (!is.null(alpha_prior)) {
    if (is.numeric(alpha_prior)) {
      priors$alpha_scale <- alpha_prior[1]
    } else if (is.list(alpha_prior)) {
      priors$alpha_scale <- alpha_prior$scale %||% alpha_prior$sd %||% 2
    }
  } else {
    priors$alpha_scale <- 2.0 # Default
  }

  # Iota priors (fit-only functional transform intercept/slope shifts)
  if (!is.null(iota_prior)) {
    if (is.numeric(iota_prior)) {
      priors$iota_scale <- iota_prior[1]
    } else if (is.list(iota_prior)) {
      priors$iota_scale <- iota_prior$scale %||% iota_prior$sd %||% 1
    }
  } else {
    priors$iota_scale <- 1.0
  }
  # LKJ prior (correlation)
  if (!is.null(lkj_prior)) {
    priors$lkj_eta <- as.numeric(lkj_prior)
  } else {
    priors$lkj_eta <- 1.0 # Default: uniform
  }

  priors
}
