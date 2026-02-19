#' Family and Transformation Utilities for Mixed-Family Joint Models
#'
#' @description
#' Helper functions for specifying mixed families per marker and transformation compositions.
#'
#' @name family_utils
NULL

# File overview:
# - Map between family names and integer codes.
# - Validate outcome ranges and parse transform compositions.

#' Parse family string to numeric code
#' @param family String family name (e.g., "gaussian", "student_t", "binomial") or vector of names
#' @return Integer family code (1-11) or vector of codes
#' @keywords internal
.parse_family <- function(family) {
  # Map human-readable family names to Stan integer codes
  family <- tolower(as.character(family))
  
  # Handle vector input
  if (length(family) > 1) {
    return(sapply(family, .parse_family))
  }
  
  switch(family,
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
      "Unknown family: ", family,
      ". Must be one of: gaussian, student_t, bernoulli, binomial, poisson, negbin2, ",
      "skew_normal, double_exponential, skew_double_exponential, beta, cumulative_logit"
    )
  )
}

#' Convert family code to string
#' @param fam Integer family code (1-11), can be vector
#' @return Character family name(s)
#' @keywords internal
.family_code_to_name <- function(fam) {
  # Map Stan integer codes to canonical family names
  fam <- as.integer(fam)
  
  # Handle vector input
  if (length(fam) > 1) {
    return(sapply(fam, .family_code_to_name))
  }
  
  switch(as.character(fam),
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
.validate_family_outcome <- function(family, y_range, name = "y") {
  # Validate observed outcome ranges by family requirements
  fam_name <- .family_code_to_name(family)
  
  if (family %in% c(1, 2, 7, 8, 9)) {
    # Gaussian, student_t, skew normal, double exponential, skew double exponential: any real
    return(TRUE)
  } else if (family == 10) {
    # Beta: must be strictly between 0 and 1
    if (any(y_range <= 0) || any(y_range >= 1)) {
      stop("Beta family requires ", name, " in (0,1), but got values in [",
           min(y_range), ", ", max(y_range), "]")
    }
  } else if (family == 11) {
    # Cumulative logit (ordinal): must be integer categories >= 1
    if (any(y_range < 1) || any(y_range != as.integer(y_range))) {
      stop("Cumulative logit family requires ", name, " to be integer categories >= 1, but got values in [",
           min(y_range), ", ", max(y_range), "]")
    }
  } else if (family == 3) {
    # Bernoulli: must be 0 or 1
    if (!all(y_range %in% c(0, 1))) {
      stop("Bernoulli family requires ", name, " in {0,1}, but got values in [", 
           min(y_range), ", ", max(y_range), "]")
    }
  } else if (family == 4) {
    # Binomial: integers in [0,n]
    if (any(y_range < 0) || any(y_range != as.integer(y_range))) {
      stop("Binomial family requires ", name, " to be non-negative integers, but got values in [", 
           min(y_range), ", ", max(y_range), "]")
    }
  } else if (family %in% c(5, 6)) {
    # Poisson, negbin2: non-negative integers
    if (any(y_range < 0) || any(y_range != as.integer(y_range))) {
      stop(fam_name, " family requires ", name, " to be non-negative integers, but got values in [", 
           min(y_range), ", ", max(y_range), "]")
    }
  }
  
  return(TRUE)
}

#' Parse transformation composition string
#' @description
#' Convert strings like "2_3" or "log_sqrt" into array of transform IDs.
#' Transform codes: 0=none, 1=id, 2=log, 3=sqrt, 4=logit, 5=inv, 6=exp
#' 
#' @param tf_string Character transformation specification. Can be:
#'   - "none" or "" or 0: identity
#'   - Numeric codes: "2", "3", "2_3", "1_2_3"
#'   - Named: "log", "sqrt", "log_sqrt", "sqrt_log"
#' @param max_length Maximum array length (default 7)
#' @return List with: codes (integer array), names (character), description (character)
#' 
#' @keywords internal
.parse_tf_composition <- function(tf_string, max_length = 7) {
  # Parse transform composition into numeric codes for Stan
  tf_string <- as.character(tf_string)
  tf_string <- tolower(trimws(tf_string))
  
  # Map named transforms
  tf_map <- list(
    "none" = c(),
    "id" = c(1L),
    "log" = c(2L),
    "sqrt" = c(3L),
    "logit" = c(4L),
    "inv" = c(5L),
    "exp" = c(6L),
    "log_sqrt" = c(2L, 3L),  # Apply sqrt first, then log: sqrt(x) then log(sqrt(x))
    "sqrt_log" = c(3L, 2L),  # Apply log first, then sqrt: log(x) then sqrt(log(x))
    "log_log" = c(2L, 2L),
    "sqrt_sqrt" = c(3L, 3L),
    "logit_log" = c(4L, 2L),
    "exp_log" = c(6L, 2L)
  )
  
  if (tf_string %in% names(tf_map)) {
    tf_codes <- tf_map[[tf_string]]
  } else if (tf_string == "0" || tf_string == "") {
    tf_codes <- c()
  } else if (grepl("^[0-6](_[0-6])*$", tf_string)) {
    # Parse numeric codes like "2_3"
    parts <- as.integer(strsplit(tf_string, "_")[[1]])
    if (any(parts < 0 | parts > 6)) {
      stop("Transform codes must be in 0-6, got: ", tf_string)
    }
    tf_codes <- parts
  } else {
    stop("Unknown transformation: ", tf_string, 
         ". Use codes (0-6) or names (none, log, sqrt, log_sqrt, etc)")
  }
  
  if (length(tf_codes) == 0) tf_codes <- 1L  # identity
  if (length(tf_codes) > max_length) {
    stop("Transformation composition too long (max ", max_length, " transforms)")
  }
  
  # Pad to max_length with 0 (will be skipped in Stan)
  tf_padded <- c(tf_codes, rep(0L, max_length - length(tf_codes)))
  
  # Create names
  tf_names <- c("none", "id", "log", "sqrt", "logit", "inv", "exp")
  tf_name_str <- paste(tf_names[tf_codes + 1], collapse = " then ")
  if (length(tf_codes) == 1 && tf_codes[1] == 1) tf_name_str <- "identity"
  
  list(
    codes = tf_padded,
    n_codes = length(tf_codes),
    description = tf_name_str
  )
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
.validate_family_list <- function(families, D, dataLong, marker_var, y_var) {
  # Validate per-marker families against observed data
  # Convert to list if single value
  if (length(families) == 1 && !is.list(families)) {
    families <- rep(list(families), D)
  } else if (!is.list(families)) {
    families <- as.list(families)
  }
  
  if (length(families) != D) {
    stop("families must have length ", D, " (number of markers), got ", length(families))
  }
  
  # Parse each family
  family_codes <- integer(D)
  for (d in seq_along(families)) {
    family_codes[d] <- .parse_family(families[[d]])
  }
  
  # Validate against actual data
  # Make sure marker_var is a factor
  if (!is.factor(dataLong[[marker_var]])) {
    dataLong[[marker_var]] <- factor(dataLong[[marker_var]])
  }
  
  marker_levels <- levels(dataLong[[marker_var]])
  if (length(marker_levels) != D) {
    stop("Number of unique markers (", length(marker_levels), ") does not match D (", D, ")")
  }
  
  for (d in seq_len(D)) {
    marker_name <- marker_levels[d]
    y_vals <- dataLong[[y_var]][dataLong[[marker_var]] == marker_name]
    if (length(y_vals) > 0) {
      .validate_family_outcome(family_codes[d], y_vals, 
                              name = paste0(y_var, "[marker=", marker_name, "]"))
    }
  }
  
  family_codes
}

#' Get required distributional parameters for a family
#' @param family Integer family code
#' @return Character vector of parameter names
#' @keywords internal
.family_distrib_params <- function(family) {
  # Enumerate distributional regression parameters per family
  switch(as.integer(family),
    "1" = c("sigma"),        # gaussian
    "2" = c("sigma", "nu"),  # student_t
    "3" = c(),               # bernoulli
    "4" = c(),               # binomial (uses trials)
    "5" = c(),               # poisson
    "6" = c("phi"),          # negbin2
    "7" = c("sigma", "alpha"), # skew_normal (scale + skew)
    "8" = c("sigma"),           # double_exponential (scale)
    "9" = c("sigma", "tau_sde"), # skew_double_exponential (scale + skew)
    "10" = c("phi_beta"),       # beta (precision)
    "11" = c(),                 # cumulative_logit (cutpoints)
    stop("Unknown family code: ", family)
  )
}

#' Generate prior predictive for outcome given family
#' @description
#' Simple prior predictive sampling for validation/simulation
#' 
#' @param n Sample size
#' @param eta Linear predictor
#' @param family Integer family code
#' @param sigma Standard deviation (if needed)
#' @param nu Degrees of freedom (if needed, student_t)
#' @param phi Dispersion (if needed, negbin2)
#' @param phi_beta Precision for beta family (if needed)
#' @param tau_sde Skew parameter for skew double exponential (if needed)
#' @param trials Number of trials (if needed, binomial)
#' 
#' @return Numeric vector of samples
#' @keywords internal
.sample_from_family <- function(
  n,
  eta,
  family,
  sigma = 1,
  nu = 4,
  phi = 1,
  phi_beta = 10,
  tau_sde = 0.5,
  trials = 1,
  skew = 0
) {
  # Prior predictive helper for a single family
  family <- as.integer(family)
  
  if (family == 1) {
    # Gaussian
    rnorm(n, eta, sigma)
  } else if (family == 2) {
    # Student-t
    rt(n, df = nu) * sigma + eta
  } else if (family == 3) {
    # Bernoulli
    rbinom(n, size = 1, prob = plogis(eta))
  } else if (family == 4) {
    # Binomial
    rbinom(n, size = trials, prob = plogis(eta))
  } else if (family == 5) {
    # Poisson
    rpois(n, lambda = exp(eta))
  } else if (family == 6) {
    # Negative binomial 2
    mu <- exp(eta)
    rnbinom(n, size = phi, mu = mu)
  } else if (family == 7) {
    # Skew-normal (using fGarch::rssnorm if available, else fallback)
    if (requireNamespace("fGarch", quietly = TRUE)) {
      fGarch::rssnorm(n, mean = eta, sd = sigma, xi = skew)
    } else {
      rnorm(n, mean = eta, sd = sigma)
    }
  } else if (family == 8) {
    # Double exponential (Laplace)
    # Sample by drawing an exponential magnitude and random sign.
    sign <- sample(c(-1, 1), size = n, replace = TRUE)
    eta + sign * rexp(n, rate = 1 / sigma)
  } else if (family == 9) {
    # Skew double exponential (asymmetric Laplace)
    # Piecewise exponential: probability mass split by tau_sde.
    u <- runif(n)
    e <- rexp(n, rate = 1)
    eta + ifelse(u < tau_sde, e * sigma / tau_sde, -e * sigma / (1 - tau_sde))
  } else if (family == 10) {
    # Beta: mean is logistic(eta), precision is phi_beta
    mu <- plogis(eta)
    shape1 <- pmax(mu * phi_beta, 1e-6)
    shape2 <- pmax((1 - mu) * phi_beta, 1e-6)
    rbeta(n, shape1 = shape1, shape2 = shape2)
  } else if (family == 11) {
    # Cumulative logit: return category labels (1..K)
    # This is a placeholder sampler; use with a valid cutpoint structure externally.
    stop("Sampling for cumulative logit requires cutpoints; use model-based sampling instead.")
  } else {
    stop("Unknown family: ", family)
  }
}

#' Build prior specification for joinme model
#'
#' @description
#' Create prior specifications for fixed effects (beta), baseline hazard coefficients (alpha),
#' and correlation structure (LKJ prior).
#'
#' @param beta_prior Named list or vector. If vector: uniform scale for all beta.
#'   If list with 'scale' or 'sd': uses normal(0, sd) for all.
#'   Example: list(scale = 2) or list(scale = c(2, 1, 1))
#' @param alpha_prior Scale for baseline hazard coefficients. Default: normal(0, 2)
#' @param lkj_prior Concentration parameter for LKJ correlation prior. Default: 1 (uniform)
#'
#' @return List with prior specifications compatible with Stan
#' @keywords internal
.build_priors <- function(beta_prior = NULL, alpha_prior = NULL, lkj_prior = NULL) {
  # Create full prior spec list consumed by joinme_standata()
  priors <- list()
  
  # Beta priors (fixed effects)
  if (!is.null(beta_prior)) {
    if (is.numeric(beta_prior)) {
      if (length(beta_prior) == 1) {
        priors$beta_scale <- rep(beta_prior, 100)  # Will be truncated to P
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
    priors$beta_scale <- rep(2.0, 100)  # Default: normal(0, 2)
  }
  
  # Alpha priors (baseline hazard)
  if (!is.null(alpha_prior)) {
    if (is.numeric(alpha_prior)) {
      priors$alpha_scale <- alpha_prior[1]
    } else if (is.list(alpha_prior)) {
      priors$alpha_scale <- alpha_prior$scale %||% alpha_prior$sd %||% 2
    }
  } else {
    priors$alpha_scale <- 2.0  # Default
  }
  
  # LKJ prior (correlation)
  if (!is.null(lkj_prior)) {
    priors$lkj_eta <- as.numeric(lkj_prior)
  } else {
    priors$lkj_eta <- 1.0  # Default: uniform
  }
  
  priors
}

#' Parse transformation composition list to Stan format
#' @param tf_composition List with keys: tf_cv_mean_comp, tf_cv_marker_comp, tf_cs_mean_comp, tf_cs_marker_comp
#' @return List with codes and n_codes for each composition
#' @keywords internal
.parse_tf_composition_list <- function(tf_composition) {
  # Expand named composition list into per-term code arrays
  tf_names <- c("tf_cv_mean_comp", "tf_cv_marker_comp", "tf_cs_mean_comp", "tf_cs_marker_comp")
  result <- list()
  
  for (name in tf_names) {
    tf_str <- tf_composition[[name]] %||% "none"
    parsed <- .parse_tf_composition(tf_str)
    result[[name]] <- list(
      codes = as.integer(parsed$codes),
      n_codes = as.integer(parsed$n_codes),
      description = parsed$description
    )
  }
  
  result
}
