#' Validate that fitting coefficient declarations are distributions
#'
#' @description
#' Fitting priors and simulation truths share a scientific component hierarchy,
#' but only the latter may contain fixed numerical coefficients. This validator
#' prevents a numerical truth from being mistaken for a probability
#' distribution in `jm_prior()`.
#'
#' @param prior_request Named fitting-prior declarations before normalisation.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.assert_joinme_prior_distributions <- function(prior_request) {
  # A prior object contains numeric hyperparameters internally, so recursion
  # stops as soon as a recognised distribution declaration is reached. Only a
  # bare number occupying a coefficient role is prohibited.
  inspect_coefficient_declaration <- function(declaration, context) {
    if (is.null(declaration) ||
        inherits(declaration, "joinme_prior_spec") ||
        inherits(declaration, "joinme_lkj_prior")) {
      return(invisible(NULL))
    }
    if (is.numeric(declaration)) {
      cli::cli_abort(c(
        x = "{.arg {context}} is a fixed numerical value, not a prior distribution.",
        i = "Use a {.fn prior_normal}, {.fn prior_student_t}, {.fn prior_laplace}, or {.fn prior_horseshoe} declaration.",
        i = "Use {.fn jm_truth} when the number is a data-generating coefficient."
      ))
    }
    if (is.list(declaration)) {
      for (field_name in names(declaration) %||% seq_along(declaration)) {
        inspect_coefficient_declaration(
          declaration[[field_name]],
          paste0(context, "$", field_name)
        )
      }
    }
    invisible(NULL)
  }

  inspect_coefficient_declaration(prior_request$intercept, "intercept")
  inspect_coefficient_declaration(prior_request$slope, "slope")
  for (component_name in c("longitudinal", "survival", "baseline", "vcov", "assoc", "functional")) {
    inspect_coefficient_declaration(prior_request[[component_name]], component_name)
  }
  for (parameter_name in names(prior_request$distributional) %||% character(0)) {
    inspect_coefficient_declaration(
      prior_request$distributional[[parameter_name]],
      parameter_name
    )
  }
  if (is.list(prior_request$marker_weights)) {
    inspect_coefficient_declaration(
      prior_request$marker_weights$intercept,
      "marker_weights$intercept"
    )
  }
  inspect_coefficient_declaration(prior_request$marker, "marker")
  if (is.list(prior_request$class)) {
    inspect_coefficient_declaration(prior_request$class$slope, "class$slope")
    inspect_coefficient_declaration(prior_request$class$family, "class$family")
  }
  invisible(NULL)
}

#' Remove fixed truths when constructing a recovery prior
#'
#' @description
#' A fixed data-generating coefficient has no probability density and therefore
#' cannot enter the fitting prior. This helper replaces such coefficients by
#' `NULL`, allowing the usual fitting defaults to apply, whilst retaining
#' structural numeric quantities such as marker offsets and Dirichlet
#' concentrations.
#'
#' @param truth_request Truth components arranged in the prior hierarchy.
#'
#' @return A named list suitable for `.joinme_priors_()`.
#' @keywords internal
#' @noRd
.remove_joinme_truth_constants <- function(truth_request) {
  remove_coefficients <- function(declaration) {
    if (is.null(declaration) ||
        inherits(declaration, "joinme_prior_spec") ||
        inherits(declaration, "joinme_lkj_prior")) {
      return(declaration)
    }
    if (is.numeric(declaration)) return(NULL)
    if (!is.list(declaration)) return(declaration)
    output <- declaration # named coefficient roles whose numerical leaves are removed recursively
    for (field_name in names(output) %||% character(0)) {
      output[[field_name]] <- remove_coefficients(output[[field_name]])
    }
    output
  }

  fitting_request <- truth_request # local fitting declaration assembled without altering the truth object
  fitting_request$intercept <- remove_coefficients(fitting_request$intercept)
  fitting_request$slope <- remove_coefficients(fitting_request$slope)
  for (component_name in c("longitudinal", "survival", "baseline", "vcov", "assoc", "functional")) {
    fitting_request[[component_name]] <- remove_coefficients(fitting_request[[component_name]])
  }
  fitting_request$distributional <- remove_coefficients(fitting_request$distributional)
  if (is.list(fitting_request$marker_weights)) {
    fitting_request$marker_weights$intercept <- remove_coefficients(
      fitting_request$marker_weights$intercept
    )
  }
  if (is.list(fitting_request$class)) {
    fitting_request$class$slope <- remove_coefficients(fitting_request$class$slope)
    fitting_request$class$family <- remove_coefficients(fitting_request$class$family)
  }
  fitting_request
}

#' Draw a correlation matrix from an LKJ distribution
#'
#' @description
#' Uses the extended-onion construction for the Cholesky factor. Row `k`
#' receives a random direction on the unit sphere and a squared radius from
#' `Beta((k - 1) / 2, eta + (dimension - k) / 2)`. The resulting lower
#' triangular factor has the `lkj_corr_cholesky(eta)` distribution used by the
#' Stan model; multiplying it by its transpose gives the corresponding LKJ
#' correlation matrix.
#'
#' @param dimension Positive correlation-matrix dimension.
#' @param eta Positive LKJ concentration.
#'
#' @return A list containing `corr` and its lower Cholesky factor `Lcorr`.
#' @keywords internal
#' @noRd
.sim_draw_lkj_correlation <- function(dimension, eta) {
  dimension <- as.integer(dimension) # number of random-effect coefficients in the covariance block
  eta <- as.numeric(eta) # LKJ concentration shared with the Stan fitting model
  if (length(dimension) != 1L || is.na(dimension) || dimension < 0L) {
    cli::cli_abort("The LKJ correlation dimension must be one non-negative integer.")
  }
  if (length(eta) != 1L || !is.finite(eta) || eta <= 0) {
    cli::cli_abort("The LKJ concentration must be one positive finite number.")
  }
  if (dimension == 0L) {
    return(list(
      corr = matrix(0, nrow = 0L, ncol = 0L),
      Lcorr = matrix(0, nrow = 0L, ncol = 0L)
    ))
  }
  if (dimension == 1L) {
    return(list(corr = matrix(1, 1L, 1L), Lcorr = matrix(1, 1L, 1L)))
  }

  # Construct the Cholesky rows independently according to the onion
  # representation. Normalising a standard-Normal vector gives a direction
  # uniformly distributed on the required sphere.
  correlation_cholesky <- matrix(0, nrow = dimension, ncol = dimension) # lower factor generated under lkj_corr_cholesky(eta)
  correlation_cholesky[1L, 1L] <- 1
  for (row_index in 2:dimension) {
    direction <- stats::rnorm(row_index - 1L) # isotropic Gaussian vector before normalisation to the unit sphere
    direction_norm <- sqrt(sum(direction^2)) # Euclidean norm used to obtain the uniform spherical direction
    while (!is.finite(direction_norm) || direction_norm <= 0) {
      direction <- stats::rnorm(row_index - 1L)
      direction_norm <- sqrt(sum(direction^2))
    }
    direction <- direction / direction_norm
    squared_radius <- stats::rbeta(
      1L,
      shape1 = 0.5 * (row_index - 1L),
      shape2 = eta + 0.5 * (dimension - row_index)
    ) # row length implied by the LKJ extended-onion law
    correlation_cholesky[row_index, seq_len(row_index - 1L)] <-
      sqrt(squared_radius) * direction
    correlation_cholesky[row_index, row_index] <- sqrt(1 - squared_radius)
  }
  correlation <- tcrossprod(correlation_cholesky) # positive-definite correlation matrix with exact unit diagonal
  diag(correlation) <- 1
  list(corr = correlation, Lcorr = correlation_cholesky)
}
