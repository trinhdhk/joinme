/**
 * @file longitudinal/generated_quantities/fit_declarations.stan
 * @brief Declare the longitudinal pointwise log likelihood.
 *
 * @details
 * The declaration file introduces the longitudinal pointwise log-likelihood
 * vector. The calculation file rebuilds each observation under its
 * marker-specific response family. Coefficients and random-effect scales need
 * no generated copies because their sampled values already use original time.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  vector[N] log_lik_long; // Role: log lik long.
