/**
 * @file assoc/generated_quantities/fit_declarations.stan
 * @brief Declare association coefficients used by the event hazard.
 *
 * @details
 * The declared outputs include channel activation and identifying sign conventions.  No later calculation is required because the transformed parameters already hold the effective coefficients.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  /**
  * @brief Association coefficients actually used in the event hazard.
   */
  real alpha_cv_total = a_cv_total; // Role: shape current value total.
  real alpha_cs_total = a_cs_total; // Role: shape current slope total.
  real alpha_cv_mean = a_cv_mean; // Role: shape current value mean.
  real alpha_cs_mean = a_cs_mean; // Role: shape current slope mean.
  vector[M_corr] alpha_corr = a_corr; // Role: shape correlation.
  vector[M_vcov] alpha_vcov = a_vcov; // Role: shape covariance.



