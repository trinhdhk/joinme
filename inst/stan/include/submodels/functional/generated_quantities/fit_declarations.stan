/**
 * @file functional/generated_quantities/fit_declarations.stan
 * @brief Declare fitted affine transformation coefficients.
 *
 * @details
 * The declared arrays expose the intercepts and slopes used within functional association transformations.  No additional calculation is needed in generated quantities.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  vector[estimate_iota_intercept_cv] iota_intercept_cv = iota_intercept_cv_eff; // Role: affine transformation intercept current value.
  vector[estimate_iota_slope_cv] iota_slope_cv = iota_slope_cv_eff; // Role: affine transformation slope current value.
  vector[estimate_iota_intercept_cs] iota_intercept_cs = iota_intercept_cs_eff; // Role: affine transformation intercept current slope.
  vector[estimate_iota_slope_cs] iota_slope_cs = iota_slope_cs_eff; // Role: affine transformation slope current slope.
  vector[M_corr * estimate_iota_intercept_corr] iota_intercept_corr = iota_intercept_corr_eff; // Role: affine transformation intercept correlation.
  vector[M_corr * estimate_iota_slope_corr] iota_slope_corr = iota_slope_corr_eff; // Role: affine transformation slope correlation.
  vector[M_vcov * estimate_iota_intercept_vcov] iota_intercept_vcov = iota_intercept_vcov_eff; // Role: affine transformation intercept covariance.
  vector[M_vcov * estimate_iota_slope_vcov] iota_slope_vcov = iota_slope_vcov_eff; // Role: affine transformation slope covariance.
  vector[estimate_iota_intercept_cv_mean] iota_intercept_cv_mean = iota_intercept_cv_mean_eff; // Role: affine transformation intercept current value mean.
  vector[estimate_iota_slope_cv_mean] iota_slope_cv_mean = iota_slope_cv_mean_eff; // Role: affine transformation slope current value mean.
  vector[estimate_iota_intercept_cv_marker] iota_intercept_cv_marker = iota_intercept_cv_marker_eff; // Role: affine transformation intercept current value marker.
  vector[estimate_iota_slope_cv_marker] iota_slope_cv_marker = iota_slope_cv_marker_eff; // Role: affine transformation slope current value marker.
  vector[estimate_iota_intercept_cs_mean] iota_intercept_cs_mean = iota_intercept_cs_mean_eff; // Role: affine transformation intercept current slope mean.
  vector[estimate_iota_slope_cs_mean] iota_slope_cs_mean = iota_slope_cs_mean_eff; // Role: affine transformation slope current slope mean.
  vector[estimate_iota_intercept_cs_marker] iota_intercept_cs_marker = iota_intercept_cs_marker_eff; // Role: affine transformation intercept current slope marker.
  vector[estimate_iota_slope_cs_marker] iota_slope_cs_marker = iota_slope_cs_marker_eff; // Role: affine transformation slope current slope marker.


