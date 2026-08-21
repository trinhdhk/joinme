/**
 * @file functional/data/fit.stan
 * @brief Data declarations for transformations of association channels.
 *
 * @details
 * Every channel has an explicit transformation mode, optional fitted affine
 * intercept and slope, bytecode instructions and constants, and optional
 * monotone spline information.  These declarations describe the function
 * applied to a longitudinal feature before its association coefficient enters
 * the log hazard.
 */
  int<lower=0> n_iota_prior; // total number of fitted affine-shift coefficients
  /* Transform configuration */
  int<lower=0, upper=7> tf_mode_cv_tot;     // transform mode for total CV
  int<lower=0, upper=7> tf_mode_cs_tot;     // transform mode for total CS
  int<lower=0, upper=7> tf_mode_corr;       // transform mode for corr term
  int<lower=0, upper=7> tf_mode_vcov;       // transform mode for vcov term
  int<lower=0, upper=7> tf_mode_cv_mean;    // transform mode for mean CV
  int<lower=0, upper=7> tf_mode_cv_marker;  // transform mode for marker CV
  int<lower=0, upper=7> tf_mode_cs_mean;    // transform mode for mean CS
  int<lower=0, upper=7> tf_mode_cs_marker;  // transform mode for marker CS

  /* Fit-only affine shift for functional transforms */
  int<lower=0> estimate_iota_intercept_cv; // Role: whether to estimate affine transformation intercept current value.
  int<lower=0> estimate_iota_slope_cv; // Role: whether to estimate affine transformation slope current value.
  int<lower=0> estimate_iota_intercept_cs; // Role: whether to estimate affine transformation intercept current slope.
  int<lower=0> estimate_iota_slope_cs; // Role: whether to estimate affine transformation slope current slope.
  int<lower=0> estimate_iota_intercept_corr; // Role: whether to estimate affine transformation intercept correlation.
  int<lower=0> estimate_iota_slope_corr; // Role: whether to estimate affine transformation slope correlation.
  int<lower=0> estimate_iota_intercept_vcov; // Role: whether to estimate affine transformation intercept covariance.
  int<lower=0> estimate_iota_slope_vcov; // Role: whether to estimate affine transformation slope covariance.
  int<lower=0> estimate_iota_intercept_cv_mean; // Role: whether to estimate affine transformation intercept current value mean.
  int<lower=0> estimate_iota_slope_cv_mean; // Role: whether to estimate affine transformation slope current value mean.
  int<lower=0> estimate_iota_intercept_cv_marker; // Role: whether to estimate affine transformation intercept current value marker.
  int<lower=0> estimate_iota_slope_cv_marker; // Role: whether to estimate affine transformation slope current value marker.
  int<lower=0> estimate_iota_intercept_cs_mean; // Role: whether to estimate affine transformation intercept current slope mean.
  int<lower=0> estimate_iota_slope_cs_mean; // Role: whether to estimate affine transformation slope current slope mean.
  int<lower=0> estimate_iota_intercept_cs_marker; // Role: whether to estimate affine transformation intercept current slope marker.
  int<lower=0> estimate_iota_slope_cs_marker; // Role: whether to estimate affine transformation slope current slope marker.

  /* Functional bytecode specifications */
  int<lower=0> n_functional_ops_cv;         // op count for total CV
  array[n_functional_ops_cv] int<lower=0, upper=27> functional_ops_cv; // bytecode stream
  array[n_functional_ops_cv] int<lower=0, upper=estimate_iota_intercept_cv> functional_iota_intercept_idx_cv; // Role: functional transformation affine transformation intercept index current value.
  array[n_functional_ops_cv] int<lower=0, upper=estimate_iota_slope_cv> functional_iota_slope_idx_cv; // Role: functional transformation affine transformation slope index current value.
  int<lower=0> n_const_cv;                  // constants used by CV bytecode
  vector[n_const_cv] const_data_cv;         // constants used by CV bytecode

  int<lower=0> n_functional_ops_cs;         // op count for total CS
  array[n_functional_ops_cs] int<lower=0, upper=27> functional_ops_cs; // bytecode stream
  array[n_functional_ops_cs] int<lower=0, upper=estimate_iota_intercept_cs> functional_iota_intercept_idx_cs; // Role: functional transformation affine transformation intercept index current slope.
  array[n_functional_ops_cs] int<lower=0, upper=estimate_iota_slope_cs> functional_iota_slope_idx_cs; // Role: functional transformation affine transformation slope index current slope.
  int<lower=0> n_const_cs;                  // constants used by CS bytecode
  vector[n_const_cs] const_data_cs;         // constants used by CS bytecode

  int<lower=0> n_functional_ops_corr;       // op count for corr
  array[n_functional_ops_corr] int<lower=0, upper=27> functional_ops_corr; // bytecode stream
  array[n_functional_ops_corr] int<lower=0, upper=estimate_iota_intercept_corr> functional_iota_intercept_idx_corr; // Role: functional transformation affine transformation intercept index correlation.
  array[n_functional_ops_corr] int<lower=0, upper=estimate_iota_slope_corr> functional_iota_slope_idx_corr; // Role: functional transformation affine transformation slope index correlation.
  int<lower=0> n_const_corr;                // constants used by corr bytecode
  vector[n_const_corr] const_data_corr;     // constants used by corr bytecode

  int<lower=0> n_functional_ops_vcov;       // op count for vcov
  array[n_functional_ops_vcov] int<lower=0, upper=27> functional_ops_vcov; // bytecode stream
  array[n_functional_ops_vcov] int<lower=0, upper=estimate_iota_intercept_vcov> functional_iota_intercept_idx_vcov; // Role: functional transformation affine transformation intercept index covariance.
  array[n_functional_ops_vcov] int<lower=0, upper=estimate_iota_slope_vcov> functional_iota_slope_idx_vcov; // Role: functional transformation affine transformation slope index covariance.
  int<lower=0> n_const_vcov;                // constants used by vcov bytecode
  vector[n_const_vcov] const_data_vcov;     // constants used by vcov bytecode

  int<lower=0> n_functional_ops_cv_mean;    // op count for mean CV
  array[n_functional_ops_cv_mean] int<lower=0, upper=27> functional_ops_cv_mean; // bytecode stream
  array[n_functional_ops_cv_mean] int<lower=0, upper=estimate_iota_intercept_cv_mean> functional_iota_intercept_idx_cv_mean; // Role: functional transformation affine transformation intercept index current value mean.
  array[n_functional_ops_cv_mean] int<lower=0, upper=estimate_iota_slope_cv_mean> functional_iota_slope_idx_cv_mean; // Role: functional transformation affine transformation slope index current value mean.
  int<lower=0> n_const_cv_mean;             // constants used by mean CV bytecode
  vector[n_const_cv_mean] const_data_cv_mean; // constants used by mean CV bytecode

  int<lower=0> n_functional_ops_cv_marker;  // op count for marker CV
  array[n_functional_ops_cv_marker] int<lower=0, upper=27> functional_ops_cv_marker; // bytecode stream
  array[n_functional_ops_cv_marker] int<lower=0, upper=estimate_iota_intercept_cv_marker> functional_iota_intercept_idx_cv_marker; // Role: functional transformation affine transformation intercept index current value marker.
  array[n_functional_ops_cv_marker] int<lower=0, upper=estimate_iota_slope_cv_marker> functional_iota_slope_idx_cv_marker; // Role: functional transformation affine transformation slope index current value marker.
  int<lower=0> n_const_cv_marker;           // constants used by marker CV bytecode
  vector[n_const_cv_marker] const_data_cv_marker; // constants used by marker CV bytecode

  int<lower=0> n_functional_ops_cs_mean;    // op count for mean CS
  array[n_functional_ops_cs_mean] int<lower=0, upper=27> functional_ops_cs_mean; // bytecode stream
  array[n_functional_ops_cs_mean] int<lower=0, upper=estimate_iota_intercept_cs_mean> functional_iota_intercept_idx_cs_mean; // Role: functional transformation affine transformation intercept index current slope mean.
  array[n_functional_ops_cs_mean] int<lower=0, upper=estimate_iota_slope_cs_mean> functional_iota_slope_idx_cs_mean; // Role: functional transformation affine transformation slope index current slope mean.
  int<lower=0> n_const_cs_mean;             // constants used by mean CS bytecode
  vector[n_const_cs_mean] const_data_cs_mean; // constants used by mean CS bytecode

  int<lower=0> n_functional_ops_cs_marker;  // op count for marker CS
  array[n_functional_ops_cs_marker] int<lower=0, upper=27> functional_ops_cs_marker; // bytecode stream
  array[n_functional_ops_cs_marker] int<lower=0, upper=estimate_iota_intercept_cs_marker> functional_iota_intercept_idx_cs_marker; // Role: functional transformation affine transformation intercept index current slope marker.
  array[n_functional_ops_cs_marker] int<lower=0, upper=estimate_iota_slope_cs_marker> functional_iota_slope_idx_cs_marker; // Role: functional transformation affine transformation slope index current slope marker.
  int<lower=0> n_const_cs_marker;           // constants used by marker CS bytecode
  vector[n_const_cs_marker] const_data_cs_marker; // constants used by marker CS bytecode

  /* Spline specifications */
  int<lower=0> n_knots_cv;        // knot count for total CV spline
  vector[n_knots_cv] knots_cv;    // knot locations for total CV spline
  int<lower=0> n_coeff_cv;        // coefficient count for total CV spline
  vector[n_coeff_cv] coeff_cv;    // coefficients for total CV spline
  int<lower=1, upper=5> spline_degree_cv; // spline degree for total CV
  int<lower=0, upper=1> estimate_spline_cv; // 1 if total CV spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cv;  // smoothness penalty strength for total CV spline
  int<lower=0> n_free_spline_cv;   // free monotone increments for total CV spline

  int<lower=0> n_knots_cs;        // knot count for total CS spline
  vector[n_knots_cs] knots_cs;    // knot locations for total CS spline
  int<lower=0> n_coeff_cs;        // coefficient count for total CS spline
  vector[n_coeff_cs] coeff_cs;    // coefficients for total CS spline
  int<lower=1, upper=5> spline_degree_cs; // spline degree for total CS
  int<lower=0, upper=1> estimate_spline_cs; // 1 if total CS spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cs;  // smoothness penalty strength for total CS spline
  int<lower=0> n_free_spline_cs;   // free monotone increments for total CS spline

  int<lower=0> n_knots_corr;      // knot count for corr spline
  vector[n_knots_corr] knots_corr; // knot locations for corr spline
  int<lower=0> n_coeff_corr;      // coefficient count for corr spline
  matrix[M_corr_tf, n_coeff_corr] coeff_corr; // coefficients for corr spline by component
  int<lower=1, upper=5> spline_degree_corr; // spline degree for corr
  int<lower=0, upper=1> estimate_spline_corr; // 1 if corr spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_corr; // smoothness penalty strength for corr spline
  int<lower=0> n_free_spline_corr;  // free monotone increments for corr spline

  int<lower=0> n_knots_vcov;      // knot count for vcov spline
  vector[n_knots_vcov] knots_vcov; // knot locations for vcov spline
  int<lower=0> n_coeff_vcov;      // coefficient count for vcov spline
  matrix[M_vcov_tf, n_coeff_vcov] coeff_vcov; // coefficients for vcov spline by component
  int<lower=1, upper=5> spline_degree_vcov; // spline degree for vcov
  int<lower=0, upper=1> estimate_spline_vcov; // 1 if vcov spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_vcov; // smoothness penalty strength for vcov spline
  int<lower=0> n_free_spline_vcov;  // free monotone increments for vcov spline

  int<lower=0> n_knots_cv_mean;   // knot count for mean CV spline
  vector[n_knots_cv_mean] knots_cv_mean; // knot locations for mean CV spline
  int<lower=0> n_coeff_cv_mean;   // coefficient count for mean CV spline
  vector[n_coeff_cv_mean] coeff_cv_mean; // coefficients for mean CV spline
  int<lower=1, upper=5> spline_degree_cv_mean; // spline degree for mean CV
  int<lower=0, upper=1> estimate_spline_cv_mean; // 1 if mean CV spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cv_mean; // smoothness penalty strength for mean CV spline
  int<lower=0> n_free_spline_cv_mean;  // free monotone increments for mean CV spline

  int<lower=0> n_knots_cv_marker; // knot count for marker CV spline
  vector[n_knots_cv_marker] knots_cv_marker; // knot locations for marker CV spline
  int<lower=0> n_coeff_cv_marker; // coefficient count for marker CV spline
  vector[n_coeff_cv_marker] coeff_cv_marker; // coefficients for marker CV spline
  int<lower=1, upper=5> spline_degree_cv_marker; // spline degree for marker CV
  int<lower=0, upper=1> estimate_spline_cv_marker; // 1 if marker CV spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cv_marker; // smoothness penalty strength for marker CV spline
  int<lower=0> n_free_spline_cv_marker;  // free monotone increments for marker CV spline

  int<lower=0> n_knots_cs_mean;   // knot count for mean CS spline
  vector[n_knots_cs_mean] knots_cs_mean; // knot locations for mean CS spline
  int<lower=0> n_coeff_cs_mean;   // coefficient count for mean CS spline
  vector[n_coeff_cs_mean] coeff_cs_mean; // coefficients for mean CS spline
  int<lower=1, upper=5> spline_degree_cs_mean; // spline degree for mean CS
  int<lower=0, upper=1> estimate_spline_cs_mean; // 1 if mean CS spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cs_mean; // smoothness penalty strength for mean CS spline
  int<lower=0> n_free_spline_cs_mean;  // free monotone increments for mean CS spline

  int<lower=0> n_knots_cs_marker; // knot count for marker CS spline
  vector[n_knots_cs_marker] knots_cs_marker; // knot locations for marker CS spline
  int<lower=0> n_coeff_cs_marker; // coefficient count for marker CS spline
  vector[n_coeff_cs_marker] coeff_cs_marker; // coefficients for marker CS spline
  int<lower=1, upper=5> spline_degree_cs_marker; // spline degree for marker CS
  int<lower=0, upper=1> estimate_spline_cs_marker; // 1 if marker CS spline coeffs are estimated in Stan
  real<lower=0> lambda_spline_cs_marker; // smoothness penalty strength for marker CS spline
  int<lower=0> n_free_spline_cs_marker;  // free monotone increments for marker CS spline


