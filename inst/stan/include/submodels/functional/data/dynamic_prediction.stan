/**
 * @file functional/data/dynamic_prediction.stan
 * @brief Declare posterior functional transformations for association channels.
 *
 * @details
 * Affine coefficients, transformation bytecode, constants, spline knots and fitted spline coefficients reproduce exactly the transformations used during fitting for each retained posterior draw.
 */
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

  array[n_draws] vector[estimate_iota_intercept_cv] iota_intercept_cv; // Role: affine transformation intercept current value.
  array[n_draws] vector[estimate_iota_slope_cv] iota_slope_cv; // Role: affine transformation slope current value.
  array[n_draws] vector[estimate_iota_intercept_cs] iota_intercept_cs; // Role: affine transformation intercept current slope.
  array[n_draws] vector[estimate_iota_slope_cs] iota_slope_cs; // Role: affine transformation slope current slope.
  array[n_draws] vector[M_corr_tf * estimate_iota_intercept_corr] iota_intercept_corr; // Role: affine transformation intercept correlation.
  array[n_draws] vector[M_corr_tf * estimate_iota_slope_corr] iota_slope_corr; // Role: affine transformation slope correlation.
  array[n_draws] vector[M_vcov_tf * estimate_iota_intercept_vcov] iota_intercept_vcov; // Role: affine transformation intercept covariance.
  array[n_draws] vector[M_vcov_tf * estimate_iota_slope_vcov] iota_slope_vcov; // Role: affine transformation slope covariance.
  array[n_draws] vector[estimate_iota_intercept_cv_mean] iota_intercept_cv_mean; // Role: affine transformation intercept current value mean.
  array[n_draws] vector[estimate_iota_slope_cv_mean] iota_slope_cv_mean; // Role: affine transformation slope current value mean.
  array[n_draws] vector[estimate_iota_intercept_cv_marker] iota_intercept_cv_marker; // Role: affine transformation intercept current value marker.
  array[n_draws] vector[estimate_iota_slope_cv_marker] iota_slope_cv_marker; // Role: affine transformation slope current value marker.
  array[n_draws] vector[estimate_iota_intercept_cs_mean] iota_intercept_cs_mean; // Role: affine transformation intercept current slope mean.
  array[n_draws] vector[estimate_iota_slope_cs_mean] iota_slope_cs_mean; // Role: affine transformation slope current slope mean.
  array[n_draws] vector[estimate_iota_intercept_cs_marker] iota_intercept_cs_marker; // Role: affine transformation intercept current slope marker.
  array[n_draws] vector[estimate_iota_slope_cs_marker] iota_slope_cs_marker; // Role: affine transformation slope current slope marker.

  int<lower=0, upper=7> tf_mode_cv_tot;        // transform mode: total CV
  int<lower=0, upper=7> tf_mode_cs_tot;        // transform mode: total CS
  int<lower=0, upper=7> tf_mode_corr;          // transform mode: corr
  int<lower=0, upper=7> tf_mode_vcov;          // transform mode: vcov
  int<lower=0, upper=7> tf_mode_cv_mean;       // transform mode: mean CV
  int<lower=0, upper=7> tf_mode_cv_marker;     // transform mode: marker CV
  int<lower=0, upper=7> tf_mode_cs_mean;       // transform mode: mean CS
  int<lower=0, upper=7> tf_mode_cs_marker;     // transform mode: marker CS

  int<lower=0> n_functional_ops_cv;            // op count for total CV
  array[n_functional_ops_cv] int<lower=0, upper=27> functional_ops_cv; // bytecode stream
  array[n_functional_ops_cv] int<lower=0, upper=estimate_iota_intercept_cv> functional_iota_intercept_idx_cv; // Role: functional transformation affine transformation intercept index current value.
  array[n_functional_ops_cv] int<lower=0, upper=estimate_iota_slope_cv> functional_iota_slope_idx_cv; // Role: functional transformation affine transformation slope index current value.
  int<lower=0> n_const_cv;                     // constants for total CV bytecode
  vector[n_const_cv] const_data_cv;            // constants for total CV bytecode

  int<lower=0> n_functional_ops_cs;            // op count for total CS
  array[n_functional_ops_cs] int<lower=0, upper=27> functional_ops_cs; // bytecode stream
  array[n_functional_ops_cs] int<lower=0, upper=estimate_iota_intercept_cs> functional_iota_intercept_idx_cs; // Role: functional transformation affine transformation intercept index current slope.
  array[n_functional_ops_cs] int<lower=0, upper=estimate_iota_slope_cs> functional_iota_slope_idx_cs; // Role: functional transformation affine transformation slope index current slope.
  int<lower=0> n_const_cs;                     // constants for total CS bytecode
  vector[n_const_cs] const_data_cs;            // constants for total CS bytecode

  int<lower=0> n_functional_ops_corr;          // op count for corr
  array[n_functional_ops_corr] int<lower=0, upper=27> functional_ops_corr; // bytecode stream
  array[n_functional_ops_corr] int<lower=0, upper=estimate_iota_intercept_corr> functional_iota_intercept_idx_corr; // Role: functional transformation affine transformation intercept index correlation.
  array[n_functional_ops_corr] int<lower=0, upper=estimate_iota_slope_corr> functional_iota_slope_idx_corr; // Role: functional transformation affine transformation slope index correlation.
  int<lower=0> n_const_corr;                   // constants for corr bytecode
  vector[n_const_corr] const_data_corr;        // constants for corr bytecode

  int<lower=0> n_functional_ops_vcov;          // op count for vcov
  array[n_functional_ops_vcov] int<lower=0, upper=27> functional_ops_vcov; // bytecode stream
  array[n_functional_ops_vcov] int<lower=0, upper=estimate_iota_intercept_vcov> functional_iota_intercept_idx_vcov; // Role: functional transformation affine transformation intercept index covariance.
  array[n_functional_ops_vcov] int<lower=0, upper=estimate_iota_slope_vcov> functional_iota_slope_idx_vcov; // Role: functional transformation affine transformation slope index covariance.
  int<lower=0> n_const_vcov;                   // constants for vcov bytecode
  vector[n_const_vcov] const_data_vcov;        // constants for vcov bytecode

  int<lower=0> n_knots_cv;                     // knots for total CV spline
  vector[n_knots_cv] knots_cv;                 // knot locations
  int<lower=0> n_coeff_cv;                     // coeff count for total CV spline
  array[n_draws] vector[n_coeff_cv] coeff_cv;  // coefficients for total CV spline by draw
  int<lower=1, upper=5> spline_degree_cv;      // spline degree for total CV

  int<lower=0> n_knots_cs;                     // knots for total CS spline
  vector[n_knots_cs] knots_cs;                 // knot locations
  int<lower=0> n_coeff_cs;                     // coeff count for total CS spline
  array[n_draws] vector[n_coeff_cs] coeff_cs;  // coefficients for total CS spline by draw
  int<lower=1, upper=5> spline_degree_cs;      // spline degree for total CS

  int<lower=0> n_knots_corr;                   // knots for corr spline
  vector[n_knots_corr] knots_corr;             // knot locations
  int<lower=0> n_coeff_corr;                   // coeff count for corr spline
  array[n_draws] matrix[M_corr_tf, n_coeff_corr] coeff_corr; // coefficients for corr spline by draw and component
  int<lower=1, upper=5> spline_degree_corr;    // spline degree for corr

  int<lower=0> n_knots_vcov;                   // knots for vcov spline
  vector[n_knots_vcov] knots_vcov;             // knot locations
  int<lower=0> n_coeff_vcov;                   // coeff count for vcov spline
  array[n_draws] matrix[M_vcov_tf, n_coeff_vcov] coeff_vcov; // coefficients for vcov spline by draw and component
  int<lower=1, upper=5> spline_degree_vcov;    // spline degree for vcov

  int<lower=0> n_functional_ops_cv_mean;        // op count for mean CV
  array[n_functional_ops_cv_mean] int<lower=0, upper=27> functional_ops_cv_mean; // bytecode stream
  array[n_functional_ops_cv_mean] int<lower=0, upper=estimate_iota_intercept_cv_mean> functional_iota_intercept_idx_cv_mean; // Role: functional transformation affine transformation intercept index current value mean.
  array[n_functional_ops_cv_mean] int<lower=0, upper=estimate_iota_slope_cv_mean> functional_iota_slope_idx_cv_mean; // Role: functional transformation affine transformation slope index current value mean.
  int<lower=0> n_const_cv_mean;                 // constants for mean CV bytecode
  vector[n_const_cv_mean] const_data_cv_mean;   // constants for mean CV bytecode

  int<lower=0> n_functional_ops_cv_marker;      // op count for marker CV
  array[n_functional_ops_cv_marker] int<lower=0, upper=27> functional_ops_cv_marker; // bytecode stream
  array[n_functional_ops_cv_marker] int<lower=0, upper=estimate_iota_intercept_cv_marker> functional_iota_intercept_idx_cv_marker; // Role: functional transformation affine transformation intercept index current value marker.
  array[n_functional_ops_cv_marker] int<lower=0, upper=estimate_iota_slope_cv_marker> functional_iota_slope_idx_cv_marker; // Role: functional transformation affine transformation slope index current value marker.
  int<lower=0> n_const_cv_marker;               // constants for marker CV bytecode
  vector[n_const_cv_marker] const_data_cv_marker; // constants for marker CV bytecode

  int<lower=0> n_functional_ops_cs_mean;        // op count for mean CS
  array[n_functional_ops_cs_mean] int<lower=0, upper=27> functional_ops_cs_mean; // bytecode stream
  array[n_functional_ops_cs_mean] int<lower=0, upper=estimate_iota_intercept_cs_mean> functional_iota_intercept_idx_cs_mean; // Role: functional transformation affine transformation intercept index current slope mean.
  array[n_functional_ops_cs_mean] int<lower=0, upper=estimate_iota_slope_cs_mean> functional_iota_slope_idx_cs_mean; // Role: functional transformation affine transformation slope index current slope mean.
  int<lower=0> n_const_cs_mean;                 // constants for mean CS bytecode
  vector[n_const_cs_mean] const_data_cs_mean;   // constants for mean CS bytecode

  int<lower=0> n_functional_ops_cs_marker;      // op count for marker CS
  array[n_functional_ops_cs_marker] int<lower=0, upper=27> functional_ops_cs_marker; // bytecode stream
  array[n_functional_ops_cs_marker] int<lower=0, upper=estimate_iota_intercept_cs_marker> functional_iota_intercept_idx_cs_marker; // Role: functional transformation affine transformation intercept index current slope marker.
  array[n_functional_ops_cs_marker] int<lower=0, upper=estimate_iota_slope_cs_marker> functional_iota_slope_idx_cs_marker; // Role: functional transformation affine transformation slope index current slope marker.
  int<lower=0> n_const_cs_marker;               // constants for marker CS bytecode
  vector[n_const_cs_marker] const_data_cs_marker; // constants for marker CS bytecode

  int<lower=0> n_knots_cv_mean;                 // knots for mean CV spline
  vector[n_knots_cv_mean] knots_cv_mean;         // knot locations
  int<lower=0> n_coeff_cv_mean;                  // coeff count for mean CV spline
  array[n_draws] vector[n_coeff_cv_mean] coeff_cv_mean; // coefficients for mean CV spline by draw
  int<lower=1, upper=5> spline_degree_cv_mean;   // spline degree for mean CV

  int<lower=0> n_knots_cv_marker;               // knots for marker CV spline
  vector[n_knots_cv_marker] knots_cv_marker;    // knot locations
  int<lower=0> n_coeff_cv_marker;               // coeff count for marker CV spline
  array[n_draws] vector[n_coeff_cv_marker] coeff_cv_marker; // coefficients for marker CV spline by draw
  int<lower=1, upper=5> spline_degree_cv_marker; // spline degree for marker CV

  int<lower=0> n_knots_cs_mean;                 // knots for mean CS spline
  vector[n_knots_cs_mean] knots_cs_mean;         // knot locations
  int<lower=0> n_coeff_cs_mean;                  // coeff count for mean CS spline
  array[n_draws] vector[n_coeff_cs_mean] coeff_cs_mean; // coefficients for mean CS spline by draw
  int<lower=1, upper=5> spline_degree_cs_mean;   // spline degree for mean CS

  int<lower=0> n_knots_cs_marker;               // knots for marker CS spline
  vector[n_knots_cs_marker] knots_cs_marker;    // knot locations
  int<lower=0> n_coeff_cs_marker;               // coeff count for marker CS spline
  array[n_draws] vector[n_coeff_cs_marker] coeff_cs_marker; // coefficients for marker CS spline by draw
  int<lower=1, upper=5> spline_degree_cs_marker; // spline degree for marker CS


