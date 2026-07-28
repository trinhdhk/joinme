  int<lower=1> n_draws; // number of posterior draws to process

  /* Observed longitudinal history */
  int<lower=1> n_obs_long;                          // number of observed rows
  array[n_obs_long] int<lower=1> idx_marker_obs;    // marker index per observed row
  int<lower=1> n_marker_types;                      // total number of markers
  vector[n_marker_types] marker_weights_cv_total;            // reference weights for cv_total
  vector[n_marker_types] marker_weights_cs_total;            // reference weights for cs_total
  vector[n_marker_types] marker_weights_cv_marker;           // reference weights for cv_marker
  vector[n_marker_types] marker_weights_cs_marker;           // reference weights for cs_marker
  matrix[n_draws, n_marker_types] marker_weights_draws_cv_total; // per-draw weights for cv_total
  matrix[n_draws, n_marker_types] marker_weights_draws_cs_total; // per-draw weights for cs_total
  matrix[n_draws, n_marker_types] marker_weights_draws_cv_marker; // per-draw weights for cv_marker
  matrix[n_draws, n_marker_types] marker_weights_draws_cs_marker; // per-draw weights for cs_marker
  vector[n_obs_long] y_real;                        // continuous outcomes
  array[n_obs_long] int y_int;                      // discrete outcomes
  array[n_obs_long] int<lower=0> trials_obs;        // binomial trials (1 otherwise)

  /* Design matrix dimensions */
  int<lower=1> n_fixed_effects;        // fixed-effect columns
  int<lower=1> n_random_id;            // id-level RE columns
  int<lower=0> n_random_marker;        // marker-only RE columns
  int<lower=0> n_random_marker_id;     // marker-by-id RE columns

  /* Observed design matrices */
  matrix[n_obs_long, n_fixed_effects] mat_fixed_obs;        // fixed effects
  matrix[n_obs_long, n_random_id] mat_id_obs;               // id RE design
  matrix[n_obs_long, n_random_marker] mat_marker_obs;       // marker RE design
  matrix[n_obs_long, n_random_marker_id] mat_marker_id_obs; // marker-by-id RE
  vector<lower=0>[n_random_marker_id] marker_id_row_scale;  // row scaling from original-time covariance regression to scaled-time marker-id coefficients

  /* Prediction grid identifiers */
  int<lower=1> n_obs_pred;                       // rows in prediction grid
  array[n_obs_pred] int<lower=1> idx_marker_pred; // marker index per pred row

  /* Distributional regression designs */
  int<lower=0> P_sigma;                          // sigma fixed-effect columns
  matrix[n_obs_long, P_sigma] X_sigma_obs;       // sigma design for observed
  matrix[n_obs_pred, P_sigma] X_sigma_pred;      // sigma design for prediction
  int<lower=0> P_nu;                             // nu fixed-effect columns
  matrix[n_obs_long, P_nu] X_nu_obs;             // nu design for observed
  matrix[n_obs_pred, P_nu] X_nu_pred;            // nu design for prediction
  int<lower=0> P_phi;                            // phi fixed-effect columns
  matrix[n_obs_long, P_phi] X_phi_obs;           // phi design for observed
  matrix[n_obs_pred, P_phi] X_phi_pred;          // phi design for prediction
  int<lower=0> P_alpha;                          // alpha fixed-effect columns
  matrix[n_obs_long, P_alpha] X_alpha_obs;       // alpha design for observed
  matrix[n_obs_pred, P_alpha] X_alpha_pred;      // alpha design for prediction
  int<lower=0> P_kappa;                       // kappa fixed-effect columns
  matrix[n_obs_long, P_kappa] X_kappa_obs; // kappa design observed
  matrix[n_obs_pred, P_kappa] X_kappa_pred; // kappa design pred
  int<lower=0> P_tau;                        // tau fixed-effect columns
  matrix[n_obs_long, P_tau] X_tau_obs;   // tau design observed
  matrix[n_obs_pred, P_tau] X_tau_pred;  // tau design pred

  /* Covariance regression covariates */
  int<lower=0> n_cov_vcov;                        // covariate count for covariance regression
  vector[n_cov_vcov] vec_cov_vcov;                // covariate vector (subject-level)

  /* Hazard covariates */
  int<lower=0> n_cov_hazard;                      // hazard covariate count
  vector[n_cov_hazard] vec_cov_hazard;            // hazard covariate vector

  /* Conditioning information */
  int<lower=1> n_basehaz_basis;                   // baseline hazard basis count
  real<lower=0> time_condition;                   // conditioning time T_cond
  int<lower=1> n_gk;                              // quadrature node count (nodes/weights hardcoded in Stan)
  matrix[n_gk, n_basehaz_basis] mat_basis_gk_cond;  // basis at quadrature nodes (0..T_cond)
  matrix[n_gk, n_fixed_effects] mat_fixed_gk_cond;  // fixed effects at quadrature nodes
  matrix[n_gk, n_random_id] mat_id_gk_cond;         // id RE design at quadrature nodes
  matrix[n_gk, n_random_marker] mat_marker_gk_cond; // marker RE design at quadrature nodes
  matrix[n_gk, n_random_marker_id] mat_marker_id_gk_cond; // marker-id RE design

  /* Forward-step design matrices for slope (finite difference) */
  matrix[n_gk, n_fixed_effects] mat_fixed_gk_cond_fwd; // Role: matrix fixed effect Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_id] mat_id_gk_cond_fwd; // Role: matrix individual Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_marker] mat_marker_gk_cond_fwd; // Role: matrix marker Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_marker_id] mat_marker_id_gk_cond_fwd; // Role: matrix marker individual Gauss-Kronrod conditioning forward difference.

  real<lower=1e-6> eps_finite_diff;               // finite-difference step

  /* Longitudinal prediction grid */
  matrix[n_obs_pred, n_fixed_effects] mat_fixed_pred;        // fixed design
  matrix[n_obs_pred, n_random_id] mat_id_pred;               // id RE design
  matrix[n_obs_pred, n_random_marker] mat_marker_pred;       // marker RE design
  matrix[n_obs_pred, n_random_marker_id] mat_marker_id_pred; // marker-id RE design
  array[n_obs_pred] int<lower=0> trials_pred;                // binomial trials

  /* Survival prediction grid */
  int<lower=1> n_times_surv;                     // number of survival time points
  vector[n_times_surv] vec_time_surv;            // times for S(t | T_cond)
  array[n_times_surv] matrix[n_gk, n_basehaz_basis] mat_basis_gk_surv; // Role: matrix basis Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv; // Role: matrix fixed effect Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv; // Role: matrix individual Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv; // Role: matrix marker Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv; // Role: matrix marker individual Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv_fwd; // Role: matrix fixed effect Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv_fwd; // Role: matrix individual Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv_fwd; // Role: matrix marker Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv_fwd; // Role: matrix marker individual Gauss-Kronrod survival forward difference.

  /* Model configuration flags and draws */
  array[n_marker_types] int<lower=1, upper=11> family_long; // family codes
  array[n_marker_types] int<lower=0, upper=5> link_long; // 0 custom VM; 1 identity, 2 log, 3 logit, 4 probit, 5 exp
  int<lower=1> max_inv_link_ops; // Role: maximum inverse link operations.
  array[n_marker_types] int<lower=0> inv_link_n_ops; // Role: inverse link number of operations.
  array[n_marker_types, max_inv_link_ops] int<lower=0, upper=27> inv_link_ops; // Role: inverse link operations.
  int<lower=1> max_inv_link_const; // Role: maximum inverse link constants.
  array[n_marker_types] int<lower=0> inv_link_n_const; // Role: inverse link number of constants.
  matrix[n_marker_types, max_inv_link_const] inv_link_const; // Role: inverse link constants.
  int<lower=0> n_family_sigma; // Role: number of family scale.
  array[n_marker_types] int<lower=0, upper=n_family_sigma> marker_to_sigma_family; // Role: marker to scale family.
  int<lower=0> n_family_nu; // Role: number of family degrees of freedom.
  array[n_marker_types] int<lower=0, upper=n_family_nu> marker_to_nu_family; // Role: marker to degrees of freedom family.
  int<lower=0> n_family_phi; // Role: number of family precision.
  array[n_marker_types] int<lower=0, upper=n_family_phi> marker_to_phi_family; // Role: marker to precision family.
  int<lower=0> n_family_alpha; // Role: number of family shape.
  array[n_marker_types] int<lower=0, upper=n_family_alpha> marker_to_alpha_family; // Role: marker to shape family.
  int<lower=0> n_family_kappa; // Role: number of family sample size.
  array[n_marker_types] int<lower=0, upper=n_family_kappa> marker_to_kappa_family; // Role: marker to sample size family.
  int<lower=0> n_family_tau; // Role: number of family quantile.
  array[n_marker_types] int<lower=0, upper=n_family_tau> marker_to_tau_family; // Role: marker to quantile family.
  // int<lower=0, upper=1> flag_resid_dim;                      // sigma dimension flag
  int<lower=0, upper=1> vcov_diag_link;                      // 0 softplus, 1 exp
  array[n_marker_types] int<lower=0, upper=1> use_tau_fixed; // marker-specific fixed-tau flags
  vector<lower=0, upper=1>[n_marker_types] tau_fixed;        // marker-specific fixed quantiles

  array[n_draws] vector[n_fixed_effects] beta_fixed;          // fixed effects per draw
  array[n_draws] vector[n_random_id] tau_id;                  // id RE scales per draw
  array[n_draws] matrix[n_random_id, n_random_id] Lcorr_id;   // id RE correlations
  array[n_draws] vector[n_random_marker] tau_marker;          // marker RE scales
  array[n_draws] matrix[n_random_marker, n_random_marker] Lcorr_marker; // marker RE corr
  array[n_draws] matrix[n_random_marker_id, n_random_marker] B_cross; // cross-corr map

  int<lower=0> num_unique_cov_entries;         // unique elements in L_i
  array[n_draws] vector[num_unique_cov_entries] alpha_vcov_reg; // intercepts
  array[n_draws] vector[num_unique_cov_entries * n_cov_vcov] beta_vcov_reg_flat; // flattened slopes
  array[n_draws] vector[num_unique_cov_entries] lambda_vcov_reg; // effective loadings

  int<lower=1> K_event;                        // number of competing risks
  array[n_draws, K_event] vector[n_basehaz_basis] bs_gamma_c; // spline coeffs
  array[n_draws, K_event] vector[n_cov_hazard] gamma_hazard;  // hazard covariates

  array[n_draws] vector[P_sigma] beta_sigma;    // sigma regression betas
  array[n_draws] vector[P_nu] beta_nu;          // nu regression betas
  array[n_draws] vector[P_phi] beta_phi;        // phi regression betas
  array[n_draws] vector[P_alpha] beta_alpha;    // alpha regression betas
  array[n_draws] vector[P_kappa] beta_kappa; // kappa regression betas
  array[n_draws] vector[P_tau] beta_tau;   // tau regression betas

  array[n_draws] vector[n_family_sigma] sigma_family;          // family-shared sigma
  array[n_draws] vector[n_family_nu] nu_family;                // family-shared nu
  array[n_draws] vector[n_family_phi] phi_family;              // family-shared phi
  array[n_draws] vector[n_family_alpha] alpha_family;          // family-shared alpha
  array[n_draws] vector[n_family_kappa] kappa_family;    // family-shared kappa
  array[n_draws] vector[n_family_tau] tau_family;      // family-shared tau

  /* Ordinal (cumulative logit) cutpoints */
  int<lower=2> K_ord;                           // number of ordinal categories
  array[n_draws] vector[K_ord - 1] cutpoints_ord; // cutpoints per draw

  array[n_draws] real coeff_assoc_cv_total;     // association coeffs: total CV
  array[n_draws] real coeff_assoc_cs_total;     // association coeffs: total CS
  array[n_draws] real coeff_assoc_cv_mean;      // association coeffs: mean CV
  array[n_draws] real coeff_assoc_cs_mean;      // association coeffs: mean CS
  array[n_draws] real coeff_assoc_cv_marker;    // association coeffs: marker CV
  array[n_draws] real coeff_assoc_cs_marker;    // association coeffs: marker CS
  array[n_draws] vector[(n_random_marker_id * (n_random_marker_id - 1)) %/% 2] coeff_assoc_corr; // corr coeffs (off-diagonal correlations)
  array[n_draws] vector[num_unique_cov_entries] coeff_assoc_vcov; // vcov coeffs for lower-triangular L entries

  int<lower=0, upper=1> flag_assoc_cv_total;    // include total CV association
  int<lower=0, upper=1> flag_assoc_cv_mean;     // include mean CV association
  int<lower=0, upper=1> flag_assoc_cv_marker;   // include marker CV association
  int<lower=0, upper=1> flag_assoc_cs_total;    // include total CS association
  int<lower=0, upper=1> flag_assoc_cs_mean;     // include mean CS association
  int<lower=0, upper=1> flag_assoc_cs_marker;   // include marker CS association
  int<lower=0, upper=1> flag_assoc_corr;        // include corr association
  int<lower=0, upper=1> flag_assoc_vcov;        // include vcov association
  int<lower=0> M_corr_tf;                       // number of corr transform components
  int<lower=0> M_vcov_tf;                       // number of vcov transform components
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
