/**
 * @file longitudinal/data/dynamic_prediction.stan
 * @brief Declare observed histories, prediction designs and fitted longitudinal draws.
 *
 * @details
 * The fragment describes marker outcomes, nested random-effect designs, distributional regressions, covariance regression and the longitudinal prediction grid.  It restores fitted longitudinal and distributional parameters draw by draw.
 */
  /* Observed longitudinal history */
  int<lower=1> n_obs_long;                          // number of observed rows
  array[n_obs_long] int<lower=1> idx_marker_obs;    // marker index per observed row
  int<lower=1> n_marker_types;                      // total number of markers
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

  /* Independent covariance regression covariates */
  int<lower=0> n_cov_vcov_sd; // covariate count for the standard-deviation regression
  vector[n_cov_vcov_sd] vec_cov_vcov_sd; // subject covariates entering standard-deviation predictors
  int<lower=0> n_cov_vcov_corr; // covariate count for the off-diagonal correlation regression
  vector[n_cov_vcov_corr] vec_cov_vcov_corr; // subject covariates entering partial-correlation predictors

  /* Longitudinal prediction grid */
  matrix[n_obs_pred, n_fixed_effects] mat_fixed_pred;        // fixed design
  matrix[n_obs_pred, n_random_id] mat_id_pred;               // id RE design
  matrix[n_obs_pred, n_random_marker] mat_marker_pred;       // marker RE design
  matrix[n_obs_pred, n_random_marker_id] mat_marker_id_pred; // marker-id RE design
  array[n_obs_pred] int<lower=0> trials_pred;                // binomial trials

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
  array[n_draws] vector[n_random_marker_id * n_cov_vcov_sd] beta_vcov_sd_flat; // row-major standard-deviation slopes
  array[n_draws] vector[((n_random_marker_id * (n_random_marker_id - 1)) %/% 2) * n_cov_vcov_corr] beta_vcov_corr_flat; // row-major off-diagonal correlation slopes
  array[n_draws] vector[num_unique_cov_entries] lambda_vcov_reg; // one nonnegative scalar latent-heterogeneity loading per packed covariance coordinate and retained fit draw

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

  /**
   * Structural flags and packed covariance indices retained from fitting.
   * These values determine how nested longitudinal random effects and their
   * subject-specific covariance factors are reconstructed for a new subject.
   */
  int flag_indep_id_re; // whether shared-subject effects are mutually independent
  int flag_indep_marker_re; // whether shared-marker effects are mutually independent
  int flag_indep_idmarker_cov; // whether subject-by-marker covariance is diagonal
  int flag_allow_marker_crosscorr; // whether marker effects may correlate across markers
  array[num_unique_cov_entries] int idx_row_cov; // row of each packed covariance entry
  array[num_unique_cov_entries] int idx_col_cov; // column of each packed covariance entry
