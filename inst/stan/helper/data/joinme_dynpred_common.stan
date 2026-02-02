  int<lower=1> n_draws; // number of posterior draws to process

  // Observed longitudinal history
  int<lower=1> n_obs_long;
  array[n_obs_long] int<lower=1> idx_marker_obs;
  int<lower=1> n_marker_types;
  vector[n_obs_long] y_real;
  array[n_obs_long] int y_int;
  array[n_obs_long] int<lower=0> trials_obs;

  // Design matrix dimensions
  int<lower=1> n_fixed_effects;
  int<lower=1> n_random_id;
  int<lower=0> n_random_marker;
  int<lower=0> n_random_marker_id;

  // Observed design matrices
  matrix[n_obs_long, n_fixed_effects] mat_fixed_obs;
  matrix[n_obs_long, n_random_id] mat_id_obs;
  matrix[n_obs_long, n_random_marker] mat_marker_obs;
  matrix[n_obs_long, n_random_marker_id] mat_marker_id_obs;

  // Prediction grid identifiers
  int<lower=1> n_obs_pred;
  array[n_obs_pred] int<lower=1> idx_marker_pred;

  // Distributional regression designs
  int<lower=0> P_sigma;
  matrix[n_obs_long, P_sigma] X_sigma_obs;
  matrix[n_obs_pred, P_sigma] X_sigma_pred;
  int<lower=0> P_nu;
  matrix[n_obs_long, P_nu] X_nu_obs;
  matrix[n_obs_pred, P_nu] X_nu_pred;
  int<lower=0> P_phi;
  matrix[n_obs_long, P_phi] X_phi_obs;
  matrix[n_obs_pred, P_phi] X_phi_pred;
  int<lower=0> P_alpha;
  matrix[n_obs_long, P_alpha] X_alpha_obs;
  matrix[n_obs_pred, P_alpha] X_alpha_pred;
  int<lower=0> P_phi_beta;
  matrix[n_obs_long, P_phi_beta] X_phi_beta_obs;
  matrix[n_obs_pred, P_phi_beta] X_phi_beta_pred;
  int<lower=0> P_tau_sde;
  matrix[n_obs_long, P_tau_sde] X_tau_sde_obs;
  matrix[n_obs_pred, P_tau_sde] X_tau_sde_pred;

  // Covariance regression covariates
  int<lower=1> n_cov_vcov;
  vector[n_cov_vcov] vec_cov_vcov;

  // Hazard covariates
  int<lower=0> n_cov_hazard;
  vector[n_cov_hazard] vec_cov_hazard;

  // Conditioning information
  int<lower=1> n_basehaz_basis;
  real<lower=0> time_condition;
  matrix[15, n_basehaz_basis] mat_basis_gk_cond;
  matrix[15, n_fixed_effects] mat_fixed_gk_cond;
  matrix[15, n_random_id] mat_id_gk_cond;
  matrix[15, n_random_marker] mat_marker_gk_cond;
  matrix[15, n_random_marker_id] mat_marker_id_gk_cond;

  matrix[15, n_fixed_effects] mat_fixed_gk_cond_fwd;
  matrix[15, n_random_id] mat_id_gk_cond_fwd;
  matrix[15, n_random_marker] mat_marker_gk_cond_fwd;
  matrix[15, n_random_marker_id] mat_marker_id_gk_cond_fwd;

  real<lower=1e-6> eps_finite_diff;

  // Longitudinal prediction grid
  matrix[n_obs_pred, n_fixed_effects] mat_fixed_pred;
  matrix[n_obs_pred, n_random_id] mat_id_pred;
  matrix[n_obs_pred, n_random_marker] mat_marker_pred;
  matrix[n_obs_pred, n_random_marker_id] mat_marker_id_pred;
  array[n_obs_pred] int<lower=0> trials_pred;

  // Survival prediction grid
  int<lower=1> n_times_surv;
  vector[n_times_surv] vec_time_surv;
  array[n_times_surv] matrix[15, n_basehaz_basis] mat_basis_gk_surv;
  array[n_times_surv] matrix[15, n_fixed_effects] mat_fixed_gk_surv;
  array[n_times_surv] matrix[15, n_random_id] mat_id_gk_surv;
  array[n_times_surv] matrix[15, n_random_marker] mat_marker_gk_surv;
  array[n_times_surv] matrix[15, n_random_marker_id] mat_marker_id_gk_surv;
  array[n_times_surv] matrix[15, n_fixed_effects] mat_fixed_gk_surv_fwd;
  array[n_times_surv] matrix[15, n_random_id] mat_id_gk_surv_fwd;
  array[n_times_surv] matrix[15, n_random_marker] mat_marker_gk_surv_fwd;
  array[n_times_surv] matrix[15, n_random_marker_id] mat_marker_id_gk_surv_fwd;

  // Model configuration flags and draws
  array[n_marker_types] int<lower=1, upper=11> family_long;
  int<lower=0, upper=1> flag_resid_dim;

  array[n_draws] vector[n_fixed_effects] beta_fixed;
  array[n_draws] vector[n_random_id] tau_id;
  array[n_draws] matrix[n_random_id, n_random_id] Lcorr_id;
  array[n_draws] vector[n_random_marker] tau_marker;
  array[n_draws] matrix[n_random_marker, n_random_marker] Lcorr_marker;
  array[n_draws] vector[n_random_marker_id] tau_marker_id;
  array[n_draws] matrix[n_random_marker_id, n_random_marker_id] Lcorr_marker_id;
  array[n_draws] matrix[n_random_marker_id, n_random_marker] B_cross;

  int<lower=0> num_unique_cov_entries;
  array[n_draws] vector[num_unique_cov_entries] alpha_vcov_reg;
  array[n_draws] vector[num_unique_cov_entries * n_cov_vcov] beta_vcov_reg_flat;
  array[n_draws] real tau_vcov_reg;
  array[n_draws] vector[num_unique_cov_entries] lambda_vcov_reg;

  int<lower=1> K_event;
  array[n_draws, K_event] real log_h0_intercept;
  array[n_draws, K_event] vector[n_basehaz_basis] bs_gamma_c;
  array[n_draws, K_event] vector[n_cov_hazard] gamma_hazard;

  array[n_draws] vector[P_sigma] beta_sigma;
  array[n_draws] vector[P_nu] beta_nu;
  array[n_draws] vector[P_phi] beta_phi;
  array[n_draws] vector[P_alpha] beta_alpha;
  array[n_draws] vector[P_phi_beta] beta_phi_beta;
  array[n_draws] vector[P_tau_sde] beta_tau_sde;

  array[n_draws] real sigma_y_shared;
  array[n_draws] vector[n_marker_types] sigma_marker_specific;
  array[n_draws] vector[n_marker_types] nu_marker;
  array[n_draws] vector[n_marker_types] phi_nb_marker;
  array[n_draws] vector[n_marker_types] alpha_skew_marker;
  array[n_draws] vector[n_marker_types] phi_beta_marker;
  array[n_draws] vector[n_marker_types] tau_sde_marker;

  // Ordinal (cumulative logit) cutpoints
  int<lower=2> K_ord;
  array[n_draws] vector[K_ord - 1] cutpoints_ord;

  array[n_draws] real coeff_assoc_cv_total;
  array[n_draws] real coeff_assoc_cs_total;
  array[n_draws] real coeff_assoc_cv_mean;
  array[n_draws] real coeff_assoc_cs_mean;
  array[n_draws] real coeff_assoc_cv_marker;
  array[n_draws] real coeff_assoc_cs_marker;
  array[n_draws] vector[n_random_marker_id] coeff_assoc_vcov_var;

  int<lower=0, upper=1> flag_assoc_cv_total;
  int<lower=0, upper=1> flag_assoc_cv_mean;
  int<lower=0, upper=1> flag_assoc_cv_marker;
  int<lower=0, upper=1> flag_assoc_cs_total;
  int<lower=0, upper=1> flag_assoc_cs_mean;
  int<lower=0, upper=1> flag_assoc_cs_marker;
  int<lower=0, upper=1> flag_assoc_vcov;

  int<lower=0, upper=3> tf_mode_cv_tot;
  int<lower=0, upper=3> tf_mode_cs_tot;
  int<lower=0, upper=3> tf_mode_vcov;
  int<lower=0, upper=3> tf_mode_cv_mean;
  int<lower=0, upper=3> tf_mode_cv_marker;
  int<lower=0, upper=3> tf_mode_cs_mean;
  int<lower=0, upper=3> tf_mode_cs_marker;

  int<lower=0> n_functional_ops_cv;
  array[n_functional_ops_cv] int<lower=0, upper=23> functional_ops_cv;
  int<lower=0> n_const_cv;
  vector[n_const_cv] const_data_cv;

  int<lower=0> n_functional_ops_cs;
  array[n_functional_ops_cs] int<lower=0, upper=23> functional_ops_cs;
  int<lower=0> n_const_cs;
  vector[n_const_cs] const_data_cs;

  int<lower=0> n_functional_ops_vcov;
  array[n_functional_ops_vcov] int<lower=0, upper=23> functional_ops_vcov;
  int<lower=0> n_const_vcov;
  vector[n_const_vcov] const_data_vcov;

  int<lower=0> n_knots_cv;
  vector[n_knots_cv] knots_cv;
  int<lower=0> n_coeff_cv;
  vector[n_coeff_cv] coeff_cv;
  int<lower=1, upper=5> spline_degree_cv;

  int<lower=0> n_knots_cs;
  vector[n_knots_cs] knots_cs;
  int<lower=0> n_coeff_cs;
  vector[n_coeff_cs] coeff_cs;
  int<lower=1, upper=5> spline_degree_cs;

  int<lower=0> n_knots_vcov;
  vector[n_knots_vcov] knots_vcov;
  int<lower=0> n_coeff_vcov;
  vector[n_coeff_vcov] coeff_vcov;
  int<lower=1, upper=5> spline_degree_vcov;

  int<lower=0> n_functional_ops_cv_mean;
  array[n_functional_ops_cv_mean] int<lower=0, upper=23> functional_ops_cv_mean;
  int<lower=0> n_const_cv_mean;
  vector[n_const_cv_mean] const_data_cv_mean;

  int<lower=0> n_functional_ops_cv_marker;
  array[n_functional_ops_cv_marker] int<lower=0, upper=23> functional_ops_cv_marker;
  int<lower=0> n_const_cv_marker;
  vector[n_const_cv_marker] const_data_cv_marker;

  int<lower=0> n_functional_ops_cs_mean;
  array[n_functional_ops_cs_mean] int<lower=0, upper=23> functional_ops_cs_mean;
  int<lower=0> n_const_cs_mean;
  vector[n_const_cs_mean] const_data_cs_mean;

  int<lower=0> n_functional_ops_cs_marker;
  array[n_functional_ops_cs_marker] int<lower=0, upper=23> functional_ops_cs_marker;
  int<lower=0> n_const_cs_marker;
  vector[n_const_cs_marker] const_data_cs_marker;

  int<lower=0> n_knots_cv_mean;
  vector[n_knots_cv_mean] knots_cv_mean;
  int<lower=0> n_coeff_cv_mean;
  vector[n_coeff_cv_mean] coeff_cv_mean;
  int<lower=1, upper=5> spline_degree_cv_mean;

  int<lower=0> n_knots_cv_marker;
  vector[n_knots_cv_marker] knots_cv_marker;
  int<lower=0> n_coeff_cv_marker;
  vector[n_coeff_cv_marker] coeff_cv_marker;
  int<lower=1, upper=5> spline_degree_cv_marker;

  int<lower=0> n_knots_cs_mean;
  vector[n_knots_cs_mean] knots_cs_mean;
  int<lower=0> n_coeff_cs_mean;
  vector[n_coeff_cs_mean] coeff_cs_mean;
  int<lower=1, upper=5> spline_degree_cs_mean;

  int<lower=0> n_knots_cs_marker;
  vector[n_knots_cs_marker] knots_cs_marker;
  int<lower=0> n_coeff_cs_marker;
  vector[n_coeff_cs_marker] coeff_cs_marker;
  int<lower=1, upper=5> spline_degree_cs_marker;
