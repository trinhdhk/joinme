  int<lower=1> n_draws; // number of posterior draws to process

  /* Observed longitudinal history */
  int<lower=1> n_obs_long;                          // number of observed rows
  array[n_obs_long] int<lower=1> idx_marker_obs;    // marker index per observed row
  int<lower=1> n_marker_types;                      // total number of markers
  vector[n_marker_types] marker_weights;            // effective weights (reference only)
  matrix[n_draws, n_marker_types] marker_weights_draws; // per-draw weights
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
  int<lower=0> P_phi_beta;                       // phi_beta fixed-effect columns
  matrix[n_obs_long, P_phi_beta] X_phi_beta_obs; // phi_beta design observed
  matrix[n_obs_pred, P_phi_beta] X_phi_beta_pred; // phi_beta design pred
  int<lower=0> P_tau_sde;                        // tau_sde fixed-effect columns
  matrix[n_obs_long, P_tau_sde] X_tau_sde_obs;   // tau_sde design observed
  matrix[n_obs_pred, P_tau_sde] X_tau_sde_pred;  // tau_sde design pred

  /* Covariance regression covariates */
  int<lower=1> n_cov_corr;                        // covariate count for corr
  vector[n_cov_corr] vec_cov_corr;                // covariate vector (subject-level)

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
  matrix[n_gk, n_fixed_effects] mat_fixed_gk_cond_fwd;
  matrix[n_gk, n_random_id] mat_id_gk_cond_fwd;
  matrix[n_gk, n_random_marker] mat_marker_gk_cond_fwd;
  matrix[n_gk, n_random_marker_id] mat_marker_id_gk_cond_fwd;

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
  array[n_times_surv] matrix[n_gk, n_basehaz_basis] mat_basis_gk_surv;
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv;
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv;
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv;
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv;
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv_fwd;
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv_fwd;
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv_fwd;
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv_fwd;

  /* Model configuration flags and draws */
  array[n_marker_types] int<lower=1, upper=11> family_long; // family codes
  array[n_marker_types] int<lower=0, upper=5> link_long; // 0 custom VM; 1 identity, 2 log, 3 logit, 4 probit, 5 exp
  int<lower=1> max_inv_link_ops;
  array[n_marker_types] int<lower=0> inv_link_n_ops;
  array[n_marker_types, max_inv_link_ops] int<lower=0, upper=26> inv_link_ops;
  int<lower=1> max_inv_link_const;
  array[n_marker_types] int<lower=0> inv_link_n_const;
  matrix[n_marker_types, max_inv_link_const] inv_link_const;
  int<lower=0> n_family_sigma;
  array[n_marker_types] int<lower=0, upper=n_family_sigma> marker_to_sigma_family;
  int<lower=0> n_family_nu;
  array[n_marker_types] int<lower=0, upper=n_family_nu> marker_to_nu_family;
  int<lower=0> n_family_phi;
  array[n_marker_types] int<lower=0, upper=n_family_phi> marker_to_phi_family;
  int<lower=0> n_family_alpha;
  array[n_marker_types] int<lower=0, upper=n_family_alpha> marker_to_alpha_family;
  int<lower=0> n_family_phi_beta;
  array[n_marker_types] int<lower=0, upper=n_family_phi_beta> marker_to_phi_beta_family;
  int<lower=0> n_family_tau_sde;
  array[n_marker_types] int<lower=0, upper=n_family_tau_sde> marker_to_tau_sde_family;
  int<lower=0, upper=1> flag_resid_dim;                      // sigma dimension flag
  int<lower=0, upper=1> corr_diag_link;                      // 0 softplus, 1 exp
  int<lower=0, upper=1> use_tau_sde_fixed;                   // 1 use fixed tau
  real<lower=0, upper=1> tau_sde_fixed;                      // fixed tau value

  array[n_draws] vector[n_fixed_effects] beta_fixed;          // fixed effects per draw
  array[n_draws] vector[n_random_id] tau_id;                  // id RE scales per draw
  array[n_draws] matrix[n_random_id, n_random_id] Lcorr_id;   // id RE correlations
  array[n_draws] vector[n_random_marker] tau_marker;          // marker RE scales
  array[n_draws] matrix[n_random_marker, n_random_marker] Lcorr_marker; // marker RE corr
  array[n_draws] vector[n_random_marker_id] tau_marker_id;    // marker-id RE scales
  array[n_draws] matrix[n_random_marker_id, n_random_marker_id] Lcorr_marker_id; // marker-id corr
  array[n_draws] matrix[n_random_marker_id, n_random_marker] B_cross; // cross-corr map

  int<lower=0> num_unique_cov_entries;         // unique elements in L_i
  array[n_draws] vector[num_unique_cov_entries] alpha_corr_reg; // intercepts
  array[n_draws] vector[num_unique_cov_entries * n_cov_corr] beta_corr_reg_flat; // flattened slopes
  array[n_draws] real tau_corr_reg;            // scale for latent corr effect
  array[n_draws] vector[num_unique_cov_entries] lambda_corr_reg; // loadings

  int<lower=1> K_event;                        // number of competing risks
  array[n_draws, K_event] vector[n_basehaz_basis] bs_gamma_c; // spline coeffs
  array[n_draws, K_event] vector[n_cov_hazard] gamma_hazard;  // hazard covariates

  array[n_draws] vector[P_sigma] beta_sigma;    // sigma regression betas
  array[n_draws] vector[P_nu] beta_nu;          // nu regression betas
  array[n_draws] vector[P_phi] beta_phi;        // phi regression betas
  array[n_draws] vector[P_alpha] beta_alpha;    // alpha regression betas
  array[n_draws] vector[P_phi_beta] beta_phi_beta; // phi_beta regression betas
  array[n_draws] vector[P_tau_sde] beta_tau_sde;   // tau_sde regression betas

  array[n_draws] vector[n_family_sigma] sigma_family;          // family-shared sigma
  array[n_draws] vector[n_family_nu] nu_family;                // family-shared nu
  array[n_draws] vector[n_family_phi] phi_family;              // family-shared phi
  array[n_draws] vector[n_family_alpha] alpha_family;          // family-shared alpha
  array[n_draws] vector[n_family_phi_beta] phi_beta_family;    // family-shared phi_beta
  array[n_draws] vector[n_family_tau_sde] tau_sde_family;      // family-shared tau_sde

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

  int<lower=0, upper=1> flag_assoc_cv_total;    // include total CV association
  int<lower=0, upper=1> flag_assoc_cv_mean;     // include mean CV association
  int<lower=0, upper=1> flag_assoc_cv_marker;   // include marker CV association
  int<lower=0, upper=1> flag_assoc_cs_total;    // include total CS association
  int<lower=0, upper=1> flag_assoc_cs_mean;     // include mean CS association
  int<lower=0, upper=1> flag_assoc_cs_marker;   // include marker CS association
  int<lower=0, upper=1> flag_assoc_corr;        // include corr association

  int<lower=0, upper=3> tf_mode_cv_tot;        // transform mode: total CV
  int<lower=0, upper=3> tf_mode_cs_tot;        // transform mode: total CS
  int<lower=0, upper=3> tf_mode_corr;          // transform mode: corr
  int<lower=0, upper=3> tf_mode_cv_mean;       // transform mode: mean CV
  int<lower=0, upper=3> tf_mode_cv_marker;     // transform mode: marker CV
  int<lower=0, upper=3> tf_mode_cs_mean;       // transform mode: mean CS
  int<lower=0, upper=3> tf_mode_cs_marker;     // transform mode: marker CS

  int<lower=0> n_functional_ops_cv;            // op count for total CV
  array[n_functional_ops_cv] int<lower=0, upper=26> functional_ops_cv; // bytecode stream
  int<lower=0> n_const_cv;                     // constants for total CV bytecode
  vector[n_const_cv] const_data_cv;            // constants for total CV bytecode

  int<lower=0> n_functional_ops_cs;            // op count for total CS
  array[n_functional_ops_cs] int<lower=0, upper=26> functional_ops_cs; // bytecode stream
  int<lower=0> n_const_cs;                     // constants for total CS bytecode
  vector[n_const_cs] const_data_cs;            // constants for total CS bytecode

  int<lower=0> n_functional_ops_corr;          // op count for corr
  array[n_functional_ops_corr] int<lower=0, upper=26> functional_ops_corr; // bytecode stream
  int<lower=0> n_const_corr;                   // constants for corr bytecode
  vector[n_const_corr] const_data_corr;        // constants for corr bytecode

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
  array[n_draws] vector[n_coeff_corr] coeff_corr; // coefficients for corr spline by draw
  int<lower=1, upper=5> spline_degree_corr;    // spline degree for corr

  int<lower=0> n_functional_ops_cv_mean;        // op count for mean CV
  array[n_functional_ops_cv_mean] int<lower=0, upper=26> functional_ops_cv_mean; // bytecode stream
  int<lower=0> n_const_cv_mean;                 // constants for mean CV bytecode
  vector[n_const_cv_mean] const_data_cv_mean;   // constants for mean CV bytecode

  int<lower=0> n_functional_ops_cv_marker;      // op count for marker CV
  array[n_functional_ops_cv_marker] int<lower=0, upper=26> functional_ops_cv_marker; // bytecode stream
  int<lower=0> n_const_cv_marker;               // constants for marker CV bytecode
  vector[n_const_cv_marker] const_data_cv_marker; // constants for marker CV bytecode

  int<lower=0> n_functional_ops_cs_mean;        // op count for mean CS
  array[n_functional_ops_cs_mean] int<lower=0, upper=26> functional_ops_cs_mean; // bytecode stream
  int<lower=0> n_const_cs_mean;                 // constants for mean CS bytecode
  vector[n_const_cs_mean] const_data_cs_mean;   // constants for mean CS bytecode

  int<lower=0> n_functional_ops_cs_marker;      // op count for marker CS
  array[n_functional_ops_cs_marker] int<lower=0, upper=26> functional_ops_cs_marker; // bytecode stream
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
