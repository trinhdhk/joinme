/**
 * @file survival/data/dynamic_prediction.stan
 * @brief Declare conditioning and future event-time prediction grids.
 *
 * @details
 * The baseline-hazard bases, event covariates and Gauss--Kronrod designs define cumulative hazards before and after the conditioning time for every competing risk.
 */
  /* Hazard covariates */
  int<lower=0> n_cov_hazard;                      // hazard covariate count

  /* Conditioning information */
  int<lower=1> n_basehaz_basis;                   // baseline hazard basis count
  real<lower=0> time_condition;                   // conditioning time T_cond divided by the fitted event-time scale
  int<lower=1> n_gk;                              // quadrature node count (nodes/weights hardcoded in Stan)
  matrix[n_gk, n_basehaz_basis] mat_basis_gk_cond;  // centred baseline basis evaluated at original-time quadrature ordinates from zero to T_cond
  matrix[n_gk, n_cov_hazard] mat_cov_hazard_gk_cond; // Cox design at conditioning quadrature nodes on original study time
  matrix[n_gk, n_fixed_effects] mat_fixed_gk_cond;  // fixed effects at quadrature nodes
  matrix[n_gk, n_random_id] mat_id_gk_cond;         // id RE design at quadrature nodes
  matrix[n_gk, n_random_marker] mat_marker_gk_cond; // marker RE design at quadrature nodes
  matrix[n_gk, n_random_marker_id] mat_marker_id_gk_cond; // marker-id RE design

  /* Forward-step design matrices for slope (finite difference) */
  matrix[n_gk, n_fixed_effects] mat_fixed_gk_cond_fwd; // Role: matrix fixed effect Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_id] mat_id_gk_cond_fwd; // Role: matrix individual Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_marker] mat_marker_gk_cond_fwd; // Role: matrix marker Gauss-Kronrod conditioning forward difference.
  matrix[n_gk, n_random_marker_id] mat_marker_id_gk_cond_fwd; // Role: matrix marker individual Gauss-Kronrod conditioning forward difference.

  real<lower=1e-6> eps_finite_diff; // finite-difference interval in original study-time units

  /* Survival prediction grid */
  int<lower=1> n_times_surv;                     // number of survival time points
  vector[n_times_surv] vec_time_surv;            // future survival times divided by the fitted event-time scale
  array[n_times_surv] matrix[n_gk, n_basehaz_basis] mat_basis_gk_surv; // centred baseline bases at original-time future quadrature ordinates
  array[n_times_surv] matrix[n_gk, n_cov_hazard] mat_cov_hazard_gk_surv; // Cox design at future quadrature nodes on original study time
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv; // Role: matrix fixed effect Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv; // Role: matrix individual Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv; // Role: matrix marker Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv; // Role: matrix marker individual Gauss-Kronrod survival.
  array[n_times_surv] matrix[n_gk, n_fixed_effects] mat_fixed_gk_surv_fwd; // Role: matrix fixed effect Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_id] mat_id_gk_surv_fwd; // Role: matrix individual Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_marker] mat_marker_gk_surv_fwd; // Role: matrix marker Gauss-Kronrod survival forward difference.
  array[n_times_surv] matrix[n_gk, n_random_marker_id] mat_marker_id_gk_surv_fwd; // Role: matrix marker individual Gauss-Kronrod survival forward difference.

  int<lower=1> K_event;                        // number of competing risks
  array[n_draws, K_event] vector[n_basehaz_basis] bs_gamma_c; // spline coeffs
  array[n_draws, K_event] vector[n_cov_hazard] gamma_hazard;  // hazard covariates
