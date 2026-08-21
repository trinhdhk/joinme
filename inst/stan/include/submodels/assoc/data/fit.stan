/**
 * @file assoc/data/fit.stan
 * @brief Data declarations linking longitudinal features to the event hazard.
 *
 * @details
 * Present and forward-shifted designs evaluate current values and finite-
 * difference slopes at event times and quadrature nodes.  Inclusion flags
 * select total, subject-mean, marker, correlation and covariance channels.
 * n_alpha_prior records the packed association coefficients followed by any
 * fitted common marker-weight locations.
 */
  /* Finite difference control and association design matrices */
  real<lower=1e-6> eps_fd;       // finite-difference step for slope approximation

  /* Mean association designs (CV_mean) */
  array[N_event] matrix[n_gk, P] X_gk_now; // Role: design matrix Gauss-Kronrod now.
  array[N_event] matrix[n_gk, P] X_gk_fwd; // Role: design matrix Gauss-Kronrod forward difference.
  array[N_event] matrix[n_gk, R_id] Z_id_gk_now; // Role: random-effect design individual Gauss-Kronrod now.
  array[N_event] matrix[n_gk, R_id] Z_id_gk_fwd; // Role: random-effect design individual Gauss-Kronrod forward difference.
  matrix[N_event, P] X_event_now; // Role: design matrix event now.
  matrix[N_event, P] X_event_fwd; // Role: design matrix event forward difference.
  matrix[N_event, R_id] Z_id_event_now; // Role: random-effect design individual event now.
  matrix[N_event, R_id] Z_id_event_fwd; // Role: random-effect design individual event forward difference.

  /* Marker association designs (CV_marker) */
  array[N_event] matrix[n_gk, R_mk] Z_mk_gk_now; // Role: random-effect design mk Gauss-Kronrod now.
  array[N_event] matrix[n_gk, R_mk] Z_mk_gk_fwd; // Role: random-effect design mk Gauss-Kronrod forward difference.
  array[N_event] matrix[n_gk, Q_idm] Z_idm_gk_now; // Role: random-effect design idm Gauss-Kronrod now.
  array[N_event] matrix[n_gk, Q_idm] Z_idm_gk_fwd; // Role: random-effect design idm Gauss-Kronrod forward difference.
  matrix[N_event, R_mk] Z_mk_event_now; // Role: random-effect design mk event now.
  matrix[N_event, R_mk] Z_mk_event_fwd; // Role: random-effect design mk event forward difference.
  matrix[N_event, Q_idm] Z_idm_event_now; // Role: random-effect design idm event now.
  matrix[N_event, Q_idm] Z_idm_event_fwd; // Role: random-effect design idm event forward difference.

  /* Association include flags */
  int<lower=0, upper=1> assoc_cv_total;   // include total current value term
  int<lower=0, upper=1> assoc_cv_mean;    // include mean current value term
  int<lower=0, upper=1> assoc_cv_marker;  // include marker current value term
  int<lower=0, upper=1> assoc_cs_total;   // include total current slope term
  int<lower=0, upper=1> assoc_cs_mean;    // include mean current slope term
  int<lower=0, upper=1> assoc_cs_marker;  // include marker current slope term
  int<lower=0, upper=1> assoc_corr;       // include corr association term
  int<lower=0, upper=1> assoc_vcov;       // include vcov association term

  int<lower=6> n_alpha_prior; // six scalar associations, corr and vcov components, then fitted marker-weight means



