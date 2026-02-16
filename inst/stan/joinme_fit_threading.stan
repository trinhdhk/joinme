/**
 * @file joinme_fit_threading.stan
 * @brief Threaded implementation of the Joint Mixed-Effects Model using reduce_sum.
 *
 * @details
 * This model implements the same joint model as `joinme_fit.stan` but uses `reduce_sum`
 * to parallelize the likelihood computation over subjects.
 *
 * It uses the shared likelihood function `partial_joinme` to ensure mathematical equivalence
 * with the serial version.
 *
 * @see joinme_fit.stan
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  
  #include helper/functions/cumhaz.stanfunctions
  
  #include helper/functions/functional_transform.stanfunctions
  
  #include helper/functions/basis_functions.stanfunctions
  
  #include helper/functions/composite_transform.stanfunctions
  
  #include helper/functions/joinme_fit_partial.stanfunctions
}

data {
  #include helper/data/joinme_fit_common.stan
  array[n_id] int<lower=1, upper=N> id_start;
  array[n_id] int<lower=1, upper=N> id_end;
  int<lower=1> grainsize;
}

transformed data {
  #include helper/transformed_data/fit_cov_index.stan

  array[n_id] int id_seq;
  for (i in 1 : n_id) {
    id_seq[i] = i;
  }
}

parameters {
  #include helper/parameters/joinme_fit_common.stan
}

transformed parameters {
  #include helper/transformed_parameters/fit_scaling_and_effects.stan
}

model {
  #include helper/model/fit_priors.stan
  
  // -------------------- Threaded likelihood
  // Each task evaluates a slice of subjects (id_seq), summing both longitudinal
  // and survival contributions. The grainsize controls chunk size; adjust based
  // on subject count and per-subject cost for best throughput.
  target += reduce_sum(
    partial_joinme,
    id_seq,
    grainsize,
    /* Longitudinal observation layout */
    id,
    marker,
    y_real,
    y_int,
    trials,
    family_long,
    flag_resid_dim,
    use_tau_sde_fixed,
    tau_sde_fixed,
    X_obs,
    Z_id_obs,
    Z_mk_obs,
    Z_idm_obs,
    R_mk,
    Q_idm,
    /* Fixed effects + distributional regression design */
    beta_scaled,
    P_sigma,
    X_sigma,
    P_nu,
    X_nu,
    P_phi,
    X_phi,
    P_alpha,
    X_alpha,
    P_phi_beta,
    X_phi_beta,
    P_tau_sde,
    X_tau_sde,
    /* Distributional regression coefficients */
    beta_sigma,
    beta_nu,
    beta_phi,
    beta_alpha,
    beta_phi_beta,
    beta_tau_sde,
    /* Distributional random effects (by submodel) */
    n_re_sigma,
    K_sigma,
    K_sigma_max,
    Z_sigma,
    J_sigma,
    tau_sigma,
    z_sigma,
    n_re_nu,
    K_nu,
    K_nu_max,
    Z_nu,
    J_nu,
    tau_nu,
    z_nu,
    n_re_phi,
    K_phi,
    K_phi_max,
    Z_phi,
    J_phi,
    tau_phi,
    z_phi,
    n_re_alpha,
    K_alpha,
    K_alpha_max,
    Z_alpha,
    J_alpha,
    tau_alpha,
    z_alpha,
    n_re_phi_beta,
    K_phi_beta,
    K_phi_beta_max,
    Z_phi_beta,
    J_phi_beta,
    tau_phi_beta,
    z_phi_beta,
    n_re_tau_sde,
    K_tau_sde,
    K_tau_sde_max,
    Z_tau_sde,
    J_tau_sde,
    tau_tau_sde,
    z_tau_sde,
    /* Realized random effects + longitudinal distributional parameters */
    u_id,
    v_marker,
    w_idscaled,
    wbar_i,
    sigma_y,
    sigma_marker,
    nu_marker,
    phi_nb_marker,
    alpha_skew_marker,
    phi_beta_marker,
    tau_sde_marker,
    K_ord,
    cutpoints_ord,
    /* Hazard covariates */
    p_w,
    W,
    gamma_w,
    /* Baseline hazard spline pieces */
    log_h0_intercept,
    bs_gamma_c,
    Kbs,
    Bs_event_c,
    Bs_gk_c,
    /* Event times + integration controls */
    S_event,
    d_event,
    K_event,
    event_type,
    eps_fd,
    /* Design matrices at GK and event grids */
    X_gk_now,
    X_gk_fwd,
    Z_id_gk_now,
    Z_id_gk_fwd,
    X_event_now,
    X_event_fwd,
    Z_id_event_now,
    Z_id_event_fwd,
    Z_mk_gk_now,
    Z_mk_gk_fwd,
    Z_idm_gk_now,
    Z_idm_gk_fwd,
    Z_mk_event_now,
    Z_mk_event_fwd,
    Z_idm_event_now,
    Z_idm_event_fwd,
    /* Association features + covariance regression */
    vbar,
    L_i,
    a_cv_total,
    a_cs_total,
    a_cv_mean,
    a_cv_marker,
    a_cs_mean,
    a_cs_marker,
    a_vcov_var,
    /* Transform modes (which path to use) */
    tf_mode_cv_tot,
    tf_mode_cs_tot,
    tf_mode_vcov,
    tf_mode_cv_mean,
    tf_mode_cv_marker,
    tf_mode_cs_mean,
    tf_mode_cs_marker,
    /* Transform opcodes + constants */
    functional_ops_cv,
    const_data_cv,
    knots_cv,
    coeff_cv,
    spline_degree_cv,
    functional_ops_cs,
    const_data_cs,
    knots_cs,
    coeff_cs,
    spline_degree_cs,
    functional_ops_vcov,
    const_data_vcov,
    knots_vcov,
    coeff_vcov,
    spline_degree_vcov,
    functional_ops_cv_mean,
    const_data_cv_mean,
    knots_cv_mean,
    coeff_cv_mean,
    spline_degree_cv_mean,
    functional_ops_cv_marker,
    const_data_cv_marker,
    knots_cv_marker,
    coeff_cv_marker,
    spline_degree_cv_marker,
    functional_ops_cs_mean,
    const_data_cs_mean,
    knots_cs_mean,
    coeff_cs_mean,
    spline_degree_cs_mean,
    functional_ops_cs_marker,
    const_data_cs_marker,
    knots_cs_marker,
    coeff_cs_marker,
    spline_degree_cs_marker,
    /* Subject ranges + time scaling */
    id_start,
    id_end,
    tmax
  );
}

generated quantities {
  #include helper/generated_quantities/fit_outputs.stan
}
