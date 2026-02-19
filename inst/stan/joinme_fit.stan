/**
 * @file joinme_fit.stan
 * @brief Serial implementation of the Joint Mixed-Effects Model.
 *
 * @details
 * This model implements a multivariate longitudinal and survival joint model.
 * It uses the same underlying likelihood function (`partial_joinme`) as the threaded version,
 * but executes serially (grainsize = n_id).
 *
 * Key components:
 * - Longitudinal submodel: Generalized linear mixed effects with flexible distributions.
 * - Survival submodel: Proportional hazards with spline-based baseline hazard.
 * - Association: Random effects and their transformations link the submodels.
 *
 * @see joinme_fit_threading.stan
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
}
transformed data {
  #include helper/transformed_data/fit_cov_index.stan

  // Build per-id row ranges so both threaded and non-threaded paths share the same slice logic.
  array[n_id] int id_start;
  array[n_id] int id_end;
  array[n_id] int id_seq;

  for (i in 1 : n_id) {
    id_start[i] = N + 1;
    id_end[i] = 0;
    id_seq[i] = i;
  }
  for (n in 1 : N) {
    int i = id[n];
    if (n < id_start[i]) id_start[i] = n;
    id_end[i] = n;
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

  // Reuse the shared partial log-likelihood so the threaded and non-threaded
  // programs stay in perfect sync. Using grainsize=n_id yields a single chunk
  // (serial behavior) while still keeping the same code path as reduce_sum.
  target += reduce_sum(
    partial_joinme,
    id_seq,
    n_id,
    /* Longitudinal observation layout */
    id,
    marker,
    y_real,
    y_int,
    trials,
    family_long,
    n_family_sigma,
    marker_to_sigma_family,
    n_family_nu,
    marker_to_nu_family,
    n_family_phi,
    marker_to_phi_family,
    n_family_alpha,
    marker_to_alpha_family,
    n_family_phi_beta,
    marker_to_phi_beta_family,
    n_family_tau_sde,
    marker_to_tau_sde_family,
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
    sigma_family,
    nu_family,
    phi_family,
    alpha_family,
    phi_beta_family,
    tau_sde_family,
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
