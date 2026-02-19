/**
 * @file joinme_dynpred_threading.stan
 * @brief Threaded dynamic prediction twin that reuses the same `partial_draw` helper as the non-threaded model.
 *
 * ### Inputs
 * - Same as joinme_dynpred.stan plus `grainsize` to size draw slices for reduce_sum.
 *
 * ### Outputs
 * - The same conditional predictions as the non-threaded path; threading only changes how work is split, not what is computed.
 *
 * ### Execution story
 * - reduce_sum splits posterior draws into slices; each slice calls `partial_draw` (pure function) so math stays identical.
 * - No random numbers appear inside tasks; determinism holds given seed and thread topology.
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/functions/functional_transform.stanfunctions
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions

  #include helper/functions/joinme_dynpred_partial.stanfunctions
}

data {
  #include helper/data/joinme_dynpred_common.stan
  int<lower=1> grainsize;

  int flag_indep_id_re;
  int flag_indep_marker_re;
  int flag_indep_marker_byid_latent_re;
  int flag_indep_idmarker_cov;
  int flag_allow_marker_crosscorr;

  array[num_unique_cov_entries] int idx_row_cov;
  array[num_unique_cov_entries] int idx_col_cov;
}

transformed data {
  array[n_draws] int draw_ids;
  for (i in 1 : n_draws) draw_ids[i] = i;
}

parameters {
  #include helper/parameters/joinme_dynpred_common.stan
}

model {
  // -------------------- Threaded draw-wise conditioning likelihood
  // Each task handles a slice of draws; this is typically more stable than
  // subject-wise splits for prediction workloads.
  target += reduce_sum(
    partial_draw,
    draw_ids,
    grainsize,
    /* Observed longitudinal layout */
    n_obs_long,
    idx_marker_obs,
    n_marker_types,
    marker_weights_draws,
    y_real,
    y_int,
    trials_obs,
    /* Longitudinal design and dimensions */
    n_fixed_effects,
    n_random_id,
    n_random_marker,
    n_random_marker_id,
    mat_fixed_obs,
    mat_id_obs,
    mat_marker_obs,
    mat_marker_id_obs,
    /* Distributional regression design */
    P_sigma,
    X_sigma_obs,
    P_nu,
    X_nu_obs,
    P_phi,
    X_phi_obs,
    P_alpha,
    X_alpha_obs,
    P_phi_beta,
    X_phi_beta_obs,
    P_tau_sde,
    X_tau_sde_obs,
    /* Covariance regression inputs */
    n_cov_vcov,
    vec_cov_vcov,
    /* Hazard baseline inputs */
    n_cov_hazard,
    vec_cov_hazard,
    n_basehaz_basis,
    time_condition,
    /* GK node design (conditioning interval) */
    mat_basis_gk_cond,
    mat_fixed_gk_cond,
    mat_id_gk_cond,
    mat_marker_gk_cond,
    mat_marker_id_gk_cond,
    mat_fixed_gk_cond_fwd,
    mat_id_gk_cond_fwd,
    mat_marker_gk_cond_fwd,
    mat_marker_id_gk_cond_fwd,
    eps_finite_diff,
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
    vcov_diag_link,
    use_tau_sde_fixed,
    tau_sde_fixed,
    num_unique_cov_entries,
    idx_row_cov,
    idx_col_cov,
    /* Draw-specific fixed + distributional coefficients */
    beta_fixed,
    beta_sigma,
    beta_nu,
    beta_phi,
    beta_alpha,
    beta_phi_beta,
    beta_tau_sde,
    /* Draw-specific random effect scales/correlations */
    tau_id,
    Lcorr_id,
    tau_marker,
    Lcorr_marker,
    tau_marker_id,
    Lcorr_marker_id,
    B_cross,
    /* Draw-specific covariance regression weights */
    alpha_vcov_reg,
    beta_vcov_reg_flat,
    tau_vcov_reg,
    lambda_vcov_reg,
    /* Draw-specific survival parameters */
    K_event,
    log_h0_intercept,
    bs_gamma_c,
    gamma_hazard,
    /* Draw-specific longitudinal distributional parameters */
    sigma_family,
    nu_family,
    phi_family,
    alpha_family,
    phi_beta_family,
    tau_sde_family,
    /* Draw-specific association weights */
    coeff_assoc_cv_total,
    coeff_assoc_cs_total,
    coeff_assoc_cv_mean,
    coeff_assoc_cs_mean,
    coeff_assoc_cv_marker,
    coeff_assoc_cs_marker,
    coeff_assoc_vcov_var,
    K_ord,
    cutpoints_ord,
    /* Association toggles */
    flag_assoc_cv_total,
    flag_assoc_cs_total,
    flag_assoc_cv_mean,
    flag_assoc_cv_marker,
    flag_assoc_cs_mean,
    flag_assoc_cs_marker,
    flag_assoc_vcov,
    /* Transform modes (per association channel) */
    tf_mode_cv_tot,
    tf_mode_cs_tot,
    tf_mode_vcov,
    tf_mode_cv_mean,
    tf_mode_cv_marker,
    tf_mode_cs_mean,
    tf_mode_cs_marker,
    /* Transform configs (opcodes, constants, splines) */
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
    /* Flags and latent draws */
    flag_indep_id_re,
    flag_indep_marker_re,
    flag_indep_marker_byid_latent_re,
    flag_allow_marker_crosscorr,
    z_u,
    z_v,
    z_w_lat,
    z_L
  );
}

generated quantities {
  #include helper/generated_quantities/dynpred_outputs.stan
}
