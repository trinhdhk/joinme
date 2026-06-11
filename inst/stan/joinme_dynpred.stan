/**
 * @file joinme_dynpred.stan
 * @brief Non-threaded dynamic prediction for joint longitudinal + survival models, now sharing one likelihood helper with the threaded twin.
 *
 * ### Inputs
 * - Posterior draws for all fitted parameters (beta, random effects scales/correlations, baseline hazard, association weights, covariance regression weights).
 * - New-subject data and design matrices for longitudinal and survival conditioning, hazard covariates, and transform controls mirroring the fit model.
 * - Flags describing independence/cross-correlation and association toggles.
 *
 * ### Outputs
 * - Inference for the new subject (random effects draws, conditional survival summaries, etc.) via generated quantities.
 *
 * ### Flow
 * - `draw_ids` and `grainsize` are built in transformed data so we can reuse the same `partial_draw` reducer path as the threaded model while staying serial here.
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  #include helper/functions/eta_corr_varonly_weighted_const.stanfunctions
  #include helper/functions/eta_vcov_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/functions/bytecode_transform.stanfunctions
  #include helper/functions/link_functions.stanfunctions
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions

  #include helper/functions/joinme_dynpred_partial.stanfunctions
}

data {
  #include helper/data/joinme_dynpred_common.stan

  int flag_indep_id_re;
  int flag_indep_marker_re;
  int flag_indep_idmarker_cov;
  int flag_allow_marker_crosscorr;

  // Mappings for covariance regression (pre-calculated indices)
  array[num_unique_cov_entries] int idx_row_cov;
  array[num_unique_cov_entries] int idx_col_cov;
}

transformed data {
  array[n_draws] int draw_ids;
  int grainsize = n_draws;
  for (i in 1 : n_draws) draw_ids[i] = i;
}
parameters {
  #include helper/parameters/joinme_dynpred_common.stan
}
model {
  target += reduce_sum(
    partial_draw,
    draw_ids,
    grainsize,
    /* Observed longitudinal layout */
    n_obs_long,
    idx_marker_obs,
    n_marker_types,
    marker_weights_draws_cv_total,
    marker_weights_draws_cs_total,
    marker_weights_draws_cv_marker,
    marker_weights_draws_cs_marker,
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
    marker_id_row_scale,
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
    n_gk,
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
    link_long,
    max_inv_link_ops,
    inv_link_n_ops,
    inv_link_ops,
    max_inv_link_const,
    inv_link_n_const,
    inv_link_const,
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
    B_cross,
    /* Draw-specific covariance regression weights */
    alpha_vcov_reg,
    beta_vcov_reg_flat,
    lambda_vcov_reg,
    /* Draw-specific survival parameters */
    K_event,
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
    coeff_assoc_corr,
    coeff_assoc_vcov,
    iota_intercept_cv,
    iota_slope_cv,
    iota_intercept_cs,
    iota_slope_cs,
    iota_intercept_corr,
    iota_slope_corr,
    iota_intercept_vcov,
    iota_slope_vcov,
    iota_intercept_cv_mean,
    iota_slope_cv_mean,
    iota_intercept_cv_marker,
    iota_slope_cv_marker,
    iota_intercept_cs_mean,
    iota_slope_cs_mean,
    iota_intercept_cs_marker,
    iota_slope_cs_marker,
    K_ord,
    cutpoints_ord,
    /* Association toggles */
    flag_assoc_cv_total,
    flag_assoc_cs_total,
    flag_assoc_cv_mean,
    flag_assoc_cv_marker,
    flag_assoc_cs_mean,
    flag_assoc_cs_marker,
    flag_assoc_corr,
    flag_assoc_vcov,
    /* Transform modes (per association channel) */
    tf_mode_cv_tot,
    tf_mode_cs_tot,
    tf_mode_corr,
    tf_mode_vcov,
    tf_mode_cv_mean,
    tf_mode_cv_marker,
    tf_mode_cs_mean,
    tf_mode_cs_marker,
    /* Transform configs (bytecode, constants, splines) */
    functional_ops_cv,
    functional_iota_intercept_idx_cv,
    functional_iota_slope_idx_cv,
    const_data_cv,
    knots_cv,
    coeff_cv,
    spline_degree_cv,
    functional_ops_cs,
    functional_iota_intercept_idx_cs,
    functional_iota_slope_idx_cs,
    const_data_cs,
    knots_cs,
    coeff_cs,
    spline_degree_cs,
    functional_ops_corr,
    functional_iota_intercept_idx_corr,
    functional_iota_slope_idx_corr,
    const_data_corr,
    knots_corr,
    coeff_corr,
    spline_degree_corr,
    functional_ops_vcov,
    functional_iota_intercept_idx_vcov,
    functional_iota_slope_idx_vcov,
    const_data_vcov,
    knots_vcov,
    coeff_vcov,
    spline_degree_vcov,
    functional_ops_cv_mean,
    functional_iota_intercept_idx_cv_mean,
    functional_iota_slope_idx_cv_mean,
    const_data_cv_mean,
    knots_cv_mean,
    coeff_cv_mean,
    spline_degree_cv_mean,
    functional_ops_cv_marker,
    functional_iota_intercept_idx_cv_marker,
    functional_iota_slope_idx_cv_marker,
    const_data_cv_marker,
    knots_cv_marker,
    coeff_cv_marker,
    spline_degree_cv_marker,
    functional_ops_cs_mean,
    functional_iota_intercept_idx_cs_mean,
    functional_iota_slope_idx_cs_mean,
    const_data_cs_mean,
    knots_cs_mean,
    coeff_cs_mean,
    spline_degree_cs_mean,
    functional_ops_cs_marker,
    functional_iota_intercept_idx_cs_marker,
    functional_iota_slope_idx_cs_marker,
    const_data_cs_marker,
    knots_cs_marker,
    coeff_cs_marker,
    spline_degree_cs_marker,
    /* Flags and latent draws */
    flag_indep_id_re,
    flag_indep_marker_re,
    flag_indep_idmarker_cov,
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
