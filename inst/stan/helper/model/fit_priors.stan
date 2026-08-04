  /**
  * @brief Prior distributions for all model parameters.
   *
   * Includes priors for:
   * - Fixed effects (beta) and distributional regression parameters
   * - Random effect variances and correlations
   * - Baseline hazard parameters
   * - Association coefficients
   * - Shrinkage scales for covariance regression
   */

  // -------------------- Basic outcome checks for discrete families
  for (n in 1:N) {
    int fam = family_long[marker[n]]; // family for this marker
    if (fam == 3) {
      if (y_int[n] != 0 && y_int[n] != 1)
        reject("Bernoulli requires y_int in {0,1}.");
    }
    if (fam == 4) {
      if (y_int[n] < 0 || y_int[n] > trials[n])
        reject("Binomial requires 0 <= y_int <= trials.");
    }
    if (fam == 5 || fam == 6) {
      if (y_int[n] < 0)
        reject("Poisson/NegBin require y_int >= 0.");
    }
    if (fam == 10) {
      if (y_real[n] <= 0 || y_real[n] >= 1)
        reject("Beta requires y_real in (0,1).");
    }
    if (fam == 11) {
      if (y_int[n] < 1 || y_int[n] > K_ord)
        reject("Cumulative logit requires y_int in {1..K_ord}.");
    }
  }
  
  // -------------------- Priors
  /* Fixed effects + distributional regressions */
  if (prior_beta_family == 1) {
    beta ~ student_t(prior_beta_df, prior_beta_mu, prior_beta_scale);
  } else if (prior_beta_family == 2) {
    beta ~ normal(prior_beta_mu, prior_beta_scale);
  } else if (prior_beta_family == 3) {
    beta ~ double_exponential(prior_beta_mu, prior_beta_scale);
  } else {
    vector[P] beta_horseshoe_scale = prior_beta_scale .* joinme_regularized_horseshoe_scale(
      horseshoe_local_beta,
      horseshoe_global_beta[1],
      horseshoe_slab_beta[1],
      prior_beta_slab_scale
    ); // regularised coefficient-specific scale after global and local shrinkage
    beta ~ normal(prior_beta_mu, beta_horseshoe_scale);
  }
  beta_sigma ~ student_t(6, 0, 1);
  beta_nu ~ student_t(6, 0, 1);
  beta_phi ~ student_t(6, 0, 1);
  beta_alpha ~ student_t(6, 0, 1);
  beta_kappa ~ student_t(6, 0, 1);
  beta_tau ~ student_t(6, 0, 1);
  /* Distributional random-effect scales + latent draws */
  for (j in 1 : n_re_sigma) { // sigma RE terms
    tau_sigma[j][1:K_sigma[j]] ~ exponential(1);
    for (g in 1 : G_sigma[j])
      target += re_weight_sigma[j][g] * std_normal_lpdf(to_vector(z_sigma[j][g, 1:K_sigma[j]]));
  }
  for (j in 1 : n_re_nu) { // nu RE terms
    tau_nu[j][1:K_nu[j]] ~ exponential(1);
    for (g in 1 : G_nu[j])
      target += re_weight_nu[j][g] * std_normal_lpdf(to_vector(z_nu[j][g, 1:K_nu[j]]));
  }
  for (j in 1 : n_re_phi) { // phi RE terms
    tau_phi[j][1:K_phi[j]] ~ exponential(1);
    for (g in 1 : G_phi[j])
      target += re_weight_phi[j][g] * std_normal_lpdf(to_vector(z_phi[j][g, 1:K_phi[j]]));
  }
  for (j in 1 : n_re_alpha) { // alpha RE terms
    tau_alpha[j][1:K_alpha[j]] ~ exponential(1);
    for (g in 1 : G_alpha[j])
      target += re_weight_alpha[j][g] * std_normal_lpdf(to_vector(z_alpha[j][g, 1:K_alpha[j]]));
  }
  for (j in 1 : n_re_kappa) { // kappa RE terms
    tau_kappa[j][1:K_kappa[j]] ~ exponential(1);
    for (g in 1 : G_kappa[j])
      target += re_weight_kappa[j][g] * std_normal_lpdf(to_vector(z_kappa[j][g, 1:K_kappa[j]]));
  }
  for (j in 1 : n_re_tau) { // tau RE terms
    tau_tau[j][1:K_tau[j]] ~ exponential(1);
    for (g in 1 : G_tau[j])
      target += re_weight_tau[j][g] * std_normal_lpdf(to_vector(z_tau[j][g, 1:K_tau[j]]));
  }
  
  /* ID-level random effects (tau_u is ORIGINAL scale; internal scaling in transformed parameters) */
  tau_u ~ exponential(1);
  Lcorr_u ~ lkj_corr_cholesky(lkj_eta);
  for (i in 1 : n_id)
    target += re_weight_id[i] * std_normal_lpdf(z_u[i]);
  
  /* Marker-only random effects (if present) */
  if (R_mk > 0) {
    tau_v ~ exponential(1);
    Lcorr_v ~ lkj_corr_cholesky(lkj_eta);
    for (d in 1 : D)
      target += re_weight_marker[d] * joinme_standard_prior_lpdf(
        z_v[d] | prior_marker_family, prior_marker_df
      );
    to_vector(B_cross) ~ std_normal();
  }
  
  /* Marker-by-id latent random effects */
  // These seeds are iid standard normal. Baseline covariance is carried by
  // alpha_L rather than by a second global covariance layer.
  if (Q_idm > 0) {
    for (i in 1 : n_id)
      for (d in 1 : D)
        target += re_weight_idm[i] * std_normal_lpdf(z_w_lat[i, d]);
  }
  
  /* Covariance regression priors */
  {
    real vcov_latent_scale = 1; // fixed half-Student-t scale for unexplained covariance-predictor heterogeneity

    // STEP 1: assign the declared family independently to the packed SD and
    // correlation coefficient blocks. Locations and scales enter through the
    // transformed-parameter mapping, whilst this density remains standardised.
    target += joinme_standard_prior_lpdf(
      vcov_sd_coefficient_raw |
      prior_vcov_sd_family,
      prior_vcov_sd_df
    );
    target += joinme_standard_prior_lpdf(
      vcov_corr_coefficient_raw |
      prior_vcov_corr_family,
      prior_vcov_corr_df
    );

    // STEP 2: add a regularised-horseshoe hierarchy only for a block which
    // selected that family. Zero-length declarations ensure all other families
    // add neither parameters nor irrelevant prior geometry.
    if (prior_vcov_sd_family == 4) {
      horseshoe_local_vcov_sd ~ student_t(prior_vcov_sd_df, 0, 1);
      horseshoe_global_vcov_sd ~ student_t(
        prior_vcov_sd_global_df, 0, prior_vcov_sd_global_scale
      );
      horseshoe_slab_vcov_sd ~ inv_gamma(
        0.5 * prior_vcov_sd_slab_df,
        0.5 * prior_vcov_sd_slab_df
      );
    }
    if (prior_vcov_corr_family == 4) {
      horseshoe_local_vcov_corr ~ student_t(prior_vcov_corr_df, 0, 1);
      horseshoe_global_vcov_corr ~ student_t(
        prior_vcov_corr_global_df, 0, prior_vcov_corr_global_scale
      );
      horseshoe_slab_vcov_corr ~ inv_gamma(
        0.5 * prior_vcov_corr_slab_df,
        0.5 * prior_vcov_corr_slab_df
      );
    }

    // STEP 3: retain the identified non-negative residual loading. This is a
    // hyperparameter, rather than an observed-covariate coefficient, so it is
    // deliberately separate from the two formula-specific coefficient priors.
    // Each lambda_L[m] is the residual standard deviation on one untransformed
    // covariance-predictor coordinate because ordinary z_L[i][m] is standard
    // Normal. The vector acts element-wise, equivalently through a diagonal
    // loading matrix. It is constrained nonnegative because its sign can
    // always be absorbed into the symmetric latent z_L coordinate.
    lambda_L ~ student_t(6, 0, vcov_latent_scale);
    for (i in 1 : n_id)
      target += re_weight_L[i] * std_normal_lpdf(z_L[i]);
  }
  
  /* Baseline hazard priors (per event type) */
  for (k_ev in 1 : K_event) { // event-specific baseline hazard
    // Intercept is encoded in the first basis column (constant 1s).
    bs_gamma_c[k_ev][1] ~ student_t(3, -5, 6);
    if (Kbs > 1)
      bs_gamma_c[k_ev][2:Kbs] ~ student_t(3, 0, 1);
    if (Kbs >= 4) // penalised spline second differences (exclude intercept)
      for (k in 4 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
    gamma_w[k_ev] ~ student_t(6, 0, 1);
  }
  
  /* Distributional parameters (marker-specific + ordinal cutpoints) */
  sigma_family ~ exponential(1);
  nu_family ~ gamma(2, 1);
  phi_family ~ exponential(1);
  alpha_family ~ normal(0, 2);
  kappa_family ~ exponential(1);
  tau_family ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  /**
   * @brief Exchangeable Dirichlet priors for piecewise-linear increments.
   *
   * @details A vector of ones gives a uniform prior over the simplex, so no
   * interval is preferred before seeing the longitudinal and event data. The
   * optional second-difference penalty below can add smoothness without
   * changing the ordering constraint. Existing I-splines retain latent
   * increment priors chosen by the alpha prior family below.
   */
  if ((tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) && n_free_spline_cv > 0)
    pwlin_simplex_cv ~ dirichlet(rep_vector(1, n_free_spline_cv));
  if ((tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) && n_free_spline_cs > 0)
    pwlin_simplex_cs ~ dirichlet(rep_vector(1, n_free_spline_cs));
  if ((tf_mode_corr == 3 || tf_mode_corr == 7) && n_free_spline_corr > 0)
    for (m in 1:M_corr)
      pwlin_simplex_corr[m] ~ dirichlet(rep_vector(1, n_free_spline_corr));
  if ((tf_mode_vcov == 3 || tf_mode_vcov == 7) && n_free_spline_vcov > 0)
    for (m in 1:M_vcov)
      pwlin_simplex_vcov[m] ~ dirichlet(rep_vector(1, n_free_spline_vcov));
  if ((tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) && n_free_spline_cv_mean > 0)
    pwlin_simplex_cv_mean ~ dirichlet(rep_vector(1, n_free_spline_cv_mean));
  if ((tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) && n_free_spline_cv_marker > 0)
    pwlin_simplex_cv_marker ~ dirichlet(rep_vector(1, n_free_spline_cv_marker));
  if ((tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) && n_free_spline_cs_mean > 0)
    pwlin_simplex_cs_mean ~ dirichlet(rep_vector(1, n_free_spline_cs_mean));
  if ((tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) && n_free_spline_cs_marker > 0)
    pwlin_simplex_cs_marker ~ dirichlet(rep_vector(1, n_free_spline_cs_marker));

  /**
   * Independent prior families for association, iota and marker weights.
   * Raw quantities are zero-centred and unit-scale; their scientific location,
   * scale and optional horseshoe transformation are applied in transformed
   * parameters. This keeps the five families independent of the class
   * component family, which is declared only by the mixture entry point.
   */
  target += joinme_standard_prior_lpdf(
    alpha_prior_raw | prior_alpha_family, prior_alpha_df
  );
  target += joinme_standard_prior_lpdf(
    z_marker_weights | prior_marker_weight_family, prior_marker_weight_df
  );
  target += joinme_standard_prior_lpdf(z_iota_intercept_cv | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cv | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_cs | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cs | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_corr | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_corr | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_vcov | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_vcov | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_cv_mean | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cv_mean | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_cv_marker | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cv_marker | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_cs_mean | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cs_mean | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_intercept_cs_marker | prior_iota_family, prior_iota_df);
  target += joinme_standard_prior_lpdf(z_iota_slope_cs_marker | prior_iota_family, prior_iota_df);

  if (prior_beta_family == 4) {
    horseshoe_local_beta ~ student_t(prior_beta_df, 0, 1);
    horseshoe_global_beta ~ student_t(prior_beta_global_df, 0, prior_beta_global_scale);
    horseshoe_slab_beta ~ inv_gamma(0.5 * prior_beta_slab_df, 0.5 * prior_beta_slab_df);
  }
  if (prior_alpha_family == 4) {
    horseshoe_local_alpha ~ student_t(prior_alpha_df, 0, 1);
    horseshoe_global_alpha ~ student_t(prior_alpha_global_df, 0, prior_alpha_global_scale);
    horseshoe_slab_alpha ~ inv_gamma(0.5 * prior_alpha_slab_df, 0.5 * prior_alpha_slab_df);
  }
  if (prior_iota_family == 4) {
    horseshoe_local_iota ~ student_t(prior_iota_df, 0, 1);
    horseshoe_global_iota ~ student_t(prior_iota_global_df, 0, prior_iota_global_scale);
    horseshoe_slab_iota ~ inv_gamma(0.5 * prior_iota_slab_df, 0.5 * prior_iota_slab_df);
  }
  if (prior_marker_family == 4) {
    horseshoe_local_marker ~ student_t(prior_marker_df, 0, 1);
    horseshoe_global_marker ~ student_t(prior_marker_global_df, 0, prior_marker_global_scale);
    horseshoe_slab_marker ~ inv_gamma(0.5 * prior_marker_slab_df, 0.5 * prior_marker_slab_df);
  }
  if (prior_marker_weight_family == 4) {
    horseshoe_local_marker_weight ~ student_t(prior_marker_weight_df, 0, 1);
    horseshoe_global_marker_weight ~ student_t(prior_marker_weight_global_df, 0, prior_marker_weight_global_scale);
    horseshoe_slab_marker_weight ~ inv_gamma(0.5 * prior_marker_weight_slab_df, 0.5 * prior_marker_weight_slab_df);
  }


  /*
   * Association-shape increments use the alpha family's standard raw law.
   * A regularised horseshoe has a standard-Normal raw law, so the simplex or
   * softmax shape coordinates remain separate from the horseshoe coefficient
   * scales and do not duplicate their global--local hierarchy.
   */
  if (prior_alpha_family == 3) {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ double_exponential(0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ double_exponential(0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ double_exponential(0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ double_exponential(0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ double_exponential(0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ double_exponential(0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ double_exponential(0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ double_exponential(0, 1);

  } else if (prior_alpha_family == 2 || prior_alpha_family == 4) {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ std_normal();
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ std_normal();
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ std_normal();
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ std_normal();
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ std_normal();
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ std_normal();
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ std_normal();
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ std_normal();

  } else {
    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ student_t(prior_alpha_df, 0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ student_t(prior_alpha_df, 0, 1);

  }

  

  // - first differences control monotone increase,
  // - second differences control wiggliness,
  // - lambda says how strongly we discourage bends.
  if (estimate_spline_cv == 1 && n_coeff_cv > 2 && lambda_spline_cv > 0) {
    for (j in 3:n_coeff_cv)
      target += -0.5 * lambda_spline_cv * square(coeff_cv_eff[j] - 2 * coeff_cv_eff[j - 1] + coeff_cv_eff[j - 2]);
  }
  if (estimate_spline_cs == 1 && n_coeff_cs > 2 && lambda_spline_cs > 0) {
    for (j in 3:n_coeff_cs)
      target += -0.5 * lambda_spline_cs * square(coeff_cs_eff[j] - 2 * coeff_cs_eff[j - 1] + coeff_cs_eff[j - 2]);
  }
  if (estimate_spline_corr == 1 && n_coeff_corr > 2 && lambda_spline_corr > 0) {
    for (m in 1:M_corr)
      for (j in 3:n_coeff_corr)
        target += -0.5 * lambda_spline_corr * square(coeff_corr_eff[m, j] - 2 * coeff_corr_eff[m, j - 1] + coeff_corr_eff[m, j - 2]);
  }
  if (estimate_spline_vcov == 1 && n_coeff_vcov > 2 && lambda_spline_vcov > 0) {
    for (m in 1:M_vcov)
      for (j in 3:n_coeff_vcov)
        target += -0.5 * lambda_spline_vcov * square(coeff_vcov_eff[m, j] - 2 * coeff_vcov_eff[m, j - 1] + coeff_vcov_eff[m, j - 2]);
  }
  if (estimate_spline_cv_mean == 1 && n_coeff_cv_mean > 2 && lambda_spline_cv_mean > 0) {
    for (j in 3:n_coeff_cv_mean)
      target += -0.5 * lambda_spline_cv_mean * square(coeff_cv_mean_eff[j] - 2 * coeff_cv_mean_eff[j - 1] + coeff_cv_mean_eff[j - 2]);
  }
  if (estimate_spline_cv_marker == 1 && n_coeff_cv_marker > 2 && lambda_spline_cv_marker > 0) {
    for (j in 3:n_coeff_cv_marker)
      target += -0.5 * lambda_spline_cv_marker * square(coeff_cv_marker_eff[j] - 2 * coeff_cv_marker_eff[j - 1] + coeff_cv_marker_eff[j - 2]);
  }
  if (estimate_spline_cs_mean == 1 && n_coeff_cs_mean > 2 && lambda_spline_cs_mean > 0) {
    for (j in 3:n_coeff_cs_mean)
      target += -0.5 * lambda_spline_cs_mean * square(coeff_cs_mean_eff[j] - 2 * coeff_cs_mean_eff[j - 1] + coeff_cs_mean_eff[j - 2]);
  }
  if (estimate_spline_cs_marker == 1 && n_coeff_cs_marker > 2 && lambda_spline_cs_marker > 0) {
    for (j in 3:n_coeff_cs_marker)
      target += -0.5 * lambda_spline_cs_marker * square(coeff_cs_marker_eff[j] - 2 * coeff_cs_marker_eff[j - 1] + coeff_cs_marker_eff[j - 2]);
  }
