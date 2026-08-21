/**
 * @file longitudinal/model/fit.stan
 * @brief Validate longitudinal outcomes and assign longitudinal priors.
 *
 * @details
 * Discrete outcomes are checked against their mathematical support.  The fragment assigns priors to distributional random effects, nested longitudinal effects, covariance regression, family-level distributional parameters, ordinal cutpoints and marker-effect horseshoe auxiliaries.
 */
  /**
  * @brief Prior distributions for all model parameters.
   *
   * Includes priors for:
   * - longitudinal population and distributional-regression coefficients;
   * - random-effect variances and correlations;
   * - baseline-hazard parameters;
   * - longitudinal--event association coefficients; and
   * - covariance-regression heterogeneity scales.
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
  
  /* Distributional random-effect scales + latent draws.
   * The same Exponential(1) law generates an omitted `sd` in
   * `jm_truth(re_params = list(dist = ...))`. One scale is drawn for the
   * simulated population; the conditional random effects are then drawn at
   * their formula-defined grouping level.
   */
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
  
  /* ID-level random effects (tau_u is ORIGINAL scale; internal scaling in transformed parameters).
   * `simulate_joinme()` uses these exact population laws when the subject
   * block in `jm_truth()$re_params` omits `sd` or `corr`: independent
   * Exponential(1) scales and an LKJ-Cholesky correlation draw with the same
   * `lkj_eta` concentration.
   */
  tau_u ~ exponential(1);
  Lcorr_u ~ lkj_corr_cholesky(lkj_eta);
  for (i in 1 : n_id)
    target += re_weight_id[i] * std_normal_lpdf(z_u[i]);
  
  /* Marker-only random effects (if present).
   * Omitted marker-block covariance truths use the same Exponential and LKJ
   * laws before marker effects are generated.
   */
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

    // Retain the identified non-negative residual loading. This is a
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
  
  /* Distributional parameters (marker-specific + ordinal cutpoints) */
  sigma_family ~ exponential(1);
  nu_family ~ gamma(2, 1);
  phi_family ~ exponential(1);
  alpha_family ~ normal(0, 2);
  kappa_family ~ exponential(1);
  tau_family ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  if (prior_marker_family == 4) {
    horseshoe_local_marker ~ student_t(prior_marker_df, 0, 1);
    horseshoe_global_marker ~ student_t(prior_marker_global_df, 0, prior_marker_global_scale);
    horseshoe_slab_marker ~ inv_gamma(0.5 * prior_marker_slab_df, 0.5 * prior_marker_slab_df);
  }

