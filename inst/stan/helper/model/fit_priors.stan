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
  beta ~ student_t(6, 0, beta_scale);
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
      target += re_weight_marker[d] * std_normal_lpdf(z_v[d]);
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
    real vcov_lp_scale = 1; // Role: covariance lp scale.
    alpha_L ~ student_t(6, 0, vcov_lp_scale);
    for (m in 1 : M_cov) 
      beta_L[m] ~ student_t(6, 0, vcov_lp_scale);
    // Each lambda_L[m] is the residual standard deviation on one untransformed
    // covariance-predictor coordinate because ordinary z_L[i][m] is standard
    // Normal. The vector acts element-wise, equivalently through a diagonal
    // loading matrix. It is constrained nonnegative because its sign can
    // always be absorbed into the symmetric latent z_L coordinate.
    lambda_L ~ student_t(6, 0, vcov_lp_scale);
    for (i in 1 : n_id)
      target += re_weight_L[i] * std_normal_lpdf(z_L[i]);
  }
  
  /* Baseline hazard priors (per event type) */
  for (k_ev in 1 : K_event) { // event-specific baseline hazard
    // Intercept is encoded in the first basis column (constant 1s).
    // bs_gamma_c[k_ev][1] ~ normal(-3, alpha_scale);
    bs_gamma_c[k_ev][1] ~ student_t(3, -5, 6);
    if (Kbs > 1)
      bs_gamma_c[k_ev][2:Kbs] ~ student_t(3, 0, alpha_scale);
    if (Kbs >= 4) // penalised spline second differences (exclude intercept)
      for (k in 4 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
    if (shrinkage == 1) {
      gamma_w[k_ev] ~ double_exponential(0, alpha_scale);
    } else if (shrinkage == 2) {
      gamma_w[k_ev] ~ normal(0, alpha_scale);
    } else {
      gamma_w[k_ev] ~ student_t(6, 0, alpha_scale);
    }
  }
  
  /* Distributional parameters (marker-specific + ordinal cutpoints) */
  sigma_family ~ exponential(1);
  nu_family ~ gamma(2, 1);
  phi_family ~ exponential(1);
  alpha_family ~ normal(0, 2);
  kappa_family ~ exponential(1);
  tau_family ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  /* Association priors */
 
  sd_alpha_cv_total ~ exponential(1);
  sd_alpha_cs_total ~ exponential(1);
  sd_alpha_cv_mean ~ exponential(1);
  sd_alpha_cs_mean ~ exponential(1);
  sd_alpha_cv_marker ~ exponential(1);
  sd_alpha_cs_marker ~ exponential(1);
  // s_corr / s_vcov are global half-normal scales for covariance-style
  // association coefficients. z_alpha_corr / z_alpha_vcov are standardized
  // coefficient latents, and the effective hazard coefficients are defined in
  // transformed parameters as s_* times those latents.
  s_corr ~ exponential(1);
  s_vcov ~ exponential(1);

  /**
   * @brief Exchangeable Dirichlet priors for piecewise-linear increments.
   *
   * @details A vector of ones gives a uniform prior over the simplex, so no
   * interval is preferred before seeing the longitudinal and event data. The
   * optional second-difference penalty below can add smoothness without
   * changing the ordering constraint. Existing I-splines retain their latent
   * increment priors in the shrinkage-family branches below.
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

  
  /* Shrinkage family switch for corr weights */
  if (shrinkage == 1) {
    z_alpha_corr ~ double_exponential(0, 1);
    z_alpha_vcov ~ double_exponential(0, 1);
    z_alpha_cv_total ~ double_exponential(0, 1);
    z_alpha_cs_total ~ double_exponential(0, 1);
    z_alpha_cv_mean ~ double_exponential(0, 1);
    z_alpha_cs_mean ~ double_exponential(0, 1);
    z_alpha_cv_marker ~ double_exponential(0, 1);
    z_alpha_cs_marker ~ double_exponential(0, 1);
    // if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1)
    z_marker_weights ~ double_exponential(0, 1);

    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ double_exponential(0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ double_exponential(0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ double_exponential(0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ double_exponential(0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ double_exponential(0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ double_exponential(0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ double_exponential(0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ double_exponential(0, 1);

    if (estimate_iota_intercept_cv > 0) z_iota_intercept_cv ~ std_normal();
    if (estimate_iota_slope_cv > 0) z_iota_slope_cv ~ std_normal();
    if (estimate_iota_intercept_cs > 0) z_iota_intercept_cs ~ std_normal();
    if (estimate_iota_slope_cs > 0) z_iota_slope_cs ~ std_normal();
    if (estimate_iota_intercept_corr > 0) z_iota_intercept_corr ~ std_normal();
    if (estimate_iota_slope_corr > 0) z_iota_slope_corr ~ std_normal();
    if (estimate_iota_intercept_vcov > 0) z_iota_intercept_vcov ~ std_normal();
    if (estimate_iota_slope_vcov > 0) z_iota_slope_vcov ~ std_normal();
    if (estimate_iota_intercept_cv_mean > 0) z_iota_intercept_cv_mean ~ std_normal();
    if (estimate_iota_slope_cv_mean > 0) z_iota_slope_cv_mean ~ std_normal();
    if (estimate_iota_intercept_cv_marker > 0) z_iota_intercept_cv_marker ~ std_normal();
    if (estimate_iota_slope_cv_marker > 0) z_iota_slope_cv_marker ~ std_normal();
    if (estimate_iota_intercept_cs_mean > 0) z_iota_intercept_cs_mean ~ std_normal();
    if (estimate_iota_slope_cs_mean > 0) z_iota_slope_cs_mean ~ std_normal();
    if (estimate_iota_intercept_cs_marker > 0) z_iota_intercept_cs_marker ~ std_normal();
    if (estimate_iota_slope_cs_marker > 0) z_iota_slope_cs_marker ~ std_normal();
  } else if (shrinkage == 2) {
    z_alpha_corr ~ std_normal();
    z_alpha_vcov ~ std_normal();
    z_alpha_cv_total ~ std_normal();
    z_alpha_cs_total ~ std_normal();
    z_alpha_cv_mean ~ std_normal();
    z_alpha_cs_mean ~ std_normal();
    z_alpha_cv_marker ~ std_normal();
    z_alpha_cs_marker ~ std_normal();
    // if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1)
    z_marker_weights ~ std_normal();

    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ std_normal();
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ std_normal();
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ std_normal();
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ std_normal();
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ std_normal();
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ std_normal();
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ std_normal();
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ std_normal();

    if (estimate_iota_intercept_cv > 0) z_iota_intercept_cv ~ std_normal();
    if (estimate_iota_slope_cv > 0) z_iota_slope_cv ~ std_normal();
    if (estimate_iota_intercept_cs > 0) z_iota_intercept_cs ~ std_normal();
    if (estimate_iota_slope_cs > 0) z_iota_slope_cs ~ std_normal();
    if (estimate_iota_intercept_corr > 0) z_iota_intercept_corr ~ std_normal();
    if (estimate_iota_slope_corr > 0) z_iota_slope_corr ~ std_normal();
    if (estimate_iota_intercept_vcov > 0) z_iota_intercept_vcov ~ std_normal();
    if (estimate_iota_slope_vcov > 0) z_iota_slope_vcov ~ std_normal();
    if (estimate_iota_intercept_cv_mean > 0) z_iota_intercept_cv_mean ~ std_normal();
    if (estimate_iota_slope_cv_mean > 0) z_iota_slope_cv_mean ~ std_normal();
    if (estimate_iota_intercept_cv_marker > 0) z_iota_intercept_cv_marker ~ std_normal();
    if (estimate_iota_slope_cv_marker > 0) z_iota_slope_cv_marker ~ std_normal();
    if (estimate_iota_intercept_cs_mean > 0) z_iota_intercept_cs_mean ~ std_normal();
    if (estimate_iota_slope_cs_mean > 0) z_iota_slope_cs_mean ~ std_normal();
    if (estimate_iota_intercept_cs_marker > 0) z_iota_intercept_cs_marker ~ std_normal();
    if (estimate_iota_slope_cs_marker > 0) z_iota_slope_cs_marker ~ std_normal();
  } else {
    z_alpha_corr ~ student_t(6, 0, 1);
    z_alpha_vcov ~ student_t(6, 0, 1);
    z_alpha_cv_total ~ student_t(6, 0, 1);
    z_alpha_cs_total ~ student_t(6, 0, 1);
    z_alpha_cv_mean ~ student_t(6, 0, 1);
    z_alpha_cs_mean ~ student_t(6, 0, 1);
    z_alpha_cv_marker ~ student_t(6, 0, 1);
    z_alpha_cs_marker ~ student_t(6, 0, 1);
    // if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1)
    z_marker_weights ~ student_t(6, 0, 1);

    /* Penalised monotone I-spline shape priors. */
    if (tf_mode_cv_tot != 3 && tf_mode_cv_tot != 7 && n_free_spline_cv > 0) z_spline_cv ~ student_t(6, 0, 1);
    if (tf_mode_cs_tot != 3 && tf_mode_cs_tot != 7 && n_free_spline_cs > 0) z_spline_cs ~ student_t(6, 0, 1);
    if (tf_mode_corr != 3 && tf_mode_corr != 7 && n_free_spline_corr > 0) to_vector(z_spline_corr) ~ student_t(6, 0, 1);
    if (tf_mode_vcov != 3 && tf_mode_vcov != 7 && n_free_spline_vcov > 0) to_vector(z_spline_vcov) ~ student_t(6, 0, 1);
    if (tf_mode_cv_mean != 3 && tf_mode_cv_mean != 7 && n_free_spline_cv_mean > 0) z_spline_cv_mean ~ student_t(6, 0, 1);
    if (tf_mode_cv_marker != 3 && tf_mode_cv_marker != 7 && n_free_spline_cv_marker > 0) z_spline_cv_marker ~ student_t(6, 0, 1);
    if (tf_mode_cs_mean != 3 && tf_mode_cs_mean != 7 && n_free_spline_cs_mean > 0) z_spline_cs_mean ~ student_t(6, 0, 1);
    if (tf_mode_cs_marker != 3 && tf_mode_cs_marker != 7 && n_free_spline_cs_marker > 0) z_spline_cs_marker ~ student_t(6, 0, 1);

    if (estimate_iota_intercept_cv > 0) z_iota_intercept_cv ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cv > 0) z_iota_slope_cv ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_cs > 0) z_iota_intercept_cs ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cs > 0) z_iota_slope_cs ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_corr > 0) z_iota_intercept_corr ~ student_t(6, 0, 1);
    if (estimate_iota_slope_corr > 0) z_iota_slope_corr ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_vcov > 0) z_iota_intercept_vcov ~ student_t(6, 0, 1);
    if (estimate_iota_slope_vcov > 0) z_iota_slope_vcov ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_cv_mean > 0) z_iota_intercept_cv_mean ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cv_mean > 0) z_iota_slope_cv_mean ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_cv_marker > 0) z_iota_intercept_cv_marker ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cv_marker > 0) z_iota_slope_cv_marker ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_cs_mean > 0) z_iota_intercept_cs_mean ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cs_mean > 0) z_iota_slope_cs_mean ~ student_t(6, 0, 1);
    if (estimate_iota_intercept_cs_marker > 0) z_iota_intercept_cs_marker ~ student_t(6, 0, 1);
    if (estimate_iota_slope_cs_marker > 0) z_iota_slope_cs_marker ~ student_t(6, 0, 1);
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
