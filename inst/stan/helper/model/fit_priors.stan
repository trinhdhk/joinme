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
  for (n in 1 : N) {
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
  beta_phi_beta ~ student_t(6, 0, 1);
  beta_tau_sde ~ student_t(6, 0, 1);
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
  for (j in 1 : n_re_phi_beta) { // phi_beta RE terms
    tau_phi_beta[j][1:K_phi_beta[j]] ~ exponential(1);
    for (g in 1 : G_phi_beta[j])
      target += re_weight_phi_beta[j][g] * std_normal_lpdf(to_vector(z_phi_beta[j][g, 1:K_phi_beta[j]]));
  }
  for (j in 1 : n_re_tau_sde) { // tau_sde RE terms
    tau_tau_sde[j][1:K_tau_sde[j]] ~ exponential(1);
    for (g in 1 : G_tau_sde[j])
      target += re_weight_tau_sde[j][g] * std_normal_lpdf(to_vector(z_tau_sde[j][g, 1:K_tau_sde[j]]));
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
  
  /* Marker-by-id latent random effects (tau_w ORIGINAL scale) */
  if (Q_idm > 0) {
    tau_w ~ exponential(1);
    Lcorr_w ~ lkj_corr_cholesky(lkj_eta);
    for (i in 1 : n_id)
      for (d in 1 : D)
        target += re_weight_idm[i] * std_normal_lpdf(z_w_lat[i, d]);
  }
  
  /* Covariance regression priors */
  {
    real corr_lp_scale = (corr_diag_link == 1) ? 0.2 : 0.3;
    alpha_L ~ student_t(6, 0, corr_lp_scale);
    for (m in 1 : M_cov) 
      beta_L[m] ~ student_t(6, 0, corr_lp_scale);
    tau_L ~ student_t(6, 0, corr_lp_scale);
    lambda_L ~ student_t(6, 0, corr_lp_scale);
    for (i in 1 : n_id)
      target += re_weight_L[i] * std_normal_lpdf(z_L[i]);
  }
  
  /* Baseline hazard priors (per event type) */
  for (k_ev in 1 : K_event) { // event-specific baseline hazard
    // Intercept is encoded in the first basis column (constant 1s).
    bs_gamma_c[k_ev][1] ~ normal(-3, alpha_scale);
    if (Kbs > 1)
      bs_gamma_c[k_ev][2:Kbs] ~ normal(0, alpha_scale);
    if (Kbs >= 4) // penalized spline second differences (exclude intercept)
      for (k in 4 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
    gamma_w[k_ev] ~ std_normal();
  }
  
  /* Distributional parameters (marker-specific + ordinal cutpoints) */
  sigma_family ~ exponential(1);
  nu_family ~ gamma(2, 1);
  phi_family ~ exponential(1);
  alpha_family ~ normal(0, 2);
  phi_beta_family ~ exponential(1);
  tau_sde_family ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  /* Association priors */
 
  sd_alpha_cv_total ~ normal(0, 0.5);
  sd_alpha_cs_total ~ normal(0, 0.5);
  sd_alpha_cv_mean ~ normal(0, 0.5);
  sd_alpha_cs_mean ~ normal(0, 0.5);
  sd_alpha_cv_marker ~ normal(0, 0.5);
  sd_alpha_cs_marker ~ normal(0, 0.5);
  s_corr ~ normal(0, 0.5);
  
  
  /* Shrinkage family switch for corr weights */
  if (shrinkage == 1) {
    alpha_corr ~ double_exponential(0, 1);
    z_alpha_cv_total ~ double_exponential(0, 1);
    z_alpha_cs_total ~ double_exponential(0, 1);
    z_alpha_cv_mean ~ double_exponential(0, 1);
    z_alpha_cs_mean ~ double_exponential(0, 1);
    z_alpha_cv_marker ~ double_exponential(0, 1);
    z_alpha_cs_marker ~ double_exponential(0, 1);
    // if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1)
    z_marker_weights ~ double_exponential(0, 1);
  } else if (shrinkage == 2) {
    alpha_corr ~ std_normal();
    z_alpha_cv_total ~ std_normal();
    z_alpha_cs_total ~ std_normal();
    z_alpha_cv_mean ~ std_normal();
    z_alpha_cs_mean ~ std_normal();
    z_alpha_cv_marker ~ std_normal();
    z_alpha_cs_marker ~ std_normal();
    // if (estimate_marker_weights == 1 && use_marker_weight_assoc == 1)
    z_marker_weights ~ std_normal();
  }
