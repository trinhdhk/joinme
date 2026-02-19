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
    tau_sigma[j][1:K_sigma[j]] ~ std_normal();
    to_vector(z_sigma[j][1:G_sigma[j], 1:K_sigma[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_nu) { // nu RE terms
    tau_nu[j][1:K_nu[j]] ~ std_normal();
    to_vector(z_nu[j][1:G_nu[j], 1:K_nu[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_phi) { // phi RE terms
    tau_phi[j][1:K_phi[j]] ~ std_normal();
    to_vector(z_phi[j][1:G_phi[j], 1:K_phi[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_alpha) { // alpha RE terms
    tau_alpha[j][1:K_alpha[j]] ~ std_normal();
    to_vector(z_alpha[j][1:G_alpha[j], 1:K_alpha[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_phi_beta) { // phi_beta RE terms
    tau_phi_beta[j][1:K_phi_beta[j]] ~ std_normal();
    to_vector(z_phi_beta[j][1:G_phi_beta[j], 1:K_phi_beta[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_tau_sde) { // tau_sde RE terms
    tau_tau_sde[j][1:K_tau_sde[j]] ~ std_normal();
    to_vector(z_tau_sde[j][1:G_tau_sde[j], 1:K_tau_sde[j]]) ~ std_normal();
  }
  
  /* ID-level random effects (tau_u is ORIGINAL scale; internal scaling in transformed parameters) */
  tau_u ~ normal(0, 0.5);
  Lcorr_u ~ lkj_corr_cholesky(lkj_eta);
  for (i in 1 : n_id) 
    z_u[i] ~ std_normal();
  
  /* Marker-only random effects (if present) */
  if (R_mk > 0) {
    tau_v ~ normal(0, 0.5);
    Lcorr_v ~ lkj_corr_cholesky(lkj_eta);
    for (d in 1 : D) 
      z_v[d] ~ std_normal();
    to_vector(B_cross) ~ normal(0, 0.2);
  }
  
  /* Marker-by-id latent random effects (tau_w ORIGINAL scale) */
  if (Q_idm > 0) {
    tau_w ~ normal(0, 0.5);
    Lcorr_w ~ lkj_corr_cholesky(lkj_eta);
    for (i in 1 : n_id)
      for (d in 1 : D)
        z_w_lat[i, d] ~ std_normal();
  }
  
  /* Covariance regression priors */
  {
    real vcov_lp_scale = (vcov_diag_link == 1) ? 0.2 : 0.3;
    alpha_L ~ normal(0, vcov_lp_scale);
    for (m in 1 : M_cov) 
      beta_L[m] ~ normal(0, vcov_lp_scale);
    tau_L ~ normal(0, vcov_lp_scale);
    lambda_L ~ normal(0, vcov_lp_scale);
    z_L ~ std_normal();
  }
  
  /* Baseline hazard priors (per event type) */
  for (k_ev in 1 : K_event) { // event-specific baseline hazard
    log_h0_intercept[k_ev] ~ normal(-3, alpha_scale);
    bs_gamma_c[k_ev] ~ normal(0, alpha_scale);
    if (Kbs >= 3) // penalized spline second differences
      for (k in 3 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
    gamma_w[k_ev] ~ normal(0, 0.5);
  }
  
  /* Distributional parameters (marker-specific + ordinal cutpoints) */
  sigma_family ~ exponential(1);
  nu_family ~ gamma(2, 1);
  phi_family ~ exponential(1);
  alpha_family ~ normal(0, 2);
  phi_beta_family ~ exponential(1);
  tau_sde_family ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  /* Association priors (mean-side) */
  alpha_cv_total ~ std_normal();
  alpha_cs_total ~ std_normal();
  alpha_cv_mean ~ std_normal();
  alpha_cs_mean ~ std_normal();
  
  /* Marker-side shrinkage scales */

  /* Marker-weight shrinkage priors */
  tau_marker_weights ~ exponential(marker_weight_scale);
  s_cv_marker ~ normal(0, 0.2);
  s_cs_marker ~ normal(0, 0.2);
  s_vcov ~ normal(0, 0.2);
  
  /* Shrinkage family switch for marker weights, marker-side alphas, and vcov weights */
  if (shrinkage == 1) {
    z_marker_weights ~ double_exponential(0, 1);
    alpha_cv_marker ~ double_exponential(0, s_cv_marker);
    alpha_cs_marker ~ double_exponential(0, s_cs_marker);
    alpha_vcov_var ~ double_exponential(0, s_vcov);
  } else {
    z_marker_weights ~ std_normal();
    alpha_cv_marker ~ normal(0, s_cv_marker);
    alpha_cs_marker ~ normal(0, s_cs_marker);
    alpha_vcov_var ~ normal(0, s_vcov);
  }
