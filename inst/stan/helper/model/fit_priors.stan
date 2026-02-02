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
    int fam = family_long[marker[n]];
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
  beta ~ normal(0, beta_scale);
  beta_sigma ~ normal(0, 1);
  beta_nu ~ normal(0, 1);
  beta_phi ~ normal(0, 1);
  beta_alpha ~ normal(0, 1);
  beta_phi_beta ~ normal(0, 1);
  beta_tau_sde ~ normal(0, 1);
  for (j in 1 : n_re_sigma) {
    tau_sigma[j][1:K_sigma[j]] ~ normal(0, 1);
    to_vector(z_sigma[j][1:G_sigma[j], 1:K_sigma[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_nu) {
    tau_nu[j][1:K_nu[j]] ~ normal(0, 1);
    to_vector(z_nu[j][1:G_nu[j], 1:K_nu[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_phi) {
    tau_phi[j][1:K_phi[j]] ~ normal(0, 1);
    to_vector(z_phi[j][1:G_phi[j], 1:K_phi[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_alpha) {
    tau_alpha[j][1:K_alpha[j]] ~ normal(0, 1);
    to_vector(z_alpha[j][1:G_alpha[j], 1:K_alpha[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_phi_beta) {
    tau_phi_beta[j][1:K_phi_beta[j]] ~ normal(0, 1);
    to_vector(z_phi_beta[j][1:G_phi_beta[j], 1:K_phi_beta[j]]) ~ std_normal();
  }
  for (j in 1 : n_re_tau_sde) {
    tau_tau_sde[j][1:K_tau_sde[j]] ~ normal(0, 1);
    to_vector(z_tau_sde[j][1:G_tau_sde[j], 1:K_tau_sde[j]]) ~ std_normal();
  }
  
  // id RE priors (tau_u is ORIGINAL scale; internal scaling happens in transformed parameters)
  tau_u ~ normal(0, 0.5);
  Lcorr_u ~ lkj_corr_cholesky(lkj_eta);
  for (i in 1 : n_id) 
    z_u[i] ~ std_normal();
  
  // marker-only priors (if present)
  if (R_mk > 0) {
    tau_v ~ normal(0, 0.5);
    Lcorr_v ~ lkj_corr_cholesky(lkj_eta);
    for (d in 1 : D) 
      z_v[d] ~ std_normal();
    to_vector(B_cross) ~ normal(0, 0.2);
  }
  
  // marker-by-id latent priors (tau_w ORIGINAL scale)
  if (Q_idm > 0) {
    tau_w ~ normal(0, 0.5);
    Lcorr_w ~ lkj_corr_cholesky(lkj_eta);
    for (i in 1 : n_id)
      for (d in 1 : D)
        z_w_lat[i, d] ~ std_normal();
  }
  
  // covariance regression priors
  alpha_L ~ normal(0, 0.3);
  for (m in 1 : M_cov) 
    beta_L[m] ~ normal(0, 0.3);
  tau_L ~ normal(0, 0.3);
  lambda_L ~ normal(0, 0.3);
  z_L ~ std_normal();
  
  // baseline hazard priors (per event type)
  for (k_ev in 1 : K_event) {
    log_h0_intercept[k_ev] ~ normal(-3, alpha_scale);
    bs_gamma_c[k_ev] ~ normal(0, alpha_scale);
    if (Kbs >= 3)
      for (k in 3 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
    gamma_w[k_ev] ~ normal(0, 0.5);
  }
  
  // distributional parameters
  sigma_y ~ std_normal();
  sigma_marker ~ std_normal();
  nu_marker ~ gamma(2, 1);
  phi_nb_marker ~ exponential(1);
  alpha_skew_marker ~ normal(0, 2);
  phi_beta_marker ~ exponential(1);
  tau_sde_marker ~ beta(2, 2);
  cutpoints_ord ~ normal(0, 2);
  
  // association priors
  alpha_cv_total ~ normal(0, 0.5);
  alpha_cs_total ~ normal(0, 0.5);
  alpha_cv_mean ~ normal(0, 0.5);
  alpha_cs_mean ~ normal(0, 0.5);
  alpha_cv_marker ~ normal(0, 0.5);
  alpha_cs_marker ~ normal(0, 0.5);
  
  // marker-side shrinkage scales
  s_cv_marker ~ normal(0, 0.2);
  s_cs_marker ~ normal(0, 0.2);
  s_vcov ~ normal(0, 0.2);
  
  // shrinkage family switch for vcov weights (and marker-side alphas if you wish to keep them shrinky)
  if (shrinkage == 1) {
    alpha_vcov_var ~ double_exponential(0, s_vcov);
  } else {
    alpha_vcov_var ~ normal(0, s_vcov);
  }
