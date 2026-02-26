  /* Fixed effects */
  vector[P] beta; // fixed-effect coefficients on original time scale

  /* id-level random effects */
  array[n_id] vector[R_id] z_u;          // latent standard normals per subject
  vector<lower=0>[R_id] tau_u;           // marginal SDs for id REs
  cholesky_factor_corr[R_id] Lcorr_u;    // Cholesky corr for id REs

  /* marker-only random effects */
  array[D] vector[R_mk] z_v;             // latent standard normals per marker
  vector<lower=0>[R_mk] tau_v;           // marginal SDs for marker REs
  cholesky_factor_corr[R_mk] Lcorr_v;    // Cholesky corr for marker REs

  /* marker-by-id latent effects */
  array[n_id, D] vector[Q_idm] z_w_lat;  // latent standard normals per id, marker
  vector<lower=0>[Q_idm] tau_w;          // marginal SDs for marker-id REs
  cholesky_factor_corr[Q_idm] Lcorr_w;   // Cholesky corr for marker-id REs

  /* Cross-correlation mapping between v_marker and z_w */
  matrix[Q_idm, R_mk] B_cross;           // cross-loadings when crosscorr enabled

  /* id-specific covariance regression for L_i */
  vector[M_cov] alpha_L;                 // intercepts for each L_i element
  array[M_cov] vector[K_cov] beta_L;     // covariate slopes for each element
  real<lower=0> tau_L;                   // scale for latent u_L
  vector[M_cov] lambda_L;                // loadings for latent u_L
  vector[n_id] z_L;                      // latent standard normals for u_L

  /* baseline hazard (cause-specific) */
  array[K_event] real log_h0_intercept;  // baseline log-hazard intercepts
  array[K_event] vector[Kbs] bs_gamma_c; // spline coefficients for baseline hazard

  /* hazard covariates (cause-specific) */
  array[K_event] vector[p_w] gamma_w;    // regression coefficients for W

  /* longitudinal distributional parameters */
  vector[P_sigma] beta_sigma;            // sigma regression coefficients
  vector[P_nu] beta_nu;                  // nu regression coefficients
  vector[P_phi] beta_phi;                // phi regression coefficients
  vector[P_alpha] beta_alpha;            // alpha regression coefficients
  vector[P_phi_beta] beta_phi_beta;      // phi_beta regression coefficients
  vector[P_tau_sde] beta_tau_sde;        // tau_sde regression coefficients
  array[n_re_sigma] vector<lower=0>[K_sigma_max] tau_sigma; // sigma RE SDs
  array[n_re_sigma] matrix[G_sigma_max, K_sigma_max] z_sigma; // sigma RE latents
  array[n_re_nu] vector<lower=0>[K_nu_max] tau_nu; // nu RE SDs
  array[n_re_nu] matrix[G_nu_max, K_nu_max] z_nu;  // nu RE latents
  array[n_re_phi] vector<lower=0>[K_phi_max] tau_phi; // phi RE SDs
  array[n_re_phi] matrix[G_phi_max, K_phi_max] z_phi;  // phi RE latents
  array[n_re_alpha] vector<lower=0>[K_alpha_max] tau_alpha; // alpha RE SDs
  array[n_re_alpha] matrix[G_alpha_max, K_alpha_max] z_alpha; // alpha RE latents
  array[n_re_phi_beta] vector<lower=0>[K_phi_beta_max] tau_phi_beta; // phi_beta RE SDs
  array[n_re_phi_beta] matrix[G_phi_beta_max, K_phi_beta_max] z_phi_beta; // phi_beta RE latents
  array[n_re_tau_sde] vector<lower=0>[K_tau_sde_max] tau_tau_sde; // tau_sde RE SDs
  array[n_re_tau_sde] matrix[G_tau_sde_max, K_tau_sde_max] z_tau_sde; // tau_sde RE latents
  vector<lower=0>[n_family_sigma] sigma_family;        // shared sigma by family
  vector<lower=2>[n_family_nu] nu_family;              // shared nu by family
  vector<lower=0>[n_family_phi] phi_family;            // shared phi by family
  vector[n_family_alpha] alpha_family;                 // shared alpha by family
  vector<lower=0>[n_family_phi_beta] phi_beta_family;  // shared phi_beta by family
  vector<lower=0, upper=1>[n_family_tau_sde] tau_sde_family; // shared tau_sde by family

  /* Ordinal cutpoints (shared across ordinal markers) */
  ordered[K_ord - 1] cutpoints_ord;     // cumulative-logit cutpoints

  /* association coefficients (non-centered for CV/CS total + marker) */
  real<lower=0> z_alpha_cv_total;               // latent for total CV association
  real<lower=0> z_alpha_cs_total;               // latent for total CS association
  real z_alpha_cv_mean;                // latent for mean CV association
  real z_alpha_cs_mean;                // latent for mean CS association
  real<lower=0> z_alpha_cv_marker;              // latent for marker CV association
  real<lower=0> z_alpha_cs_marker;              // latent for marker CS association
  vector[M_corr] alpha_corr; // corr association coefficients (off-diagonal correlations)

  /* marker-weight perturbations (global across ids) */
  vector[D * estimate_marker_weights * use_marker_weight_assoc] z_marker_weights; // latent signed perturbations (active only when needed)

  /* association scales */
  real<lower=0> sd_alpha_cv_total;     // scale for total CV association
  real<lower=0> sd_alpha_cs_total;     // scale for total CS association
  real<lower=0> sd_alpha_cv_mean;      // scale for mean CV association
  real<lower=0> sd_alpha_cs_mean;      // scale for mean CS association
  real<lower=0> sd_alpha_cv_marker;    // scale for marker CV association
  real<lower=0> sd_alpha_cs_marker;    // scale for marker CS association
  real<lower=0> s_corr;                // scale for corr association
