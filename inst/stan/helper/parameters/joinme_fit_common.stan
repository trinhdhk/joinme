  // Fixed effects
  vector[P] beta;

  // id-level random effects
  array[n_id] vector[R_id] z_u;
  vector<lower=0>[R_id] tau_u;
  cholesky_factor_corr[R_id] Lcorr_u;

  // marker-only random effects
  array[D] vector[R_mk] z_v;
  vector<lower=0>[R_mk] tau_v;
  cholesky_factor_corr[R_mk] Lcorr_v;

  // marker-by-id latent effects
  array[n_id, D] vector[Q_idm] z_w_lat;
  vector<lower=0>[Q_idm] tau_w;
  cholesky_factor_corr[Q_idm] Lcorr_w;

  // Cross-correlation mapping between v_marker and z_w
  matrix[Q_idm, R_mk] B_cross;

  // id-specific covariance regression for L_i
  vector[M_cov] alpha_L;
  array[M_cov] vector[K_cov] beta_L;
  real<lower=0> tau_L;
  vector[M_cov] lambda_L;
  vector[n_id] z_L;

  // baseline hazard (cause-specific)
  array[K_event] real log_h0_intercept;
  array[K_event] vector[Kbs] bs_gamma_c;

  // hazard covariates (cause-specific)
  array[K_event] vector[p_w] gamma_w;

  // longitudinal distributional parameters
  vector[P_sigma] beta_sigma;
  vector[P_nu] beta_nu;
  vector[P_phi] beta_phi;
  vector[P_alpha] beta_alpha;
  vector[P_phi_beta] beta_phi_beta;
  vector[P_tau_sde] beta_tau_sde;
  array[n_re_sigma] vector<lower=0>[K_sigma_max] tau_sigma;
  array[n_re_sigma] matrix[G_sigma_max, K_sigma_max] z_sigma;
  array[n_re_nu] vector<lower=0>[K_nu_max] tau_nu;
  array[n_re_nu] matrix[G_nu_max, K_nu_max] z_nu;
  array[n_re_phi] vector<lower=0>[K_phi_max] tau_phi;
  array[n_re_phi] matrix[G_phi_max, K_phi_max] z_phi;
  array[n_re_alpha] vector<lower=0>[K_alpha_max] tau_alpha;
  array[n_re_alpha] matrix[G_alpha_max, K_alpha_max] z_alpha;
  array[n_re_phi_beta] vector<lower=0>[K_phi_beta_max] tau_phi_beta;
  array[n_re_phi_beta] matrix[G_phi_beta_max, K_phi_beta_max] z_phi_beta;
  array[n_re_tau_sde] vector<lower=0>[K_tau_sde_max] tau_tau_sde;
  array[n_re_tau_sde] matrix[G_tau_sde_max, K_tau_sde_max] z_tau_sde;
  real<lower=0> sigma_y;
  vector<lower=0>[D] sigma_marker;
  vector<lower=2>[D] nu_marker;
  vector<lower=0>[D] phi_nb_marker;
  vector[D] alpha_skew_marker;
  vector<lower=0>[D] phi_beta_marker;
  vector<lower=0, upper=1>[D] tau_sde_marker;

  // Ordinal cutpoints (shared across ordinal markers)
  ordered[K_ord - 1] cutpoints_ord;

  // association coefficients
  real alpha_cv_total;
  real alpha_cs_total;
  real alpha_cv_mean;
  real alpha_cs_mean;
  real alpha_cv_marker;
  real alpha_cs_marker;
  vector[Q_idm] alpha_vcov_var;

  // marker-weight shrinkage (global across ids)
  vector[D] z_marker_weights;
  real<lower=0> tau_marker_weights;

  // marker-side shrinkage scales
  real<lower=0> s_cv_marker;
  real<lower=0> s_cs_marker;
  real<lower=0> s_vcov;
