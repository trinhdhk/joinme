/**
 * @file longitudinal/parameters/fit.stan
 * @brief Parameters for longitudinal population, distributional and nested random effects.
 *
 * @details
 * The fragment contains non-centred population coefficients, subject, marker and subject-by-marker effects, covariance-regression coefficients and latents, distributional regressions, ordinal cutpoints and the marker-effect horseshoe hierarchy.
 */
  /* Fixed effects */
  vector[P] z_beta; // standardised longitudinal population coefficients, transformed by the common prior programme

  /* id-level random effects */
  array[n_id] vector[R_id] z_u;          // latent standard normals per subject
  vector<lower=0>[R_id] tau_u;           // marginal SDs for id REs
  cholesky_factor_corr[R_id] Lcorr_u;    // Cholesky corr for id REs

  /* marker-only random effects */
  array[D] vector[R_mk] z_v;             // latent standard normals per marker
  vector<lower=0>[R_mk] tau_v;           // marginal SDs for marker REs
  cholesky_factor_corr[R_mk] Lcorr_v;    // Cholesky corr for marker REs

  /* marker-by-id latent effects */
  array[n_id, D] vector[Q_idm] z_w_lat;  // iid standard normals per id, marker; baseline covariance now lives in L_i

  /* Cross-correlation mapping between v_marker and z_w */
  matrix[Q_idm, R_mk] B_cross;           // cross-loadings when crosscorr enabled

  /* id-specific covariance regression for L_i */
  vector[P_vcov_sd] vcov_sd_coefficient_raw; // standardised SD-regression intercepts followed by row-major slopes
  vector[P_vcov_corr] vcov_corr_coefficient_raw; // standardised off-diagonal correlation intercepts followed by row-major slopes
  vector<lower=0>[M_cov] lambda_L;       // one nonnegative scalar loading per packed covariance coordinate; collectively a diagonal, not full, loading matrix
  array[n_id] vector[M_cov] z_L;         // latent standard normals for covariance regression, one per subject and L_i element

  /* longitudinal distributional parameters */
  vector[P_sigma] z_beta_sigma; // standardised sigma-regression coefficients
  vector[P_nu] z_beta_nu; // standardised degrees-of-freedom-regression coefficients
  vector[P_phi] z_beta_phi; // standardised dispersion-regression coefficients
  vector[P_alpha] z_beta_alpha; // standardised skewness-regression coefficients
  vector[P_kappa] z_beta_kappa; // standardised beta-sample-size-regression coefficients
  vector[P_tau] z_beta_tau; // standardised quantile-regression coefficients
  array[n_re_sigma] vector<lower=0>[K_sigma_max] tau_sigma; // sigma RE SDs
  array[n_re_sigma] matrix[G_sigma_max, K_sigma_max] z_sigma; // sigma RE latents
  array[n_re_nu] vector<lower=0>[K_nu_max] tau_nu; // nu RE SDs
  array[n_re_nu] matrix[G_nu_max, K_nu_max] z_nu;  // nu RE latents
  array[n_re_phi] vector<lower=0>[K_phi_max] tau_phi; // phi RE SDs
  array[n_re_phi] matrix[G_phi_max, K_phi_max] z_phi;  // phi RE latents
  array[n_re_alpha] vector<lower=0>[K_alpha_max] tau_alpha; // alpha RE SDs
  array[n_re_alpha] matrix[G_alpha_max, K_alpha_max] z_alpha; // alpha RE latents
  array[n_re_kappa] vector<lower=0>[K_kappa_max] tau_kappa; // kappa RE SDs
  array[n_re_kappa] matrix[G_kappa_max, K_kappa_max] z_kappa; // kappa RE latents
  array[n_re_tau] vector<lower=0>[K_tau_max] tau_tau; // tau RE SDs
  array[n_re_tau] matrix[G_tau_max, K_tau_max] z_tau; // tau RE latents
  vector<lower=0>[n_family_sigma] sigma_family;        // shared sigma by family
  vector<lower=2>[n_family_nu] nu_family;              // shared nu by family
  vector<lower=0>[n_family_phi] phi_family;            // shared phi by family
  vector[n_family_alpha] alpha_family;                 // shared alpha by family
  vector<lower=0>[n_family_kappa] kappa_family;  // shared kappa by family
  vector<lower=0, upper=1>[n_family_tau] tau_family; // shared tau by family

  /* Ordinal cutpoints (shared across ordinal markers) */
  ordered[K_ord - 1] cutpoints_ord;     // cumulative-logit cutpoints

  vector<lower=0>[prior_marker_family == 4 ? D * R_mk : 0] horseshoe_local_marker; // local scales for marker random effects
  vector<lower=0>[prior_marker_family == 4 ? 1 : 0] horseshoe_global_marker; // global scale for marker random effects
  vector<lower=0>[prior_marker_family == 4 ? 1 : 0] horseshoe_slab_marker; // slab auxiliary for marker random effects


