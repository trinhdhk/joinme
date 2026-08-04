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
  array[n_id, D] vector[Q_idm] z_w_lat;  // iid standard normals per id, marker; baseline covariance now lives in L_i

  /* Cross-correlation mapping between v_marker and z_w */
  matrix[Q_idm, R_mk] B_cross;           // cross-loadings when crosscorr enabled

  /* id-specific covariance regression for L_i */
  vector[P_vcov_sd] vcov_sd_coefficient_raw; // standardised SD-regression intercepts followed by row-major slopes
  vector[P_vcov_corr] vcov_corr_coefficient_raw; // standardised off-diagonal correlation intercepts followed by row-major slopes
  vector<lower=0>[prior_vcov_sd_family == 4 ? P_vcov_sd : 0] horseshoe_local_vcov_sd; // local regularised-horseshoe scales for SD coefficients
  vector<lower=0>[prior_vcov_sd_family == 4 ? 1 : 0] horseshoe_global_vcov_sd; // global regularised-horseshoe scale for SD coefficients
  vector<lower=0>[prior_vcov_sd_family == 4 ? 1 : 0] horseshoe_slab_vcov_sd; // finite-slab multiplier for SD coefficients
  vector<lower=0>[prior_vcov_corr_family == 4 ? P_vcov_corr : 0] horseshoe_local_vcov_corr; // local regularised-horseshoe scales for correlation coefficients
  vector<lower=0>[prior_vcov_corr_family == 4 ? 1 : 0] horseshoe_global_vcov_corr; // global regularised-horseshoe scale for correlation coefficients
  vector<lower=0>[prior_vcov_corr_family == 4 ? 1 : 0] horseshoe_slab_vcov_corr; // finite-slab multiplier for correlation coefficients
  vector<lower=0>[M_cov] lambda_L;       // one nonnegative scalar loading per packed covariance coordinate; collectively a diagonal, not full, loading matrix
  array[n_id] vector[M_cov] z_L;         // latent standard normals for covariance regression, one per subject and L_i element

  /* baseline hazard (cause-specific) */
  array[K_event] vector[Kbs] bs_gamma_c; // spline coefficients for baseline hazard (includes intercept basis)

  /* hazard covariates (cause-specific) */
  array[K_event] vector[p_w] gamma_w;    // regression coefficients for W

  /* longitudinal distributional parameters */
  vector[P_sigma] beta_sigma;            // sigma regression coefficients
  vector[P_nu] beta_nu;                  // nu regression coefficients
  vector[P_phi] beta_phi;                // phi regression coefficients
  vector[P_alpha] beta_alpha;            // alpha regression coefficients
  vector[P_kappa] beta_kappa;      // kappa regression coefficients
  vector[P_tau] beta_tau;        // tau regression coefficients
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

  /* association coefficients (non-centred for totals/mean/marker) */
  real z_alpha_cv_total;               // latent for total CV association
  real z_alpha_cs_total;               // latent for total CS association
  real z_alpha_cv_mean;               // latent for mean CV association
  real z_alpha_cs_mean;               // latent for mean CS association
  real z_alpha_cv_marker;              // latent for marker CV association
  real z_alpha_cs_marker;              // latent for marker CS association
  vector[M_corr] z_alpha_corr;        // standardized corr association coefficient latents
  vector[M_vcov] z_alpha_vcov;        // standardized vcov association coefficient latents

  /* marker-weight perturbations (global across ids) */
  vector[D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets] z_marker_weights; // latent signed perturbations for shared or term-specific weight structures

  /**
   * @brief Latent increments for Stan-estimated monotone I-splines.
   *
   * @details These parameters preserve the established I-spline
   * parameterisation and its shrinkage-family prior.  Piecewise-linear modes
   * omit them because those modes use the explicit simplexes declared below.
   */
  vector[(tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) ? 0 : n_free_spline_cv] z_spline_cv; // Role: standardised latent value spline current value.
  vector[(tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) ? 0 : n_free_spline_cs] z_spline_cs; // Role: standardised latent value spline current slope.
  matrix[M_corr, (tf_mode_corr == 3 || tf_mode_corr == 7) ? 0 : n_free_spline_corr] z_spline_corr; // Role: standardised latent value spline correlation.
  matrix[M_vcov, (tf_mode_vcov == 3 || tf_mode_vcov == 7) ? 0 : n_free_spline_vcov] z_spline_vcov; // Role: standardised latent value spline covariance.
  vector[(tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) ? 0 : n_free_spline_cv_mean] z_spline_cv_mean; // Role: standardised latent value spline current value mean.
  vector[(tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) ? 0 : n_free_spline_cv_marker] z_spline_cv_marker; // Role: standardised latent value spline current value marker.
  vector[(tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) ? 0 : n_free_spline_cs_mean] z_spline_cs_mean; // Role: standardised latent value spline current slope mean.
  vector[(tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) ? 0 : n_free_spline_cs_marker] z_spline_cs_marker; // Role: standardised latent value spline current slope marker.

  /**
   * @brief Simplex increments for ordered piecewise-linear associations.
   *
   * @details For K ordinates, a K-1 simplex allocates the unit association span
   * between adjacent knots. Cumulative sums in transformed parameters recover
   * the ordered ordinates. An inactive channel has a one-element simplex,
   * which has no free parameter and permits one compiled model to cover every
   * transform choice without changing existing I-spline parameters.
   */
  simplex[(tf_mode_cv_tot == 3 || tf_mode_cv_tot == 7) ? n_free_spline_cv : 1] pwlin_simplex_cv; // Role: piecewise-linear simplex current value.
  simplex[(tf_mode_cs_tot == 3 || tf_mode_cs_tot == 7) ? n_free_spline_cs : 1] pwlin_simplex_cs; // Role: piecewise-linear simplex current slope.
  array[M_corr] simplex[(tf_mode_corr == 3 || tf_mode_corr == 7) ? n_free_spline_corr : 1] pwlin_simplex_corr; // Role: piecewise-linear simplex correlation.
  array[M_vcov] simplex[(tf_mode_vcov == 3 || tf_mode_vcov == 7) ? n_free_spline_vcov : 1] pwlin_simplex_vcov; // Role: piecewise-linear simplex covariance.
  simplex[(tf_mode_cv_mean == 3 || tf_mode_cv_mean == 7) ? n_free_spline_cv_mean : 1] pwlin_simplex_cv_mean; // Role: piecewise-linear simplex current value mean.
  simplex[(tf_mode_cv_marker == 3 || tf_mode_cv_marker == 7) ? n_free_spline_cv_marker : 1] pwlin_simplex_cv_marker; // Role: piecewise-linear simplex current value marker.
  simplex[(tf_mode_cs_mean == 3 || tf_mode_cs_mean == 7) ? n_free_spline_cs_mean : 1] pwlin_simplex_cs_mean; // Role: piecewise-linear simplex current slope mean.
  simplex[(tf_mode_cs_marker == 3 || tf_mode_cs_marker == 7) ? n_free_spline_cs_marker : 1] pwlin_simplex_cs_marker; // Role: piecewise-linear simplex current slope marker.

  /* fit-only affine-shift parameters for functional transforms */
  vector[estimate_iota_intercept_cv] z_iota_intercept_cv; // Role: standardised latent value affine transformation intercept current value.
  vector[estimate_iota_slope_cv] z_iota_slope_cv; // Role: standardised latent value affine transformation slope current value.
  vector[estimate_iota_intercept_cs] z_iota_intercept_cs; // Role: standardised latent value affine transformation intercept current slope.
  vector[estimate_iota_slope_cs] z_iota_slope_cs; // Role: standardised latent value affine transformation slope current slope.
  vector[M_corr * estimate_iota_intercept_corr] z_iota_intercept_corr; // Role: standardised latent value affine transformation intercept correlation.
  vector[M_corr * estimate_iota_slope_corr] z_iota_slope_corr; // Role: standardised latent value affine transformation slope correlation.
  vector[M_vcov * estimate_iota_intercept_vcov] z_iota_intercept_vcov; // Role: standardised latent value affine transformation intercept covariance.
  vector[M_vcov * estimate_iota_slope_vcov] z_iota_slope_vcov; // Role: standardised latent value affine transformation slope covariance.
  vector[estimate_iota_intercept_cv_mean] z_iota_intercept_cv_mean; // Role: standardised latent value affine transformation intercept current value mean.
  vector[estimate_iota_slope_cv_mean] z_iota_slope_cv_mean; // Role: standardised latent value affine transformation slope current value mean.
  vector[estimate_iota_intercept_cv_marker] z_iota_intercept_cv_marker; // Role: standardised latent value affine transformation intercept current value marker.
  vector[estimate_iota_slope_cv_marker] z_iota_slope_cv_marker; // Role: standardised latent value affine transformation slope current value marker.
  vector[estimate_iota_intercept_cs_mean] z_iota_intercept_cs_mean; // Role: standardised latent value affine transformation intercept current slope mean.
  vector[estimate_iota_slope_cs_mean] z_iota_slope_cs_mean; // Role: standardised latent value affine transformation slope current slope mean.
  vector[estimate_iota_intercept_cs_marker] z_iota_intercept_cs_marker; // Role: standardised latent value affine transformation intercept current slope marker.
  vector[estimate_iota_slope_cs_marker] z_iota_slope_cs_marker; // Role: standardised latent value affine transformation slope current slope marker.

  /**
   * Regularised-horseshoe auxiliaries.
   *
   * Each vector has positive length only when its block selects family code 4.
   * Zero-length declarations keep all other prior families free of unused
   * parameters. The slab auxiliary has an inverse-gamma prior and represents
   * the random squared width of the finite regularising slab.
   */
  vector<lower=0>[prior_beta_family == 4 ? P : 0] horseshoe_local_beta; // local shrinkage scales for beta
  vector<lower=0>[prior_beta_family == 4 ? 1 : 0] horseshoe_global_beta; // global shrinkage scale for beta
  vector<lower=0>[prior_beta_family == 4 ? 1 : 0] horseshoe_slab_beta; // slab-variance auxiliary for beta

  vector<lower=0>[prior_alpha_family == 4 ? 6 + M_corr + M_vcov : 0] horseshoe_local_alpha; // local scales for association coefficients
  vector<lower=0>[prior_alpha_family == 4 ? 1 : 0] horseshoe_global_alpha; // global scale for association coefficients
  vector<lower=0>[prior_alpha_family == 4 ? 1 : 0] horseshoe_slab_alpha; // slab auxiliary for association coefficients

  vector<lower=0>[prior_iota_family == 4 ? n_iota_prior : 0] horseshoe_local_iota; // local scales for fitted affine shifts
  vector<lower=0>[prior_iota_family == 4 ? 1 : 0] horseshoe_global_iota; // global scale for fitted affine shifts
  vector<lower=0>[prior_iota_family == 4 ? 1 : 0] horseshoe_slab_iota; // slab auxiliary for fitted affine shifts

  vector<lower=0>[prior_marker_family == 4 ? D * R_mk : 0] horseshoe_local_marker; // local scales for marker random effects
  vector<lower=0>[prior_marker_family == 4 ? 1 : 0] horseshoe_global_marker; // global scale for marker random effects
  vector<lower=0>[prior_marker_family == 4 ? 1 : 0] horseshoe_slab_marker; // slab auxiliary for marker random effects

  vector<lower=0>[prior_marker_weight_family == 4 ? D * estimate_marker_weights * use_marker_weight_assoc * n_marker_weight_sets : 0] horseshoe_local_marker_weight; // local scales for marker weights
  vector<lower=0>[prior_marker_weight_family == 4 ? 1 : 0] horseshoe_global_marker_weight; // global scale for marker weights
  vector<lower=0>[prior_marker_weight_family == 4 ? 1 : 0] horseshoe_slab_marker_weight; // slab auxiliary for marker weights
