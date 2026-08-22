/**
 * @file longitudinal/data/fit.stan
 * @brief Data declarations for the multivariate longitudinal submodel.
 *
 * @details
 * This fragment declares the observed marker outcomes, marker-specific
 * families and inverse links, population and nested random-effect designs,
 * distributional regressions, covariance regression, ordinal categories,
 * correlation priors and survival time-scale metadata. Subject and event indexing are
 * introduced here because the longitudinal history supplies the shared subject
 * axis on which the event process and dynamic prediction are conditioned.
 *
 * Include this fragment before the other fit-data submodels: later fragments
 * use n_id, N, D, P and the random-effect dimensions declared here.
 */
  /* -------------------- Longitudinal (long format) */
  int<lower=1> n_id;             // number of subjects (ids)
  int<lower=1> N_event;          // number of event-process rows (intervals)
  int<lower=1> N;                // total longitudinal rows across all subjects
  array[N] int<lower=1, upper=n_id> id; // subject index for each row
  array[N_event] int<lower=1, upper=n_id> event_id; // subject index for each event interval row
  array[n_id] int<lower=1, upper=N_event> event_start_idx; // first interval row per subject
  array[n_id] int<lower=1, upper=N_event> event_end_idx;   // last interval row per subject
  array[N] int<lower=1> marker;  // marker index per row (1..D)
  int<lower=1> D;                // total number of markers
  vector<lower=0>[n_id] subject_weights; // subject-level weights for longitudinal/survival likelihood
  /* Outcomes by family type */
  vector[N] y_real;              // continuous outcomes (gaussian, student_t, beta, skew families)
  array[N] int y_int;            // discrete outcomes (bernoulli, binomial, poisson, negbin2, ordinal)
  array[N] int<lower=0> trials;  // binomial trials (1 for non-binomial families)

  /* Longitudinal family per marker */
  // 1 gaussian, 2 student_t, 3 bernoulli, 4 binomial, 5 poisson, 6 negbin2,
  // 7 skew_normal, 8 double_exponential, 9 skew_double_exponential, 10 beta,
  // 11 cumulative_logit
  array[D] int<lower=1, upper=11> family_long; // Role: family long.
  array[D] int<lower=0, upper=5> link_long; // 0 custom VM; 1 identity, 2 log, 3 logit, 4 probit, 5 exp
  int<lower=1> max_inv_link_ops; // Role: maximum inverse link operations.
  array[D] int<lower=0> inv_link_n_ops; // Role: inverse link number of operations.
  array[D, max_inv_link_ops] int<lower=0, upper=27> inv_link_ops; // Role: inverse link operations.
  int<lower=1> max_inv_link_const; // Role: maximum inverse link constants.
  array[D] int<lower=0> inv_link_n_const; // Role: inverse link number of constants.
  matrix[D, max_inv_link_const] inv_link_const; // Role: inverse link constants.

  /* Family-level distributional parameter indexing */
  int<lower=0> n_family_sigma;                    // number of families using sigma
  array[D] int<lower=0, upper=n_family_sigma> marker_to_sigma_family; // Role: marker to scale family.
  int<lower=0> n_family_nu;                       // number of families using nu
  array[D] int<lower=0, upper=n_family_nu> marker_to_nu_family; // Role: marker to degrees of freedom family.
  int<lower=0> n_family_phi;                      // number of families using phi
  array[D] int<lower=0, upper=n_family_phi> marker_to_phi_family; // Role: marker to precision family.
  int<lower=0> n_family_alpha;                    // number of families using alpha
  array[D] int<lower=0, upper=n_family_alpha> marker_to_alpha_family; // Role: marker to shape family.
  int<lower=0> n_family_kappa;                 // number of families using kappa
  array[D] int<lower=0, upper=n_family_kappa> marker_to_kappa_family; // Role: marker to sample size family.
  int<lower=0> n_family_tau;                  // number of families using tau
  array[D] int<lower=0, upper=n_family_tau> marker_to_tau_family; // Role: marker to quantile family.

  /* Distributional regression designs */
  int<lower=0> P_sigma;                // fixed-effect column count for sigma regression
  matrix[N, P_sigma] X_sigma;          // fixed-effect design for sigma
  int<lower=0> n_re_sigma;             // number of random-effect terms for sigma
  array[n_re_sigma] int<lower=1> K_sigma; // columns per random-effect term
  array[n_re_sigma] int<lower=1> G_sigma; // grouping levels per random-effect term
  int<lower=0> K_sigma_max;            // max columns across sigma RE terms
  int<lower=0> G_sigma_max;            // max groups across sigma RE terms
  array[n_re_sigma] matrix[N, K_sigma_max] Z_sigma; // random-effect design matrices
  array[n_re_sigma, N] int<lower=1> J_sigma;        // group index per row
  array[n_re_sigma] vector<lower=0>[G_sigma_max] re_weight_sigma; // group weights for sigma RE priors

  int<lower=0> P_nu;                   // fixed-effect columns for nu regression
  matrix[N, P_nu] X_nu;                // fixed-effect design for nu
  int<lower=0> n_re_nu;                // number of random-effect terms for nu
  array[n_re_nu] int<lower=1> K_nu;    // columns per nu RE term
  array[n_re_nu] int<lower=1> G_nu;    // grouping levels per nu RE term
  int<lower=0> K_nu_max;               // max columns across nu RE terms
  int<lower=0> G_nu_max;               // max groups across nu RE terms
  array[n_re_nu] matrix[N, K_nu_max] Z_nu; // random-effect designs
  array[n_re_nu, N] int<lower=1> J_nu;      // group index per row
  array[n_re_nu] vector<lower=0>[G_nu_max] re_weight_nu; // group weights for nu RE priors

  int<lower=0> P_phi;                  // fixed-effect columns for phi regression
  matrix[N, P_phi] X_phi;              // fixed-effect design for phi
  int<lower=0> n_re_phi;               // number of random-effect terms for phi
  array[n_re_phi] int<lower=1> K_phi;  // columns per phi RE term
  array[n_re_phi] int<lower=1> G_phi;  // grouping levels per phi RE term
  int<lower=0> K_phi_max;              // max columns across phi RE terms
  int<lower=0> G_phi_max;              // max groups across phi RE terms
  array[n_re_phi] matrix[N, K_phi_max] Z_phi; // random-effect designs
  array[n_re_phi, N] int<lower=1> J_phi;      // group index per row
  array[n_re_phi] vector<lower=0>[G_phi_max] re_weight_phi; // group weights for phi RE priors

  int<lower=0> P_alpha;                // fixed-effect columns for alpha (skew) regression
  matrix[N, P_alpha] X_alpha;          // fixed-effect design for alpha
  int<lower=0> n_re_alpha;             // number of random-effect terms for alpha
  array[n_re_alpha] int<lower=1> K_alpha; // columns per alpha RE term
  array[n_re_alpha] int<lower=1> G_alpha; // grouping levels per alpha RE term
  int<lower=0> K_alpha_max;            // max columns across alpha RE terms
  int<lower=0> G_alpha_max;            // max groups across alpha RE terms
  array[n_re_alpha] matrix[N, K_alpha_max] Z_alpha; // random-effect designs
  array[n_re_alpha, N] int<lower=1> J_alpha;        // group index per row
  array[n_re_alpha] vector<lower=0>[G_alpha_max] re_weight_alpha; // group weights for alpha RE priors

  // Beta sample-size (kappa) regression designs
  int<lower=0> P_kappa;             // fixed-effect columns for beta sample size
  matrix[N, P_kappa] X_kappa;    // fixed-effect design for kappa
  int<lower=0> n_re_kappa;          // number of random-effect terms for kappa
  array[n_re_kappa] int<lower=1> K_kappa; // columns per kappa RE term
  array[n_re_kappa] int<lower=1> G_kappa; // grouping levels per kappa RE term
  int<lower=0> K_kappa_max;         // max columns across kappa RE terms
  int<lower=0> G_kappa_max;         // max groups across kappa RE terms
  array[n_re_kappa] matrix[N, K_kappa_max] Z_kappa; // random-effect designs
  array[n_re_kappa, N] int<lower=1> J_kappa;           // group index per row
  array[n_re_kappa] vector<lower=0>[G_kappa_max] re_weight_kappa; // group weights for kappa RE priors

  // Skew-double-exponential skew (tau) regression designs
  int<lower=0> P_tau;              // fixed-effect columns for tau regression
  matrix[N, P_tau] X_tau;      // fixed-effect design for tau
  int<lower=0> n_re_tau;           // number of random-effect terms for tau
  array[n_re_tau] int<lower=1> K_tau; // columns per tau RE term
  array[n_re_tau] int<lower=1> G_tau; // grouping levels per tau RE term
  int<lower=0> K_tau_max;          // max columns across tau RE terms
  int<lower=0> G_tau_max;          // max groups across tau RE terms
  array[n_re_tau] matrix[N, K_tau_max] Z_tau; // random-effect designs
  array[n_re_tau, N] int<lower=1> J_tau;          // group index per row
  array[n_re_tau] vector<lower=0>[G_tau_max] re_weight_tau; // group weights for tau RE priors

  /* Fixed effects design */
  int<lower=1> P;                // columns in fixed-effect design
  matrix[N, P] X_obs;            // fixed-effect design matrix for longitudinal outcomes

  /* Random effects designs */
  int<lower=1> R_id;             // number of id-level random effects
  matrix[N, R_id] Z_id_obs;      // id-level random-effect design
  vector<lower=0>[n_id] re_weight_id; // group weights for id random-effects priors
  int<lower=0> R_mk;             // number of marker-only random effects
  matrix[N, R_mk] Z_mk_obs;      // marker-only random-effect design
  vector<lower=0>[D] re_weight_marker; // group weights for marker random-effects priors
  int<lower=0> Q_idm;            // number of marker-by-id random effects
  matrix[N, Q_idm] Z_idm_obs;    // marker-by-id random-effect design
  vector<lower=0>[n_id] re_weight_idm; // group weights for marker-by-id latent RE priors
  vector<lower=0>[n_id] re_weight_L;   // group weights for covariance latent priors

  /* Residual SD option for continuous families */
  // int<lower=0, upper=1> flag_resid_dim; // 1 -> sigma_marker[d], 0 -> sigma_y

  /* Independence/structure flags */
  int<lower=0, upper=1> indep_id_re; // Role: independence of individual random effect.
  int<lower=0, upper=1> indep_marker_re; // Role: independence of marker random effect.
  int<lower=0, upper=1> indep_marker_byid_latent_re; // Role: independence of marker byid latent random effect.
  int<lower=0, upper=1> indep_idmarker_cov; // Role: independence of idmarker covariance.
  int<lower=0> M_corr_tf;                    // number of corr transform components
  int<lower=0> M_vcov_tf;                    // number of vcov transform components
  int<lower=0, upper=1> allow_marker_crosscorr; // 1 allow cross via B_cross
  int<lower=0, upper=1> vcov_diag_link;         // 0 softplus, 1 exp for covariance diagonals
  array[D] int<lower=0, upper=1> use_tau_fixed; // marker d uses its family-fixed tau
  vector<lower=0, upper=1>[D] tau_fixed;        // marker-specific fixed quantiles

  /* Independent covariance regressions for the SD and correlation components of L_i */
  int<lower=0> K_cov_sd; // number of observed predictors in the standard-deviation regression
  matrix[n_id, K_cov_sd] Xcov_sd; // subject-aligned predictors for covariance standard deviations
  int<lower=0> K_cov_corr; // number of observed predictors in the off-diagonal partial-correlation regression
  matrix[n_id, K_cov_corr] Xcov_corr; // subject-aligned predictors for off-diagonal partial correlations

  /* Ordinal (cumulative logit) categories */
  int<lower=2> K_ord;            // number of ordinal categories (shared)

  /* Correlation prior */
  real<lower=0> lkj_eta;       // LKJ concentration shared by random-effect correlation matrices
  int<lower=1, upper=4> prior_marker_family; // family of standardised marker random effects
  real<lower=0> prior_marker_df; // fixed Student-t degrees of freedom for marker effects
  real<lower=0> prior_marker_global_df; // unit-default horseshoe global degrees of freedom for marker effects
  real<lower=0> prior_marker_global_scale; // configured horseshoe global regularisation scale for marker effects
  real<lower=0> prior_marker_slab_df; // fixed horseshoe slab degrees of freedom for marker effects
  real<lower=0> prior_marker_slab_scale; // configured finite-slab scale for marker effects

  int<lower=0> P_vcov_sd; // packed count: Q_idm SD intercepts, then row-major SD slopes
  int<lower=0> P_vcov_corr; // packed count: off-diagonal intercepts, then row-major correlation slopes

  /* Survival integration uses tmax to map [0,1] quadrature increments back to
   * original time when calculating current slopes. Longitudinal designs and
   * their coefficients already use original time and need no index metadata. */
  real<lower=0> tmax; // maximum event time converting scaled integration differences to original-time derivatives

