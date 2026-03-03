  /* -------------------- Longitudinal (long format) */
  int<lower=1> n_id;             // number of subjects (ids)
  int<lower=1> N;                // total longitudinal rows across all subjects
  array[N] int<lower=1, upper=n_id> id; // subject index for each row
  array[N] int<lower=1> marker;  // marker index per row (1..D)
  int<lower=1> D;                // total number of markers
  vector<lower=0>[n_id] subject_weights; // subject-level weights for longitudinal/survival likelihood
  vector[D] marker_weights;      // base weights (prior offsets) for marker-averaged CV/CS terms
  int<lower=0, upper=1> estimate_marker_weights; // 1 enables signed shrinkage perturbations around base weights
  int<lower=0, upper=1> use_marker_weight_assoc; // 1 when marker-weighted assoc terms are active

  /* Outcomes by family type */
  vector[N] y_real;              // continuous outcomes (gaussian, student_t, beta, skew families)
  array[N] int y_int;            // discrete outcomes (bernoulli, binomial, poisson, negbin2, ordinal)
  array[N] int<lower=0> trials;  // binomial trials (1 for non-binomial families)

  /* Longitudinal family per marker */
  // 1 gaussian, 2 student_t, 3 bernoulli, 4 binomial, 5 poisson, 6 negbin2,
  // 7 skew_normal, 8 double_exponential, 9 skew_double_exponential, 10 beta,
  // 11 cumulative_logit
  array[D] int<lower=1, upper=11> family_long;
  array[D] int<lower=0, upper=5> link_long; // 0 custom VM; 1 identity, 2 log, 3 logit, 4 probit, 5 exp
  int<lower=1> max_inv_link_ops;
  array[D] int<lower=0> inv_link_n_ops;
  array[D, max_inv_link_ops] int<lower=0, upper=26> inv_link_ops;
  int<lower=1> max_inv_link_const;
  array[D] int<lower=0> inv_link_n_const;
  matrix[D, max_inv_link_const] inv_link_const;

  /* Family-level distributional parameter indexing */
  int<lower=0> n_family_sigma;                    // number of families using sigma
  array[D] int<lower=0, upper=n_family_sigma> marker_to_sigma_family;
  int<lower=0> n_family_nu;                       // number of families using nu
  array[D] int<lower=0, upper=n_family_nu> marker_to_nu_family;
  int<lower=0> n_family_phi;                      // number of families using phi
  array[D] int<lower=0, upper=n_family_phi> marker_to_phi_family;
  int<lower=0> n_family_alpha;                    // number of families using alpha
  array[D] int<lower=0, upper=n_family_alpha> marker_to_alpha_family;
  int<lower=0> n_family_phi_beta;                 // number of families using phi_beta
  array[D] int<lower=0, upper=n_family_phi_beta> marker_to_phi_beta_family;
  int<lower=0> n_family_tau_sde;                  // number of families using tau_sde
  array[D] int<lower=0, upper=n_family_tau_sde> marker_to_tau_sde_family;

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

  // Beta precision (phi_beta) regression designs
  int<lower=0> P_phi_beta;             // fixed-effect columns for beta precision
  matrix[N, P_phi_beta] X_phi_beta;    // fixed-effect design for phi_beta
  int<lower=0> n_re_phi_beta;          // number of random-effect terms for phi_beta
  array[n_re_phi_beta] int<lower=1> K_phi_beta; // columns per phi_beta RE term
  array[n_re_phi_beta] int<lower=1> G_phi_beta; // grouping levels per phi_beta RE term
  int<lower=0> K_phi_beta_max;         // max columns across phi_beta RE terms
  int<lower=0> G_phi_beta_max;         // max groups across phi_beta RE terms
  array[n_re_phi_beta] matrix[N, K_phi_beta_max] Z_phi_beta; // random-effect designs
  array[n_re_phi_beta, N] int<lower=1> J_phi_beta;           // group index per row
  array[n_re_phi_beta] vector<lower=0>[G_phi_beta_max] re_weight_phi_beta; // group weights for phi_beta RE priors

  // Skew-double-exponential skew (tau_sde) regression designs
  int<lower=0> P_tau_sde;              // fixed-effect columns for tau_sde regression
  matrix[N, P_tau_sde] X_tau_sde;      // fixed-effect design for tau_sde
  int<lower=0> n_re_tau_sde;           // number of random-effect terms for tau_sde
  array[n_re_tau_sde] int<lower=1> K_tau_sde; // columns per tau_sde RE term
  array[n_re_tau_sde] int<lower=1> G_tau_sde; // grouping levels per tau_sde RE term
  int<lower=0> K_tau_sde_max;          // max columns across tau_sde RE terms
  int<lower=0> G_tau_sde_max;          // max groups across tau_sde RE terms
  array[n_re_tau_sde] matrix[N, K_tau_sde_max] Z_tau_sde; // random-effect designs
  array[n_re_tau_sde, N] int<lower=1> J_tau_sde;          // group index per row
  array[n_re_tau_sde] vector<lower=0>[G_tau_sde_max] re_weight_tau_sde; // group weights for tau_sde RE priors

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
  int<lower=0, upper=1> flag_resid_dim; // 1 -> sigma_marker[d], 0 -> sigma_y

  /* Independence/structure flags */
  int<lower=0, upper=1> indep_id_re;
  int<lower=0, upper=1> indep_marker_re;
  int<lower=0, upper=1> indep_marker_byid_latent_re;
  int<lower=0, upper=1> indep_idmarker_cov;
  int<lower=0, upper=1> allow_marker_crosscorr; // 1 allow cross via B_cross
  int<lower=0, upper=1> corr_diag_link;         // 0 softplus, 1 exp for corr diag
  int<lower=0, upper=1> use_tau_sde_fixed;      // 1 use fixed tau for SDE
  real<lower=0, upper=1> tau_sde_fixed;         // fixed tau value when enabled

  /* Hazard covariates */
  int<lower=0> p_w;              // number of baseline hazard covariates
  matrix[n_id, p_w] W;           // hazard covariate design per subject

  /* Covariance regression covariates for L_i */
  int<lower=1> K_cov;            // number of covariates in covariance regression
  matrix[n_id, K_cov] Xcov;      // covariate matrix for covariance regression

  /* Baseline hazard spline (centered basis) */
  int<lower=1> Kbs;              // number of baseline hazard basis functions
  matrix[n_id, Kbs] Bs_event_c;  // event-time basis per subject (centered)
  int<lower=1> n_gk;             // quadrature node count (nodes/weights hardcoded in Stan)
  array[n_id] matrix[n_gk, Kbs] Bs_gk_c; // quadrature-node basis per subject (centered)
  real<lower=0> tau_spline;      // spline penalty scale for baseline hazard

  /* Survival outcomes */
  vector<lower=0, upper=1>[n_id] S_event; // scaled event/censoring time in [0,1]
  array[n_id] int<lower=0, upper=1> d_event; // event indicator (1 event, 0 censored)
  int<lower=1> K_event;                    // number of competing risks
  array[n_id] int<lower=1, upper=K_event> event_type; // event type per subject

  /* Ordinal (cumulative logit) categories */
  int<lower=2> K_ord;            // number of ordinal categories (shared)

  /* Finite difference control and association design matrices */
  real<lower=1e-6> eps_fd;       // finite-difference step for slope approximation

  /* Mean association designs (CV_mean) */
  array[n_id] matrix[n_gk, P] X_gk_now;
  array[n_id] matrix[n_gk, P] X_gk_fwd;
  array[n_id] matrix[n_gk, R_id] Z_id_gk_now;
  array[n_id] matrix[n_gk, R_id] Z_id_gk_fwd;
  matrix[n_id, P] X_event_now;
  matrix[n_id, P] X_event_fwd;
  matrix[n_id, R_id] Z_id_event_now;
  matrix[n_id, R_id] Z_id_event_fwd;

  /* Marker association designs (CV_marker) */
  array[n_id] matrix[n_gk, R_mk] Z_mk_gk_now;
  array[n_id] matrix[n_gk, R_mk] Z_mk_gk_fwd;
  array[n_id] matrix[n_gk, Q_idm] Z_idm_gk_now;
  array[n_id] matrix[n_gk, Q_idm] Z_idm_gk_fwd;
  matrix[n_id, R_mk] Z_mk_event_now;
  matrix[n_id, R_mk] Z_mk_event_fwd;
  matrix[n_id, Q_idm] Z_idm_event_now;
  matrix[n_id, Q_idm] Z_idm_event_fwd;

  /* Association include flags */
  int<lower=0, upper=1> assoc_cv_total;   // include total current value term
  int<lower=0, upper=1> assoc_cv_mean;    // include mean current value term
  int<lower=0, upper=1> assoc_cv_marker;  // include marker current value term
  int<lower=0, upper=1> assoc_cs_total;   // include total current slope term
  int<lower=0, upper=1> assoc_cs_mean;    // include mean current slope term
  int<lower=0, upper=1> assoc_cs_marker;  // include marker current slope term
  int<lower=0, upper=1> assoc_corr;       // include corr association term

  /* Prior scales */
  vector[P] beta_scale;        // per-coefficient scale for fixed effects
  real<lower=0> alpha_scale;   // scale for association + baseline hazard priors
  real<lower=0> lkj_eta;       // LKJ concentration for correlation priors

  /* Marker-side association shrinkage */
  int<lower=1, upper=2> shrinkage; // 1 Laplace, 2 Gaussian

  /* Time scaling metadata */
  real<lower=0> tmax;                         // original time scale max
  int<lower=0> n_time_beta;                   // count of time columns in X_obs
  array[n_time_beta] int<lower=1, upper=P> idx_time_beta; // indices of time columns in beta
  int<lower=0> n_time_uid;                    // count of time cols in id RE design
  array[n_time_uid] int<lower=1, upper=R_id> idx_time_uid; // indices in tau_u to scale
  int<lower=0> n_time_vmk;                    // count of time cols in marker RE design
  array[n_time_vmk] int<lower=1> idx_time_vmk; // indices in tau_v to scale
  int<lower=0> n_time_widm;                   // count of time cols in marker-by-id design
  array[n_time_widm] int<lower=0, upper=Q_idm> idx_time_widm; // indices in tau_w to scale

  /* Transform configuration */
  int<lower=0, upper=3> tf_mode_cv_tot;     // transform mode for total CV
  int<lower=0, upper=3> tf_mode_cs_tot;     // transform mode for total CS
  int<lower=0, upper=3> tf_mode_corr;       // transform mode for corr term
  int<lower=0, upper=3> tf_mode_cv_mean;    // transform mode for mean CV
  int<lower=0, upper=3> tf_mode_cv_marker;  // transform mode for marker CV
  int<lower=0, upper=3> tf_mode_cs_mean;    // transform mode for mean CS
  int<lower=0, upper=3> tf_mode_cs_marker;  // transform mode for marker CS

  /* Functional bytecode specifications */
  int<lower=0> n_functional_ops_cv;         // op count for total CV
  array[n_functional_ops_cv] int<lower=0, upper=26> functional_ops_cv; // bytecode stream
  int<lower=0> n_const_cv;                  // constants used by CV bytecode
  vector[n_const_cv] const_data_cv;         // constants used by CV bytecode

  int<lower=0> n_functional_ops_cs;         // op count for total CS
  array[n_functional_ops_cs] int<lower=0, upper=26> functional_ops_cs; // bytecode stream
  int<lower=0> n_const_cs;                  // constants used by CS bytecode
  vector[n_const_cs] const_data_cs;         // constants used by CS bytecode

  int<lower=0> n_functional_ops_corr;       // op count for corr
  array[n_functional_ops_corr] int<lower=0, upper=26> functional_ops_corr; // bytecode stream
  int<lower=0> n_const_corr;                // constants used by corr bytecode
  vector[n_const_corr] const_data_corr;     // constants used by corr bytecode

  int<lower=0> n_functional_ops_cv_mean;    // op count for mean CV
  array[n_functional_ops_cv_mean] int<lower=0, upper=26> functional_ops_cv_mean; // bytecode stream
  int<lower=0> n_const_cv_mean;             // constants used by mean CV bytecode
  vector[n_const_cv_mean] const_data_cv_mean; // constants used by mean CV bytecode

  int<lower=0> n_functional_ops_cv_marker;  // op count for marker CV
  array[n_functional_ops_cv_marker] int<lower=0, upper=26> functional_ops_cv_marker; // bytecode stream
  int<lower=0> n_const_cv_marker;           // constants used by marker CV bytecode
  vector[n_const_cv_marker] const_data_cv_marker; // constants used by marker CV bytecode

  int<lower=0> n_functional_ops_cs_mean;    // op count for mean CS
  array[n_functional_ops_cs_mean] int<lower=0, upper=26> functional_ops_cs_mean; // bytecode stream
  int<lower=0> n_const_cs_mean;             // constants used by mean CS bytecode
  vector[n_const_cs_mean] const_data_cs_mean; // constants used by mean CS bytecode

  int<lower=0> n_functional_ops_cs_marker;  // op count for marker CS
  array[n_functional_ops_cs_marker] int<lower=0, upper=26> functional_ops_cs_marker; // bytecode stream
  int<lower=0> n_const_cs_marker;           // constants used by marker CS bytecode
  vector[n_const_cs_marker] const_data_cs_marker; // constants used by marker CS bytecode

  /* Spline specifications */
  int<lower=0> n_knots_cv;        // knot count for total CV spline
  vector[n_knots_cv] knots_cv;    // knot locations for total CV spline
  int<lower=0> n_coeff_cv;        // coefficient count for total CV spline
  vector[n_coeff_cv] coeff_cv;    // coefficients for total CV spline
  int<lower=1, upper=5> spline_degree_cv; // spline degree for total CV

  int<lower=0> n_knots_cs;        // knot count for total CS spline
  vector[n_knots_cs] knots_cs;    // knot locations for total CS spline
  int<lower=0> n_coeff_cs;        // coefficient count for total CS spline
  vector[n_coeff_cs] coeff_cs;    // coefficients for total CS spline
  int<lower=1, upper=5> spline_degree_cs; // spline degree for total CS

  int<lower=0> n_knots_corr;      // knot count for corr spline
  vector[n_knots_corr] knots_corr; // knot locations for corr spline
  int<lower=0> n_coeff_corr;      // coefficient count for corr spline
  vector[n_coeff_corr] coeff_corr; // coefficients for corr spline
  int<lower=1, upper=5> spline_degree_corr; // spline degree for corr

  int<lower=0> n_knots_cv_mean;   // knot count for mean CV spline
  vector[n_knots_cv_mean] knots_cv_mean; // knot locations for mean CV spline
  int<lower=0> n_coeff_cv_mean;   // coefficient count for mean CV spline
  vector[n_coeff_cv_mean] coeff_cv_mean; // coefficients for mean CV spline
  int<lower=1, upper=5> spline_degree_cv_mean; // spline degree for mean CV

  int<lower=0> n_knots_cv_marker; // knot count for marker CV spline
  vector[n_knots_cv_marker] knots_cv_marker; // knot locations for marker CV spline
  int<lower=0> n_coeff_cv_marker; // coefficient count for marker CV spline
  vector[n_coeff_cv_marker] coeff_cv_marker; // coefficients for marker CV spline
  int<lower=1, upper=5> spline_degree_cv_marker; // spline degree for marker CV

  int<lower=0> n_knots_cs_mean;   // knot count for mean CS spline
  vector[n_knots_cs_mean] knots_cs_mean; // knot locations for mean CS spline
  int<lower=0> n_coeff_cs_mean;   // coefficient count for mean CS spline
  vector[n_coeff_cs_mean] coeff_cs_mean; // coefficients for mean CS spline
  int<lower=1, upper=5> spline_degree_cs_mean; // spline degree for mean CS

  int<lower=0> n_knots_cs_marker; // knot count for marker CS spline
  vector[n_knots_cs_marker] knots_cs_marker; // knot locations for marker CS spline
  int<lower=0> n_coeff_cs_marker; // coefficient count for marker CS spline
  vector[n_coeff_cs_marker] coeff_cs_marker; // coefficients for marker CS spline
  int<lower=1, upper=5> spline_degree_cs_marker; // spline degree for marker CS
