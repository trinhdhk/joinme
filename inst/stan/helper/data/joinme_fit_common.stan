  // -------------------- Longitudinal (long format)
  int<lower=1> n_id;             // number of subjects
  int<lower=1> N;                // total longitudinal observations
  array[N] int<lower=1, upper=n_id> id;
  array[N] int<lower=1> marker;  // marker index 1..D
  int<lower=1> D;                // number of markers
  vector[D] marker_weights;      // normalized weights for marker-level averaging (CV + CS)
  int<lower=0, upper=1> estimate_marker_weights; // 1 to estimate weights
  real<lower=0> marker_weight_scale; // prior scale for weight shrinkage

  // Outcomes by family type
  vector[N] y_real;              // gaussian/student_t
  array[N] int y_int;            // bernoulli/binomial/poisson/negbin2
  array[N] int<lower=0> trials;  // binomial trials (1 otherwise)

  // Longitudinal family per marker
  // 1 gaussian, 2 student_t, 3 bernoulli, 4 binomial, 5 poisson, 6 negbin2,
  // 7 skew_normal, 8 double_exponential, 9 skew_double_exponential, 10 beta,
  // 11 cumulative_logit
  array[D] int<lower=1, upper=11> family_long;

  // Distributional regression designs
  int<lower=0> P_sigma;
  matrix[N, P_sigma] X_sigma;
  int<lower=0> n_re_sigma;
  array[n_re_sigma] int<lower=1> K_sigma;
  array[n_re_sigma] int<lower=1> G_sigma;
  int<lower=0> K_sigma_max;
  int<lower=0> G_sigma_max;
  array[n_re_sigma] matrix[N, K_sigma_max] Z_sigma;
  array[n_re_sigma, N] int<lower=1> J_sigma;

  int<lower=0> P_nu;
  matrix[N, P_nu] X_nu;
  int<lower=0> n_re_nu;
  array[n_re_nu] int<lower=1> K_nu;
  array[n_re_nu] int<lower=1> G_nu;
  int<lower=0> K_nu_max;
  int<lower=0> G_nu_max;
  array[n_re_nu] matrix[N, K_nu_max] Z_nu;
  array[n_re_nu, N] int<lower=1> J_nu;

  int<lower=0> P_phi;
  matrix[N, P_phi] X_phi;
  int<lower=0> n_re_phi;
  array[n_re_phi] int<lower=1> K_phi;
  array[n_re_phi] int<lower=1> G_phi;
  int<lower=0> K_phi_max;
  int<lower=0> G_phi_max;
  array[n_re_phi] matrix[N, K_phi_max] Z_phi;
  array[n_re_phi, N] int<lower=1> J_phi;

  int<lower=0> P_alpha;
  matrix[N, P_alpha] X_alpha;
  int<lower=0> n_re_alpha;
  array[n_re_alpha] int<lower=1> K_alpha;
  array[n_re_alpha] int<lower=1> G_alpha;
  int<lower=0> K_alpha_max;
  int<lower=0> G_alpha_max;
  array[n_re_alpha] matrix[N, K_alpha_max] Z_alpha;
  array[n_re_alpha, N] int<lower=1> J_alpha;

  // Beta precision (phi_beta) regression designs
  int<lower=0> P_phi_beta;
  matrix[N, P_phi_beta] X_phi_beta;
  int<lower=0> n_re_phi_beta;
  array[n_re_phi_beta] int<lower=1> K_phi_beta;
  array[n_re_phi_beta] int<lower=1> G_phi_beta;
  int<lower=0> K_phi_beta_max;
  int<lower=0> G_phi_beta_max;
  array[n_re_phi_beta] matrix[N, K_phi_beta_max] Z_phi_beta;
  array[n_re_phi_beta, N] int<lower=1> J_phi_beta;

  // Skew-double-exponential skew (tau_sde) regression designs
  int<lower=0> P_tau_sde;
  matrix[N, P_tau_sde] X_tau_sde;
  int<lower=0> n_re_tau_sde;
  array[n_re_tau_sde] int<lower=1> K_tau_sde;
  array[n_re_tau_sde] int<lower=1> G_tau_sde;
  int<lower=0> K_tau_sde_max;
  int<lower=0> G_tau_sde_max;
  array[n_re_tau_sde] matrix[N, K_tau_sde_max] Z_tau_sde;
  array[n_re_tau_sde, N] int<lower=1> J_tau_sde;

  // Fixed effects design
  int<lower=1> P;
  matrix[N, P] X_obs;

  // Random effects designs
  int<lower=1> R_id;
  matrix[N, R_id] Z_id_obs;
  int<lower=0> R_mk;
  matrix[N, R_mk] Z_mk_obs;
  int<lower=0> Q_idm;
  matrix[N, Q_idm] Z_idm_obs;

  // Residual SD option for continuous families
  int<lower=0, upper=1> flag_resid_dim; // 1 -> sigma_marker[d], 0 -> sigma_y

  // Independence/structure flags
  int<lower=0, upper=1> indep_id_re;
  int<lower=0, upper=1> indep_marker_re;
  int<lower=0, upper=1> indep_marker_byid_latent_re;
  int<lower=0, upper=1> indep_idmarker_cov;
  int<lower=0, upper=1> allow_marker_crosscorr; // 1 allow cross via B_cross

  // Hazard covariates
  int<lower=0> p_w;
  matrix[n_id, p_w] W;

  // Covariance regression covariates for L_i
  int<lower=1> K_cov;
  matrix[n_id, K_cov] Xcov;

  // Baseline hazard spline (centered basis)
  int<lower=1> Kbs;
  matrix[n_id, Kbs] Bs_event_c;
  array[n_id] matrix[15, Kbs] Bs_gk_c;
  real<lower=0> tau_spline;

  // Survival outcomes
  vector<lower=0, upper=1>[n_id] S_event; // scaled survival times in [0,1]
  array[n_id] int<lower=0, upper=1> d_event;
  int<lower=1> K_event;                    // number of event types / risks
  array[n_id] int<lower=1, upper=K_event> event_type;

  // Ordinal (cumulative logit) categories
  int<lower=2> K_ord;

  // Finite difference control and association design matrices
  real<lower=1e-6> eps_fd;

  // Mean association designs (CV_mean)
  array[n_id] matrix[15, P] X_gk_now;
  array[n_id] matrix[15, P] X_gk_fwd;
  array[n_id] matrix[15, R_id] Z_id_gk_now;
  array[n_id] matrix[15, R_id] Z_id_gk_fwd;
  matrix[n_id, P] X_event_now;
  matrix[n_id, P] X_event_fwd;
  matrix[n_id, R_id] Z_id_event_now;
  matrix[n_id, R_id] Z_id_event_fwd;

  // Marker association designs (CV_marker)
  array[n_id] matrix[15, R_mk] Z_mk_gk_now;
  array[n_id] matrix[15, R_mk] Z_mk_gk_fwd;
  array[n_id] matrix[15, Q_idm] Z_idm_gk_now;
  array[n_id] matrix[15, Q_idm] Z_idm_gk_fwd;
  matrix[n_id, R_mk] Z_mk_event_now;
  matrix[n_id, R_mk] Z_mk_event_fwd;
  matrix[n_id, Q_idm] Z_idm_event_now;
  matrix[n_id, Q_idm] Z_idm_event_fwd;

  // Association include flags
  int<lower=0, upper=1> assoc_cv_total;
  int<lower=0, upper=1> assoc_cv_mean;
  int<lower=0, upper=1> assoc_cv_marker;
  int<lower=0, upper=1> assoc_cs_total;
  int<lower=0, upper=1> assoc_cs_mean;
  int<lower=0, upper=1> assoc_cs_marker;
  int<lower=0, upper=1> assoc_vcov;

  // Prior scales
  vector[P] beta_scale;
  real<lower=0> alpha_scale;
  real<lower=0> lkj_eta;

  // Marker-side association shrinkage
  int<lower=1, upper=2> shrinkage; // 1 Laplace, 2 Gaussian

  // Time scaling metadata
  real<lower=0> tmax;
  int<lower=0> n_time_beta;
  array[n_time_beta] int<lower=1, upper=P> idx_time_beta;
  int<lower=0> n_time_uid;
  array[n_time_uid] int<lower=1, upper=R_id> idx_time_uid;
  int<lower=0> n_time_vmk;
  array[n_time_vmk] int<lower=1> idx_time_vmk;
  int<lower=0> n_time_widm;
  array[n_time_widm] int<lower=0, upper=Q_idm> idx_time_widm;

  // Transform configuration
  int<lower=0, upper=3> tf_mode_cv_tot;
  int<lower=0, upper=3> tf_mode_cs_tot;
  int<lower=0, upper=3> tf_mode_vcov;
  int<lower=0, upper=3> tf_mode_cv_mean;
  int<lower=0, upper=3> tf_mode_cv_marker;
  int<lower=0, upper=3> tf_mode_cs_mean;
  int<lower=0, upper=3> tf_mode_cs_marker;

  // Functional opcode specifications
  int<lower=0> n_functional_ops_cv;
  array[n_functional_ops_cv] int<lower=0, upper=23> functional_ops_cv;
  int<lower=0> n_const_cv;
  vector[n_const_cv] const_data_cv;

  int<lower=0> n_functional_ops_cs;
  array[n_functional_ops_cs] int<lower=0, upper=23> functional_ops_cs;
  int<lower=0> n_const_cs;
  vector[n_const_cs] const_data_cs;

  int<lower=0> n_functional_ops_vcov;
  array[n_functional_ops_vcov] int<lower=0, upper=23> functional_ops_vcov;
  int<lower=0> n_const_vcov;
  vector[n_const_vcov] const_data_vcov;

  int<lower=0> n_functional_ops_cv_mean;
  array[n_functional_ops_cv_mean] int<lower=0, upper=23> functional_ops_cv_mean;
  int<lower=0> n_const_cv_mean;
  vector[n_const_cv_mean] const_data_cv_mean;

  int<lower=0> n_functional_ops_cv_marker;
  array[n_functional_ops_cv_marker] int<lower=0, upper=23> functional_ops_cv_marker;
  int<lower=0> n_const_cv_marker;
  vector[n_const_cv_marker] const_data_cv_marker;

  int<lower=0> n_functional_ops_cs_mean;
  array[n_functional_ops_cs_mean] int<lower=0, upper=23> functional_ops_cs_mean;
  int<lower=0> n_const_cs_mean;
  vector[n_const_cs_mean] const_data_cs_mean;

  int<lower=0> n_functional_ops_cs_marker;
  array[n_functional_ops_cs_marker] int<lower=0, upper=23> functional_ops_cs_marker;
  int<lower=0> n_const_cs_marker;
  vector[n_const_cs_marker] const_data_cs_marker;

  // Spline specifications
  int<lower=0> n_knots_cv;
  vector[n_knots_cv] knots_cv;
  int<lower=0> n_coeff_cv;
  vector[n_coeff_cv] coeff_cv;
  int<lower=1, upper=5> spline_degree_cv;

  int<lower=0> n_knots_cs;
  vector[n_knots_cs] knots_cs;
  int<lower=0> n_coeff_cs;
  vector[n_coeff_cs] coeff_cs;
  int<lower=1, upper=5> spline_degree_cs;

  int<lower=0> n_knots_vcov;
  vector[n_knots_vcov] knots_vcov;
  int<lower=0> n_coeff_vcov;
  vector[n_coeff_vcov] coeff_vcov;
  int<lower=1, upper=5> spline_degree_vcov;

  int<lower=0> n_knots_cv_mean;
  vector[n_knots_cv_mean] knots_cv_mean;
  int<lower=0> n_coeff_cv_mean;
  vector[n_coeff_cv_mean] coeff_cv_mean;
  int<lower=1, upper=5> spline_degree_cv_mean;

  int<lower=0> n_knots_cv_marker;
  vector[n_knots_cv_marker] knots_cv_marker;
  int<lower=0> n_coeff_cv_marker;
  vector[n_coeff_cv_marker] coeff_cv_marker;
  int<lower=1, upper=5> spline_degree_cv_marker;

  int<lower=0> n_knots_cs_mean;
  vector[n_knots_cs_mean] knots_cs_mean;
  int<lower=0> n_coeff_cs_mean;
  vector[n_coeff_cs_mean] coeff_cs_mean;
  int<lower=1, upper=5> spline_degree_cs_mean;

  int<lower=0> n_knots_cs_marker;
  vector[n_knots_cs_marker] knots_cs_marker;
  int<lower=0> n_coeff_cs_marker;
  vector[n_coeff_cs_marker] coeff_cs_marker;
  int<lower=1, upper=5> spline_degree_cs_marker;
