/* Outputs per draw and per row/time */
matrix[n_draws, n_obs_long] y_fit_linpred;
matrix[n_draws, n_obs_long] y_fit_epred;
matrix[n_draws, n_obs_pred] y_pred_linpred;
matrix[n_draws, n_obs_pred] y_pred_epred;
matrix[n_draws, n_obs_pred] y_pred;
matrix[n_draws, n_times_surv] surv_prob;
matrix[n_draws, n_times_surv] cumhaz_cond;

for (k in 1 : n_draws) {

  // -------------------------------------------------------------
  // 1. Reconstruct Subject Random Effects (Scale and Correlate)
  // -------------------------------------------------------------

  // Shared ID Random Effects: u_id
  matrix[n_random_id, n_random_id] L_u;
  if (flag_indep_id_re == 1)
    L_u = diag_matrix(tau_id[k]);
  else
    L_u = diag_pre_multiply(tau_id[k], Lcorr_id[k]);
  vector[n_random_id] u_id = L_u * z_u[k]; // realized id random effects

  // Marker-Specific Random Effects (optional): v_marker
  matrix[n_random_marker, n_random_marker] L_v;
  if (n_random_marker > 0) {
    if (flag_indep_marker_re == 1)
      L_v = diag_matrix(tau_marker[k]);
    else
      L_v = diag_pre_multiply(tau_marker[k], Lcorr_marker[k]);
  }
  array[n_marker_types] vector[n_random_marker] v_marker; // realized marker REs
  if (n_random_marker > 0) {
    for (d in 1 : n_marker_types)
      v_marker[d] = L_v * z_v[k, d];
  }

  // Marker-by-ID Latent Effects: z_w
  // These latent seeds are iid standard normal. The full subject-specific
  // covariance now lives in Li, so there is no separate marker-by-id baseline
  // covariance matrix to reconstruct here.
  array[n_marker_types] vector[n_random_marker_id] z_w; // latent marker-id effects
  for (d in 1 : n_marker_types) {
    vector[n_random_marker_id] cross = rep_vector(0.0, n_random_marker_id);
    if (flag_allow_marker_crosscorr == 1 && n_random_marker > 0)
      cross = B_cross[k] * v_marker[d];
    z_w[d] = cross + z_w_lat[k, d];
  }

  // Use covariance regression to construct L_i directly from component-specific z_L[k][m]
  matrix[n_random_marker_id, n_random_marker_id] Li = rep_matrix(0.0,
                                                                 n_random_marker_id,
                                                                 n_random_marker_id);
  for (m in 1 : num_unique_cov_entries) {
    int r_ = idx_row_cov[m];
    int c_ = idx_col_cov[m];
    // unpack flattened beta
    vector[n_cov_vcov] bL_m = beta_vcov_reg_flat[k][((m - 1) * n_cov_vcov
                             + 1) : (m * n_cov_vcov)];
        real lp = alpha_vcov_reg[k][m] + dot_product(bL_m, vec_cov_vcov)
          + lambda_vcov_reg[k][m] * z_L[k][m];
    Li[r_, c_] = (r_ == c_) ? log1p_exp(lp) : lp;
  }

  matrix[n_random_marker_id, n_random_marker_id] Li_used = diag_pre_multiply(marker_id_row_scale, Li);

  // Scale latent effects: w_idscaled
  array[n_marker_types] vector[n_random_marker_id] w_idscaled; // scaled marker-id REs
  for (d in 1 : n_marker_types)
    w_idscaled[d] = Li_used * z_w[d];

  // -------------------------------------------------------------
  // 2. Prepare for Association (Averages)
  // -------------------------------------------------------------
  vector[n_marker_types] marker_weights_k = to_vector(marker_weights_draws[k, ]); // draw-specific weights
  // Aggregation is by marker count D to match fit-time construction.
  // Marker weights are applied in transformed marker aggregation, not here.
  vector[n_random_marker] vbar; // unweighted mean marker REs
  if (n_random_marker > 0) {
    for (r in 1 : n_random_marker) {
      real acc = 0;
      for (d in 1 : n_marker_types)
        acc += v_marker[d][r];
      vbar[r] = acc / n_marker_types;
    }
  }
  vector[n_random_marker_id] zbar; // unweighted mean marker-id latents
  for (q in 1 : n_random_marker_id) {
    real acc = 0;
    for (d in 1 : n_marker_types)
      acc += z_w[d][q];
    zbar[q] = acc / n_marker_types;
  }
  vector[n_random_marker_id] wbar_i = Li_used * zbar; // scaled mean marker-id effects
    array[n_marker_types] vector[n_random_marker_id] w_draw;
    for (d in 1 : n_marker_types)
      w_draw[d] = Li_used * z_w[d];

  // Global Association Scales
  real a_cv_total = flag_assoc_cv_total * coeff_assoc_cv_total[k];
  real a_cs_total = flag_assoc_cs_total * coeff_assoc_cs_total[k];
  real a_cv_mean = flag_assoc_cv_mean * coeff_assoc_cv_mean[k];
  real a_cs_mean = flag_assoc_cs_mean * coeff_assoc_cs_mean[k];
  real a_cv_marker = flag_assoc_cv_marker * coeff_assoc_cv_marker[k];
  real a_cs_marker = flag_assoc_cs_marker * coeff_assoc_cs_marker[k];
  vector[(n_random_marker_id * (n_random_marker_id - 1)) %/% 2] a_corr = flag_assoc_corr
                      * coeff_assoc_corr[k];
  vector[num_unique_cov_entries] a_vcov = flag_assoc_vcov * coeff_assoc_vcov[k];

  // -------------------------------------------------------------
  // 3. Fitted values for observed history (scale-aware)
  // -------------------------------------------------------------
  for (n in 1 : n_obs_long) {
    int d = idx_marker_obs[n];
     real eta_long = dot_product(mat_fixed_obs[n], beta_fixed[k])
                    + dot_product(mat_id_obs[n], u_id)
                    + ((n_random_marker > 0)
                       ? dot_product(mat_marker_obs[n], v_marker[d]) : 0.0)
                    + dot_product(mat_marker_id_obs[n], w_idscaled[d]);
     // eta_long: linear predictor for observed history

    y_fit_linpred[k, n] = eta_long;
    {
      int link_d = canonical_link_code_from_program(d, inv_link_n_ops, inv_link_ops, inv_link_n_const);
      real mu_long = inv_link_bytecode(eta_long, d, inv_link_n_ops, inv_link_ops, inv_link_n_const, inv_link_const);

      if (family_long[d] == 1) {
        y_fit_epred[k, n] = mu_long;
      } else if (family_long[d] == 2) {
        y_fit_epred[k, n] = mu_long;
      } else if (family_long[d] == 3) {
        y_fit_epred[k, n] = fmin(fmax(mu_long, 1e-12), 1 - 1e-12);
      } else if (family_long[d] == 4) {
        y_fit_epred[k, n] = trials_obs[n] * fmin(fmax(mu_long, 1e-12), 1 - 1e-12);
      } else if (family_long[d] == 5) {
        y_fit_epred[k, n] = fmax(mu_long, 1e-12);
      } else if (family_long[d] == 6) {
        y_fit_epred[k, n] = fmax(mu_long, 1e-12);
      } else if (family_long[d] == 7) {
        y_fit_epred[k, n] = mu_long;
      } else if (family_long[d] == 8) {
        y_fit_epred[k, n] = mu_long;
      } else if (family_long[d] == 9) {
        y_fit_epred[k, n] = mu_long;
      } else if (family_long[d] == 10) {
        y_fit_epred[k, n] = fmin(fmax(mu_long, 1e-12), 1 - 1e-12);
      } else {
      // Cumulative logit: report expected category
      vector[K_ord] p;
      p[1] = inv_logit(cutpoints_ord[k][1] - eta_long);
      for (c in 2 : (K_ord - 1))
        p[c] = inv_logit(cutpoints_ord[k][c] - eta_long)
               - inv_logit(cutpoints_ord[k][c - 1] - eta_long);
      p[K_ord] = 1 - inv_logit(cutpoints_ord[k][K_ord - 1] - eta_long);
      real acc = 0;
      for (c in 1 : K_ord) acc += c * p[c];
      y_fit_epred[k, n] = acc;
      }
    }
  }

  // -------------------------------------------------------------
  // 4. Longitudinal Prediction
  // -------------------------------------------------------------
  for (n in 1 : n_obs_pred) {
    int d = idx_marker_pred[n];
     real eta_long = dot_product(mat_fixed_pred[n], beta_fixed[k])
                    + dot_product(mat_id_pred[n], u_id)
                    + ((n_random_marker > 0)
                       ? dot_product(mat_marker_pred[n], v_marker[d]) : 0.0)
                    + dot_product(mat_marker_id_pred[n], w_idscaled[d]);
     // eta_long: linear predictor for prediction grid

    y_pred_linpred[k, n] = eta_long;
    {
      int link_d = canonical_link_code_from_program(d, inv_link_n_ops, inv_link_ops, inv_link_n_const);
      real mu_long = inv_link_bytecode(eta_long, d, inv_link_n_ops, inv_link_ops, inv_link_n_const, inv_link_const);

      if (family_long[d] == 1) {
      // Gaussian
      real sig = (P_sigma > 0) ? exp(fmin(dot_product(X_sigma_pred[n], beta_sigma[k]), 20))
             : sigma_family[k][marker_to_sigma_family[d]];
      // sig: residual scale for Gaussian prediction
      y_pred_epred[k, n] = mu_long;
      y_pred[k, n] = normal_rng(mu_long, sig);
      } else if (family_long[d] == 2) {
      // Student-t
      real sig = (P_sigma > 0) ? exp(fmin(dot_product(X_sigma_pred[n], beta_sigma[k]), 20))
             : sigma_family[k][marker_to_sigma_family[d]];
      real nu = (P_nu > 0) ? (2 + exp(fmin(dot_product(X_nu_pred[n], beta_nu[k]), 20)))
            : nu_family[k][marker_to_nu_family[d]];
      // nu: degrees of freedom for Student-t prediction
      y_pred_epred[k, n] = mu_long;
      y_pred[k, n] = student_t_rng(nu, mu_long, sig);
      } else if (family_long[d] == 3) {
      // Bernoulli (logit link)
      real p = fmin(fmax(mu_long, 1e-12), 1 - 1e-12); // Bernoulli probability
      y_pred_epred[k, n] = p;
      y_pred[k, n] = bernoulli_rng(p);
      } else if (family_long[d] == 4) {
      // Binomial (logit link)
      real p = fmin(fmax(mu_long, 1e-12), 1 - 1e-12); // Binomial success probability
      y_pred_epred[k, n] = trials_pred[n] * p;
      y_pred[k, n] = binomial_rng(trials_pred[n], p);
      } else if (family_long[d] == 5) {
      // Poisson (log link)
      real mu = fmax(mu_long, 1e-12); // Poisson mean
      y_pred_epred[k, n] = mu;
      y_pred[k, n] = poisson_rng(mu);
      } else if (family_long[d] == 6) {
      // Negative binomial 2 (log link)
      real mu = fmax(mu_long, 1e-12); // NegBin mean
      y_pred_epred[k, n] = mu;
      {
        real phi = (P_phi > 0) ? exp(dot_product(X_phi_pred[n], beta_phi[k]))
             : phi_family[k][marker_to_phi_family[d]];
        // phi: dispersion for NegBin2
        y_pred[k, n] = neg_binomial_2_rng(mu, phi);
      }
      } else if (family_long[d] == 7) {
      // Skew-normal
      real sig = (P_sigma > 0) ? exp(fmin(dot_product(X_sigma_pred[n], beta_sigma[k]), 20))
             : sigma_family[k][marker_to_sigma_family[d]];
      real alpha = (P_alpha > 0) ? dot_product(X_alpha_pred[n], beta_alpha[k])
          : alpha_family[k][marker_to_alpha_family[d]];
      // sig: residual scale, alpha: skew parameter
      y_pred_epred[k, n] = mu_long;
      y_pred[k, n] = skew_normal_rng(mu_long, sig, alpha);
      } else if (family_long[d] == 8) {
      // Double exponential (Laplace)
      real sig = (P_sigma > 0) ? exp(fmin(dot_product(X_sigma_pred[n], beta_sigma[k]), 20))
             : sigma_family[k][marker_to_sigma_family[d]];
      // sig: Laplace scale
      y_pred_epred[k, n] = mu_long;
      y_pred[k, n] = double_exponential_rng(mu_long, sig);
      } else if (family_long[d] == 9) {
      // Skew double exponential (asymmetric Laplace)
        real sig = (P_sigma > 0) ? exp(fmin(dot_product(X_sigma_pred[n], beta_sigma[k]), 20))
             : sigma_family[k][marker_to_sigma_family[d]];
      real tau_sde = (P_tau_sde > 0) ? inv_logit(dot_product(X_tau_sde_pred[n], beta_tau_sde[k]))
               : tau_sde_family[k][marker_to_tau_sde_family[d]];
      // sig: scale, tau_sde: skewness in (0,1)
      y_pred_epred[k, n] = mu_long;
      y_pred[k, n] = skew_double_exponential_rng(mu_long, sig, tau_sde);
      } else if (family_long[d] == 10) {
      // Beta
      real phi_beta = (P_phi_beta > 0) ? exp(dot_product(X_phi_beta_pred[n], beta_phi_beta[k]))
                      : phi_beta_family[k][marker_to_phi_beta_family[d]];
      real mu = fmin(fmax(mu_long, 1e-12), 1 - 1e-12);
      real shape1 = fmax(mu * phi_beta, 1e-6);
      real shape2 = fmax((1 - mu) * phi_beta, 1e-6);
      // mu: mean, phi_beta: precision, shape1/shape2: beta shapes
      y_pred_epred[k, n] = mu;
      y_pred[k, n] = beta_rng(shape1, shape2);
      } else {
      // Cumulative logit (ordered logistic)
      vector[K_ord] p;
      p[1] = inv_logit(cutpoints_ord[k][1] - eta_long);
      for (c in 2 : (K_ord - 1))
        p[c] = inv_logit(cutpoints_ord[k][c] - eta_long)
               - inv_logit(cutpoints_ord[k][c - 1] - eta_long);
      p[K_ord] = 1 - inv_logit(cutpoints_ord[k][K_ord - 1] - eta_long);
      real acc = 0;
      for (c in 1 : K_ord) acc += c * p[c];
      // p: category probabilities, acc: expected category
      y_pred_epred[k, n] = acc;
      y_pred[k, n] = ordered_logistic_rng(eta_long, cutpoints_ord[k]);
      }
    }
  }

  // -------------------------------------------------------------
  // 5. Survival Prediction S(t | T_cond)
  // -------------------------------------------------------------
  // Calculate Cumulative Hazard at Conditioning Time: H(T_cond)
  real H_cond; // cumulative hazard at conditioning time
  {
    vector[n_gk] cvm, cvk, cvm_f, cvk_f; // mean/marker CV at nodes and forward shift
    for (j in 1 : n_gk) {
      // Mean trajectory association (Current Value)
      cvm[j] = dot_product(mat_fixed_gk_cond[j], beta_fixed[k])
               + dot_product(mat_id_gk_cond[j], u_id);
      cvm_f[j] = dot_product(mat_fixed_gk_cond_fwd[j], beta_fixed[k])
                 + dot_product(mat_id_gk_cond_fwd[j], u_id);

      // Marker-Specific trajectory association
      real mk_part = 0, mk_part_f = 0;
      if (n_random_marker > 0) {
        mk_part = dot_product(mat_marker_gk_cond[j], vbar);
        mk_part_f = dot_product(mat_marker_gk_cond_fwd[j], vbar);
      }
      cvk[j] = mk_part + dot_product(mat_marker_id_gk_cond[j], wbar_i);
      cvk_f[j] = mk_part_f
                 + dot_product(mat_marker_id_gk_cond_fwd[j], wbar_i);
    }

    vector[n_gk] csm_raw = eta_fd(cvm, cvm_f, eps_finite_diff); // mean slope
    vector[n_gk] csk_raw; // marker slope (weighted across markers)
    int M_corr_local = num_elements(a_corr);
    vector[M_corr_local] corr_terms_raw = eta_corr_varonly_weighted_const(Li);
    int M_vcov_local = num_elements(a_vcov);
    vector[M_vcov_local] vcov_terms_raw = eta_vcov_weighted_const(Li, flag_indep_idmarker_cov);

    vector[n_gk] cv_tot = cvm + cvk; // total current value

    vector[n_gk] cv_tot_tf;
    vector[n_gk] cv_mean_tf = apply_transform_vector(
      cvm,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
        coeff_cv_mean[k],
      spline_degree_cv_mean
    );
    vector[n_gk] cv_marker_tf;
    vector[n_gk] cs_tot_tf;
    vector[n_gk] cs_marker_tf;
    for (j in 1 : n_gk) {
      real acc_cv_tot_tf = 0;
      real acc_cv_marker_tf = 0;
      real acc_csk_raw = 0;
      real acc_cs_tot_tf = 0;
      real acc_cs_marker_tf = 0;
      for (d in 1 : n_marker_types) {
        real mk_part_d = 0;
        real mk_part_f_d = 0;
        if (n_random_marker > 0) {
          mk_part_d = dot_product(mat_marker_gk_cond[j], v_marker[d]);
          mk_part_f_d = dot_product(mat_marker_gk_cond_fwd[j], v_marker[d]);
        }
        real cvk_d = mk_part_d + dot_product(mat_marker_id_gk_cond[j], w_draw[d]);
        real cvk_f_d = mk_part_f_d + dot_product(mat_marker_id_gk_cond_fwd[j], w_draw[d]);
        real cv_tot_d = cvm[j] + cvk_d;

        acc_cv_tot_tf += marker_weights_k[d]
          * apply_transform_scalar(cv_tot_d, tf_mode_cv_tot,
                                   functional_ops_cv, const_data_cv,
                                   knots_cv, coeff_cv[k], spline_degree_cv);
        acc_cv_marker_tf += marker_weights_k[d]
          * apply_transform_scalar(cvk_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker[k], spline_degree_cv_marker);
        {
          real csk_raw_d = (cvk_f_d - cvk_d) / eps_finite_diff;
          real cs_tot_raw_d = csm_raw[j] + csk_raw_d;

          acc_csk_raw += marker_weights_k[d] * csk_raw_d;
          acc_cs_tot_tf += marker_weights_k[d]
            * apply_transform_scalar(cs_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, const_data_cs,
                                     knots_cs, coeff_cs[k], spline_degree_cs);
          acc_cs_marker_tf += marker_weights_k[d]
            * apply_transform_scalar(csk_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker[k], spline_degree_cs_marker);
        }
      }
      cv_tot_tf[j] = acc_cv_tot_tf / n_marker_types;
      cv_marker_tf[j] = acc_cv_marker_tf / n_marker_types;
      csk_raw[j] = acc_csk_raw / n_marker_types;
      cs_tot_tf[j] = acc_cs_tot_tf / n_marker_types;
      cs_marker_tf[j] = acc_cs_marker_tf / n_marker_types;
    }

    vector[n_gk] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
        coeff_cs_mean[k],
      spline_degree_cs_mean
    );
    vector[M_corr_local] corr_terms_tf = apply_transform_vector_by_component(
      corr_terms_raw,
      tf_mode_corr,
      functional_ops_corr,
      const_data_corr,
      knots_corr,
        coeff_corr[k],
      spline_degree_corr
    );
      vector[M_corr_local] corr_terms_ref = apply_transform_vector_by_component(
        rep_vector(0, M_corr_local),
        tf_mode_corr,
        functional_ops_corr,
        const_data_corr,
        knots_corr,
        coeff_corr[k],
        spline_degree_corr
      );
      corr_terms_tf = corr_terms_tf - corr_terms_ref;
    real corr_assoc_scalar = dot_product(a_corr, corr_terms_tf);
    vector[n_gk] corr_assoc = rep_vector(corr_assoc_scalar, n_gk);
    vector[M_vcov_local] vcov_terms_tf = apply_transform_vector_by_component(
      vcov_terms_raw,
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
        coeff_vcov[k],
      spline_degree_vcov
    );
    vector[M_vcov_local] vcov_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_vcov_local),
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov[k],
      spline_degree_vcov
    );
    vcov_terms_tf = vcov_terms_tf - vcov_terms_ref;
    real vcov_assoc_scalar = dot_product(a_vcov, vcov_terms_tf);
    vector[n_gk] vcov_assoc = rep_vector(vcov_assoc_scalar, n_gk);

    vector[n_gk] eta_assoc_nodes = a_cv_total * cv_tot_tf
                                 + a_cv_mean * cv_mean_tf
                                 + a_cv_marker * cv_marker_tf
                                 + a_cs_total * cs_tot_tf
                                 + a_cs_mean * cs_mean_tf
                                 + a_cs_marker * cs_marker_tf
                                 + corr_assoc
                                 + vcov_assoc;
    // eta_assoc_nodes: association predictor at GK nodes

    vector[n_gk] log_h_total; // total log-hazard at quadrature nodes
    for (j in 1 : n_gk) {
      vector[K_event] log_h_cause; // cause-specific log hazards
      for (k_ev in 1 : K_event) {
        real eta_w = 0;
        if (n_cov_hazard > 0)
          eta_w = dot_product(vec_cov_hazard, gamma_hazard[k, k_ev]);
        log_h_cause[k_ev] = dot_product(mat_basis_gk_cond[j], bs_gamma_c[k, k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    H_cond = cumhaz(time_condition, n_gk, log_h_total, rep_vector(0.0, n_gk));
    // H_cond: cumulative hazard at time_condition
  }

  // Calculate H(t_surv) for each prediction time and compute S(t_surv | T_cond)
  for (s in 1 : n_times_surv) {
    real H_t; // cumulative hazard at prediction time s

    vector[n_gk] cvm, cvk, cvm_f, cvk_f; // CV features at quadrature nodes for time s
    for (j in 1 : n_gk) {
      cvm[j] = dot_product(mat_fixed_gk_surv[s][j], beta_fixed[k])
               + dot_product(mat_id_gk_surv[s][j], u_id);
      cvm_f[j] = dot_product(mat_fixed_gk_surv_fwd[s][j], beta_fixed[k])
                 + dot_product(mat_id_gk_surv_fwd[s][j], u_id);

      real mk_part = 0, mk_part_f = 0; // marker-only contribution
      if (n_random_marker > 0) {
        mk_part = dot_product(mat_marker_gk_surv[s][j], vbar);
        mk_part_f = dot_product(mat_marker_gk_surv_fwd[s][j], vbar);
      }
      cvk[j] = mk_part + dot_product(mat_marker_id_gk_surv[s][j], wbar_i);
      cvk_f[j] = mk_part_f
                 + dot_product(mat_marker_id_gk_surv_fwd[s][j], wbar_i);
    }

    vector[n_gk] csm_raw = eta_fd(cvm, cvm_f, eps_finite_diff); // mean slope
    vector[n_gk] csk_raw; // marker slope (weighted across markers)
    int M_corr_local2 = num_elements(a_corr);
    vector[M_corr_local2] corr_terms_raw = eta_corr_varonly_weighted_const(Li);
    int M_vcov_local2 = num_elements(a_vcov);
    vector[M_vcov_local2] vcov_terms_raw = eta_vcov_weighted_const(Li, flag_indep_idmarker_cov);

    vector[n_gk] cv_tot = cvm + cvk; // total current value

    vector[n_gk] cv_tot_tf;
    vector[n_gk] cv_mean_tf = apply_transform_vector(
      cvm,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
        coeff_cv_mean[k],
      spline_degree_cv_mean
    );
    vector[n_gk] cv_marker_tf;
    vector[n_gk] cs_tot_tf;
    vector[n_gk] cs_marker_tf;
    for (j in 1 : n_gk) {
      real acc_cv_tot_tf = 0;
      real acc_cv_marker_tf = 0;
      real acc_csk_raw = 0;
      real acc_cs_tot_tf = 0;
      real acc_cs_marker_tf = 0;
      for (d in 1 : n_marker_types) {
        real mk_part_d = 0;
        real mk_part_f_d = 0;
        if (n_random_marker > 0) {
          mk_part_d = dot_product(mat_marker_gk_surv[s][j], v_marker[d]);
          mk_part_f_d = dot_product(mat_marker_gk_surv_fwd[s][j], v_marker[d]);
        }
        real cvk_d = mk_part_d + dot_product(mat_marker_id_gk_surv[s][j], w_draw[d]);
        real cvk_f_d = mk_part_f_d + dot_product(mat_marker_id_gk_surv_fwd[s][j], w_draw[d]);
        real cv_tot_d = cvm[j] + cvk_d;

        acc_cv_tot_tf += marker_weights_k[d]
          * apply_transform_scalar(cv_tot_d, tf_mode_cv_tot,
                                   functional_ops_cv, const_data_cv,
                                   knots_cv, coeff_cv[k], spline_degree_cv);
        acc_cv_marker_tf += marker_weights_k[d]
          * apply_transform_scalar(cvk_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker[k], spline_degree_cv_marker);
        {
          real csk_raw_d = (cvk_f_d - cvk_d) / eps_finite_diff;
          real cs_tot_raw_d = csm_raw[j] + csk_raw_d;

          acc_csk_raw += marker_weights_k[d] * csk_raw_d;
          acc_cs_tot_tf += marker_weights_k[d]
            * apply_transform_scalar(cs_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, const_data_cs,
                                     knots_cs, coeff_cs[k], spline_degree_cs);
          acc_cs_marker_tf += marker_weights_k[d]
            * apply_transform_scalar(csk_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker[k], spline_degree_cs_marker);
        }
      }
      cv_tot_tf[j] = acc_cv_tot_tf / n_marker_types;
      cv_marker_tf[j] = acc_cv_marker_tf / n_marker_types;
      csk_raw[j] = acc_csk_raw / n_marker_types;
      cs_tot_tf[j] = acc_cs_tot_tf / n_marker_types;
      cs_marker_tf[j] = acc_cs_marker_tf / n_marker_types;
    }

    vector[n_gk] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
        coeff_cs_mean[k],
      spline_degree_cs_mean
    );
    vector[M_corr_local2] corr_terms_tf = apply_transform_vector_by_component(
      corr_terms_raw,
      tf_mode_corr,
      functional_ops_corr,
      const_data_corr,
      knots_corr,
        coeff_corr[k],
      spline_degree_corr
    );
      vector[M_corr_local2] corr_terms_ref = apply_transform_vector_by_component(
        rep_vector(0, M_corr_local2),
        tf_mode_corr,
        functional_ops_corr,
        const_data_corr,
        knots_corr,
        coeff_corr[k],
        spline_degree_corr
      );
      corr_terms_tf = corr_terms_tf - corr_terms_ref;
    real corr_assoc_scalar = dot_product(a_corr, corr_terms_tf);
    vector[n_gk] corr_assoc = rep_vector(corr_assoc_scalar, n_gk);
    vector[M_vcov_local2] vcov_terms_tf = apply_transform_vector_by_component(
      vcov_terms_raw,
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
        coeff_vcov[k],
      spline_degree_vcov
    );
    vector[M_vcov_local2] vcov_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_vcov_local2),
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov[k],
      spline_degree_vcov
    );
    vcov_terms_tf = vcov_terms_tf - vcov_terms_ref;
    real vcov_assoc_scalar = dot_product(a_vcov, vcov_terms_tf);
    vector[n_gk] vcov_assoc = rep_vector(vcov_assoc_scalar, n_gk);

    vector[n_gk] eta_assoc_nodes = a_cv_total * cv_tot_tf
                                 + a_cv_mean * cv_mean_tf
                                 + a_cv_marker * cv_marker_tf
                                 + a_cs_total * cs_tot_tf
                                 + a_cs_mean * cs_mean_tf
                                 + a_cs_marker * cs_marker_tf
                                 + corr_assoc
                                 + vcov_assoc;
    // eta_assoc_nodes: association predictor at GK nodes

    vector[n_gk] log_h_total; // total log-hazard at quadrature nodes
    for (j in 1 : n_gk) {
      vector[K_event] log_h_cause; // cause-specific log hazards
      for (k_ev in 1 : K_event) {
        real eta_w = 0;
        if (n_cov_hazard > 0)
          eta_w = dot_product(vec_cov_hazard, gamma_hazard[k, k_ev]);
        log_h_cause[k_ev] = dot_product(mat_basis_gk_surv[s][j], bs_gamma_c[k, k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    H_t = cumhaz(vec_time_surv[s], n_gk, log_h_total, rep_vector(0.0, n_gk));
    surv_prob[k, s] = exp(H_cond - H_t); // S(t | T_cond)
    cumhaz_cond[k, s] = H_t - H_cond;    // H(t) - H(T_cond)
  }
}
