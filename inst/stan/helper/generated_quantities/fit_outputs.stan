  /**
  * @brief Convenience outputs for debugging and interpretation.
   * beta is on ORIGINAL time scale by construction.
   * beta_scaled is the coefficient vector actually used with scaled-time design matrices.
   */
  vector[P] beta_used_in_likelihood = beta_scaled;
  
  /**
  * @brief Association coefficients actually used for totals.
   */
  real alpha_cv_total_used = a_cv_total;
  real alpha_cs_total_used = a_cs_total;

  /**
  * @brief Per-observation log-likelihoods for model assessment (loo/waic).
   *
   * We compute:
   * - log_lik_long[n]: longitudinal log-likelihood contribution for row n.
   * - log_lik_surv[i]: survival log-likelihood contribution for subject i.
   */
  vector[N] log_lik_long;
  vector[n_id] log_lik_surv;
  vector[n_id] cumhaz_event;
  vector[n_id] surv_prob_event;

  /* -------------------- Longitudinal log-likelihood per observation */
  for (n in 1 : N) {
    int i = id[n];
    int d = marker[n];
     real eta_long = dot_product(X_obs[n], beta_scaled)
                    + dot_product(Z_id_obs[n], u_id[i])
                    + ((R_mk > 0) ? dot_product(Z_mk_obs[n], v_marker[d])
                       : 0.0)
                    + dot_product(Z_idm_obs[n], w_idscaled[i, d]);
     // eta_long: linear predictor for longitudinal outcome
    int link_d = canonical_link_code_from_program(d, inv_link_n_ops, inv_link_ops, inv_link_n_const);
    real mu_long = inv_link_bytecode(eta_long, d, inv_link_n_ops, inv_link_ops, inv_link_n_const, inv_link_const);

    if (family_long[d] == 1) {
      real eta_sigma = 0; // linear predictor for log sigma
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          // b: random-effect coefficients for sigma term j
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for Gaussian
      log_lik_long[n] = normal_lpdf(y_real[n] | mu_long, sig);
    } else if (family_long[d] == 2) {
      real eta_sigma = 0; // linear predictor for log sigma
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          // b: random-effect coefficients for sigma term j
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for Student-t
      real eta_nu = 0; // linear predictor for df (nu)
      if (P_nu > 0) eta_nu += dot_product(X_nu[n], beta_nu);
      if (n_re_nu > 0) {
        for (j in 1 : n_re_nu) {
          vector[K_nu[j]] b = (to_vector(z_nu[j][J_nu[j, n], 1:K_nu[j]])
                               .* tau_nu[j][1:K_nu[j]]);
          // b: random-effect coefficients for nu term j
          eta_nu += dot_product(Z_nu[j][n, 1:K_nu[j]], b);
        }
      }
      real nu = (P_nu > 0 || n_re_nu > 0) ? (2 + exp(eta_nu)) : nu_family[marker_to_nu_family[d]];
      // nu: degrees of freedom for Student-t
      log_lik_long[n] = student_t_lpdf(y_real[n] | nu, mu_long, sig);
    } else if (family_long[d] == 3) {
      if (link_d == 3)
        log_lik_long[n] = bernoulli_logit_lpmf(y_int[n] | eta_long);
      else
        log_lik_long[n] = bernoulli_lpmf(y_int[n] | fmin(fmax(mu_long, 1e-12), 1 - 1e-12));
    } else if (family_long[d] == 4) {
      if (link_d == 3)
        log_lik_long[n] = binomial_logit_lpmf(y_int[n] | trials[n], eta_long);
      else
        log_lik_long[n] = binomial_lpmf(y_int[n] | trials[n], fmin(fmax(mu_long, 1e-12), 1 - 1e-12));
    } else if (family_long[d] == 5) {
      if (link_d == 5)
        log_lik_long[n] = poisson_log_lpmf(y_int[n] | eta_long);
      else
        log_lik_long[n] = poisson_lpmf(y_int[n] | fmax(mu_long, 1e-12));
    } else if (family_long[d] == 6) {
      real eta_phi = 0; // linear predictor for log phi (dispersion)
      if (P_phi > 0) eta_phi += dot_product(X_phi[n], beta_phi);
      if (n_re_phi > 0) {
        for (j in 1 : n_re_phi) {
          vector[K_phi[j]] b = (to_vector(z_phi[j][J_phi[j, n], 1:K_phi[j]])
                                .* tau_phi[j][1:K_phi[j]]);
          // b: random-effect coefficients for phi term j
          eta_phi += dot_product(Z_phi[j][n, 1:K_phi[j]], b);
        }
      }
      real phi = (P_phi > 0 || n_re_phi > 0) ? exp(eta_phi) : phi_family[marker_to_phi_family[d]];
      // phi: negbin2 dispersion parameter
      if (link_d == 5)
        log_lik_long[n] = neg_binomial_2_log_lpmf(y_int[n] | eta_long, phi);
      else
        log_lik_long[n] = neg_binomial_2_lpmf(y_int[n] | fmax(mu_long, 1e-12), phi);
    } else if (family_long[d] == 7) {
      real eta_sigma = 0; // linear predictor for log sigma
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          // b: random-effect coefficients for sigma term j
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for skew normal
      real eta_alpha = 0; // linear predictor for skew alpha
      if (P_alpha > 0) eta_alpha += dot_product(X_alpha[n], beta_alpha);
      if (n_re_alpha > 0) {
        for (j in 1 : n_re_alpha) {
          vector[K_alpha[j]] b = (to_vector(z_alpha[j][J_alpha[j, n], 1:K_alpha[j]])
                                  .* tau_alpha[j][1:K_alpha[j]]);
          // b: random-effect coefficients for alpha term j
          eta_alpha += dot_product(Z_alpha[j][n, 1:K_alpha[j]], b);
        }
      }
      real alpha = (P_alpha > 0 || n_re_alpha > 0) ? eta_alpha : alpha_family[marker_to_alpha_family[d]];
      // alpha: skewness parameter
      log_lik_long[n] = skew_normal_lpdf(y_real[n] | mu_long, sig, alpha);
    } else if (family_long[d] == 8) {
      real eta_sigma = 0; // linear predictor for log sigma
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          // b: random-effect coefficients for sigma term j
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for Laplace
      log_lik_long[n] = double_exponential_lpdf(y_real[n] | mu_long, sig);
    } else if (family_long[d] == 9) {
      real eta_sigma = 0; // linear predictor for log sigma
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          // b: random-effect coefficients for sigma term j
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for skew Laplace
      real eta_tau = 0; // linear predictor for tau_sde
      if (P_tau_sde > 0) eta_tau += dot_product(X_tau_sde[n], beta_tau_sde);
      if (n_re_tau_sde > 0) {
        for (j in 1 : n_re_tau_sde) {
          vector[K_tau_sde[j]] b = (to_vector(z_tau_sde[j][J_tau_sde[j, n], 1:K_tau_sde[j]])
                                    .* tau_tau_sde[j][1:K_tau_sde[j]]);
          // b: random-effect coefficients for tau_sde term j
          eta_tau += dot_product(Z_tau_sde[j][n, 1:K_tau_sde[j]], b);
        }
      }
      real tau_sde = (P_tau_sde > 0 || n_re_tau_sde > 0) ? inv_logit(eta_tau) : tau_sde_family[marker_to_tau_sde_family[d]];
      // tau_sde: skewness parameter in (0,1)
      log_lik_long[n] = skew_double_exponential_lpdf(y_real[n] | mu_long, sig, tau_sde);
    } else if (family_long[d] == 10) {
      real eta_phi_beta = 0; // linear predictor for log phi_beta
      if (P_phi_beta > 0) eta_phi_beta += dot_product(X_phi_beta[n], beta_phi_beta);
      if (n_re_phi_beta > 0) {
        for (j in 1 : n_re_phi_beta) {
          vector[K_phi_beta[j]] b = (to_vector(z_phi_beta[j][J_phi_beta[j, n], 1:K_phi_beta[j]])
                                     .* tau_phi_beta[j][1:K_phi_beta[j]]);
          // b: random-effect coefficients for phi_beta term j
          eta_phi_beta += dot_product(Z_phi_beta[j][n, 1:K_phi_beta[j]], b);
        }
      }
      real phi_beta = (P_phi_beta > 0 || n_re_phi_beta > 0) ? exp(eta_phi_beta) : phi_beta_family[marker_to_phi_beta_family[d]];
      real mu = mu_long; // beta mean on (0,1)
      real shape1 = phi_beta; // beta shape1
      real shape2 = (1 - mu) * phi_beta; // beta shape2
      log_lik_long[n] = beta_lpdf(y_real[n] | shape1, shape2);
    } else {
      log_lik_long[n] = ordered_logistic_lpmf(y_int[n] | eta_long, cutpoints_ord);
    }
  }

  /* -------------------- Survival log-likelihood per subject */
  for (i in 1 : n_id) {
    vector[n_gk] cvm_now;
    vector[n_gk] cvm_fwd;

    for (j in 1 : n_gk) {
      cvm_now[j] = dot_product(X_gk_now[i][j], beta_scaled)
                   + dot_product(Z_id_gk_now[i][j], u_id[i]);
      cvm_fwd[j] = dot_product(X_gk_fwd[i][j], beta_scaled)
                   + dot_product(Z_id_gk_fwd[i][j], u_id[i]);
      // cvm_now/cvm_fwd: mean CV at GK node and forward shift

    }

    vector[n_gk] csm_raw = eta_fd(cvm_now, cvm_fwd, eps_fd) / tmax; // mean slope
    vector[n_gk] csk_raw; // marker slope (weighted across markers)
    int M_corr_local = num_elements(a_corr);
    vector[M_corr_local] corr_terms_raw = eta_corr_varonly_weighted_const(L_i[i]); // raw off-diagonal corr terms
    int M_vcov_local = num_elements(a_vcov);
    vector[M_vcov_local] vcov_terms_raw = eta_vcov_weighted_const(L_i[i], (M_vcov_local == Q_idm)); // raw lower-triangular L entries

    int D_mkrs = size(v_marker);
    vector[n_gk] cv_tot_tf;
    vector[n_gk] cv_mean_tf = apply_transform_vector(
      cvm_now,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
      coeff_cv_mean_eff,
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
      for (d in 1 : D_mkrs) {
        real mk_part_now_d = 0;
        real mk_part_fwd_d = 0;
        if (R_mk > 0) {
          mk_part_now_d = dot_product(Z_mk_gk_now[i][j], v_marker[d]);
          mk_part_fwd_d = dot_product(Z_mk_gk_fwd[i][j], v_marker[d]);
        }
        real cvk_now_d = mk_part_now_d + dot_product(Z_idm_gk_now[i][j], w_idscaled[i, d]);
        real cvk_fwd_d = mk_part_fwd_d + dot_product(Z_idm_gk_fwd[i][j], w_idscaled[i, d]);
        real cv_tot_now_d = cvm_now[j] + cvk_now_d;
        acc_cv_tot_tf += marker_weights_eff[d]
          * apply_transform_scalar(cv_tot_now_d, tf_mode_cv_tot,
                                   functional_ops_cv, const_data_cv,
                                   knots_cv, coeff_cv_eff, spline_degree_cv);
        acc_cv_marker_tf += marker_weights_eff[d]
          * apply_transform_scalar(cvk_now_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker_eff, spline_degree_cv_marker);
        {
          real csk_raw_d = ((cvk_fwd_d - cvk_now_d) / eps_fd) / tmax;
          real cs_tot_raw_d = csm_raw[j] + csk_raw_d;

          acc_csk_raw += marker_weights_eff[d] * csk_raw_d;
          acc_cs_tot_tf += marker_weights_eff[d]
            * apply_transform_scalar(cs_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, const_data_cs,
                                     knots_cs, coeff_cs_eff, spline_degree_cs);
          acc_cs_marker_tf += marker_weights_eff[d]
            * apply_transform_scalar(csk_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker_eff, spline_degree_cs_marker);
        }
      }
      cv_tot_tf[j] = acc_cv_tot_tf / D_mkrs;
      cv_marker_tf[j] = acc_cv_marker_tf / D_mkrs;
      csk_raw[j] = acc_csk_raw / D_mkrs;
      cs_tot_tf[j] = acc_cs_tot_tf / D_mkrs;
      cs_marker_tf[j] = acc_cs_marker_tf / D_mkrs;
    }

    vector[n_gk] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
      coeff_cs_mean_eff,
      spline_degree_cs_mean
    );
    vector[M_corr_local] corr_terms_tf = apply_transform_vector_by_component(
      corr_terms_raw,
      tf_mode_corr,
      functional_ops_corr,
      const_data_corr,
      knots_corr,
      coeff_corr_eff,
      spline_degree_corr
    );
    vector[M_corr_local] corr_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_corr_local),
      tf_mode_corr,
      functional_ops_corr,
      const_data_corr,
      knots_corr,
      coeff_corr_eff,
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
      coeff_vcov_eff,
      spline_degree_vcov
    );
    vector[M_vcov_local] vcov_terms_ref = apply_transform_vector_by_component(
      rep_vector(0, M_vcov_local),
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov_eff,
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
    // eta_assoc_nodes: association predictor across GK nodes

    // Event contribution
    real event_term = 0; // log-hazard at event time (if event)
    if (d_event[i] == 1) {
      int k_ev = event_type[i];
      real eta_w_ev = 0; // hazard covariate contribution
      if (p_w > 0)
        eta_w_ev = dot_product(to_vector(W[i]'), gamma_w[k_ev]);
      real log_h0_S = dot_product(Bs_event_c[i], bs_gamma_c[k_ev]);
      // log_h0_S: baseline log-hazard at event time

      real cvm_S = dot_product(X_event_now[i], beta_scaled)
                   + dot_product(Z_id_event_now[i], u_id[i]);
      real cvm_S_fwd = dot_product(X_event_fwd[i], beta_scaled)
                       + dot_product(Z_id_event_fwd[i], u_id[i]);
      // cvm_S/cvm_S_fwd: mean CV at event time and forward shift
      real mk_part_S = 0;
      real mk_part_S_fwd = 0;
      // mk_part_S*: marker-only CV at event time
      if (R_mk > 0) {
        mk_part_S = dot_product(Z_mk_event_now[i], vbar);
        mk_part_S_fwd = dot_product(Z_mk_event_fwd[i], vbar);
      }
      real cvk_S = mk_part_S + dot_product(Z_idm_event_now[i], wbar_i[i]);
      real cvk_S_fwd = mk_part_S_fwd
               + dot_product(Z_idm_event_fwd[i], wbar_i[i]);
      // cvk_S/cvk_S_fwd: marker-id CV at event time and forward shift

      real csm_S_raw = ((cvm_S_fwd - cvm_S) / eps_fd) / tmax; // mean slope at event
      real csk_S_raw = 0; // marker slope at event (weighted across markers)

      real corr_S_assoc = corr_assoc_scalar; // same transformed-and-weighted constant at event time

      real cv_S_tot_tf = 0;
      real cv_S_mean_tf = apply_transform_scalar(
        cvm_S,
        tf_mode_cv_mean,
        functional_ops_cv_mean,
        const_data_cv_mean,
        knots_cv_mean,
        coeff_cv_mean_eff,
        spline_degree_cv_mean
      );
      real cv_S_marker_tf = 0;
      real cs_S_tot_tf = 0;
      real cs_S_marker_tf = 0;
      for (d in 1 : D_mkrs) {
        real mk_part_S_d = 0;
        real mk_part_S_fwd_d = 0;
        if (R_mk > 0) {
          mk_part_S_d = dot_product(Z_mk_event_now[i], v_marker[d]);
          mk_part_S_fwd_d = dot_product(Z_mk_event_fwd[i], v_marker[d]);
        }
        real cvk_S_d = mk_part_S_d + dot_product(Z_idm_event_now[i], w_idscaled[i, d]);
        real cvk_S_fwd_d = mk_part_S_fwd_d + dot_product(Z_idm_event_fwd[i], w_idscaled[i, d]);
        real cv_S_tot_d = cvm_S + cvk_S_d;

        cv_S_tot_tf += marker_weights_eff[d]
          * apply_transform_scalar(cv_S_tot_d, tf_mode_cv_tot,
                                   functional_ops_cv, const_data_cv,
                                   knots_cv, coeff_cv_eff, spline_degree_cv);
        cv_S_marker_tf += marker_weights_eff[d]
          * apply_transform_scalar(cvk_S_d, tf_mode_cv_marker,
                                   functional_ops_cv_marker, const_data_cv_marker,
                                   knots_cv_marker, coeff_cv_marker_eff, spline_degree_cv_marker);
        {
          real csk_S_raw_d = ((cvk_S_fwd_d - cvk_S_d) / eps_fd) / tmax;
          real cs_S_tot_raw_d = csm_S_raw + csk_S_raw_d;

          csk_S_raw += marker_weights_eff[d] * csk_S_raw_d;
          cs_S_tot_tf += marker_weights_eff[d]
            * apply_transform_scalar(cs_S_tot_raw_d, tf_mode_cs_tot,
                                     functional_ops_cs, const_data_cs,
                                     knots_cs, coeff_cs_eff, spline_degree_cs);
          cs_S_marker_tf += marker_weights_eff[d]
            * apply_transform_scalar(csk_S_raw_d, tf_mode_cs_marker,
                                     functional_ops_cs_marker, const_data_cs_marker,
                                     knots_cs_marker, coeff_cs_marker_eff, spline_degree_cs_marker);
        }
      }
      cv_S_tot_tf /= D_mkrs;
      cv_S_marker_tf /= D_mkrs;
      csk_S_raw /= D_mkrs;
      cs_S_tot_tf /= D_mkrs;
      cs_S_marker_tf /= D_mkrs;

      real cs_S_mean_tf = apply_transform_scalar(
        csm_S_raw,
        tf_mode_cs_mean,
        functional_ops_cs_mean,
        const_data_cs_mean,
        knots_cs_mean,
        coeff_cs_mean_eff,
        spline_degree_cs_mean
      );
      real eta_assoc_S = a_cv_total * cv_S_tot_tf
                         + a_cv_mean * cv_S_mean_tf
                         + a_cv_marker * cv_S_marker_tf
                         + a_cs_total * cs_S_tot_tf
                         + a_cs_mean * cs_S_mean_tf
                         + a_cs_marker * cs_S_marker_tf
                         + corr_S_assoc;
      // eta_assoc_S: association predictor at event time

      event_term = log_h0_S + eta_w_ev + eta_assoc_S;
    }

    vector[n_gk] log_h_total; // total log-hazard at quadrature nodes
    for (j in 1 : n_gk) {
      vector[K_event] log_h_cause; // cause-specific log hazards
      for (k_ev in 1 : K_event) {
        real eta_w = 0;
        if (p_w > 0)
          eta_w = dot_product(to_vector(W[i]'), gamma_w[k_ev]);
        log_h_cause[k_ev] = dot_product(Bs_gk_c[i][j], bs_gamma_c[k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    {
      real H_event = cumhaz(S_event[i], n_gk, log_h_total, rep_vector(0.0, n_gk));
      cumhaz_event[i] = H_event;
      surv_prob_event[i] = exp(-H_event);
      log_lik_surv[i] = event_term - H_event;
    }
    // log_lik_surv: event term minus cumulative hazard
  }
