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

  // -------------------- Longitudinal log-likelihood per observation
  for (n in 1 : N) {
    int i = id[n];
    int d = marker[n];
    real eta_long = dot_product(X_obs[n], beta_scaled)
                    + dot_product(Z_id_obs[n], u_id[i])
                    + ((R_mk > 0) ? dot_product(Z_mk_obs[n], v_marker[d])
                       : 0.0)
                    + dot_product(Z_idm_obs[n], w_idscaled[i, d]);

    if (family_long[d] == 1) {
      real eta_sigma = 0;
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : ((flag_resid_dim == 1) ? sigma_marker[d] : sigma_y);
      log_lik_long[n] = normal_lpdf(y_real[n] | eta_long, sig);
    } else if (family_long[d] == 2) {
      real eta_sigma = 0;
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : ((flag_resid_dim == 1) ? sigma_marker[d] : sigma_y);
      real eta_nu = 0;
      if (P_nu > 0) eta_nu += dot_product(X_nu[n], beta_nu);
      if (n_re_nu > 0) {
        for (j in 1 : n_re_nu) {
          vector[K_nu[j]] b = (to_vector(z_nu[j][J_nu[j, n], 1:K_nu[j]])
                               .* tau_nu[j][1:K_nu[j]]);
          eta_nu += dot_product(Z_nu[j][n, 1:K_nu[j]], b);
        }
      }
      real nu = (P_nu > 0 || n_re_nu > 0) ? (2 + exp(eta_nu)) : nu_marker[d];
      log_lik_long[n] = student_t_lpdf(y_real[n] | nu, eta_long, sig);
    } else if (family_long[d] == 3) {
      log_lik_long[n] = bernoulli_logit_lpmf(y_int[n] | eta_long);
    } else if (family_long[d] == 4) {
      log_lik_long[n] = binomial_logit_lpmf(y_int[n] | trials[n], eta_long);
    } else if (family_long[d] == 5) {
      log_lik_long[n] = poisson_log_lpmf(y_int[n] | eta_long);
    } else if (family_long[d] == 6) {
      real eta_phi = 0;
      if (P_phi > 0) eta_phi += dot_product(X_phi[n], beta_phi);
      if (n_re_phi > 0) {
        for (j in 1 : n_re_phi) {
          vector[K_phi[j]] b = (to_vector(z_phi[j][J_phi[j, n], 1:K_phi[j]])
                                .* tau_phi[j][1:K_phi[j]]);
          eta_phi += dot_product(Z_phi[j][n, 1:K_phi[j]], b);
        }
      }
      real phi = (P_phi > 0 || n_re_phi > 0) ? exp(eta_phi) : phi_nb_marker[d];
      log_lik_long[n] = neg_binomial_2_log_lpmf(y_int[n] | eta_long, phi);
    } else if (family_long[d] == 7) {
      real eta_sigma = 0;
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : ((flag_resid_dim == 1) ? sigma_marker[d] : sigma_y);
      real eta_alpha = 0;
      if (P_alpha > 0) eta_alpha += dot_product(X_alpha[n], beta_alpha);
      if (n_re_alpha > 0) {
        for (j in 1 : n_re_alpha) {
          vector[K_alpha[j]] b = (to_vector(z_alpha[j][J_alpha[j, n], 1:K_alpha[j]])
                                  .* tau_alpha[j][1:K_alpha[j]]);
          eta_alpha += dot_product(Z_alpha[j][n, 1:K_alpha[j]], b);
        }
      }
      real alpha = (P_alpha > 0 || n_re_alpha > 0) ? eta_alpha : alpha_skew_marker[d];
      log_lik_long[n] = skew_normal_lpdf(y_real[n] | eta_long, sig, alpha);
    } else if (family_long[d] == 8) {
      real eta_sigma = 0;
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : ((flag_resid_dim == 1) ? sigma_marker[d] : sigma_y);
      log_lik_long[n] = double_exponential_lpdf(y_real[n] | eta_long, sig);
    } else if (family_long[d] == 9) {
      real eta_sigma = 0;
      if (P_sigma > 0) eta_sigma += dot_product(X_sigma[n], beta_sigma);
      if (n_re_sigma > 0) {
        for (j in 1 : n_re_sigma) {
          vector[K_sigma[j]] b = (to_vector(z_sigma[j][J_sigma[j, n], 1:K_sigma[j]])
                                  .* tau_sigma[j][1:K_sigma[j]]);
          eta_sigma += dot_product(Z_sigma[j][n, 1:K_sigma[j]], b);
        }
      }
      real sig = (P_sigma > 0 || n_re_sigma > 0)
                 ? exp(eta_sigma)
                 : ((flag_resid_dim == 1) ? sigma_marker[d] : sigma_y);
      real eta_tau = 0;
      if (P_tau_sde > 0) eta_tau += dot_product(X_tau_sde[n], beta_tau_sde);
      if (n_re_tau_sde > 0) {
        for (j in 1 : n_re_tau_sde) {
          vector[K_tau_sde[j]] b = (to_vector(z_tau_sde[j][J_tau_sde[j, n], 1:K_tau_sde[j]])
                                    .* tau_tau_sde[j][1:K_tau_sde[j]]);
          eta_tau += dot_product(Z_tau_sde[j][n, 1:K_tau_sde[j]], b);
        }
      }
      real tau_sde = (P_tau_sde > 0 || n_re_tau_sde > 0) ? inv_logit(eta_tau) : tau_sde_marker[d];
      log_lik_long[n] = skew_double_exponential_lpdf(y_real[n] | eta_long, sig, tau_sde);
    } else if (family_long[d] == 10) {
      real eta_phi_beta = 0;
      if (P_phi_beta > 0) eta_phi_beta += dot_product(X_phi_beta[n], beta_phi_beta);
      if (n_re_phi_beta > 0) {
        for (j in 1 : n_re_phi_beta) {
          vector[K_phi_beta[j]] b = (to_vector(z_phi_beta[j][J_phi_beta[j, n], 1:K_phi_beta[j]])
                                     .* tau_phi_beta[j][1:K_phi_beta[j]]);
          eta_phi_beta += dot_product(Z_phi_beta[j][n, 1:K_phi_beta[j]], b);
        }
      }
      real phi_beta = (P_phi_beta > 0 || n_re_phi_beta > 0) ? exp(eta_phi_beta) : phi_beta_marker[d];
      real mu = inv_logit(eta_long);
      real shape1 = fmax(mu * phi_beta, 1e-6);
      real shape2 = fmax((1 - mu) * phi_beta, 1e-6);
      log_lik_long[n] = beta_lpdf(y_real[n] | shape1, shape2);
    } else {
      log_lik_long[n] = ordered_logistic_lpmf(y_int[n] | eta_long, cutpoints_ord);
    }
  }

  // -------------------- Survival log-likelihood per subject
  for (i in 1 : n_id) {
    vector[15] cvm_now;
    vector[15] cvm_fwd;
    vector[15] cvk_now;
    vector[15] cvk_fwd;

    for (j in 1 : 15) {
      cvm_now[j] = dot_product(X_gk_now[i][j], beta_scaled)
                   + dot_product(Z_id_gk_now[i][j], u_id[i]);
      cvm_fwd[j] = dot_product(X_gk_fwd[i][j], beta_scaled)
                   + dot_product(Z_id_gk_fwd[i][j], u_id[i]);

      real mk_part_now = 0;
      real mk_part_fwd = 0;
      if (R_mk > 0) {
        mk_part_now = dot_product(Z_mk_gk_now[i][j], vbar);
        mk_part_fwd = dot_product(Z_mk_gk_fwd[i][j], vbar);
      }

      cvk_now[j] = mk_part_now + dot_product(Z_idm_gk_now[i][j], wbar_i[i]);
      cvk_fwd[j] = mk_part_fwd + dot_product(Z_idm_gk_fwd[i][j], wbar_i[i]);
    }

    vector[15] csm_raw = eta_fd(cvm_now, cvm_fwd, eps_fd) / tmax;
    vector[15] csk_raw = eta_fd(cvk_now, cvk_fwd, eps_fd) / tmax;
    vector[15] vcov_raw = eta_vcov_varonly_weighted_const(15, L_i[i], a_vcov_var);

    vector[15] cv_tot = cvm_now + cvk_now;
    vector[15] cs_tot = csm_raw + csk_raw;

    vector[15] cv_tot_tf = apply_transform_vector(
      cv_tot,
      tf_mode_cv_tot,
      functional_ops_cv,
      const_data_cv,
      knots_cv,
      coeff_cv,
      spline_degree_cv
    );
    vector[15] cv_mean_tf = apply_transform_vector(
      cvm_now,
      tf_mode_cv_mean,
      functional_ops_cv_mean,
      const_data_cv_mean,
      knots_cv_mean,
      coeff_cv_mean,
      spline_degree_cv_mean
    );
    vector[15] cv_marker_tf = apply_transform_vector(
      cvk_now,
      tf_mode_cv_marker,
      functional_ops_cv_marker,
      const_data_cv_marker,
      knots_cv_marker,
      coeff_cv_marker,
      spline_degree_cv_marker
    );

    vector[15] cs_tot_tf = apply_transform_vector(
      cs_tot,
      tf_mode_cs_tot,
      functional_ops_cs,
      const_data_cs,
      knots_cs,
      coeff_cs,
      spline_degree_cs
    );
    vector[15] cs_mean_tf = apply_transform_vector(
      csm_raw,
      tf_mode_cs_mean,
      functional_ops_cs_mean,
      const_data_cs_mean,
      knots_cs_mean,
      coeff_cs_mean,
      spline_degree_cs_mean
    );
    vector[15] cs_marker_tf = apply_transform_vector(
      csk_raw,
      tf_mode_cs_marker,
      functional_ops_cs_marker,
      const_data_cs_marker,
      knots_cs_marker,
      coeff_cs_marker,
      spline_degree_cs_marker
    );
    vector[15] vcov_tf = apply_transform_vector(
      vcov_raw,
      tf_mode_vcov,
      functional_ops_vcov,
      const_data_vcov,
      knots_vcov,
      coeff_vcov,
      spline_degree_vcov
    );

    vector[15] eta_assoc_nodes = a_cv_total * cv_tot_tf
                                 + a_cv_mean * cv_mean_tf
                                 + a_cv_marker * cv_marker_tf
                                 + a_cs_total * cs_tot_tf
                                 + a_cs_mean * cs_mean_tf
                                 + a_cs_marker * cs_marker_tf
                                 + vcov_tf;

    // Event contribution
    real event_term = 0;
    if (d_event[i] == 1) {
      int k_ev = event_type[i];
      real eta_w_ev = 0;
      if (p_w > 0)
        eta_w_ev = dot_product(to_vector(W[i]'), gamma_w[k_ev]);
      real log_h0_S = log_h0_intercept[k_ev]
                      + dot_product(Bs_event_c[i], bs_gamma_c[k_ev]);

      real cvm_S = dot_product(X_event_now[i], beta_scaled)
                   + dot_product(Z_id_event_now[i], u_id[i]);
      real cvm_S_fwd = dot_product(X_event_fwd[i], beta_scaled)
                       + dot_product(Z_id_event_fwd[i], u_id[i]);
      real mk_part_S = 0;
      real mk_part_S_fwd = 0;
      if (R_mk > 0) {
        mk_part_S = dot_product(Z_mk_event_now[i], vbar);
        mk_part_S_fwd = dot_product(Z_mk_event_fwd[i], vbar);
      }
      real cvk_S = mk_part_S + dot_product(Z_idm_event_now[i], wbar_i[i]);
      real cvk_S_fwd = mk_part_S_fwd
                       + dot_product(Z_idm_event_fwd[i], wbar_i[i]);

      real csm_S_raw = ((cvm_S_fwd - cvm_S) / eps_fd) / tmax;
      real csk_S_raw = ((cvk_S_fwd - cvk_S) / eps_fd) / tmax;

      real vcov_S_raw = 0;
      for (q in 1 : Q_idm) {
        real var_q = 0;
        for (k in 1 : q)
          var_q += square(L_i[i][q, k]);
        vcov_S_raw += a_vcov_var[q] * log(var_q);
      }

      real cv_S_tot = cvm_S + cvk_S;
      real cs_S_tot = csm_S_raw + csk_S_raw;
      real cv_S_tot_tf = apply_transform_scalar(
        cv_S_tot,
        tf_mode_cv_tot,
        functional_ops_cv,
        const_data_cv,
        knots_cv,
        coeff_cv,
        spline_degree_cv
      );
      real cv_S_mean_tf = apply_transform_scalar(
        cvm_S,
        tf_mode_cv_mean,
        functional_ops_cv_mean,
        const_data_cv_mean,
        knots_cv_mean,
        coeff_cv_mean,
        spline_degree_cv_mean
      );
      real cv_S_marker_tf = apply_transform_scalar(
        cvk_S,
        tf_mode_cv_marker,
        functional_ops_cv_marker,
        const_data_cv_marker,
        knots_cv_marker,
        coeff_cv_marker,
        spline_degree_cv_marker
      );

      real cs_S_tot_tf = apply_transform_scalar(
        cs_S_tot,
        tf_mode_cs_tot,
        functional_ops_cs,
        const_data_cs,
        knots_cs,
        coeff_cs,
        spline_degree_cs
      );
      real cs_S_mean_tf = apply_transform_scalar(
        csm_S_raw,
        tf_mode_cs_mean,
        functional_ops_cs_mean,
        const_data_cs_mean,
        knots_cs_mean,
        coeff_cs_mean,
        spline_degree_cs_mean
      );
      real cs_S_marker_tf = apply_transform_scalar(
        csk_S_raw,
        tf_mode_cs_marker,
        functional_ops_cs_marker,
        const_data_cs_marker,
        knots_cs_marker,
        coeff_cs_marker,
        spline_degree_cs_marker
      );
      real vcov_S_tf = apply_transform_scalar(
        vcov_S_raw,
        tf_mode_vcov,
        functional_ops_vcov,
        const_data_vcov,
        knots_vcov,
        coeff_vcov,
        spline_degree_vcov
      );

      real eta_assoc_S = a_cv_total * cv_S_tot_tf
                         + a_cv_mean * cv_S_mean_tf
                         + a_cv_marker * cv_S_marker_tf
                         + a_cs_total * cs_S_tot_tf
                         + a_cs_mean * cs_S_mean_tf
                         + a_cs_marker * cs_S_marker_tf
                         + vcov_S_tf;

      event_term = log_h0_S + eta_w_ev + eta_assoc_S;
    }

    vector[15] log_h_total;
    for (j in 1 : 15) {
      vector[K_event] log_h_cause;
      for (k_ev in 1 : K_event) {
        real eta_w = 0;
        if (p_w > 0)
          eta_w = dot_product(to_vector(W[i]'), gamma_w[k_ev]);
        log_h_cause[k_ev] = log_h0_intercept[k_ev]
                            + dot_product(Bs_gk_c[i][j], bs_gamma_c[k_ev])
                            + eta_w
                            + eta_assoc_nodes[j];
      }
      log_h_total[j] = log_sum_exp(log_h_cause);
    }

    log_lik_surv[i] = event_term - cumhaz(S_event[i], log_h_total, rep_vector(0.0, 15));
  }
