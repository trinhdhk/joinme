/**
 * @file joinme_fit_threading.stan
 * @brief Threaded joint mixed-effects model for multivariate long-format longitudinal data + survival.
 *
 * ## Threading strategy
 * This threaded variant uses `reduce_sum` to parallelize the log-likelihood
 * across subjects. Each task computes the longitudinal + survival contribution
 * for a slice of subjects, while the prior terms remain on the main thread.
 * This preserves deterministic behavior for priors and avoids redundant
 * prior evaluation across worker threads.
 *
 * The model structure mirrors the non-threaded version:
 * - Same longitudinal and survival submodels.
 * - Same association terms and transform pipeline (functional, spline, pwlin).
 * - Same internal time scaling with rescaling to original time inside Stan.
 *
 * Additional data fields:
 * - `id_start`, `id_end`: per-subject observation ranges (longitudinal rows).
 * - `grainsize`: chunk size for `reduce_sum` (smaller = more tasks, larger = less overhead).
 *
 * ## Reproducibility notes
 * RNG is not used inside `reduce_sum` tasks; all randomness is handled by the
 * sampler. Results are deterministic given the same seed and threading config.
 */

functions {
  #include helper/functions/apply_tf.stanfunctions
  #include helper/functions/apply_tf_vec.stanfunctions
  
  #include helper/functions/eta_fd.stanfunctions
  
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  
  #include helper/functions/cumhaz.stanfunctions
  
  #include helper/functions/functional_transform.stanfunctions
  
  #include helper/functions/basis_functions.stanfunctions
  
  // Inline composite_transform to avoid include-related stanc3 ICE on Windows
  real apply_transform_scalar(
      real x,
      int mode,
      array[] int opcodes,
      vector const_data,
      vector spline_knots,
      vector spline_coeff,
      int spline_degree
  ) {
    if (mode == 0) {
      return x;
    } else if (mode == 1) {
      return eval_functional_scalar(x, opcodes, const_data);
    } else if (mode == 2) {
      return ispline_eval(x, spline_knots, spline_degree, spline_coeff);
    } else if (mode == 3) {
      return pwlin_eval(spline_knots, spline_coeff, x);
    } else {
      reject("Unknown transformation mode: ", mode);
    }
  }

  vector apply_transform_vector(
      vector x,
      int mode,
      array[] int opcodes,
      vector const_data,
      vector spline_knots,
      vector spline_coeff,
      int spline_degree
  ) {
    if (mode == 0) {
      return x;
    } else if (mode == 1) {
      return eval_functional_vector(x, opcodes, const_data);
    } else if (mode == 2) {
      return ispline_eval_vec(x, spline_knots, spline_degree, spline_coeff);
    } else if (mode == 3) {
      return pwlin_eval_vec(spline_knots, spline_coeff, x);
    } else {
      reject("Unknown transformation mode: ", mode);
    }
  }

  /**
   * @brief Per-subject log-likelihood slice for reduce_sum.
   *
   * The function receives a slice of subject indices (via id_seq[start:end])
   * and accumulates the longitudinal + survival log-likelihood for those
   * subjects only. It does NOT include priors (handled in the main thread)
   * and uses precomputed id_start/id_end to iterate each subject's rows.
   *
   * All inputs are read-only; the returned `lp` is the sum of contributions
   * for the slice and is added to the global target via reduce_sum.
   */
  real partial_joinme(
      array[] int id_seq,
      int start,
      int end,
      array[] int id,
      array[] int marker,
      vector y_real,
      array[] int y_int,
      array[] int trials,
      array[] int family_long,
      int flag_resid_dim,
      matrix X_obs,
      matrix Z_id_obs,
      matrix Z_mk_obs,
      matrix Z_idm_obs,
      int R_mk,
      int Q_idm,
      vector beta_scaled,
      int P_sigma,
      matrix X_sigma,
      int P_nu,
      matrix X_nu,
      int P_phi,
      matrix X_phi,
      int P_alpha,
      matrix X_alpha,
      int P_phi_beta,
      matrix X_phi_beta,
      int P_tau_sde,
      matrix X_tau_sde,
      vector beta_sigma,
      vector beta_nu,
      vector beta_phi,
      vector beta_alpha,
      vector beta_phi_beta,
      vector beta_tau_sde,
      int n_re_sigma,
      array[] int K_sigma,
      int K_sigma_max,
      array[] matrix Z_sigma,
      array[,] int J_sigma,
      array[] vector tau_sigma,
      array[] matrix z_sigma,
      int n_re_nu,
      array[] int K_nu,
      int K_nu_max,
      array[] matrix Z_nu,
      array[,] int J_nu,
      array[] vector tau_nu,
      array[] matrix z_nu,
      int n_re_phi,
      array[] int K_phi,
      int K_phi_max,
      array[] matrix Z_phi,
      array[,] int J_phi,
      array[] vector tau_phi,
      array[] matrix z_phi,
      int n_re_alpha,
      array[] int K_alpha,
      int K_alpha_max,
      array[] matrix Z_alpha,
      array[,] int J_alpha,
      array[] vector tau_alpha,
      array[] matrix z_alpha,
      int n_re_phi_beta,
      array[] int K_phi_beta,
      int K_phi_beta_max,
      array[] matrix Z_phi_beta,
      array[,] int J_phi_beta,
      array[] vector tau_phi_beta,
      array[] matrix z_phi_beta,
      int n_re_tau_sde,
      array[] int K_tau_sde,
      int K_tau_sde_max,
      array[] matrix Z_tau_sde,
      array[,] int J_tau_sde,
      array[] vector tau_tau_sde,
      array[] matrix z_tau_sde,
      array[] vector u_id,
      array[] vector v_marker,
      array[,] vector w_idscaled,
      array[] vector wbar_i,
      real sigma_y,
      vector sigma_marker,
      vector nu_marker,
      vector phi_nb_marker,
      vector alpha_skew_marker,
      vector phi_beta_marker,
      vector tau_sde_marker,
      int K_ord,
      vector cutpoints_ord,
      int p_w,
      matrix W,
      array[] vector gamma_w,
      array[] real log_h0_intercept,
      array[] vector bs_gamma_c,
      int Kbs,
      matrix Bs_event_c,
      array[] matrix Bs_gk_c,
      vector S_event,
      array[] int d_event,
      int K_event,
      array[] int event_type,
      real eps_fd,
      array[] matrix X_gk_now,
      array[] matrix X_gk_fwd,
      array[] matrix Z_id_gk_now,
      array[] matrix Z_id_gk_fwd,
      matrix X_event_now,
      matrix X_event_fwd,
      matrix Z_id_event_now,
      matrix Z_id_event_fwd,
      array[] matrix Z_mk_gk_now,
      array[] matrix Z_mk_gk_fwd,
      array[] matrix Z_idm_gk_now,
      array[] matrix Z_idm_gk_fwd,
      matrix Z_mk_event_now,
      matrix Z_mk_event_fwd,
      matrix Z_idm_event_now,
      matrix Z_idm_event_fwd,
      vector vbar,
      array[] matrix L_i,
      real a_cv_total,
      real a_cs_total,
      real a_cv_mean,
      real a_cv_marker,
      real a_cs_mean,
      real a_cs_marker,
      vector a_vcov_var,
      int tf_mode_cv_tot,
      int tf_mode_cs_tot,
      int tf_mode_vcov,
      int tf_mode_cv_mean,
      int tf_mode_cv_marker,
      int tf_mode_cs_mean,
      int tf_mode_cs_marker,
      array[] int functional_ops_cv,
      vector const_data_cv,
      vector knots_cv,
      vector coeff_cv,
      int spline_degree_cv,
      array[] int functional_ops_cs,
      vector const_data_cs,
      vector knots_cs,
      vector coeff_cs,
      int spline_degree_cs,
      array[] int functional_ops_vcov,
      vector const_data_vcov,
      vector knots_vcov,
      vector coeff_vcov,
      int spline_degree_vcov,
      array[] int functional_ops_cv_mean,
      vector const_data_cv_mean,
      vector knots_cv_mean,
      vector coeff_cv_mean,
      int spline_degree_cv_mean,
      array[] int functional_ops_cv_marker,
      vector const_data_cv_marker,
      vector knots_cv_marker,
      vector coeff_cv_marker,
      int spline_degree_cv_marker,
      array[] int functional_ops_cs_mean,
      vector const_data_cs_mean,
      vector knots_cs_mean,
      vector coeff_cs_mean,
      int spline_degree_cs_mean,
      array[] int functional_ops_cs_marker,
      vector const_data_cs_marker,
      vector knots_cs_marker,
      vector coeff_cs_marker,
      int spline_degree_cs_marker,
      array[] int id_start,
      array[] int id_end,
      real tmax
  ) {
    real lp = 0;

    for (ii in start : end) {
      int idx_local = ii - start + 1;
      int i = id_seq[idx_local];

      vector[Q_idm] wbar_i_i = wbar_i[i];

      // -------------------- Longitudinal likelihood per subject
      for (n in id_start[i] : id_end[i]) {
        int d = marker[n];
          vector[Q_idm] w_idscaled_i = w_idscaled[i, d];
        real eta_long = dot_product(X_obs[n], beta_scaled)
                        + dot_product(Z_id_obs[n], u_id[i])
                        + ((R_mk > 0) ? dot_product(Z_mk_obs[n], v_marker[d])
                           : 0.0)
                    + dot_product(Z_idm_obs[n], w_idscaled_i);

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
          lp += normal_lpdf(y_real[n] | eta_long, sig);
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
          lp += student_t_lpdf(y_real[n] | nu, eta_long, sig);
        } else if (family_long[d] == 3) {
          lp += bernoulli_logit_lpmf(y_int[n] | eta_long);
        } else if (family_long[d] == 4) {
          lp += binomial_logit_lpmf(y_int[n] | trials[n], eta_long);
        } else if (family_long[d] == 5) {
          lp += poisson_log_lpmf(y_int[n] | eta_long);
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
          lp += neg_binomial_2_log_lpmf(y_int[n] | eta_long, phi);
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
          lp += skew_normal_lpdf(y_real[n] | eta_long, sig, alpha);
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
          lp += double_exponential_lpdf(y_real[n] | eta_long, sig);
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
          lp += skew_double_exponential_lpdf(y_real[n] | eta_long, sig, tau_sde);
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
          lp += beta_lpdf(y_real[n] | shape1, shape2);
        } else {
          lp += ordered_logistic_lpmf(y_int[n] | eta_long, cutpoints_ord);
        }
      }

      // -------------------- Survival likelihood per subject
      // Cause-specific baseline hazards and covariates are combined below.

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

        cvk_now[j] = mk_part_now + dot_product(Z_idm_gk_now[i][j], wbar_i_i);
        cvk_fwd[j] = mk_part_fwd + dot_product(Z_idm_gk_fwd[i][j], wbar_i_i);
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

        real cvk_S = mk_part_S + dot_product(Z_idm_event_now[i], wbar_i_i);
        real cvk_S_fwd = mk_part_S_fwd + dot_product(Z_idm_event_fwd[i], wbar_i_i);

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
        lp += log_h0_S + eta_w_ev + eta_assoc_S;
      }

      // Total cumulative hazard across causes
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

      lp += -cumhaz(S_event[i], log_h_total, rep_vector(0.0, 15));
    }

    return lp;
  }
}

data {
  #include helper/data/joinme_fit_common.stan
  array[n_id] int<lower=1, upper=N> id_start;
  array[n_id] int<lower=1, upper=N> id_end;
  int<lower=1> grainsize;
}

transformed data {
  #include helper/transformed_data/fit_cov_index.stan

  array[n_id] int id_seq;
  for (i in 1 : n_id) {
    id_seq[i] = i;
  }
}

parameters {
  #include helper/parameters/joinme_fit_common.stan
}

transformed parameters {
  #include helper/transformed_parameters/fit_scaling_and_effects.stan
}

model {
  #include helper/model/fit_priors.stan
  
  // -------------------- Threaded likelihood
  // Each task evaluates a slice of subjects (id_seq), summing both longitudinal
  // and survival contributions. The grainsize controls chunk size; adjust based
  // on subject count and per-subject cost for best throughput.
  target += reduce_sum(
    partial_joinme,
    id_seq,
    grainsize,
    id,
    marker,
    y_real,
    y_int,
    trials,
    family_long,
    flag_resid_dim,
    X_obs,
    Z_id_obs,
    Z_mk_obs,
    Z_idm_obs,
    R_mk,
    Q_idm,
    beta_scaled,
    P_sigma,
    X_sigma,
    P_nu,
    X_nu,
    P_phi,
    X_phi,
    P_alpha,
    X_alpha,
    P_phi_beta,
    X_phi_beta,
    P_tau_sde,
    X_tau_sde,
    beta_sigma,
    beta_nu,
    beta_phi,
    beta_alpha,
    beta_phi_beta,
    beta_tau_sde,
    n_re_sigma,
    K_sigma,
    K_sigma_max,
    Z_sigma,
    J_sigma,
    tau_sigma,
    z_sigma,
    n_re_nu,
    K_nu,
    K_nu_max,
    Z_nu,
    J_nu,
    tau_nu,
    z_nu,
    n_re_phi,
    K_phi,
    K_phi_max,
    Z_phi,
    J_phi,
    tau_phi,
    z_phi,
    n_re_alpha,
    K_alpha,
    K_alpha_max,
    Z_alpha,
    J_alpha,
    tau_alpha,
    z_alpha,
    n_re_phi_beta,
    K_phi_beta,
    K_phi_beta_max,
    Z_phi_beta,
    J_phi_beta,
    tau_phi_beta,
    z_phi_beta,
    n_re_tau_sde,
    K_tau_sde,
    K_tau_sde_max,
    Z_tau_sde,
    J_tau_sde,
    tau_tau_sde,
    z_tau_sde,
    u_id,
    v_marker,
    w_idscaled,
    wbar_i,
    sigma_y,
    sigma_marker,
    nu_marker,
    phi_nb_marker,
    alpha_skew_marker,
    phi_beta_marker,
    tau_sde_marker,
    K_ord,
    cutpoints_ord,
    p_w,
    W,
    gamma_w,
    log_h0_intercept,
    bs_gamma_c,
    Kbs,
    Bs_event_c,
    Bs_gk_c,
    S_event,
    d_event,
    K_event,
    event_type,
    eps_fd,
    X_gk_now,
    X_gk_fwd,
    Z_id_gk_now,
    Z_id_gk_fwd,
    X_event_now,
    X_event_fwd,
    Z_id_event_now,
    Z_id_event_fwd,
    Z_mk_gk_now,
    Z_mk_gk_fwd,
    Z_idm_gk_now,
    Z_idm_gk_fwd,
    Z_mk_event_now,
    Z_mk_event_fwd,
    Z_idm_event_now,
    Z_idm_event_fwd,
    vbar,
    L_i,
    a_cv_total,
    a_cs_total,
    a_cv_mean,
    a_cv_marker,
    a_cs_mean,
    a_cs_marker,
    a_vcov_var,
    tf_mode_cv_tot,
    tf_mode_cs_tot,
    tf_mode_vcov,
    tf_mode_cv_mean,
    tf_mode_cv_marker,
    tf_mode_cs_mean,
    tf_mode_cs_marker,
    functional_ops_cv,
    const_data_cv,
    knots_cv,
    coeff_cv,
    spline_degree_cv,
    functional_ops_cs,
    const_data_cs,
    knots_cs,
    coeff_cs,
    spline_degree_cs,
    functional_ops_vcov,
    const_data_vcov,
    knots_vcov,
    coeff_vcov,
    spline_degree_vcov,
    functional_ops_cv_mean,
    const_data_cv_mean,
    knots_cv_mean,
    coeff_cv_mean,
    spline_degree_cv_mean,
    functional_ops_cv_marker,
    const_data_cv_marker,
    knots_cv_marker,
    coeff_cv_marker,
    spline_degree_cv_marker,
    functional_ops_cs_mean,
    const_data_cs_mean,
    knots_cs_mean,
    coeff_cs_mean,
    spline_degree_cs_mean,
    functional_ops_cs_marker,
    const_data_cs_marker,
    knots_cs_marker,
    coeff_cs_marker,
    spline_degree_cs_marker,
    id_start,
    id_end,
    tmax
  );
}

generated quantities {
  #include helper/generated_quantities/fit_outputs.stan
}
