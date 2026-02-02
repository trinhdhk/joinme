  /**
   * @brief Longitudinal likelihood for each observation (very explicit).
   *
   * We want the logic to be simple enough to explain to a child:
   *   - First, build one prediction number (eta_long) for a row.
   *   - Then pick the correct distribution for the marker.
   *   - Finally, add the probability of the observed value.
   *
   * Family codes:
   *  1 = Gaussian (Normal)
   *  2 = Student‑t
   *  3 = Bernoulli (logit)
   *  4 = Binomial (logit)
   *  5 = Poisson (log)
   *  6 = Negative binomial 2 (log)
   *  7 = Skew normal
   *  8 = Double exponential (Laplace)
   *  9 = Skew double exponential (asymmetric Laplace)
   * 10 = Beta (mean via logit link, precision via phi_beta)
   * 11 = Cumulative logit (ordered logistic)
   */

  for (n in 1 : N) {
    // Identify which subject and which marker this row belongs to
    int i = id[n];
    int d = marker[n];

    // Build the linear predictor (the “center” of the distribution)
    real eta_long = dot_product(X_obs[n], beta_scaled)
                    + dot_product(Z_id_obs[n], u_id[i])
                    + ((R_mk > 0)
                       ? dot_product(Z_mk_obs[n], v_marker[d])
                       : 0.0)
                    + dot_product(Z_idm_obs[n], w_idscaled[i, d]);

    // Choose the likelihood based on family code
    if (family_long[d] == 1) {
      // Gaussian: y ~ Normal(eta_long, sigma)
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
      y_real[n] ~ normal(eta_long, sig);
    } else if (family_long[d] == 2) {
      // Student‑t: y ~ Student_t(nu, eta_long, sigma)
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
      y_real[n] ~ student_t(nu, eta_long, sig);
    } else if (family_long[d] == 3) {
      // Bernoulli (logit): P(y=1) = inv_logit(eta_long)
      y_int[n] ~ bernoulli_logit(eta_long);
    } else if (family_long[d] == 4) {
      // Binomial (logit): P(y) = Binomial(trials, inv_logit(eta_long))
      y_int[n] ~ binomial_logit(trials[n], eta_long);
    } else if (family_long[d] == 5) {
      // Poisson (log): E[y] = exp(eta_long)
      y_int[n] ~ poisson_log(eta_long);
    } else if (family_long[d] == 6) {
      // NegBin2 (log): E[y] = exp(eta_long), dispersion=phi
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
      y_int[n] ~ neg_binomial_2_log(eta_long, phi);
    } else if (family_long[d] == 7) {
      // Skew normal: location=eta_long, scale=sigma, skew=alpha
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
      y_real[n] ~ skew_normal(eta_long, sig, alpha);
    } else if (family_long[d] == 8) {
      // Double exponential (Laplace): location=eta_long, scale=sigma
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
      y_real[n] ~ double_exponential(eta_long, sig);
    } else if (family_long[d] == 9) {
      // Skew double exponential: location=eta_long, scale=sigma, skew=tau_sde
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
      real tau_sde = (P_tau_sde > 0 || n_re_tau_sde > 0)
                     ? inv_logit(eta_tau)
                     : tau_sde_marker[d];
      y_real[n] ~ skew_double_exponential(eta_long, sig, tau_sde);
    } else if (family_long[d] == 10) {
      // Beta: mean=inv_logit(eta_long), precision=phi_beta
      real eta_phi_beta = 0;
      if (P_phi_beta > 0) eta_phi_beta += dot_product(X_phi_beta[n], beta_phi_beta);
      if (n_re_phi_beta > 0) {
        for (j in 1 : n_re_phi_beta) {
          vector[K_phi_beta[j]] b = (to_vector(z_phi_beta[j][J_phi_beta[j, n], 1:K_phi_beta[j]])
                                     .* tau_phi_beta[j][1:K_phi_beta[j]]);
          eta_phi_beta += dot_product(Z_phi_beta[j][n, 1:K_phi_beta[j]], b);
        }
      }
      real phi_beta = (P_phi_beta > 0 || n_re_phi_beta > 0)
                      ? exp(eta_phi_beta)
                      : phi_beta_marker[d];
      real mu = inv_logit(eta_long);
      real shape1 = fmax(mu * phi_beta, 1e-6);
      real shape2 = fmax((1 - mu) * phi_beta, 1e-6);
      y_real[n] ~ beta(shape1, shape2);
    } else {
      // Cumulative logit (ordered logistic): y in {1..K_ord}
      y_int[n] ~ ordered_logistic(eta_long, cutpoints_ord);
    }
  }
