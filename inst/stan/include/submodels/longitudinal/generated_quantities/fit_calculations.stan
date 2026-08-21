/**
 * @file longitudinal/generated_quantities/fit_calculations.stan
 * @brief Calculate longitudinal effective scales and pointwise likelihood.
 *
 * @details
 * The declaration file introduces effective coefficient summaries and the longitudinal pointwise log-likelihood vector.  The calculation file rebuilds each observation under its marker-specific response family.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  /* -------------------- Longitudinal log-likelihood per observation */
  for (n in 1 : N) {
    int i = id[n]; // Role: individual loop position.
    int d = marker[n]; // Role: marker loop position.
     real eta_long = dot_product(X_obs[n], beta)
                    + dot_product(Z_id_obs[n], u_id[i])
                    + ((R_mk > 0) ? dot_product(Z_mk_obs[n], v_marker[d])
                       : 0.0)
                    + dot_product(Z_idm_obs[n], w_idm[i, d]);
     // eta_long: linear predictor for longitudinal outcome
    int link_d = make_canonical_link(d, inv_link_n_ops, inv_link_ops, inv_link_n_const); // Role: link d.
    real mu_long = inv_link_bytecode(eta_long, d, inv_link_n_ops, inv_link_ops, inv_link_n_const, inv_link_const); // Role: mu long.

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
                 ? exp(fmin(eta_sigma, 20))
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
                 ? exp(fmin(eta_sigma, 20))
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
      real nu = (P_nu > 0 || n_re_nu > 0) ? (2 + exp(fmin(eta_nu, 20))) : nu_family[marker_to_nu_family[d]]; // Role: degrees of freedom.
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
      real phi = (P_phi > 0 || n_re_phi > 0) ? exp(eta_phi) : phi_family[marker_to_phi_family[d]]; // Role: precision.
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
                 ? exp(fmin(eta_sigma, 20))
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
      real alpha = (P_alpha > 0 || n_re_alpha > 0) ? eta_alpha : alpha_family[marker_to_alpha_family[d]]; // Role: shape.
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
                 ? exp(fmin(eta_sigma, 20))
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
                 ? exp(fmin(eta_sigma, 20))
                 : sigma_family[marker_to_sigma_family[d]];
      // sig: residual scale for skew Laplace
      real eta_tau = 0; // linear predictor for tau
      if (P_tau > 0) eta_tau += dot_product(X_tau[n], beta_tau);
      if (n_re_tau > 0) {
        for (j in 1 : n_re_tau) {
          vector[K_tau[j]] b = (to_vector(z_tau[j][J_tau[j, n], 1:K_tau[j]])
                                    .* tau_tau[j][1:K_tau[j]]);
          // b: random-effect coefficients for tau term j
          eta_tau += dot_product(Z_tau[j][n, 1:K_tau[j]], b);
        }
      }
      // Step 1: preserve a fixed quantile/asymmetry parameter in pointwise
      // likelihood quantities exactly as it is preserved in the model block.
      // Step 2: otherwise evaluate the fitted regression or family-level value.
      real tau = use_tau_fixed[d] == 1
                 ? tau_fixed[d]
                 : ((P_tau > 0 || n_re_tau > 0)
                    ? inv_logit(eta_tau)
                    : tau_family[marker_to_tau_family[d]]);
      // tau: skewness parameter in (0,1)
      log_lik_long[n] = skew_double_exponential_lpdf(y_real[n] | mu_long, sig, tau);
    } else if (family_long[d] == 10) {
      real eta_kappa = 0; // linear predictor for log kappa
      if (P_kappa > 0) eta_kappa += dot_product(X_kappa[n], beta_kappa);
      if (n_re_kappa > 0) {
        for (j in 1 : n_re_kappa) {
          vector[K_kappa[j]] b = (to_vector(z_kappa[j][J_kappa[j, n], 1:K_kappa[j]])
                                     .* tau_kappa[j][1:K_kappa[j]]);
          // b: random-effect coefficients for kappa term j
          eta_kappa += dot_product(Z_kappa[j][n, 1:K_kappa[j]], b);
        }
      }
      real kappa = (P_kappa > 0 || n_re_kappa > 0) ? exp(eta_kappa) : kappa_family[marker_to_kappa_family[d]]; // Role: sample size.
      // Step 1: keep the response-scale mean strictly inside its mathematical
      // support, matching the fitted likelihood calculation.
      real mu = fmin(fmax(mu_long, 1e-12), 1 - 1e-12); // Role: mu.

      // Step 2: construct both Beta shapes from the same mean/sample-size
      // parameterisation used for estimation and posterior prediction.
      real shape1 = fmax(mu * kappa, 1e-6); // Role: shape1.
      real shape2 = fmax((1 - mu) * kappa, 1e-6); // Role: shape2.
      log_lik_long[n] = beta_lpdf(y_real[n] | shape1, shape2);
    } else {
      log_lik_long[n] = ordered_logistic_lpmf(y_int[n] | eta_long, cutpoints_ord);
    }
  }


