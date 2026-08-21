/**
 * @file survival/model/fit.stan
 * @brief Assign priors to cause-specific baseline hazards.
 *
 * @details
 * The first spline basis coefficient supplies the baseline log-hazard level.  Remaining coefficients receive regularising priors, with second-difference penalties controlling baseline-hazard curvature.
 */
  /* Baseline hazard priors (per event type) */
  for (k_ev in 1 : K_event) { // event-specific baseline hazard
    // Intercept is encoded in the first basis column (constant 1s).
    bs_gamma_c[k_ev][1] ~ student_t(3, -5, 6);
    if (Kbs > 1)
      bs_gamma_c[k_ev][2:Kbs] ~ student_t(3, 0, 1);
    if (Kbs >= 4) // penalised spline second differences (exclude intercept)
      for (k in 4 : Kbs)
        target += normal_lpdf(
                              bs_gamma_c[k_ev][k] - 2 * bs_gamma_c[k_ev][k - 1]
                              + bs_gamma_c[k_ev][k - 2] | 0, tau_spline);
  }
  


