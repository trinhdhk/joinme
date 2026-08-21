/**
 * @file survival/parameters/fit.stan
 * @brief Parameters for the baseline hazard and event covariate effects.
 *
 * @details
 * Each cause has its own baseline-hazard spline coefficients and event-covariate slope vector.  Their scientific priors are assigned in the survival model fragment.
 */
  /* baseline hazard (cause-specific) */
  array[K_event] vector[Kbs] bs_gamma_c; // spline coefficients for baseline hazard (includes intercept basis)

  /* hazard covariates (cause-specific) */
  array[K_event] vector[p_w] z_gamma_w; // standardised event-covariate slopes, transformed by the common prior programme


