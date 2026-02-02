/**
 * @file joinme_fit.stan
 * @brief Joint mixed-effects model for multivariate long-format longitudinal data + survival,
 *        with internal time scaling and corrected association transforms.
 *
 * ## What this model does
 * - Longitudinal submodel:
 *   - Handles D biomarkers in long format (irregular measurement schedules).
 *   - Fixed effects: X_obs * beta
 *   - Random effects:
 *     (1) id-level RE u_id[i] from Z_id_obs
 *     (2) optional marker-only RE v_marker[d] from Z_mk_obs
 *     (3) optional marker-by-id coefficients w_idscaled[i,d] from Z_idm_obs
 *         where w_idscaled[i,d] = L_i[i] * z_w[i,d] (id-specific latent effects)
 *
 * - Survival submodel:
 *   - Hazard: h_i(t) = h0(t) * exp( W_i * gamma_w + eta_assoc_i(t) )
 *   - Baseline log-hazard uses a centered spline basis.
 *   - Association uses:
 *     - current value total: CV_tot = CV_mean + mean(CV_marker)
 *     - current slope total: CS_tot = CS_mean + mean(CS_marker)
 *     - vcov feature (variance-only from id-specific covariance)
 *
 * @note This program requires tmax and index lists telling Stan which coefficients correspond to time.
 */

functions {
  #include helper/functions/eta_fd.stanfunctions
  
  #include helper/functions/eta_vcov_varonly_weighted_const.stanfunctions
  
  #include helper/functions/cumhaz.stanfunctions
  
  #include helper/functions/functional_transform.stanfunctions
  
  #include helper/functions/basis_functions.stanfunctions
  
  #include helper/functions/composite_transform.stanfunctions

}
data {
  #include helper/data/joinme_fit_common.stan
}
transformed data {
  #include helper/transformed_data/fit_cov_index.stan
}
parameters {
  #include helper/parameters/joinme_fit_common.stan
}
transformed parameters {
  #include helper/transformed_parameters/fit_scaling_and_effects.stan
}
model {
  #include helper/model/fit_priors.stan
  #include helper/model/fit_longitudinal_likelihood.stan
  #include helper/model/fit_survival_likelihood.stan
}
generated quantities {
  #include helper/generated_quantities/fit_outputs.stan
}
