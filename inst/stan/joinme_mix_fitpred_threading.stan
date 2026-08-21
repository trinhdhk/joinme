/**
 * @file joinme_mix_fitpred_threading.stan
 * @brief Prediction for fitted subjects in a latent-class joint model.
 *
 * @details
 * This programme evaluates longitudinal and survival
 * predictions from population parameters and trained random effects retained
 * in paired posterior draws from `joinme_mix()`.  It does not sample a new
 * latent effect or a new class allocation.  The fitted subject and marker
 * allocation probabilities are therefore retained draw by draw by the R
 * prediction method rather than recalculated in Stan. 
 */

functions {
  #include include/etc/functions/eta_fd.stanfunctions
  #include include/etc/functions/eta_chol_corr.stanfunctions
  #include include/etc/functions/eta_vcov_weighted_const.stanfunctions
  #include include/etc/functions/cumhaz.stanfunctions
  #include include/etc/bytecode/interpreter.stanfunctions
  #include include/etc/functions/link_functions.stanfunctions
  #include include/etc/functions/basis_functions.stanfunctions
  #include include/etc/functions/composite_transform.stanfunctions
}

data {
  #include include/etc/data/dynamic_prediction_draw_data.stan
  #include include/submodels/longitudinal/data/dynamic_prediction.stan
  #include include/etc/data/fitted_random_effect_draw_data.stan
  #include include/etc/data/fixed_prediction_latent_placeholders.stan
  #include include/submodels/marker_weight/data/dynamic_prediction.stan
  #include include/submodels/survival/data/dynamic_prediction.stan
  #include include/submodels/assoc/data/dynamic_prediction.stan
  #include include/submodels/functional/data/dynamic_prediction.stan
}

transformed data {
  #include include/submodels/longitudinal/transformed_data/dynamic_prediction.stan
  #include include/submodels/marker_weight/transformed_data/dynamic_prediction.stan
  #include include/submodels/survival/transformed_data/dynamic_prediction.stan
  #include include/submodels/assoc/transformed_data/dynamic_prediction.stan
  #include include/submodels/functional/transformed_data/dynamic_prediction.stan
  #include include/etc/transformed_data/dynamic_prediction_draw_indices.stan
}

generated quantities {
  #include include/submodels/longitudinal/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/survival/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/assoc/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/marker_weight/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/functional/generated_quantities/dynamic_prediction_declarations.stan

  #include include/submodels/longitudinal/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/survival/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/assoc/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/marker_weight/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/functional/generated_quantities/dynamic_prediction_calculations.stan
  #include include/etc/generated_quantities/dynpred_output_calculations.stan
}
