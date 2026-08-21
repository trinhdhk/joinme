/**
 * @file joinme_mix_dynpred_threading.stan
 * @brief Thread-capable latent-class dynamic prediction.
 *
 * ### Inputs
 * - Includes `grainsize` to size draw slices for reduce_sum.
 *
 * ### Outputs
 * - The same conditional predictions regardless of thread count; threading only changes how work is split, not what is computed.
 *
 * ### Execution story
 * - reduce_sum splits posterior draws into slices; each slice calls `partial_draw` (pure function) so math stays identical.
 * - No random numbers appear inside tasks; determinism holds given seed and thread topology.
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
  #include include/submodels/latent_class/functions/component_density.stanfunctions

  #include include/etc/functions/joinme_dynpred_partial.stanfunctions
}

data {
  #include include/etc/data/dynamic_prediction_draw_data.stan
  #include include/submodels/longitudinal/data/dynamic_prediction.stan
  #include include/etc/data/fitted_random_effect_draw_data.stan
  #include include/submodels/marker_weight/data/dynamic_prediction.stan
  #include include/submodels/survival/data/dynamic_prediction.stan
  #include include/submodels/assoc/data/dynamic_prediction.stan
  #include include/submodels/functional/data/dynamic_prediction.stan
  #include include/submodels/latent_class/data/dynamic_prediction.stan
}

transformed data {
  #include include/submodels/longitudinal/transformed_data/dynamic_prediction.stan
  #include include/submodels/marker_weight/transformed_data/dynamic_prediction.stan
  #include include/submodels/survival/transformed_data/dynamic_prediction.stan
  #include include/submodels/assoc/transformed_data/dynamic_prediction.stan
  #include include/submodels/functional/transformed_data/dynamic_prediction.stan
  #include include/submodels/latent_class/transformed_data/dynamic_prediction.stan
  #include include/etc/transformed_data/dynamic_prediction_draw_indices.stan
}

parameters {
  #include include/submodels/longitudinal/parameters/dynamic_prediction.stan
  #include include/submodels/marker_weight/parameters/dynamic_prediction.stan
  #include include/submodels/survival/parameters/dynamic_prediction.stan
  #include include/submodels/assoc/parameters/dynamic_prediction.stan
  #include include/submodels/functional/parameters/dynamic_prediction.stan
  #include include/submodels/latent_class/parameters/dynamic_prediction.stan
  #include include/etc/parameters/dynamic_prediction_subject_effects.stan
}

transformed parameters {
  #include include/submodels/longitudinal/transformed_parameters/dynamic_prediction.stan
  #include include/submodels/marker_weight/transformed_parameters/dynamic_prediction.stan
  #include include/submodels/survival/transformed_parameters/dynamic_prediction.stan
  #include include/submodels/assoc/transformed_parameters/dynamic_prediction.stan
  #include include/submodels/functional/transformed_parameters/dynamic_prediction.stan
  #include include/submodels/latent_class/transformed_parameters/dynamic_prediction.stan
}

model {
  #include include/submodels/longitudinal/model/dynamic_prediction.stan
  #include include/submodels/marker_weight/model/dynamic_prediction.stan
  #include include/submodels/survival/model/dynamic_prediction.stan
  #include include/submodels/assoc/model/dynamic_prediction.stan
  #include include/submodels/functional/model/dynamic_prediction.stan
  #include include/etc/model/dynpred_threaded_likelihood.stan
  #include include/submodels/latent_class/model/dynamic_prediction.stan
}

generated quantities {
  #include include/submodels/longitudinal/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/survival/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/assoc/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/marker_weight/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/functional/generated_quantities/dynamic_prediction_declarations.stan
  #include include/submodels/latent_class/generated_quantities/dynamic_prediction_declarations.stan

  #include include/submodels/latent_class/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/longitudinal/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/survival/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/assoc/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/marker_weight/generated_quantities/dynamic_prediction_calculations.stan
  #include include/submodels/functional/generated_quantities/dynamic_prediction_calculations.stan
  #include include/etc/generated_quantities/dynpred_output_calculations.stan
}
