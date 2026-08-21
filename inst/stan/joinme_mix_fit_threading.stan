/**
 * @file joinme_mix_fit_threading.stan
 * @brief Threaded latent-progress mixture implementation using reduce_sum.
 *
 * @details
 * This is the thread-capable joinme fitting program used for both serial and
 * parallel execution. It uses `reduce_sum` to parallelize the likelihood
 * computation over subjects when more than one thread is requested.
 *
 * It uses the shared likelihood function `partial_joinme` so the same kernel can
 * run with one thread or many threads without changing the target density.
 *
 * @see joinme::joinme_mix
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
  #include include/etc/functions/prior_families.stanfunctions
  #include include/submodels/latent_class/functions/component_density.stanfunctions
  #include include/submodels/latent_class/functions/class_probability.stanfunctions
  #include include/etc/functions/joinme_fit_partial.stanfunctions
}

data {
  // The master programme assembles scientific data before class additions.
  #include include/submodels/longitudinal/data/fit.stan
  #include include/submodels/marker_weight/data/fit.stan
  #include include/submodels/survival/data/fit.stan
  #include include/submodels/assoc/data/fit.stan
  #include include/submodels/functional/data/fit.stan
  #include include/etc/data/regression_prior_data.stan
  #include include/submodels/latent_class/data/fit.stan
  array[n_id] int<lower=1, upper=N> id_start; // first longitudinal row belonging to each individual
  array[n_id] int<lower=1, upper=N> id_end; // final longitudinal row belonging to each individual
  int<lower=1> grainsize; // number of individuals assigned to each reduce-sum slice
}

transformed data {
  #include include/submodels/longitudinal/transformed_data/fit.stan
  #include include/submodels/marker_weight/transformed_data/fit.stan
  #include include/submodels/survival/transformed_data/fit.stan
  #include include/submodels/assoc/transformed_data/fit.stan
  #include include/submodels/functional/transformed_data/fit.stan
  #include include/submodels/latent_class/transformed_data/fit.stan

  array[n_id] int id_seq; // consecutive individual indices passed to reduce_sum
  for (i in 1 : n_id) {
    id_seq[i] = i;
  }
}

parameters {
  #include include/submodels/longitudinal/parameters/fit.stan
  #include include/submodels/survival/parameters/fit.stan
  #include include/submodels/assoc/parameters/fit.stan
  #include include/submodels/marker_weight/parameters/fit.stan
  #include include/submodels/functional/parameters/fit.stan
  #include include/etc/parameters/regression_prior_parameters.stan
  #include include/submodels/latent_class/parameters/fit.stan
}

transformed parameters {
  #include include/etc/transformed_parameters/regression_coefficients.stan
  #include include/submodels/longitudinal/transformed_parameters/fit.stan
  #include include/submodels/marker_weight/transformed_parameters/fit.stan
  #include include/submodels/functional/transformed_parameters/fit.stan
  #include include/submodels/assoc/transformed_parameters/fit.stan
  #include include/submodels/survival/transformed_parameters/fit.stan
  #include include/submodels/latent_class/transformed_parameters/fit.stan
}

model {
  #include include/submodels/longitudinal/model/fit.stan
  #include include/etc/model/regression_priors.stan
  #include include/submodels/survival/model/fit.stan
  #include include/submodels/assoc/model/fit.stan
  #include include/submodels/marker_weight/model/fit.stan
  #include include/submodels/functional/model/fit.stan
  #include include/submodels/latent_class/model/fit.stan
  #include include/etc/model/fit_threaded_likelihood.stan
}

generated quantities {
  #include include/submodels/longitudinal/generated_quantities/fit_declarations.stan
  #include include/submodels/survival/generated_quantities/fit_declarations.stan
  #include include/submodels/assoc/generated_quantities/fit_declarations.stan
  #include include/submodels/marker_weight/generated_quantities/fit_declarations.stan
  #include include/submodels/functional/generated_quantities/fit_declarations.stan

  #include include/submodels/longitudinal/generated_quantities/fit_calculations.stan
  #include include/submodels/survival/generated_quantities/fit_calculations.stan
  #include include/submodels/assoc/generated_quantities/fit_calculations.stan
  #include include/submodels/marker_weight/generated_quantities/fit_calculations.stan
  #include include/submodels/functional/generated_quantities/fit_calculations.stan
  #include include/submodels/latent_class/generated_quantities/fit.stan
}
