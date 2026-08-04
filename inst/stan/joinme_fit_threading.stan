/**
 * @file joinme_fit_threading.stan
 * @brief Threaded implementation of the Joint Mixed-Effects Model using reduce_sum.
 *
 * @details
 * This is the ordinary, non-mixture fitting programme used for both serial and
 * parallel execution. It uses `reduce_sum` to parallelize the likelihood
 * computation over subjects when more than one thread is requested.
 *
 * It uses the shared likelihood function `partial_joinme` so the same kernel can
 * run with one thread or many threads without changing the target density.
 *
 * @see joinme::joinme
 */

functions {
  #include helper/functions/eta_fd.stanfunctions 
  #include helper/functions/eta_chol_corr.stanfunctions
  #include helper/functions/eta_vcov_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/bytecode/interpreter.stanfunctions
  #include helper/functions/link_functions.stanfunctions 
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions
  #include helper/functions/prior_families.stanfunctions
  #include helper/functions/joinme_fit_partial.stanfunctions
}

data {
  #include helper/data/fit_data.stan
  array[n_id] int<lower=1, upper=N> id_start; // first longitudinal row belonging to each individual
  array[n_id] int<lower=1, upper=N> id_end; // final longitudinal row belonging to each individual
  int<lower=1> grainsize; // number of individuals assigned to each reduce-sum slice
}

transformed data {
  #include helper/transformed_data/fit_cov_index.stan

  array[n_id] int id_seq; // consecutive individual indices passed to reduce_sum
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

  #include helper/model/fit_threaded_likelihood.stan

}

generated quantities {
  #include helper/generated_quantities/fit_outputs.stan
}
