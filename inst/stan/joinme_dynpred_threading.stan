/**
 * @file joinme_dynpred_threading.stan
 * @brief Ordinary thread-capable dynamic prediction without latent classes.
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
  #include helper/functions/eta_fd.stanfunctions
  #include helper/functions/eta_chol_corr.stanfunctions
  #include helper/functions/eta_vcov_weighted_const.stanfunctions
  #include helper/functions/cumhaz.stanfunctions
  #include helper/bytecode/interpreter.stanfunctions
  #include helper/functions/link_functions.stanfunctions
  #include helper/functions/basis_functions.stanfunctions
  #include helper/functions/composite_transform.stanfunctions

  #include helper/functions/joinme_dynpred_partial.stanfunctions
}

data {
  #include helper/data/dynpred_data.stan
  int<lower=1> grainsize; // number of fitted draws assigned to each reduce-sum slice

  int flag_indep_id_re; // whether shared-individual effects are mutually independent
  int flag_indep_marker_re; // whether shared-marker effects are mutually independent
  int flag_indep_idmarker_cov; // whether marker-by-individual covariance is diagonal
  int flag_allow_marker_crosscorr; // whether marker-specific effects may correlate across markers

  array[num_unique_cov_entries] int idx_row_cov; // row index for each unique covariance entry
  array[num_unique_cov_entries] int idx_col_cov; // column index for each unique covariance entry
}

transformed data {
  array[n_draws] int draw_ids; // consecutive retained-fit draw indices passed to reduce_sum
  for (i in 1 : n_draws) draw_ids[i] = i;
}

parameters {
  #include helper/parameters/joinme_dynpred_common.stan
}

model {
  #include helper/model/dynpred_threaded_likelihood.stan

}

generated quantities {
  #include helper/generated_quantities/dynpred_output_declarations.stan
  #include helper/generated_quantities/dynpred_output_calculations.stan
}
