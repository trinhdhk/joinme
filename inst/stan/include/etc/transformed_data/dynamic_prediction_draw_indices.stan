/**
 * @file include/etc/transformed_data/dynamic_prediction_draw_indices.stan
 * @brief Construct the retained-draw sequence used by reduce_sum.
 *
 * @details Consecutive indices let the threaded function retrieve every fitted
 * quantity and new-subject latent vector from the same draw position.
 */
array[n_draws] int draw_ids; // consecutive retained-fit draw indices
for (draw_index in 1 : n_draws)
  draw_ids[draw_index] = draw_index;
