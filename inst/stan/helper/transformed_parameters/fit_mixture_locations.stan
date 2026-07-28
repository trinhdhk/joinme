/**
 * @file fit_mixture_locations.stan
 * @brief Present baseline probabilities and locations in established layouts.
 *
 * @details Only the selected random-intercept coordinate is ordered by
 * location. Other locations, including slopes, remain unconstrained. When
 * probability ordering is requested, ordered negative log-odds followed by a
 * zero reference log-odds give strictly increasing baseline probabilities.
 */
simplex[n_clusters] mix_probability; // baseline class probabilities in reporting layout
matrix[n_clusters, K_mix] mix_location; // class locations in reporting layout
if (mix_ordering == 2) {
  vector[n_clusters] baseline_logit; // ordered additive-log-ratio coordinates plus reference
  for (group in 1 : (n_clusters - 1)) {
    baseline_logit[group] =
      -mix_probability_ordered_gap[1][n_clusters - group];
  }
  baseline_logit[n_clusters] = 0;
  mix_probability = softmax(baseline_logit);
} else {
  mix_probability = mix_probability_unordered[1];
}

{
  int next_unordered_coordinate = 1; // next compact unrestricted-location column
  for (coordinate in 1 : K_mix) {
    if (coordinate == mix_ordered_location_coordinate) {
      mix_location[, coordinate] = mix_location_ordered[1];
    } else {
      mix_location[, coordinate] =
        mix_location_unordered[, next_unordered_coordinate];
      next_unordered_coordinate += 1;
    }
  }
}
