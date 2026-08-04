/**
 * @file fit_mixture_locations.stan
 * @brief Present baseline probabilities and locations in established layouts.
 *
 * @details Only the selected random-intercept coordinate is ordered by
 * location. Other locations, including slopes, remain unconstrained. When
 * probability ordering is requested, ordered negative log-odds followed by a
 * zero reference log-odds give strictly increasing baseline probabilities.
 */
simplex[n_classes] mix_probability; // baseline class probabilities in reporting layout
matrix[n_classes, K_mix] mix_location; // ordered and unrestricted class locations in reporting layout
vector[P_class_subject] mix_class_coefficient_subject = rep_vector(0, P_class_subject); // subject-domain formulaClass coefficients on their multinomial-logit scale
vector[P_class_marker] mix_class_coefficient_marker = rep_vector(0, P_class_marker); // marker-domain formulaClass coefficients on their multinomial-logit scale
if (mix_ordering == 2) {
  vector[n_classes] baseline_logit; // ordered additive-log-ratio coordinates plus reference
  for (group in 1 : (n_classes - 1)) {
    baseline_logit[group] =
      -mix_probability_ordered_gap[1][n_classes - group];
  }
  baseline_logit[n_classes] = 0;
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

{
  vector[P_class_subject + P_class_marker] mix_class_coefficient_effect = joinme_prior_transform(
    mix_class_coefficient_raw,
    prior_class_regression_family,
    prior_class_regression_mu,
    prior_class_regression_scale,
    horseshoe_local_class_regression,
    horseshoe_global_class_regression,
    horseshoe_slab_class_regression,
    prior_class_regression_slab_scale
  ); // complete class-regression vector after its independently selected prior transformation

  if (P_class_subject > 0) {
    mix_class_coefficient_subject = mix_class_coefficient_effect[1 : P_class_subject];
  }
  if (P_class_marker > 0) {
    mix_class_coefficient_marker = mix_class_coefficient_effect[
      P_class_subject + 1 : P_class_subject + P_class_marker
    ];
  }
}
