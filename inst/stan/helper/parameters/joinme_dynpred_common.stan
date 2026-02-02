  // Subject-specific random effects (standard normal to be scaled)
  array[n_draws] vector[n_random_id] z_u;
  array[n_draws, n_marker_types] vector[n_random_marker] z_v;
  array[n_draws, n_marker_types] vector[n_random_marker_id] z_w_lat;
  vector[n_draws] z_L;
