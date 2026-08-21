/**
 * @file include/etc/data/dynamic_prediction_draw_data.stan
 * @brief Declare the retained posterior-draw axis.
 *
 * @details
 * All dynamic-prediction submodels are evaluated for the same retained draws,
 * so the master programme introduces this common draw axis once.
 */
  int<lower=1> n_draws; // number of posterior draws to process
  int<lower=1> grainsize; // retained posterior draws assigned to each reduce-sum slice
