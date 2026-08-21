/**
 * @file survival/generated_quantities/fit_declarations.stan
 * @brief Declare subject-level event likelihood, cumulative hazard and survival.
 *
 * @details
 * The declaration file introduces one output per subject.  The calculation file reconstructs transformed longitudinal features at quadrature nodes and event endpoints before integrating the cause-specific hazard.
 * All generated-quantity declarations are assembled before any calculation,
 * as required by the Stan language.
 */
  vector[n_id] log_lik_surv; // subject-level full event log likelihood, summed across all risk rows
  vector[n_id] cumhaz_event; // all-cause cumulative hazard across the subject's represented risk intervals
  vector[n_id] surv_prob_event; // exponential of minus cumhaz_event at the final represented upper limit

