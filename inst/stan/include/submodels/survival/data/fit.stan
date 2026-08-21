/**
 * @file survival/data/fit.stan
 * @brief Data declarations for the cause-specific event-time submodel.
 *
 * @details
 * This fragment declares event covariates, the centred baseline-hazard basis,
 * quadrature bases, censoring representation and competing-risk labels.  The
 * event rows share the subject and interval indexing introduced by the
 * longitudinal fragment.
 */
  /* Hazard covariates */
  int<lower=0> p_w;              // number of ordinary Cox regression covariates after removal of the formula intercept
  matrix[N_event, p_w] W;        // Cox design at each interval endpoint, evaluated on original study time

  /* Baseline hazard spline (centred basis) */
  int<lower=1> Kbs;              // number of baseline hazard basis functions
  matrix[N_event, Kbs] Bs_event_c;  // centred baseline-hazard basis at original-time interval endpoints
  int<lower=1> n_gk;             // quadrature node count (nodes/weights hardcoded in Stan)
  array[N_event] matrix[n_gk, Kbs] Bs_gk_c; // centred baseline-hazard basis at original-time quadrature ordinates
  array[N_event] matrix[n_gk, p_w] W_gk; // Cox covariate design at quadrature nodes, evaluated on original study time
  real<lower=0> tau_spline;      // spline penalty scale for baseline hazard

  /* Survival outcomes */
  vector<lower=0, upper=1>[N_event] S_entry; // scaled lower limits of the risk rows; an added interval2 survival row begins at zero
  vector<lower=0, upper=1>[N_event] S_event; // scaled upper limits of the risk rows; every upper limit is strictly above its lower limit
  array[N_event] int<lower=0, upper=1> d_event; // one only when an exact event supplies a log-hazard contribution at the upper limit
  array[N_event] int<lower=0, upper=3> event_censor_type; // 0 survival through the interval, 1 exact event at its end, 2 failure in (0,R], 3 failure in (L,R]
  int<lower=1> K_event;                    // number of competing risks
  array[N_event] int<lower=1, upper=K_event> event_type; // event type per interval endpoint
