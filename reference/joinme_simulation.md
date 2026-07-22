# Simulation for JoiNMe

Simulator for the JoiNMe joint model aligned with the v11 Stan
semantics.

Key features:

- Multivariate longitudinal outcomes in long format with irregular
  times.

- Survival hazard depends on CV_total(t) = CV_mean(t) + CV_marker(t).

- Survival times drawn by inverse transform sampling using:

  - numeric integration (integrate) for cumulative hazard

  - numeric root solving (uniroot) for event time.
