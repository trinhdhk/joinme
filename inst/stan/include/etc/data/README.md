# Data-stage helpers

Most declarations are owned by scientific submodels. This directory contains
only data shared across several stories: `regression_prior_data.stan` defines
the common coefficient-prior layout and `dynamic_prediction_draw_data.stan`
defines the retained-draw axis. Latent-class data are documented under
`../../submodels/latent_class/data/`. Master entry points include these files
explicitly.

`fitted_random_effect_draw_data.stan` declares realised posterior effects for
a subject already represented in the fit. `fixed_prediction_latent_placeholders.stan`
provides dimensionally neutral arrays solely to let the parameter-free fitted-
effect entry point reuse the common prediction equations. These placeholders
do not define a probability distribution and never enter a likelihood.
