# Model-stage helpers

Submodel-specific prior contributions live under `../../submodels/`.
`regression_priors.stan` handles only the common coefficient-prior programme.
This directory also contains the joint threaded likelihood calls, which
necessarily operate across several scientific submodels. Latent-class density
replacement lives in `../../submodels/latent_class/model/`. Include order is stated
directly in each master file.
