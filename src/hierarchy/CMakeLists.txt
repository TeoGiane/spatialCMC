# Add sources to bnp-downscaling
target_sources(spatialCMC PUBLIC
  poisson_likelihood.h
  poisson_likelihood.cc
  gamma_prior_model.h
  gamma_prior_model.cc
  poisson_gamma_updater.h
  poisson_gamma_updater.cc
  poisson_gamma_hierarchy.h
  poisson_gamma_hierarchy.cc
)
