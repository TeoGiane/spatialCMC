# Add sources directories
add_subdirectory(likelihoods)
add_subdirectory(priors)
add_subdirectory(updaters)

# Add sources files
target_sources(spatialCMC PUBLIC
  poisson_gamma_hierarchy.h
)
