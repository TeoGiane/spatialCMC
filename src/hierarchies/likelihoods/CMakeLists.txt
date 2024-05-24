# Add sources files
target_sources(spatialCMC PUBLIC
  empty_likelihood.h
  empty_likelihood.cc
  poisson_likelihood.h
  poisson_likelihood.cc
  poisson_regression_likelihood.h
  poisson_regression_likelihood.cc
)

# Add sources directories
add_subdirectory(states)