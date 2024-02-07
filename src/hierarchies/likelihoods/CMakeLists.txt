# Add sources files
target_sources(spatialCMC PUBLIC
  poisson_likelihood.h
  poisson_likelihood.cc
)

# Add sources directories
add_subdirectory(states)