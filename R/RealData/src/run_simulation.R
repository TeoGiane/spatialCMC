# Required packages
library("RspatialCMC")
library("sf")

# Load Dataset (remove the NAs)
full_data <- read_sf("input/full_dataset.shp")
full_data <- full_data[which(full_data$ET != 0),]
full_data <- full_data[full_data$NAME_REG == "Lombardia",]

# Choose if MCMC or CMC
use_cmc <- TRUE

# Set hierarchy parameters
hier_prior =
  "
  fixed_values {
    shape: 250
    rate: 50
  }
  "

# Set mixing parameters
mix_prior =
  "
  fixed_value {
    totalmass: 1.0
    lambda: 0.35
  }
  "

# Set algorithm parameters
algo_params =
  "
  algo_id: 'Neal2'
  rng_seed: 10092022
  iterations: 5000
  burnin: 1000
  init_num_clusters: 5
  "

# Run SpatialCMC sampler (either MCMC or CMC)
if(use_cmc) {
  shard_alloc <- rep(0,nrow(full_data))
  # shard_alloc <- full_data$COD_REG-1
  fit <- run_cmc(data = as.numeric(full_data$T_20), geometry = st_geometry(full_data), shard_allocation = shard_alloc,
                 algo_params = algo_params, "PoissonGamma", hier_prior,"sPP", mix_prior, covariates = as.matrix(full_data$ET))
} else {
  fit <- run_mcmc(data = as.numeric(full_data$T_20), geometry = st_geometry(full_data),
                  algo_params = algo_params, "PoissonGamma", hier_prior,"sPP", mix_prior, covariates = as.matrix(full_data$ET))
}

# Save chain to output
save(fit, file = "chain.dat")
