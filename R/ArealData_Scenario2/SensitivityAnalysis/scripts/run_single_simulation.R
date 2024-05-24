# Import libraries
library("ggplot2")
library("sf")

# Import RspatialCMC
library("RspatialCMC")

# Set current directory
setwd(sprintf("%s/ArealData_Scenario2",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Build spatialCMC if needed
# build_spatialCMC()

# Import data
load("input/clean_data.dat")

# Set hierarchy parameters
hier_params =
  "
  fixed_values {
    mean: 0
    var_scaling: 0.1
    shape: 2
    scale: 2
  }
  "

# Set mixing parameters
mix_params =
  "
  fixed_value {
    totalmass: 2.0
    lambda: 0.1
  }
  "

# Set algorithm parameters
algo_params =
  "
  algo_id: 'Neal3'
  rng_seed: 10092022
  iterations: 4000
  burnin: 1000
  init_num_clusters: 5
  "

# Create simulations folder if needed
# out_dirname <- "simulations"
# if(!dir.exists(sprintf("output/%s", out_dirname))) {
#   dir.create(sprintf("output/%s", out_dirname))
# }

# Generate commont timestamp
# timestamp <- sprintf("%s", format(Sys.time(), "%Y%m%d-%H%M"))

# Run CMC algorithm
fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), as.numeric(sf_mun$PROVINCIA)-1,
                   "NNIG", hier_params, "sPP", mix_params, algo_params)

timestamp <- sprintf("%s", format(Sys.time(), "%Y%m%d-%H%M"))
save(list=c("fit_cmc", "hier_params", "mix_params", "algo_params"),
     file=sprintf("output/fit_cmc-%s.dat", timestamp))


# Save CMC output
# save(list=c("fit_cmc", "hier_params", "mix_params", "algo_params"),
#      file=sprintf("output/%s/fit_cmc-%s.dat", out_dirname, timestamp))

# Run MCMC algorithm
# fit_mcmc <- run_mcmc(sf_mun$DATA, st_geometry(sf_mun),
#                      "NNIG", hier_params, mix_params, algo_params)

# Save MCMC output
# save(list=c("fit_mcmc", "hier_params", "mix_params", "algo_params"),
#      file=sprintf("output/%s/fit_mcmc-%s.dat", out_dirname, timestamp))
