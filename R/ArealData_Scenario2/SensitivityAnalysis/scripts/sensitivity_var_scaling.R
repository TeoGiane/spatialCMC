# Import libraries
library("ggplot2")
library("sf")

# Import RspatialCMC
library("RspatialCMC")

# Set current directory
setwd(sprintf("%s/ArealData_Scenario2",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Build spatialCMC if needed
build_spatialCMC()

# Import data
load("input/clean_data.dat")

# Set hierarchy parameters
hier_params_template =
  "
  fixed_values {
    mean: 0
    var_scaling: %g
    shape: 3
    scale: 2
  }
  "

# Set mixing parameters
mix_params =
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
  iterations: 4000
  burnin: 1000
  init_num_clusters: 5
  "

# Create output directory if necessary
out_dir <- "output/sensitivity_var_scaling"
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Select values for var_scaling parameter
var_scalings <- seq(0.1, 1.0, by=0.1)

# Run simulation study
for (curr_vs in var_scalings) {
  # generate hier_params
  hier_params <- sprintf(hier_params_template, curr_vs)
  # Run CMC algorithm
  fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), as.numeric(sf_mun$PROVINCIA)-1,
                     "NNIG", hier_params, mix_params, algo_params)
  # Save CMC output
  save(list=c("fit_cmc", "hier_params", "mix_params", "algo_params"),
       file=sprintf("%s/fit_cmc-vs%g.dat", out_dir, curr_vs))
}
