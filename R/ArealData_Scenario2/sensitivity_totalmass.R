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
hier_params =
  "
  fixed_values {
    mean: 0
    var_scaling: 0.1
    shape: 3
    scale: 2
  }
  "

# Set mixing parameters
mix_params_template =
  "
  fixed_value {
    totalmass: %g
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
out_dir <- "output/sensitivity_totalmass"
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
}

# Select values for var_scaling parameter
total_masses <- c(0.1, 1.0)

# Run simulation study
for (curr_tm in total_masses) {
  # generate hier_params
  mix_params <- sprintf(mix_params_template, curr_tm)
  # Run CMC algorithm
  fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), as.numeric(sf_mun$PROVINCIA)-1,
                     "NNIG", hier_params, mix_params, algo_params)
  # Save CMC output
  save(list=c("fit_cmc", "hier_params", "mix_params", "algo_params"),
       file=sprintf("%s/fit_cmc-tm%g.dat", out_dir, curr_tm))
}
