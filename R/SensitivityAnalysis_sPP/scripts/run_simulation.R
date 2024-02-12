# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_simulation", hide.opts = TRUE,
                         description = "Run RspatialCMC::sample_eppf with specified input parameters")
opt_parser <- add_argument(opt_parser, arg = "geometry", default = NULL,
                           help = "Path to .dat file where the geometry associated to the data is stored")
opt_parser <- add_argument(opt_parser, arg = "--alpha", short = "-a", type = "double", default = 1.0,
                           help = "Value of 'alpha' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--lambda", short = "-l", type = "double", default = 0.35,
                           help = "Value for the 'lambda' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--output-file", type = "character", default = "./output/chain.dat",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)


# Run simulation ----------------------------------------------------------

library("RspatialCMC")

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(dirname(sub("--file=", "", args[grep("--file=", args)])))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(basedir)
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check geometry file
geometry_file <- file.path(getwd(), extra_args$geometry)
if(!file.exists(geometry_file)){
  stop(sprintf("%s does not exist", geometry_file))
}
geometry_file <- normalizePath(geometry_file)
cat(sprintf("Geometry file: %s\n", geometry_file)) # Log

# Load geometry file
load(geometry_file)

# Deduce number of data
n <- length(grid_geometry)
cat(sprintf("N° of data: %g\n", n)) # Log

# Define mixing parameters template
mix_params_template =
"
fixed_value {
  totalmass: %g
  lambda: %g
}
"

# Get mixing parameter
totalmass <- extra_args$alpha
lambda <- extra_args$lambda

# Set mixing parameters
mix_params <- sprintf(mix_params_template, totalmass, lambda)
cat(sprintf("Mixing parameters:%s", mix_params))

# Set algorithm parameters
algo_params =
"
algo_id: 'Neal2'
rng_seed: 10092022
iterations: 4000
burnin: 1000
init_num_clusters: 5
"

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log

# Run Sampler
fit <- sample_eppf(n, geometry = grid_geometry,
                   mix_type = "sPP", mix_params = mix_params,
                   algo_params = algo_params)

# Save to output file
if (exists("fit")) {
  save(fit, file = out_file)
}