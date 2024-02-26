# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_cmc", hide.opts = TRUE,
                         description = "Run RspatialCMC::run_cmc with specified input parameters")
opt_parser <- add_argument(opt_parser, arg = "input-file", default = NULL,
                           help = "Path to .dat file where the input data, is sf format, are stored")
opt_parser <- add_argument(opt_parser, arg = "--seed", short = "-s", default = 1,
                           help = "Seed for reproducibility")
opt_parser <- add_argument(opt_parser, arg = "--hier-params-file", short = "-h",
                           help = "Path to .asciipb file where the hierarchy parameters are stored")
opt_parser <- add_argument(opt_parser, arg = "--alpha", short = "-a", type = "double", default = 1.0,
                           help = "Value of 'alpha' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--lambda", short = "-l", type = "double", default = 0.35,
                           help = "Value for the 'lambda' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--output-file", type = "character", default = "./output/chain.dat",
                           help = "Relative path to the output file")
extra_args <- parse_args(opt_parser)


# Run Simulation ----------------------------------------------------------

# Required libraries
library("RspatialCMC")
library("ggplot2")
library("sf")

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(dirname(sub("--file=", "", args[grep("--file=", args)])))
setwd(normalizePath(basedir))
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check input file
input_file <- normalizePath(extra_args$input_file, mustWork = TRUE)
cat(sprintf("Input file: %s\n", input_file)) # Log

# Load input file
load(input_file)

# Define hier type
scenario <- as.character(regmatches(input_file, gregexpr('Scenario\\d*', input_file)))
if (scenario == "Scenario1" || scenario == "Scenario5") {
  hier_type <- "PoissonGamma"
} else if (scenario == "Scenario2" || scenario == "Scenario4") {
  hier_type <- "NNIG"
}
cat(sprintf("Hier Type: %s\n", hier_type)) # Log

# Check hier_params file
hier_params_file <- normalizePath(extra_args$hier_params_file, mustWork = TRUE)
cat(sprintf("Hier params file: %s\n", hier_params_file)) # Log

# Define mixing parameters template
mix_params_template =
"
fixed_value {
  totalmass: %g
  lambda: %g
}
"

# Set mixing parameters
mix_params <- sprintf(mix_params_template, extra_args$alpha, extra_args$lambda)
cat(sprintf("Mixing parameters:%s", mix_params))

# Set algorithm parameters
algo_params =
"
algo_id: 'Neal2'
rng_seed: 10092022
iterations: 50
burnin: 10
init_num_clusters: 5
"

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log

# Run SpatialCMC sampler
tcb <- addTaskCallback(function(...) {set.seed(extra_args$seed); TRUE})
fit <- run_cmc(sf_mun$data, st_geometry(sf_mun), sf_mun$province_idx,
               hier_type, hier_params_file, "sPP", mix_params, algo_params)
success <- removeTaskCallback(tcb)

# Save to output file
if (exists("fit")) {
  save(fit, file = out_file)
}
