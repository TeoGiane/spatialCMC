# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "run_sample_eppf", hide.opts = TRUE,
                         description = "Run RspatialCMC::sample_eppf with specified number of data points and parameters")
opt_parser <- add_argument(opt_parser, arg = "--num-splits", short = "-n", type = "double", default = 30,
                           help = "Number of splits for each dimension of the unit squared geometry")
opt_parser <- add_argument(opt_parser, arg = "--alpha", short = "-a", type = "double", default = 1.0,
                           help = "Value of 'alpha' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--lambda", short = "-l", type = "double", default = 0.35,
                           help = "Value for the 'lambda' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character", default = "./output",
                           help = "Relative path to the output directory")
extra_args <- parse_args(opt_parser)


# Run simulation ----------------------------------------------------------

# Load libraries
library("sf")
library("RspatialCMC")

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(dirname(sub("--file=", "", args[grep("--file=", args)])))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(basedir)
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Generate the geometry from the number of splits
box <- rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0))
square <- st_polygon(list(box))
geometry <- st_make_grid(square, n = rep(extra_args$num_splits, 2), offset = c(0,0))
n <- length(geometry)

# Parse and define common parameters
level = 0.05; alpha = extra_args$alpha; lambda = extra_args$lambda

# Set algorithm parameters
algo_params =
"
algo_id: 'Neal2'
rng_seed: 10092022
iterations: 4000
burnin: 1000
init_num_clusters: %g
"
algo_params <- sprintf(algo_params, min(n, 5))

# Set mixing parameters
mix_params =
"
fixed_value {
  totalmass: %g
  lambda: %g
}
"
mix_params <- sprintf(mix_params, alpha, lambda)

# Create directory for output if does not exist
out_dir <- file.path(getwd(), extra_args$output_dir)
if(!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(out_dir))) # Log

# Run Sampler
fit <- sample_eppf(n, geometry = geometry,
                   mix_type = "sPP", mix_params = mix_params,
                   algo_params = algo_params)

# Retrieve Nclust chain
chain <- sapply(fit, function(x){read(bayesmix.AlgorithmState,x)})
cluster_allocs <- t(sapply(chain, function(x){x$cluster_allocs}))
Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})

# Store Nclust in .csv file
filename <- sprintf("NumClust_alpha-%g_lambda-%g-%g.csv", alpha, lambda, extra_args$num_splits)
write.table(Nclust, file = file.path(out_dir, filename), row.names = FALSE, col.names = FALSE, sep = ",")
