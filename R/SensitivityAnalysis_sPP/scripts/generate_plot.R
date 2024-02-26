# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_plot", hide.opts = TRUE,
                         description = "Plot prior number of clusters with specified input parameters")
opt_parser <- add_argument(opt_parser, arg = "chain-file",
                           help = "Relative path to the .dat file where the MCMC chain is stored")
opt_parser <- add_argument(opt_parser, arg = "--output-file", short = "-o",
                           help = "Relative path to the output pdf plot file")
extra_args <- parse_args(opt_parser)


# Generate barplot --------------------------------------------------------

library("RspatialCMC")

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(dirname(sub("--file=", "", args[grep("--file=", args)])))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(basedir)
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check geometry file
chain_file <- file.path(getwd(), extra_args$chain_file)
if(!file.exists(chain_file)){
  stop(sprintf("%s does not exist", chain_file))
}
chain_file <- normalizePath(chain_file)
cat(sprintf("Chain file: %s\n", chain_file)) # Log

# Load chain file
load(chain_file)

# Get params from chain_file
alpha <- regmatches(chain_file, gregexpr('alpha_\\d*\\.*\\d*', chain_file))
alpha <- as.numeric(gsub("alpha_", "", alpha))
lambda <- regmatches(chain_file, gregexpr('lambda_\\d*\\.*\\d*', chain_file))
lambda <- as.numeric(gsub("lambda_", "", lambda))
cat(sprintf("Parameters: alpha = %g, lambda = %g\n", alpha, lambda)) # Log

# Create directory for output if does not exist
out_file <- file.path(getwd(), extra_args$output_file)
if(!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE)
}
cat(sprintf("Output directory: %s\n", normalizePath(dirname(out_file)))) # Log

# Compute N_clust chain
import_protobuf_messages()
chains <- sapply(fit, function(x){read(bayesmix.AlgorithmState, x)})
clust_allocs <- t(sapply(chains, function(x){x$cluster_allocs}))
N_clust <- apply(clust_allocs, 1, function(x){length(unique(x))})

# Plot
pdf(file = out_file, height = 4, width = 4)
barplot(table(N_clust)/length(N_clust), ylim = c(0,1))
title(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)))
dev.off()
