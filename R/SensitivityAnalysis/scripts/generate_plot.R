# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_plot", hide.opts = TRUE,
                         description = "Plot posterior number of clusters and bets cluster estimate with specified input parameters")
opt_parser <- add_argument(opt_parser, arg = "input-file", default = NULL,
                           help = "Path to .dat file where the input data, is sf format, are stored")
opt_parser <- add_argument(opt_parser, arg = "--chain-file", short = "-c",
                           help = "Relative path to the .dat file where the MCMC chain is stored")
opt_parser <- add_argument(opt_parser, arg = "--output-file", short = "-o",
                           help = "Relative path to the output pdf plot file")
extra_args <- parse_args(opt_parser)


# Auxiliary Functions -----------------------------------------------------

# Extract cluster_allocation matrix from the MCMC chain
get_cluster_allocs <- function(chain) {
  t(sapply(chain, function(state){state$cluster_allocs}))
}

# Extract unique_values list from the MCMC chain
get_unique_values <- function(chain) {
  extract_unique_values <- function(cluster_state) {
    sapply(cluster_state, function(x){
      unp_x <- read(spatialcmc.PoissonState, x$custom_state$value)
      return(unp_x$rate)
    })
  }
  lapply(chain, function(state){extract_unique_values(state$cluster_states)})
}


# Generate plot -----------------------------------------------------------

# Required libraries
library("RspatialCMC")
library("ggplot2")
library("sf")

# Import protobuf messages
import_protobuf_messages()

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

# Check chain file
chain_file <- normalizePath(extra_args$chain_file, mustWork = TRUE)
cat(sprintf("Chain file: %s\n", chain_file)) # Log

# Load chain file
load(chain_file)

# Get params from chain filename
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

# Deserialize chain
chain <- sapply(fit, function(x){read(bayesmix.AlgorithmState,x)})

# Get quantity of interest
cluster_allocs <- get_cluster_allocs(chain)
unique_values <- get_unique_values(chain)

# Compute findings from the approximated posterior distribution
Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
sf_mun$best_clust <- as.factor(salso::salso(cluster_allocs, loss = "VI"))

# Plot - Posterior number of clusters
plt_nclust <- ggplot(data = data.frame(prop.table(table(Nclust))), aes(x=Nclust,y=Freq)) +
  geom_bar(stat = "identity", color=NA, linewidth=0, fill='white') +
  geom_bar(stat = "identity", color='steelblue', alpha=0.4, linewidth=0.7, fill='steelblue') +
  xlab("N° of Clusters") + ylab("Post. Prob.")

# Plot - Best cluster on the geometry
plt_best_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=best_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  guides(fill = guide_legend(title = element_blank(), title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom")) +
  theme_void() + theme(legend.position = "bottom")

# Show posterior findings
titletext <- grid::textGrob(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)),
                            gp=grid::gpar(fontsize=16))
pdf(file = out_file, height = 4, width = 7)
gridExtra::grid.arrange(plt_nclust, plt_best_clust, ncol=2, top = titletext)
dev.off()
