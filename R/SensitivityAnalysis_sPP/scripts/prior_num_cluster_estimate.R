# Load libraries
library("sf")
library("RspatialCMC")
library("ggplot2")

# Set working directory
setwd(sprintf("%s/SensitivityAnalysis_sPP",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Create (common) Geometry
box <- rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0))
square <- st_polygon(list(box))
# Common parameters
level = 0.05; Npoints <- 5; alpha <- 1; lambda <- 1e-3
# Prepare final buffer
out <- data.frame("N" = c(1,numeric(Npoints-1)),
                  "Mean" = c(1,numeric(Npoints-1)),
                  "q_0.025" = c(1,numeric(Npoints-1)),
                  "q_0.5" = c(1,numeric(Npoints-1)),
                  "q_0.975" = c(1,numeric(Npoints-1)))

# Simulation
for (geom_dim in 2:Npoints) {
  # Make grid from polygon
  geometry <- st_make_grid(square, n = rep(geom_dim, 2), offset = c(0,0))
  n <- length(geometry)
  # Define mixing parameters template
  mix_params =
    "
    fixed_value {
      totalmass: %g
      lambda: %g
    }
    "
  mix_params <- sprintf(mix_params, alpha, lambda)
  # Set algorithm parameters
  algo_params =
    "
    algo_id: 'Neal2'
    rng_seed: 10092022
    iterations: 4000
    burnin: 1000
    init_num_clusters: %g
    "
  algo_params <- sprintf(algo_params, min(geom_dim, 5))
  # Run Sampler
  fit <- sample_eppf(n, geometry = geometry, mix_type = "sPP",
                     mix_params = mix_params, algo_params = algo_params)
  # Get posterior mean and CI 
  chain <- sapply(fit, function(x){read(bayesmix.AlgorithmState,x)})
  cluster_allocs <- t(sapply(chain, function(x){x$cluster_allocs}))
  Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
  # Fill buffer
  out[geom_dim,] <- c(n, mean(Nclust), quantile(Nclust, c(level/2, 0.5, 1-level/2)))
}

# Store .csv file
filename <- sprintf("NumClust_alpha-%g_lambda-%g.csv", alpha, lambda)
write.csv(out, file = filename, row.names = FALSE)

# Generate plot
# plt <- ggplot(data=out) +
#   geom_ribbon(aes(x=N, ymin=q_0.025, ymax=q_0.975), color='blue', alpha=0.4, linewidth=0) +
#   geom_line(aes(x=N, y=Mean)) + geom_line(aes(x=N, y=q_0.5), color='blue') +
#   ylab("Num. Clust.") + ggtitle(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)))
# # Save plot
# plt
# pdf(file = "prior_num_clust_")
