# # ---- Third test on SBNP Clustering for Massive datasets ---- # #

# Required libraries
library("RspatialCMC")
library("ggplot2")
library("parallel")
library("sf")

# Build spatialCMC (if necessary)
build_spatialCMC()

###########################################################################
# Auxiliary functions -----------------------------------------------------

# # Generate the true cluster
generate_clust_allocs <- function(Nmun) {
  clust_allocs <- matrix(1, sqrt(Nmun), sqrt(Nmun))
  for (i in 1:((sqrt(Nmun)/2)-2)) {
    clust_allocs[1:i,i] <- 2
    clust_allocs[1:i,(sqrt(Nmun)-i+1)] <- 2
    clust_allocs[rev((sqrt(Nmun)-i+1):(sqrt(Nmun))),i] <- 3
    clust_allocs[rev((sqrt(Nmun)-i+1):(sqrt(Nmun))),(sqrt(Nmun)-i+1)] <- 3
  }
  clust_allocs[1:(sqrt(Nmun)/2-2),(sqrt(Nmun)/2-2):(sqrt(Nmun)/2+2)] <- 2
  clust_allocs[1:(sqrt(Nmun)/5),(sqrt(Nmun)/2-2):(sqrt(Nmun)/2+2)] <- 3
  clust_allocs[rev((sqrt(Nmun)/2+3):sqrt(Nmun)),(sqrt(Nmun)/2-2):(sqrt(Nmun)/2+2)] <- 3
  clust_allocs[rev((4*sqrt(Nmun)/5+1):sqrt(Nmun)),(sqrt(Nmun)/2-2):(sqrt(Nmun)/2+2)] <- 2
  clust_allocs <- as.numeric(clust_allocs)
  return(clust_allocs)
}

# # Extract cluster_allocation matrix from the MCMC chain
get_cluster_allocs <- function(chain) {
  t(sapply(chain, function(state){state$cluster_allocs}))
}

# # Extract unique_values list from the MCMC chain
get_unique_values <- function(chain) {
  extract_unique_values <- function(cluster_state) {
    sapply(cluster_state, function(x){
      unp_x <- read(spatialcmc.PoissonState, x$custom_state$value)
      return(unp_x$rate)
    })
  }
  lapply(chain, function(state){extract_unique_values(state$cluster_states)})
}

# Predict response in each area
compute_y_pred <- function(chain) {
  Niter <- length(chain)
  clus_allocs <- get_cluster_allocs(chain) + 1
  unique_values <- get_unique_values(chain)
  Ndata <- dim(clus_allocs)[2]
  out <- matrix(nrow = Niter, ncol = Ndata)
  for (i in 1:Niter) {
    out[i,] <- rpois(Ndata, unique_values[[i]][clus_allocs[i,]])
  }
  return(out)
}

###########################################################################

###########################################################################
# Generate Data and Set Parameters ----------------------------------------

# Unit square polygon
box <- rbind(c(0,0), c(0,1), c(1,1), c(1,0), c(0,0))
square <- st_polygon(list(box))

# Generate municipalities and province geometries
Nprov <- 3; Nmun <- 900
geom_mun <- st_make_grid(square, n = rep(sqrt(Nmun), 2), offset = c(0,0))
geom_prov <- st_make_grid(square, n = c(Nprov, 1), offset = c(0,0))
# geom_prov <- st_make_grid(square, n = c(1, Nprov), offset = c(0,0))

# Generate indicator for province assignment
prov_allocs <- rep(rep(c(0,1,2), each = sqrt(Nmun)/Nprov), 30)
# prov_allocs <- rep(c(0,1,2), each = 300)
province_idx <- lapply(unique(prov_allocs), function(x) which(prov_allocs == x))

# Generate cluster_allocs
Ndata <- length(geom_mun)
gamma <- c(10,40,80)
clust_allocs <- generate_clust_allocs(Nmun)

# Generate data in municipality
set.seed(1996)
data <- rpois(Ndata, gamma[clust_allocs])

# Generate sf object
df_mun <- data.frame("clus_allocs" = clust_allocs,
                     "province_idx" = prov_allocs,
                     "data" = data)
sf_mun <- st_sf(df_mun, geometry = geom_mun)

###########################################################################

###########################################################################
# SpatialCMC run ----------------------------------------------------------

# Set hierarchy parameters
hier_params =
  "
  fixed_values {
    shape: 3.25
    rate: 0.25
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

# Run SpatialCMC sampler
fit <- run_cmc(sf_mun$data, st_geometry(sf_mun), sf_mun$province_idx,
               "PoissonGamma", hier_params, "sPP", mix_params, algo_params)

# fit_mcmc <- run_mcmc(sf_mun$data, st_geometry(sf_mun),
#                      "PoissonGamma", hier_params, mix_params, algo_params)

###########################################################################

###########################################################################
# Posterior inference -----------------------------------------------------

# Deserialize chain
chain <- sapply(fit, function(x){read(bayesmix.AlgorithmState,x)})

# Get quantity of interest
cluster_allocs <- get_cluster_allocs(chain)
unique_values <- get_unique_values(chain)

# Compute findings from the approximated posterior distribution
Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
psm <- salso::psm(cluster_allocs)
sf_mun$best_clust <- as.factor(salso::salso(cluster_allocs, loss = "VI"))
sf_mun$true_clust <- as.factor(clust_allocs)

# Plot - Posterior number of clusters
plt_nclust <- ggplot(data = data.frame(prop.table(table(Nclust))), aes(x=Nclust,y=Freq)) +
  geom_bar(stat = "identity", color=NA, linewidth=0, fill='white') +
  geom_bar(stat = "identity", color='steelblue', alpha=0.4, linewidth=0.7, fill='steelblue') +
  xlab("N° of Clusters") + ylab("Post. Prob.")

# Plot - Posterior similarity matrix
plt_psm <- ggplot(data = reshape2::melt(psm, c("x", "y"))) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient("Post. Prob. of Co-Clustering", low='steelblue', high='darkorange',
                      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5,
                                             direction = "horizontal", barwidth = unit(3,"in"))) +
  geom_rect(xmin=0.5, ymin=0.5, xmax=Ndata+0.5, ymax=Ndata+0.5, fill=NA, color='gray25', linewidth=0.7) +
  theme_void() + theme(legend.position = "bottom") + coord_equal()

# Plot - True cluster on the geometry
plt_true_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=true_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  # geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  scale_fill_manual(values = c("1" = "steelblue", "2" = "darkorange")) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")

# Plot - Best cluster on the geometry
plt_best_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=best_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  scale_fill_manual(values = c("1" = "darkorange", "2" = "steelblue")) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")

# Show posterior findings
gridExtra::grid.arrange(grobs = list(plt_true_clust, plt_best_clust), ncol=2)

# Plot - absolute difference between true and predictev values
y_pred <- compute_y_pred(chain)
# tmp <- apply(y_pred, 2, median)
sf_mun$std_diff <- abs((colMeans(y_pred) - sf_mun$data)) / apply(y_pred, 2, sd)
plt_diff <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=std_diff), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  scale_fill_gradient(low="gray80", high = "orange") +
  guides(fill = guide_colorbar(title = "Diff.", title.position = "bottom",
                               title.hjust=0.5, barwidth = unit(3,"in"))) +
  theme_void() + theme(legend.position = "bottom")
plt_diff

###########################################################################
