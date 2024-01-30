# # ---- First test on SBNP Clustering for Massive datasets ---- # #

# Required libraries
library("RspatialCMC")
library("ggplot2")
library("parallel")
library("sf")

# Build spatialCMC (if necessary)
build_spatialCMC()

###########################################################################
# Auxiliary functions -----------------------------------------------------

# # Extract cluster_allocation matrix from the MCMC chain
get_cluster_allocs <- function(chain) {
  t(sapply(chain, function(state){state$cluster_allocs}))
}

# # Extract unique_values list from the MCMC chain
get_unique_values <- function(chain) {
  extract_unique_values <- function(cluster_state) {
    sapply(cluster_state, function(x){x$general_state$data})
  }
  lapply(chain, function(state){extract_unique_values(state$cluster_states)})
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

# Generate indicator for province assignment
prov_allocs <- rep(rep(c(0,1,2), each = sqrt(Nmun)/Nprov), 30)
province_idx <- lapply(unique(prov_allocs), function(x) which(prov_allocs == x))

# Generate data in municipality
Ndata <- length(geom_mun)
gamma <- c(10,20)
clust_allocs <- c(rep(c(rep(1,15),rep(2,15)),15),
                  rep(c(rep(2,15), rep(1,15)),15))
set.seed(1996)
data <- rpois(Ndata, gamma[clust_allocs])

# Generate sf object
df_mun <- data.frame("province_idx" = prov_allocs,
                     "data" = data)
sf_mun <- st_sf(df_mun, geometry = geom_mun)

# Generate common timestamp (for scenario and output matching)
timestamp <- format(Sys.time(), "%Y%m%d-%H%M")

# Save generated scenario scenario
filename <- sprintf("%s/input/scenario1_%s.dat", getwd(), timestamp)
save.image(file = filename)

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
               "PoissonGamma", hier_params, mix_params, algo_params)

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

# Plot - Best cluster on the geometry
plt_true_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=true_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  scale_fill_manual(values = c("1" = "steelblue", "2" = "darkorange")) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")

plt_best_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=best_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = geom_prov, color='darkred', fill=NA, linewidth=2) +
  scale_fill_manual(values = c("1" = "steelblue", "2" = "darkorange")) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")

# Show posterior findings
gridExtra::grid.arrange(grobs = list(plt_true_clust, plt_best_clust), ncol=2)
# plt_psm

###########################################################################
