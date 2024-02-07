# Import libraries
library("dplyr")
library("ggplot2")
library("sf")

# Load RspatialCMC
library("RspatialCMC")

# Set current directory
setwd(sprintf("%s/ArealData_Scenario5",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

###########################################################################
# Auxiliary functions -----------------------------------------------------

# Extract cluster_allocation matrix from the MCMC chain --> MIGLIORARE E INTEGRARE NEL PACCHETTO
get_cluster_allocs <- function(chain) {
  t(sapply(chain, function(state){state$cluster_allocs}))
}

# Extract unique_values list from the MCMC chain --> MIGLIORARE E UNTEGRARE NEL PACCHETTO
get_unique_values <- function(chain) {
  extract_unique_values <- function(cluster_state) {
    sapply(cluster_state, function(x){c(x$uni_ls_state$mean, x$uni_ls_state$var)})
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

# Get code for Lombardia region from region shapefile
lombardia <- read_sf("input/shp/IT-regioni_2020/IT-regioni_2020.shp") %>%
  filter(DEN_REG == "Lombardia") %>% select(COD_REG) %>%
  st_drop_geometry() %>% as.numeric()

# Associate code and name of provinces from province shapefile
sf_prov <- read_sf("input/shp/IT-province_2020/IT-province_2020.shp") %>%
  filter(COD_REG == lombardia)

# Map code to province names
sigle_prov <- sf_prov %>%
  select(COD_PROV, SIGLA) %>%
  st_drop_geometry()

# Generate sf_mun dataset
sf_mun <- read_sf("input/shp/IT-comuni_2020/IT-comuni_2020.shp") %>%
  filter(COD_REG == lombardia) %>%
  merge(sigle_prov, by = "COD_PROV") %>%
  mutate(PROVINCIA=as.factor(SIGLA)) %>%
  select(COMUNE, PROVINCIA)

# Select means and standard deviations for Gaussian kernels
gammas <- c(1, 5, 10, 20)

# Generate cluster allocs
bergamo_coords <- sf_mun %>% filter(COMUNE == "Bergamo") %>%
  select(geometry) %>% st_centroid() %>% st_coordinates()
all_coords <- sf_mun %>% select(geometry) %>% st_centroid() %>% st_coordinates()
assign2clust <- function(coord){
  if(coord[1] >= bergamo_coords[1] & coord[2] >= bergamo_coords[2]){
    return(1L)
  }
  if(coord[1] >= bergamo_coords[1] & coord[2] < bergamo_coords[2]){
    return(2L)
  }
  if(coord[1] < bergamo_coords[1] & coord[2] >= bergamo_coords[2]){
    return(3L)
  }
  if(coord[1] < bergamo_coords[1] & coord[2] < bergamo_coords[2]){
    return(4L)
  }
}
sf_mun$GROUP <- as.factor(apply(all_coords, 1, assign2clust))

# Generate data
set.seed(1996)
sample_data <- function(clust_allocs){
  return(rpois(1, gammas[clust_allocs]))
}
sf_mun$DATA <- sapply(sf_mun$GROUP, sample_data)

###########################################################################

###########################################################################
# spatialCMC run ----------------------------------------------------------

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
fit <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), as.numeric(sf_mun$PROVINCIA)-1,
               "PoissonGamma", hier_params, mix_params, algo_params)

fit_mcmc <- run_mcmc(sf_mun$DATA, st_geometry(sf_mun),
                     "PoissonGamma", hier_params, mix_params, algo_params)

###########################################################################

###########################################################################
# Posterior inference -----------------------------------------------------

# Deserialize chain
chain <- sapply(fit_mcmc, function(x){read(bayesmix.AlgorithmState,x)})

# Get quantity of interest
cluster_allocs <- get_cluster_allocs(chain)
unique_values <- get_unique_values(chain)

# Compute findings from the approximated posterior distribution
Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
psm <- salso::psm(cluster_allocs)
sf_mun$best_clust <- as.factor(salso::salso(cluster_allocs, loss = "VI"))

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Plot data on the map
plt_data <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=DATA)) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(direction = "horizontal", title.hjust = 0.5,
                                             title.position = "bottom", barwidth = unit(3,"in"))) +
  geom_sf(data = sf_prov, col="darkred", linewidth=0.8, fill=NA) +
  theme_void() + theme(legend.position = "bottom")
plt_data

# Plot true cluster allocation on the map
plt_group <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=GROUP)) +
  guides(fill=guide_legend(direction = "horizontal", title.hjust = 0.5,
                           title.position = "bottom", label.position = "bottom")) +
  geom_sf(data = sf_prov, col="darkred", linewidth=0.8, fill=NA) +
  theme_void() + theme(legend.position = "bottom")
plt_group

# Plot - Posterior number of clusters
plt_nclust <- ggplot(data = data.frame(prop.table(table(Nclust))), aes(x=Nclust,y=Freq)) +
  geom_bar(stat = "identity", color=NA, linewidth=0, fill='white') +
  geom_bar(stat = "identity", color='steelblue', alpha=0.4, linewidth=0.7, fill='steelblue') +
  xlab("N° of Clusters") + ylab("Post. Prob.")
plt_nclust

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
  geom_sf(data = sf_mun, aes(fill=GROUP), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = sf_prov, color='darkred', fill=NA, linewidth=1.2) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")
plt_true_clust

# Plot - Best cluster on the geometry
plt_best_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=best_clust), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = sf_prov, color='darkred', fill=NA, linewidth=1.2) +
  guides(fill = guide_legend(title = "Cluster", title.position = "bottom", title.hjust=0.5,
                             label.position = "bottom", keywidth = unit(1,"cm"))) +
  theme_void() + theme(legend.position = "none")
plt_best_clust

# Plot - absolute difference between true and predictev values
y_pred <- compute_y_pred(chain)
sf_mun$std_diff <- abs(colMeans(y_pred) - sf_mun$DATA) / apply(y_pred, 2, sd)
plt_diff <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=std_diff), color='gray25', linewidth=0.5, alpha=0.75) +
  geom_sf(data = sf_prov, color='darkred', fill=NA, linewidth=1.2) +
  scale_fill_gradient(low="gray80", high = "orange") +
  guides(fill = guide_colorbar(title = "Diff.", title.position = "bottom",
                               title.hjust=0.5, barwidth = unit(3,"in"))) +
  theme_void() + theme(legend.position = "bottom")
plt_diff

###########################################################################
