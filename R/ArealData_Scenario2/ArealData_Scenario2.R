# Import libraries
library("ggplot2")
library("sf")
library("RspatialCMC")

# # Auxiliary functions # #
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

# Set current directory
setwd(sprintf("%s/ArealData_Scenario2",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Build spatialCMC if needed
# build_spatialCMC()

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
mix_params =
  "
  fixed_value {
    totalmass: 1.0
    lambda: 0.0
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

# Run CMC algorithm
fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), (as.numeric(sf_mun$PROVINCIA)-1),
                   "NNIG", hier_params, "sPP", mix_params, algo_params)
# fit_mcmc <- run_mcmc(sf_mun$DATA, st_geometry(sf_mun),
#                      "NNIG", hier_params, "sPP", mix_params, algo_params)

# fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), sf_mun$province_idx,
#                    "NNIG", hier_params, "sPP", mix_params, algo_params)

# Set timestamp and save CMC output
timestamp <- sprintf("%s", format(Sys.time(), "%Y%m%d-%H%M"))
save(list=c("fit_cmc", "hier_params", "mix_params", "algo_params"),
     file=sprintf("output/fit_cmc-%s.dat", timestamp))

# Deserialize chain
Ndata <- nrow(sf_mun)
chain <- sapply(fit_cmc, function(x){read(bayesmix.AlgorithmState,x)})

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
plt_nclust

# Plot - Best cluster VI on the geometry
plt_best_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=best_clust), color='gray25', show.legend = F) +
  geom_sf(data = sf_prov, col='darkred', fill=NA, linewidth=1) +
  theme_void()
plt_best_clust

# Plot - True cluster on the geometry
plt_true_clust <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=GROUP), color='gray25', show.legend = F) +
  geom_sf(data = sf_prov, col='darkred', fill=NA, linewidth=1) +
  theme_void()
plt_true_clust

# titletext <- grid::textGrob(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)),
#                             gp=grid::gpar(fontsize=16))
# pdf(file = sprintf("output/plt_cmc-%s.pdf", timestamp), height = 4, width = 7)
titletext <- grid::textGrob("2log(BF) > 0", gp=grid::gpar(fontsize=16))
gridExtra::grid.arrange(plt_nclust, plt_best_clust, ncol=2)
# dev.off()

# Esempio con block random shard partitioning --- NO
Nshards <- 6 #; Nrng <- 10

# set.seed(1996); idx <- sample(1:nrow(sf_mun), size = Nrng)

# points <- st_centroid(st_geometry(sf_mun[idx,]))
points <- st_centroid(st_geometry(sf_prov))
tessellation <- st_collection_extract(st_voronoi(st_union(points)), "POLYGON")
cropped_tessellation <- st_intersection(tessellation, st_union(sf_mun))
geom_microshard <- list()
geom_microshard <- lapply(1:length(cropped_tessellation),
                          function(i){
                            idx <- unlist(st_contains(cropped_tessellation[i], st_centroid(st_geometry(sf_mun))))
                            geom_microshard[[i]] <- st_union(st_geometry(sf_mun)[idx])
                          })
geom_microshard <- st_sfc(do.call(rbind, geom_microshard))

# Generate final shards
shard_idx <- sample(1:Nshards, size = length(geom_microshard), replace = T)
geom_shards <- lapply(1:Nshards, function(i){st_union(geom_microshard[shard_idx == i])})
geom_shards <- st_sfc(do.call(rbind, geom_shards))
st_crs(geom_shards) <- st_crs(sf_mun)

# Compute prov_allocs
get_prov_allocs <- function(geom_mun, geom_prov){
  out <- numeric(length(geom_mun))
  for (i in 1:length(geom_prov)) {
    idx <- unlist(st_contains(geom_prov[i], st_centroid(geom_mun)))
    out[idx] <- (i-1)
  }
  return(out)
}
prov_allocs <- get_prov_allocs(st_geometry(sf_mun), geom_shards)

# Compute province_idx
province_idx <- lapply(unique(prov_allocs), function(x) which(prov_allocs == x))
sf_mun$province_idx <- prov_allocs


# Riconduciamoci al caso griglia regolare per la partizione in shard -- NO
Nshards <- 5
regular_grid <- st_make_grid(sf_mun, n = c(10,10))
geom_microshard <- list()
geom_microshard <- lapply(1:length(regular_grid),
                          function(i){
                            idx <- unlist(st_contains(regular_grid[i], st_centroid(st_geometry(sf_mun))))
                            geom_microshard[[i]] <- st_union(st_geometry(sf_mun)[idx])
                          })
geom_microshard <- st_sfc(do.call(rbind, geom_microshard))

# Generate final shards
shard_idx <- sample(1:Nshards, size = length(geom_microshard), replace = T)
geom_shards <- lapply(1:Nshards, function(i){st_union(geom_microshard[shard_idx == i])})
geom_shards <- st_sfc(do.call(rbind, geom_shards))
st_crs(geom_shards) <- st_crs(sf_mun)

# Compute prov_allocs
get_prov_allocs <- function(geom_mun, geom_prov){
  out <- numeric(length(geom_mun))
  for (i in 1:length(geom_prov)) {
    idx <- unlist(st_contains(geom_prov[i], st_centroid(geom_mun)))
    out[idx] <- (i-1)
  }
  return(out)
}
prov_allocs <- get_prov_allocs(st_geometry(sf_mun), geom_shards)

# Compute province_idx
province_idx <- lapply(unique(prov_allocs), function(x) which(prov_allocs == x))
sf_mun$province_idx <- prov_allocs
