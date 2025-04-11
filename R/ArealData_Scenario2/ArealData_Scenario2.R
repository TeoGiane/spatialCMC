# Import libraries
library("ggmap")
library("ggplot2")
library("RspatialCMC")
library("sf")

# # Auxiliary functions # #

# Fix for ggmap
sf_ggmap <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")

  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector,
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")),
                       c("ymin", "xmin", "ymax", "xmax"))

  # Coonvert the bbox to an sf polygon, transform it to 3857, and convert back to a bbox
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))

  # Overwrite the bbox of the ggmap object with the transformed coordinates
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]

  # Return
  return(map)
}

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
    out[i,] <- rnorm(Ndata, unique_values[[i]][1,clus_allocs[i,]],
                            sqrt(unique_values[[i]][2,clus_allocs[i,]]))
  }
  return(out)
}

# Set current directory
setwd(sprintf("%s/ArealData_Scenario2",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Build spatialCMC if needed
# build_spatialCMC()

# Import data
load("input/clean_data.dat")

# Import map from Stadia maps
lombardy_centre <- st_coordinates(st_transform(st_centroid(st_union(sf_prov)), 4326))
lombardy_bbox <- unname(c(lombardy_centre[,"X"]-1.5,
                          lombardy_centre[,"Y"]-1.1,
                          lombardy_centre[,"X"]+1.8,
                          lombardy_centre[,"Y"]+1.2))
lombardy_map <- sf_ggmap(get_map(lombardy_bbox, maptype = "stamen_terrain", source = "stadia", crop = T))

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

# Run CMC algorithm
fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), (as.numeric(sf_mun$PROVINCIA)-1),
                   "NNIG", hier_params, "sPP", mix_params, algo_params)
fit_mcmc <- run_mcmc(sf_mun$DATA, st_geometry(sf_mun),
                     "NNIG", hier_params, "sPP", mix_params, algo_params)

# fit_cmc <- run_cmc(sf_mun$DATA, st_geometry(sf_mun), sf_mun$province_idx,
#                    "NNIG", hier_params, "sPP", mix_params, algo_params)

# Set timestamp and save CMC output
timestamp <- sprintf("%s", format(Sys.time(), "%Y%m%d-%H%M"))
save(list=c("fit_mcmc", "hier_params", "mix_params", "algo_params"),
     file=sprintf("output/fit_cmc-%s.dat", timestamp))

# Deserialize chain
Ndata <- nrow(sf_mun)
chain <- sapply(fit_mcmc, function(x){read(bayesmix.AlgorithmState,x)})

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
sf_mun_3857 <- st_transform(sf_mun, 3857)
sf_prov_3857 <- st_transform(sf_prov, 3857)
plt_best_clust <- ggmap(lombardy_map) +
  geom_sf(data = sf_mun_3857, aes(fill=best_clust), alpha = 0.75, color='gray25', inherit.aes = F) +
  geom_sf(data = sf_prov_3857, col='darkred', fill=NA, linewidth=1, inherit.aes = F) +
  scale_fill_manual(values = c("1"='salmon',"2"='steelblue',"3"='lightgreen',"4"='darkorange')) +
  guides(fill = guide_legend("Cluster", position = "bottom", direction = "horizontal",
                             title.position = "bottom", label.position = "bottom", title.hjust = 0.5)) +
  theme_void()
plt_best_clust

# Plot - True cluster on the geometry
sf_mun_3857 <- st_transform(sf_mun, 3857)
sf_prov_3857 <- st_transform(sf_prov, 3857)
plt_true_clust <- ggmap(lombardy_map) +
  geom_sf(data = sf_mun_3857, aes(fill=GROUP), alpha = 0.75, color='gray25', inherit.aes = F) +
  geom_sf(data = sf_prov_3857, col='darkred', fill=NA, linewidth=1, inherit.aes = F) +
  scale_fill_manual(values = c("1"='steelblue',"2"='darkorange',"3"='salmon',"4"='lightgreen')) +
  guides(fill = guide_legend("Cluster", position = "bottom", direction = "horizontal",
                             title.position = "bottom", label.position = "bottom", title.hjust = 0.5)) +
  theme_void()
plt_true_clust

# Plot - Data on the geometry
sf_mun_3857 <- st_transform(sf_mun, 3857)
sf_prov_3857 <- st_transform(sf_prov, 3857)
plt_data <- ggmap(lombardy_map) +
  geom_sf(data = sf_mun_3857, aes(fill=DATA), alpha = 0.75, color='gray25', inherit.aes = F) +
  geom_sf(data = sf_prov_3857, col='darkred', fill=NA, linewidth=1, inherit.aes = F) +
  scale_fill_gradient(low='steelblue', high = 'darkorange') +
  guides(fill = guide_colorbar("Data", position = "bottom", direction = "horizontal",
                               title.position = "bottom", title.hjust = 0.5, barwidth = unit(3,"in"))) +
  theme_void()
plt_data

# Plot - absolute difference between true and predictev values
y_pred <- compute_y_pred(chain)
sf_mun$std_diff <- abs((colMeans(y_pred) - sf_mun$DATA)) / apply(y_pred, 2, sd)
sf_mun_3857 <- st_transform(sf_mun, 3857)
sf_prov_3857 <- st_transform(sf_prov, 3857)
plt_diff <- ggmap(lombardy_map) +
  geom_sf(data = sf_mun_3857, aes(fill=std_diff), alpha=0.75, color='gray25', inherit.aes = F) +
  geom_sf(data = sf_prov_3857, color='darkred', fill=NA, linewidth=1, inherit.aes = F) +
  scale_fill_gradient(low="gray80", high = "orange") +
  guides(fill = guide_colorbar(title = "Std. Diff.", position = "bottom", direction = "horizontal",
                               title.position = "bottom", title.hjust=0.5, barwidth = unit(3,"in"))) +
  theme_void()
plt_diff

# titletext <- grid::textGrob(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)),
#                             gp=grid::gpar(fontsize=16))
# pdf(file = sprintf("output/plt_cmc-%s.pdf", timestamp), height = 4, width = 7)
# titletext <- grid::textGrob("2log(BF) > 0", gp=grid::gpar(fontsize=16))
# gridExtra::grid.arrange(plt_nclust, plt_best_clust, ncol=2)
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
