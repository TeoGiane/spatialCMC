library(ggplot2)
library("RspatialCMC")
library("sf")

show_plot = TRUE

full_data <- read_sf("input/full_dataset.shp")
full_data <- full_data[which(full_data$ET != 0),]

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

load("chain.dat")

# Deserialize chain
chain <- sapply(fit, function(x){read(bayesmix.AlgorithmState,x)})

# Get quantity of interest
cluster_allocs <- get_cluster_allocs(chain)
unique_values <- get_unique_values(chain)

# Compute findings from the approximated posterior distribution
Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
full_data$BEST_CLUST <- as.factor(salso::salso(cluster_allocs, loss = "VI"))

# Plot - Posterior number of clusters
plt_nclust <- ggplot(data = data.frame(prop.table(table(Nclust))), aes(x=Nclust,y=Freq)) +
  geom_bar(stat = "identity", color=NA, linewidth=0, fill='white') +
  geom_bar(stat = "identity", color='steelblue', alpha=0.4, linewidth=0.7, fill='steelblue') +
  xlab("N° of Clusters") + ylab("Post. Prob.")
if(show_plot){
  x11(height = 4, width = 4); plt_nclust
} else {
  pdf("plt_nclust.pdf", height = 4, width = 4); print(plt_nclust); dev.off()
}

# Plot - Best cluster on the map
ggplot() +
  geom_sf(data = full_data, aes(fill=BEST_CLUST))
