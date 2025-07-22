# # optparse library to select the folder of the simulation study
# library("optparse")
# 
# # Option Parser
# option_list <- list(
#   make_option(c("--sim_study"), type="character", default=NULL,
#               help="name of the simulation study. Is a folder inside the 'output/' directory",
#               metavar="character")
# )
# opt_parser <- OptionParser(option_list=option_list)
# args <- parse_args(opt_parser)
# 
# # Import libraries
# library("ggplot2")
# library("sf")
# 
# # Import RspatialCMC
# library("RspatialCMC")
# 
# # Load protobuf messages
# RspatialCMC::import_protobuf_messages()

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

# # Set current directory
# setwd(sprintf("%s/ArealData_Scenario2",
#               dirname(Sys.getenv("RSPATIALCMC_HOME"))))

# Import data
# load("input/clean_data.dat")

# pdf(file = "lombardy_geometry.pdf", height = 4, width = 4)
# ggplot() + geom_sf(data = sf_mun, color='gray25', fill=NA, linewidth=0.75) + theme_void()
# dev.off()

# Select simulation study folder
# sim_study <- sprintf("output/%s", args$sim_study)

# Generate plots for every file in folder
# list_files <- list.files(sim_study, pattern = "*.dat", full.names = T)
# for (curr_file in list_files) {
#   
#   # Load output chain and parameter
#   load(curr_file)

  # Deserialize chain
  Ndata <- nrow(sf_mun)
  chain <- sapply(fit_cmc, function(x){read(bayesmix.AlgorithmState,x)})
  
  # Get quantity of interest
  cluster_allocs <- get_cluster_allocs(chain)
  unique_values <- get_unique_values(chain)
  
  # Compute findings from the approximated posterior distribution
  Nclust <- apply(cluster_allocs, 1, function(x){length(unique(x))})
  # psm <- salso::psm(cluster_allocs)
  sf_mun$best_clust_VI <- as.factor(salso::salso(cluster_allocs, loss = "VI"))
  # sf_mun$best_clust_binder <- as.factor(salso::salso(cluster_allocs, loss = "binder"))
  
  # Plot - Posterior number of clusters
  plt_nclust <- ggplot(data = data.frame(prop.table(table(Nclust))), aes(x=Nclust,y=Freq)) +
    geom_bar(stat = "identity", color=NA, linewidth=0, fill='white') +
    geom_bar(stat = "identity", color='steelblue', alpha=0.4, linewidth=0.7, fill='steelblue') +
    xlab("N° of Clusters") + ylab("Post. Prob.")
  plt_nclust
  
  # Plot - Best cluster VI on the geometry
  # levels(sf_mun$best_clust_VI) <- list("1"="4", "2"="3", "3"="1", "4"="2")
  plt_best_clust_VI <- ggplot() +
    geom_sf(data = sf_mun, aes(fill=best_clust_VI), color='gray25', show.legend = F) +
    geom_sf(data = sf_prov, col='darkred', fill=NA, linewidth=1) +
    theme_void()
  plt_best_clust_VI

  # Plot - Best cluster Binder on the geometry
  # plt_best_clust_binder <- ggplot() +
  #   geom_sf(data = sf_mun, aes(fill=best_clust_binder), color='gray25', show.legend = F) +
  #   geom_sf(data = sf_prov, col='darkred', fill=NA, linewidth=1) +
  #   theme_void()
  # plt_best_clust_binder
  
  # Plot - True cluster on the geometry
  plt_true_clust <- ggplot() +
    geom_sf(data = sf_mun, aes(fill=GROUP), color='gray25', show.legend = F) +
    geom_sf(data = sf_prov, col='darkred', fill=NA, linewidth=1) +
    theme_void()
  plt_true_clust
  
  titletext <- grid::textGrob(bquote(alpha~"="~.(alpha)~","~lambda~"="~.(lambda)),
                              gp=grid::gpar(fontsize=16))
  pdf(file = sprintf("output/plt_cmc-%s.pdf", timestamp), height = 4, width = 7)
  print(gridExtra::grid.arrange(plt_nclust, plt_best_clust_VI, ncol=2, top = titletext))
  dev.off()
  
  # # Save plot in proper folder
  # pdf(sprintf("%s-plt_nclust.pdf", fs::path_ext_remove(curr_file)), height = 3, width = 3)
  # print(plt_nclust)
  # dev.off()
  # pdf(sprintf("%s-best_clust_VI.pdf", fs::path_ext_remove(curr_file)), height = 3, width = 3)
  # print(plt_best_clust_VI)
  # dev.off()
  # pdf(sprintf("%s-best_clust_binder.pdf", fs::path_ext_remove(curr_file)), height = 3, width = 3)
  # print(plt_best_clust_binder)
  # dev.off()

# }
