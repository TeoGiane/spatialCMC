# Import libraries
library("dplyr")
library("ggplot2")
library("sf")

# Load RspatialCMC
library("RspatialCMC")

# Set current directory
setwd(sprintf("%s/ArealData_Scenario2",
              dirname(Sys.getenv("RSPATIALCMC_HOME"))))

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
means <- c(-5,-2,2,5); std_devs <- c(1,1,1,1)

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
sample_data <- function(clust_allocs){
  return(rnorm(1, means[clust_allocs], std_devs[clust_allocs]))
}
sf_mun$DATA <- sapply(sf_mun$GROUP, sample_data)

# Save data
save(list = c("sf_mun", "sf_prov"), file = "input/clean_data.dat")

# Plot data on the map
plt_data <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=DATA)) +
  scale_fill_gradient(low = "steelblue", high = "darkorange",
                      guide = guide_colorbar(direction = "horizontal", title.hjust = 0.5,
                                             title.position = "bottom", barwidth = unit(3,"in"))) +
  geom_sf(data = sf_prov, col="darkred", linewidth=0.8, fill=NA) +
  theme_void() + theme(legend.position = "bottom")

# Plot true cluster allocation on the map
plt_group <- ggplot() +
  geom_sf(data = sf_mun, aes(fill=GROUP)) +
  guides(fill=guide_legend(direction = "horizontal", title.hjust = 0.5,
                           title.position = "bottom", label.position = "bottom")) +
  geom_sf(data = sf_prov, col="darkred", linewidth=0.8, fill=NA) +
  theme_void() + theme(legend.position = "bottom")

# Save plots in pdf
pdf("output/plt_Scenario2.pdf", height = 4, width = 8)
gridExtra::grid.arrange(plt_data, plt_group, ncol=2)
dev.off()
