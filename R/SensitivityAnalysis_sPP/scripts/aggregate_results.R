# Command line input options via argparser --------------------------------
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "aggregate_results", hide.opts = TRUE,
                         description = "Aggregate eppf simulations with same parameters")
opt_parser <- add_argument(opt_parser, arg = "input-dir",
                           help = "Relative path to the directory where .csv files are stored")
opt_parser <- add_argument(opt_parser, arg = "--alpha", short = "-a", type = "double", default = 1.0,
                           help = "Value of 'alpha' parameter of the sPP mixing measure")
opt_parser <- add_argument(opt_parser, arg = "--lambda", short = "-l", type = "double", default = 0.35,
                           help = "Value for the 'lambda' parameter of the sPP mixing measure")
extra_args <- parse_args(opt_parser)

# Summarise results in single .csv file -----------------------------------

# Find parent folder of current file and set working directory
args <- commandArgs()
basedir <- dirname(dirname(sub("--file=", "", args[grep("--file=", args)])))
basedir <- normalizePath(file.path(getwd(), basedir))
setwd(basedir)
cat(sprintf("Current Directory: %s\n", getwd())) # Log

# Check input directory exists
input_dir <- file.path(getwd(), extra_args$input_dir)
if(!dir.exists(input_dir)){
  stop(sprintf("%s does not exist", input_dir))
}
input_dir <- normalizePath(input_dir)
cat(sprintf("Search directory: %s\n", input_dir)) # Log

# Get parameters to set pattern to search
selected_files <- sprintf("NumClust_alpha-%g_lambda-%g-[0-9]+.csv", extra_args$alpha, extra_args$lambda)
files_to_parse <- list.files(input_dir, pattern = selected_files, full.names = TRUE)

# Set CI level
level = 0.05

# Prepare and fill buffer
out <- data.frame("N" = 1, "Mean" = 1, "q_0.025" = 1, "q_0.5" = 1, "q_0.975" = 1)
for (file in files_to_parse) {
  # Deduce number of splits from file
  num_splits <- regmatches(file, gregexpr('[0-9]+.csv', file))
  num_splits <- as.numeric(gsub(".csv", "", num_splits))
  # Parse file
  Nclust <- as.numeric(unlist(read.csv(file, header = F)))
  # Add new line to out data frame
  out <- rbind(out, c(num_splits^2, mean(Nclust), quantile(Nclust, c(level/2, 0.5, 1-level/2))))
}

# Sort ascending according to N
out <- out[order(out$N),]

# Remove all used files
success <- file.remove(files_to_parse)

# Save summary file 
filename <- sprintf("%s/NumClust_alpha-%g_lambda-%g.csv", input_dir, extra_args$alpha, extra_args$lambda)
write.csv(out, file = filename, row.names = FALSE)

# # Generate plot
# library("ggplot2")
# library("dplyr")
# # Get 0.001
# Nclust_df <- read.csv("output/NumClust_alpha-1_lambda-0.001.csv")
# Nclust_df <- Nclust_df %>% 
#   select(c("N", "Mean")) %>%
#   left_join(data.frame("N"=Nclust_df$N, "alpha"=rep(1,length(Nclust_df)), "lambda"=rep(1e-3,length(Nclust_df))), by = "N")
# # Add 0.01
# tmp <- read.csv("output/NumClust_alpha-1_lambda-0.01.csv") %>% 
#   select(c("N", "Mean")) %>%
#   left_join(data.frame("N"=Nclust_df$N, "alpha"=rep(1,length(Nclust_df)), "lambda"=rep(1e-2,length(Nclust_df))), by = "N")
# Nclust_df <- rbind(Nclust_df, tmp); rm(tmp)
# # Add 0.1
# tmp <- read.csv("output/NumClust_alpha-1_lambda-0.1.csv") %>% 
#   select(c("N", "Mean")) %>%
#   left_join(data.frame("N"=Nclust_df$N, "alpha"=rep(1,length(Nclust_df)), "lambda"=rep(1e-1,length(Nclust_df))), by = "N")
# Nclust_df <- rbind(Nclust_df, tmp); rm(tmp)
# # Proper cast
# Nclust_df$alpha <- as.factor(Nclust_df$alpha)
# Nclust_df$lambda <- as.factor(Nclust_df$lambda)
# # plot
# plt <- ggplot(data = Nclust_df) +
#   geom_line(aes(x=N, y=Mean, colour=lambda)) + # , linetype=alpha)) +
#   xlab("Num. Obs.") + ylab("Num. Clust") + 
#   guides(color=guide_legend(bquote(lambda), position = "bottom", direction = "horizontal"))

# # Barplot number of neighbours (LOMBARDIA)
# load("../ArealData_Scenario2/input/clean_data.dat")
# nb_mat_lombardy <- spdep::nb2mat(spdep::poly2nb(st_geometry(sf_mun)), style = "B", zero.policy = TRUE)
# NB <- apply(nb_mat_lombardy, 1, sum)
# title <- sprintf("../ArealData_Scenario2/output/NB_%s.pdf", "Lombardia")
# pdf(file = title, height = 4, width = 4); barplot(table(NB), main = "Lombardia"); dev.off()
# # Barplot number of neighbours (PROVINCE)
# for (prov in levels(sf_mun$PROVINCIA)) {
#   adj_mat <- spdep::nb2mat(spdep::poly2nb(st_geometry(sf_mun[which(sf_mun$PROVINCIA == prov),])), style = "B", zero.policy = TRUE)
#   NB <- apply(adj_mat, 1, sum)
#   title <- sprintf("../ArealData_Scenario2/output/NB_%s.pdf", prov)
#   pdf(file = title, height = 4, width = 4); barplot(table(NB), main = prov); dev.off()
# }

