#!/bin/sh

# Install RspatialCMC package
# R CMD INSTALL RspatialCMC/ --clean
Rscript -e 'devtools::install("RspatialCMC/", quick = TRUE)'
Rscript -e 'RspatialCMC::build_spatialCMC()'
