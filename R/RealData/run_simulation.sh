#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=96
#PBS -l walltime=24:00:00
#PBS -N spatialCMC_RealData-no-covariates
#PBS -o log/spatialCMC_RealData-no-covariates.out
#PBS -e log/spatialCMC_RealData-no-covariates.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Load required modules
source /opt/mox/spack/share/spack/setup-env.sh
spack load r@4.4.0
spack load sqlite@3.43.2
spack load libtiff@4.5.1
spack load curl@8.4.0 /eedle4k
spack load protobuf@3.21.12

# Rscript execution
Rscript --vanilla src/run_simulation.R