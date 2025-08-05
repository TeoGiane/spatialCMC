#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=96
#PBS -l walltime=24:00:00
#PBS -N BGS_sim_study
#PBS -o log/run_BGS_sim_study.out
#PBS -e log/run_BGS_sim_study.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Load required modules
source /opt/mox/spack/share/spack/setup-env.sh
spack load r@4.4.0