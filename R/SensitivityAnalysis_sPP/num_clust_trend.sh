#!/bin/bash

# Set values for ALPHA, LAMBDA and NUM_SPLITS
ALPHA=1
LAMBDA=0.01
MAX_NUM_SPLITS=60

# Make log folders if not present
mkdir -p log

# Execute sample_eppf_Nvarying.R in parallel via GNU parallel
parallel -j 6 \
    'Rscript --vanilla scripts/sample_eppf_Nvarying.R -n {1} -a {2} -l {3} -o ./output &> ./log/NumClust_alpha-{2}_lambda-{3}-{1}.log' \
    ::: $(seq 2 ${MAX_NUM_SPLITS}) ::: $ALPHA ::: $LAMBDA
    
# Execute aggregate_results.R
parallel -j 1 \
    'Rscript --vanilla scripts/aggregate_results.R -a {1} -l {2} ./output &> ./log/NumClust_alpha-{1}_lambda-{2}.log' \
    ::: $ALPHA ::: $LAMBDA
