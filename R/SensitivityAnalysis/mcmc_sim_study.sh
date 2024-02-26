#!/bin/bash

# Set bash folder
BASH_FOLDER=$(dirname $BASH_SOURCE)

# Set values for ALPHA, LAMBDA and SCENARIO
ALPHA=(0.1 1 5 10)
LAMBDA=(0) # (0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 1 5)
SCENARIO=(Scenario1) # (Scenario1 Scenario2 Scenario4 Scenario5)

# Make log folders if not present
for scenario in ${SCENARIO[@]}; do
  mkdir -p $BASH_FOLDER/log/$scenario
done

# Execute run_mcmc.R in parallel via GNU parallel
parallel -j 6 \
  'Rscript --vanilla {4}/scripts/run_mcmc.R input/{1}.dat \
    --seed {#} \
    --hier-params-file input/{1}_params.asciipb \
    --alpha {2} \
    --lambda {3} \
    --output-file output/{1}/mcmc/chain_alpha_{2}_lambda_{3}.dat \
    &> {4}/log/{1}/run_mcmc_alpha_{2}_lambda_{3}.log' \
  ::: ${SCENARIO[@]} ::: ${ALPHA[@]} ::: ${LAMBDA[@]} ::: ${BASH_FOLDER}

# Execute generate_plot.R in parallel via GNU parallel
parallel -j 6 \
  'Rscript --vanilla {4}/scripts/generate_plot.R input/{1}.dat \
    --chain-file output/{1}/mcmc/chain_alpha_{2}_lambda_{3}.dat \
    --output-file plots/{1}/mcmc/plot_alpha_{2}_lambda_{3}.pdf \
    &> {4}/log/{1}/generate_plot_alpha_{2}_lambda_{3}.log' \
  ::: ${SCENARIO[@]} ::: ${ALPHA[@]} ::: ${LAMBDA[@]} ::: ${BASH_FOLDER}
