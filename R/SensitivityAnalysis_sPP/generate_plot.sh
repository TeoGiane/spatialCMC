#!/bin/bash

# Set values for ALPHA, LAMBDA and GEOM
ALPHA=(0.1 1 5 10)
LAMBDA=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 1 5)
GEOM=(grid_geometry lombardy_geometry)

# Make log folders if not present
for geom in ${GEOM[@]}; do
  mkdir -p log/$geom
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 6 \
    'Rscript --vanilla scripts/generate_plot.R ./output/{1}/chain_alpha_{2}_lambda_{3}.dat -o ./plots/{1}/plot_alpha_{2}_lambda_{3}.pdf &> ./log/{1}/generate_plot_alpha_{2}_lambda_{3}.log' \
    ::: ${GEOM[@]} ::: ${ALPHA[@]} ::: ${LAMBDA[@]}
