# A Consensus MC Algorithm for Large Areal Data Clustering

**Authors**: Matteo Gianella, Alessandra Guglielmi, Fernando Quintana

## Summary

This repository, at this stage, is a collection of ideas towards the design of a Consensus MC Algorithm for cluster large areal datasets.We rely on a "Divide and Conquer" approach to split the dataset into "shards", on which the same clustering algorithm is executed, in parallel.

The MCMC output in each shard is then merged in a unique "global" Markov Chain, which is of course an approximation of the true one. 
