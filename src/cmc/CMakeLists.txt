# Add sources files
target_sources(spatialCMC PUBLIC
	shard.h
	shard.cc
	shard_merger.h
	shard_merger.cc
	shard_partition.h
	shard_partition.cc
	spatialcmc_utils.h
	spatialcmc_utils.cc
	local_cluster.h
	local_cluster_imp.h
	local_cluster_merger.h
	local_cluster_merger_imp.h
	regression_coefficients_merger.h
	regression_coefficients_merger.cc
)
