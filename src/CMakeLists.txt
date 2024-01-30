# Add sources directories
add_subdirectory(hierarchy)
add_subdirectory(mixing)

# Add sources files
target_sources(spatialCMC PUBLIC
	shard.h
	shard.cc
	shard_merger.h
	shard_merger.cc
	shard_partition.h
	shard_partition.cc
	utils.h
	utils.cc
)
