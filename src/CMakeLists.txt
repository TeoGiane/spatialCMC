# Add sources directories
add_subdirectory(hierarchy)
add_subdirectory(mixing)

# Add sources files
target_sources(spatialCMC PUBLIC
	shard.h
	shard.cc
	shard_merger.h
	shard_merger.cc
	spatial_partition.h
	spatial_partition.cc
	utils.h
	utils.cc
)
