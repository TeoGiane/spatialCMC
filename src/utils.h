#ifndef SPATIALCMC_UTILS_H
#define SPATIALCMC_UTILS_H

// STL
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include<vector>

// GEOS C API
#include<geos_c.h>

// bayesmix
#include<src/includes.h>

// Define empty string for ArgParse
#define EMPTYSTR std::string("\"\"")

namespace spatialcmc {

	// Utility that checks if a the path provided is actually an existing file
	void check_file_is_readable(const std::string & filename);

	// Utility that parse a .csv file representing the shard allocation of each data
	std::vector<unsigned int> read_shard_allocation_file(const std::string & filename);

} // namespace spatialcmc

#endif // SPATIALCMC_UTILS_H