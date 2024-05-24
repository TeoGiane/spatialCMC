#ifndef SPATIALCMC_UTILS_H
#define SPATIALCMC_UTILS_H

// STL
#include <fstream>
#include <iostream>
#include <stdarg.h>
#include <vector>

// Protobuf
#include <google/protobuf/any.pb.h>

// bayesmix
#include<src/includes.h>

// Define empty string for ArgParse
#define EMPTYSTR std::string("\"\"")

namespace spatialcmc {

	// Utility that checks if a the path provided is actually an existing file
	void check_file_is_readable(const std::string & filename);

	// Utility that parse a .csv file representing the shard allocation of each data
	std::vector<unsigned int> read_shard_allocation_file(const std::string & filename);

  // Unpack google::protobuf::Any message
	template <typename T>
	T unpack_to(const google::protobuf::Any & _packed_msg) {
		T unpacked_msg;
		if (_packed_msg.Is<T>()) {
			_packed_msg.UnpackTo(&unpacked_msg);
		} else {
			throw std::runtime_error("Unable to unpack google::protobuf::Any type.");
		}
		return unpacked_msg;
	};

} // namespace spatialcmc

#endif // SPATIALCMC_UTILS_H