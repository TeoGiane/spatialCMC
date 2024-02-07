#include "spatialcmc_utils.h"

namespace spatialcmc {

	void check_file_is_readable(const std::string & filename) {
		std::ifstream ifstr(filename);
		if (ifstr.fail())
			throw std::runtime_error("ERROR: File '" + filename + "' can not be read");
		return;
	}

	// void geos_msg_handler(const char* fmt, ...) {
	// 	va_list ap;
	// 	va_start(ap, fmt);
	// 	vprintf (fmt, ap);
	// 	va_end(ap);
	// }

	// void read_geometry_file(const std::string & filename, GEOSGeometry* & output_geom) {
	
	// 	// Initialize local context
	// 	GEOSContextHandle_t ctx = GEOS_init_r();
	// 	GEOSContext_setNoticeHandler_r(ctx, geos_msg_handler);
	// 	GEOSContext_setErrorHandler_r(ctx, geos_msg_handler);
		
	// 	// Build GEOS reader and final buffer
	// 	GEOSWKTReader * reader = GEOSWKTReader_create_r(ctx);
	// 	std::vector<GEOSGeometry*> geom_vect;

	// 	// Parse geometry file
	// 	std::ifstream geom_file(filename);
	// 	for (std::string curr_geom; std::getline(geom_file, curr_geom); ) {
	// 		geom_vect.push_back(GEOSWKTReader_read_r(ctx, reader, curr_geom.c_str()));
	// 	}
	// 	geom_file.close();

	// 	// Create collection from the vector
	// 	output_geom = GEOSGeom_createCollection_r(ctx, GEOS_GEOMETRYCOLLECTION, geom_vect.data(), geom_vect.size());

	// 	// Clean up everything we allocated
	// 	GEOSWKTReader_destroy_r(ctx, reader);
	// 	GEOS_finish_r(ctx);
	// }

	std::vector<unsigned int> read_shard_allocation_file(const std::string & filename) {
  	auto shard_allocs_matrix = bayesmix::read_eigen_matrix(filename);
  	std::vector<unsigned int> shard_allocs(shard_allocs_matrix.data(), shard_allocs_matrix.data() + shard_allocs_matrix.size());
  	return shard_allocs;
	}

} // namespace spatialcmc