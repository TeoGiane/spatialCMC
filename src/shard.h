#ifndef SPATIALCMC_SHARD_H_
#define SPATIALCMC_SHARD_H_

#include <memory>

// #include <geos_c.h>
// #include <geos/geom/util/GeometryCombiner.h>

#include <src/includes.h>
#include <stan/math.hpp>

#include "spatialcmc_utils.h"

#include "mcmc_chain.pb.h"

class Shard {
 private:

	// GEOS context handler
	// GEOSContextHandle_t ctx;

	// Just for printing stuff
	// GEOSWKTWriter* writer = nullptr;
	// char* geom_string = nullptr;

	// What is required for the proper initialization of the algorithm
	bayesmix::AlgorithmParams sampler_parameters;

	// What is needed to run mcmc in shard (per ora fisso, opportuna estensione bayesmix)
	std::shared_ptr<AbstractHierarchy> hierarchy;
	std::shared_ptr<AbstractMixing> mixing;
	std::shared_ptr<BaseAlgorithm> sampler;

	// Data we need to store in each shard
	Eigen::MatrixXd data;
	// GEOSGeometry* geometry = nullptr;
	Eigen::MatrixXd adj_matrix;
	Eigen::MatrixXd hier_covariates;

	// Chain output
	spatialcmc::MCMCChain chain;

	// Verbose parameter
	bool verbose = true;

 public:

	// Protected utilities
	// void compute_adjacency_matrix();

	// Constructor and destructor
	Shard() = default;
	// {
	// 	// Initialize GEOS
	// 	ctx = GEOS_init_r();
	// 	GEOSContext_setNoticeHandler_r(ctx, spatialcmc::geos_msg_handler);
	// 	GEOSContext_setErrorHandler_r(ctx, spatialcmc::geos_msg_handler);
	// 	// Initialize writer
	// 	writer = GEOSWKTWriter_create_r(ctx);
	// };

	// Shard & operator=(const Shard & s) {
	// 	sampler_parameters = s.sampler_parameters;
	// 	hierarchy = s.hierarchy;
	// 	mixing = s.mixing;
	// 	sampler = s.sampler;
	// 	data = s.data;
	// 	geometry = GEOSGeom_clone_r(ctx, s.geometry);
	// 	adj_matrix = s.adj_matrix;
	// 	hier_covariates = s.hier_covariates;
	// 	chain = s.chain;
	// 	verbose = s.verbose;
	// 	return *this;
	// };

	// Shard(const Shard & s) : Shard() {
	// 	*this = s;
	// };

	~Shard() = default;
	// {
	// 	// Destroy GEOS variables
	// 	GEOSGeom_destroy_r(ctx, geometry);
	// 	GEOSWKTWriter_destroy_r(ctx, writer);
	// 	// Terminate geos
	// 	GEOS_finish_r(ctx);
	// };

	// Setters
	void set_data(const Eigen::MatrixXd & _data);

	// void set_geometry(GEOSGeometry* & _geometry) {
	// 	GEOSGeometry* old_geom = geometry;
	// 	geometry = GEOSGeom_clone_r(ctx, _geometry);
	// 	GEOSGeom_destroy_r(ctx, old_geom);
	// };

	void set_adjacency_matrix(const Eigen::MatrixXd & _adj_matrix) {
		adj_matrix = _adj_matrix;
		return;
	};

	void set_hier_covariates(const Eigen::MatrixXd & _hier_covariates);
	void set_hierarchy(const std::shared_ptr<AbstractHierarchy> _hierarchy, const std::string & _hierarchy_prior);
	void set_mixing(const std::shared_ptr<AbstractMixing> _mixing, const std::string & _mixing_prior);
	void set_sampler_parameters(const std::string & _sampler_parameters_proto);
	void set_verbose(bool _verbose);

	// DEBUG
	std::string print_hier_prior() const;
	std::string print_mix_prior() const;

	// Getters
	Eigen::MatrixXd get_data() const;
	// GEOSGeometry* get_geometry() const { return GEOSGeom_clone_r(ctx, geometry); };
	// void* get_geometry_addess() const { return geometry; }; // Just for debug purposes
	std::shared_ptr<AbstractHierarchy> get_hierarchy() const;
	Eigen::MatrixXd get_adjacency_matrix() const { return adj_matrix; };

	// Utilities
	void initialize_sampler();
	spatialcmc::MCMCChain sample();

	void print() {
		// Read geometry
		// geom_string = GEOSWKTWriter_write_r(ctx, writer, geometry);
		// Print
		std::cout << "Shard debug info" << std::endl;
		std::cout << "data: " << data.transpose() << std::endl;
		std::cout << "adj_matrix: (" << adj_matrix.rows() << "," << adj_matrix.cols() << "):" << std::endl;
		std::cout << adj_matrix << std::endl;
		// std::cout << "geometry: " << geom_string << std::endl;
		std::cout << "hier_id: " << hierarchy->get_id() << std::endl;
		std::cout << "hier_prior:\n" << print_hier_prior() << std::endl;
		std::cout << "mix_prior:\n" << print_mix_prior() << std::endl;
		std::cout << "mix_id: " << mixing->get_id() << std::endl;
		std::cout << "Algorithm Params:\n" << sampler_parameters.DebugString() << std::endl;
		// Destroy created variables
		// GEOSFree_r(ctx, geom_string);
	};
};

#endif // SPATIALCMC_SHARD_H_
