#include "shard.h"

void Shard::set_data(const Eigen::MatrixXd & _data) {
	data = _data;
	return;
};

void Shard::set_hier_covariates(const Eigen::MatrixXd & _hier_covariates) {
	hier_covariates = _hier_covariates;
	return;
};

void Shard::set_hierarchy(const std::shared_ptr<AbstractHierarchy> _hierarchy, const std::string & _hierarchy_prior) {
	hierarchy = _hierarchy;
	bayesmix::read_proto_from_file(_hierarchy_prior, hierarchy->get_mutable_prior());
	return;
};

void Shard::set_mixing(const std::shared_ptr<AbstractMixing> _mixing, const std::string & _mixing_prior) {
	mixing = _mixing;
	bayesmix::read_proto_from_file(_mixing_prior, mixing->get_mutable_prior());
	return;
};

void Shard::set_sampler_parameters(const std::string & _sampler_parameters_proto) {
	bayesmix::read_proto_from_file(_sampler_parameters_proto, &sampler_parameters);
	return;
};

void Shard::set_verbose(bool _verbose) {
	verbose = _verbose;
	return;
};

// DEBUG
std::string Shard::print_hier_prior() const {
	return hierarchy->get_mutable_prior()->DebugString();
};

std::string Shard::print_mix_prior() const {
	return mixing->get_mutable_prior()->DebugString();
};

// Getters
Eigen::MatrixXd Shard::get_data() const {
	return data;
};

// geos::geom::Geometry::Ptr Shard::get_geometry() const {
// 	return std::move(geometry->clone()); //->releaseGeometries());
// };

std::shared_ptr<AbstractHierarchy> Shard::get_hierarchy() const {
	return hierarchy;
}

// Utilities
void Shard::initialize_sampler() {
	if(verbose){
		std::cout << "Initializing sampler in shard... ";
	}
	// Build appropriate algorithm object
	auto & algorithm_factory = AlgorithmFactory::Instance();
	sampler = algorithm_factory.create_object(sampler_parameters.algo_id());
	// Set up the algorithm
	sampler->read_params_from_proto(sampler_parameters);
	sampler->set_data(data);
	sampler->set_hierarchy(hierarchy);
	sampler->set_mixing(mixing);
	if(hierarchy->is_dependent()) {
		sampler->set_hier_covariates(hier_covariates);
	}
	if(mixing->is_dependent()) {
		// compute_adjacency_matrix();
		sampler->set_mix_covariates(adj_matrix);
	}
	sampler->set_verbose(verbose);
	if(verbose){
		std::cout << "Done" << std::endl;
	}
	return;
};

spatialcmc::MCMCChain Shard::sample() {
	// Run sampler
	initialize_sampler();
	BaseCollector * coll = new MemoryCollector();
	sampler->run(coll);
	// Fill the chain container
	bayesmix::AlgorithmState state;
	for (int i = 0; i < coll->get_size(); i++) {
		coll->get_next_state(&state);
		chain.add_state()->CopyFrom(state);
	}
	// Deallocate collector
	delete coll;
	// Return serialized MCMC chain
	return chain;
};

// void Shard::compute_adjacency_matrix() {
// 	// Resize matrix
// 	int N_geoms = GEOSGetNumGeometries_r(ctx, geometry);
// 	adj_matrix.resize(N_geoms, N_geoms);
// 	// Fill it with binary response
// 	for (size_t i = 0; i < adj_matrix.rows()-1; i++) {
// 		const GEOSGeometry* geom_i = GEOSGetGeometryN_r(ctx, geometry, i);
//     for (size_t j = i+1; j < adj_matrix.cols(); j++) {
// 			const GEOSGeometry* geom_j = GEOSGetGeometryN_r(ctx, geometry, j);
// 			adj_matrix(i,j) = GEOSIntersects_r(ctx, geom_i, geom_j) ? 1 : 0;
// 			// GEOSGeom_destroy_r(ctx, geom_j);
//     }
// 		// GEOSGeom_destroy_r(ctx, geom_i);
//   }
//   adj_matrix = adj_matrix.selfadjointView<Eigen::Upper>();
// 	return;
// };