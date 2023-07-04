#include "spatial_partition.h"

// SpatialPartition::SpatialPartition() {
// 	// Initialize GEOS
// 	ctx = GEOS_init_r();
// 	GEOSContext_setNoticeHandler_r(ctx, spatialcmc::geos_msg_handler);
// 	GEOSContext_setErrorHandler_r(ctx, spatialcmc::geos_msg_handler);
// 	// Initialize writer
// 	writer = GEOSWKTWriter_create_r(ctx);
// }

// SpatialPartition::SpatialPartition(const SpatialPartition & sp) : SpatialPartition() {
// 	*this = sp;
// }

// SpatialPartition & SpatialPartition::operator=(const SpatialPartition & sp) {
// 	data = sp.data;
// 	num_shard = sp.num_shard;
// 	num_cluster = sp.num_cluster;
// 	global_clust_idx = sp.global_clust_idx;
// 	prior_params = sp.prior_params;
// 	geometry = GEOSGeom_clone_r(ctx, sp.geometry);
// 	hierarchy = sp.hierarchy;
// 	return *this;
// }

// SpatialPartition::~SpatialPartition() {
// 	// Destroy GEOS variables
// 	GEOSGeom_destroy_r(ctx, geometry);
// 	GEOSWKTWriter_destroy_r(ctx, writer);
// 	// Terminate geos
// 	GEOS_finish_r(ctx);
// };

void SpatialPartition::set_shard_name(const size_t & _shard) {
	num_shard = _shard;
}

void SpatialPartition::set_cluster_name(const size_t & _global_cluster) {
	num_cluster = _global_cluster;
}

void SpatialPartition::set_data(const Eigen::MatrixXd & _data) {
	data = _data;
}

void SpatialPartition::set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx) {
	global_clust_idx = _global_clust_idx;
}

void SpatialPartition::set_prior_params(const bayesmix::AlgorithmState::HierarchyHypers & _prior_params) {
	prior_params = _prior_params;
}

// void SpatialPartition::set_geometry(GEOSGeometry* & _geometry) {
// 	GEOSGeometry* old_geom = geometry;
// 	geometry = GEOSGeom_clone_r(ctx, _geometry);
// 	GEOSGeom_destroy_r(ctx, old_geom);
// }

void SpatialPartition::set_hierarchy(std::shared_ptr<AbstractHierarchy> _hierarchy) { 
	// Deep clone of the input
	hierarchy = _hierarchy->deep_clone();		
	// Add data to hierarchy
	if(data.rows() > 0) {
		for (size_t i = 0; i < data.rows(); i++) {
			hierarchy->add_datum(i, data.row(i), false);
		}
	} else {
		throw std::runtime_error("You first need to set data!");
	}
}

unsigned int SpatialPartition::get_card() const {
	return hierarchy->get_card();
}

size_t SpatialPartition::get_shard_name() const {
	return num_shard;
}

std::deque<unsigned int> SpatialPartition::get_global_cluster_idx() const {
	return global_clust_idx;
}

bayesmix::AlgorithmState::HierarchyHypers SpatialPartition::get_prior_params() const {
	return prior_params;
}

size_t SpatialPartition::get_num_cluster() const {
	return num_cluster;
}

// GEOSGeometry* SpatialPartition::get_geometry() const {
// 	return GEOSGeom_clone_r(ctx, geometry);
// }

void SpatialPartition::merge(const SpatialPartition & rhs) {
	// New cluster name
	num_cluster = std::min(num_cluster, rhs.num_cluster);
	// Update data in hierarchy here
	for (size_t i = 0; i < rhs.data.rows(); i++){
		hierarchy->add_datum(data.rows() + i, rhs.data.row(i), true);
	}
	// Merge data
	data.conservativeResize(data.rows() + rhs.data.rows(), data.cols());
	data.bottomRightCorner(rhs.data.size(), data.cols()) = rhs.data;
	// Merge global cluster indexes
	global_clust_idx.insert(global_clust_idx.end(), rhs.global_clust_idx.begin(), rhs.global_clust_idx.end());
	// Merge Geometries
	// GEOSGeometry* geom_union = GEOSUnion_r(ctx, geometry, rhs.geometry);
	// geometry = GEOSMakeValid_r(ctx, GEOSGeom_clone_r(ctx, geom_union));
	// GEOSGeom_destroy_r(ctx, geom_union);
	// Return
	return;
}

Eigen::VectorXd SpatialPartition::sample(size_t n, bool prior) {
	// Generate buffer
	Eigen::VectorXd out(n);
	// if(hierarchy->get_state_proto()->has_uni_ls_state()){
	// 	out.resize(n, 2);
	// } else if (hierarchy->get_state_proto()->has_general_state()){
	// 	out.resize(n, hierarchy->get_state_proto()->general_state().size());
	// }

	// Sample
	for (size_t i = 0; i < n; i++) {
		if(prior) {
			hierarchy->sample_prior();
			// hierarchy->get_state_proto()->PrintDebugString();
		} else {
			hierarchy->sample_full_cond();
			// hierarchy->get_state_proto()->PrintDebugString();
		}
		out(i) = get_merging_parameter(*(hierarchy->get_state_proto()));
	}
	return out;
}

// Questa cosa va migliorata perché può dipendere dalla gerarchia
double SpatialPartition::get_merging_parameter(const ClustState & state) {
	if(state.has_uni_ls_state()){
		return state.uni_ls_state().mean() / std::sqrt(state.uni_ls_state().var());
	} else if (state.has_general_state()) {
		return state.general_state().data(0);
	} else {
		throw std::runtime_error("get_merging_parameter() not implemented for this state.");
	}
}

void SpatialPartition::print() {
	// Read geometry
	// geom_string = GEOSWKTWriter_write_r(ctx, writer, geometry);
	// Print debug info
	std::cout << "Spatial Partition: " << num_cluster << " coming from shard " << num_shard << std::endl;
	std::cout << "Data in this partition: " << data.transpose() << std::endl;
	std::cout << "Data size: " << data.size() << std::endl;
	std::cout << "Global data indexes: "; for (auto && elem : global_clust_idx) { std::cout << elem << " "; }; std::cout << std::endl;
	std::cout << "Global data index size: " << global_clust_idx.size() << std::endl;
	std::cout << "Unique Value: " << sample(1, false).row(0) << std::endl;
	std::cout << "HyperParams: " << prior_params.DebugString() << std::endl;
	// std::cout << "Geometry: " << geom_string << std::endl;
	std::cout << std::endl;
	// Destroy created variables
	// GEOSFree_r(ctx, geom_string);
	// Return
	return;
}