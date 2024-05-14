#include "shard_partition.h"

void ShardPartition::set_shard_name(const size_t & _shard) {
	num_shard = _shard;
}

void ShardPartition::set_cluster_name(const size_t & _global_cluster) {
	num_cluster = _global_cluster;
}

void ShardPartition::set_data(const Eigen::MatrixXd & _data) {
	data = _data;
}

void ShardPartition::set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx) {
	global_clust_idx = _global_clust_idx;
}

void ShardPartition::set_prior_params(const bayesmix::AlgorithmState::HierarchyHypers & _prior_params) {
	prior_params = _prior_params;
}

void ShardPartition::set_hierarchy(std::shared_ptr<AbstractHierarchy> _hierarchy) { 
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

unsigned int ShardPartition::get_card() const {
	return hierarchy->get_card();
}

size_t ShardPartition::get_shard_name() const {
	return num_shard;
}

std::deque<unsigned int> ShardPartition::get_global_cluster_idx() const {
	return global_clust_idx;
}

bayesmix::AlgorithmState::HierarchyHypers ShardPartition::get_prior_params() const {
	return prior_params;
}

size_t ShardPartition::get_num_cluster() const {
	return num_cluster;
}

void ShardPartition::merge(const ShardPartition & rhs) {
	// New cluster name
	num_shard = std::min(num_shard, rhs.num_shard);
	num_cluster = std::min(num_cluster, rhs.num_cluster);
	// Update data in hierarchy here
	for (size_t i = 0; i < rhs.data.rows(); i++){
		hierarchy->add_datum(data.rows() + i, rhs.data.row(i));
	}
	// Merge data
	data.conservativeResize(data.rows() + rhs.data.rows(), data.cols());
	data.bottomRightCorner(rhs.data.size(), data.cols()) = rhs.data;
	// Merge global cluster indexes
	global_clust_idx.insert(global_clust_idx.end(), rhs.global_clust_idx.begin(), rhs.global_clust_idx.end());
}

Eigen::VectorXd ShardPartition::sample_qoi(size_t n, bool prior) {
	// Generate buffer
	Eigen::VectorXd out(n);
	// Sample
	for (size_t i = 0; i < n; i++) {
		if(prior) {
			hierarchy->sample_prior();
		} else {
			hierarchy->sample_full_cond();
		}
		// Compute quantity of interest from the current sampled value
		out(i) = qoi_from_state(*(hierarchy->get_state_proto()));
	}
	// Return
	return out;
}

double ShardPartition::qoi_from_state(const ClustState & state) const {
	if(state.has_uni_ls_state()) {
		return state.uni_ls_state().mean() / std::sqrt(state.uni_ls_state().var());
	} else if (state.has_custom_state()) {
		// Unpack custom state
		// std::cout << "Here!!" << std::endl;
		spatialcmc::PoissonState unp_state; state.custom_state().UnpackTo(&unp_state);
    	// auto unp_state = spatialcmc::unpack_protobuf_any<spatialcmc::PoissonState>(state.custom_state());
		return unp_state.rate() / std::sqrt(unp_state.rate());
	} else {
		throw std::runtime_error("qoi_from_state() not implemented for this state.");
	}
}

void ShardPartition::print() {
	std::cout << "Partition: " << num_cluster << " coming from shard " << num_shard << std::endl;
	std::cout << "Data in this partition: " << data.transpose() << std::endl;
	std::cout << "Data size: " << data.size() << std::endl;
	std::cout << "Global data indexes: "; for (auto && elem : global_clust_idx) { std::cout << elem << " "; }; std::cout << std::endl;
	std::cout << "Global data index size: " << global_clust_idx.size() << std::endl;
	std::cout << "Unique Value:\n " << sample_full_cond()->DebugString() << std::endl;
	std::cout << "HyperParams:\n " << prior_params.DebugString() << std::endl;
	std::cout << std::endl;
}