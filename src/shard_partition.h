#ifndef SPATIALCMC_SHARD_PARTITION_H
#define SPATIALCMC_SHARD_PARTITION_H

// STL
#include <iostream>
#include <deque>

// bayesmix
#include <src/includes.h>

// stan-math
#include <stan/math.hpp>

// spatialcmc
#include "poisson_state.pb.h"
#include "spatialcmc_utils.h"

class ShardPartition {

 private:

	// Class members
	Eigen::MatrixXd data;
	size_t num_shard;
	size_t num_cluster;
	std::deque<unsigned int> global_clust_idx;
	bayesmix::AlgorithmState::HierarchyHypers prior_params;
	std::shared_ptr<AbstractHierarchy> hierarchy;

	// Type aliases
	using ClustState = bayesmix::AlgorithmState::ClusterState;

	public:

	// Constructor
	ShardPartition() = default;

	// Destructor
	~ShardPartition() = default;

	// Setters
	void set_shard_name(const size_t & _shard);
	void set_cluster_name(const size_t & _global_cluster);
	void set_data(const Eigen::MatrixXd & _data);
	void set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx);
	void set_prior_params(const bayesmix::AlgorithmState::HierarchyHypers & _prior_params);
	void set_hierarchy(std::shared_ptr<AbstractHierarchy> _hierarchy);

	// Getters
	unsigned int get_card() const;
	size_t get_shard_name() const;
	std::deque<unsigned int> get_global_cluster_idx() const;
	bayesmix::AlgorithmState::HierarchyHypers get_prior_params() const;
	size_t get_num_cluster() const;

	// Merge current spatial partition with the input spatial partition
	void merge(const ShardPartition & rhs);

	// Returns a sample from the full conditional distribution of the hierarchy
	std::shared_ptr<ClustState> sample_full_cond() {
		hierarchy->sample_full_cond();
		return hierarchy->get_state_proto();
	}

	// Print utility for debug
	void print();

	// Generate n samples of the quantity of interest either from the prior
	// or the full conditional of the hierarchy.
	Eigen::VectorXd sample_qoi(size_t n, bool prior);

	private:

	// Compute quantity of interest from sampled state of the hierarchy
	double qoi_from_state(const ClustState & state) const ;

};

#endif // SPATIALCMC_SHARD_PARTITION_H