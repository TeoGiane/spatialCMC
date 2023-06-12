#ifndef SPATIALCMC_SPATIAL_PARTITION_H
#define SPATIALCMC_SPATIAL_PARTITION_H

// STL
#include <iostream>
#include <deque>

// GEOS C API
#include <geos_c.h>

// bayesmix
#include <src/includes.h>

// stan-math
#include <stan/math.hpp>

// spatialcmc
#include "utils.h"

class SpatialPartition {
 private:

	// GEOS context handler
	// GEOSContextHandle_t ctx;

	// Just for printing stuff
	// GEOSWKTWriter* writer = nullptr;
	// char* geom_string = nullptr;

	// Class members
	Eigen::MatrixXd data;
	size_t num_shard;
	size_t num_cluster;
	std::deque<unsigned int> global_clust_idx;
	bayesmix::AlgorithmState::HierarchyHypers prior_params;
	// GEOSGeometry* geometry = nullptr;
	std::shared_ptr<AbstractHierarchy> hierarchy;

	// Type aliases
	using ClustState = bayesmix::AlgorithmState::ClusterState;

 public:

	// Constructor
	SpatialPartition() = default;
	
	// Copy-constructor
	// SpatialPartition(const SpatialPartition & sp);

	// Copy-assignment
	// SpatialPartition & operator=(const SpatialPartition & sp);

	// Destructor
	~SpatialPartition() = default;

	// Setters
	void set_shard_name(const size_t & _shard);
	void set_cluster_name(const size_t & _global_cluster);
	void set_data(const Eigen::MatrixXd & _data);
	void set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx);
	void set_prior_params(const bayesmix::AlgorithmState::HierarchyHypers & _prior_params);
	// void set_geometry(GEOSGeometry* & _geometry);
	void set_hierarchy(std::shared_ptr<AbstractHierarchy> _hierarchy);

	// Getters
	unsigned int get_card() const;
	size_t get_shard_name() const;
	std::deque<unsigned int> get_global_cluster_idx() const;
	bayesmix::AlgorithmState::HierarchyHypers get_prior_params() const;
	size_t get_num_cluster() const;
	// GEOSGeometry* get_geometry() const;

	// Merge current spatial partition with the input spatial partition
	void merge(const SpatialPartition & rhs);

	// Generate n samples from the prior or the full conditional of the hierarchy
	Eigen::VectorXd sample(size_t n, bool prior);

	// Questa cosa va migliorata perché può dipendere dalla gerarchia
	double get_merging_parameter(const ClustState & state);

	// Print utility for debug
	void print();
};

#endif // SPATIALCMC_SPATIAL_PARTITION_H