#ifndef SPATIALCMC_SHARD_MERGER_H
#define SPATIALCMC_SHARD_MERGER_H

// STL
#include <deque>

// GEOS C API
#include <geos_c.h>

// bayesmix
#include "src/includes.h"

// spatialcmc
#include "shard.h"
#include "shard_partition.h"
#include "utils.h"

class ShardMerger {
 private:

	// Class members
	std::vector<Shard> shards;
	std::vector<spatialcmc::MCMCChain> mcmc_chains;
	std::vector<std::deque<int>> data_idx_in_shard;
	Eigen::MatrixXd global_adj_matrix;
	
	unsigned int global_card = 0;
	unsigned int num_shards = 0;
	unsigned int num_iter = 0;

 public:

	// Constructor 
	ShardMerger(const std::vector<Shard> & _shards,
							const std::vector<spatialcmc::MCMCChain> & _mcmc_chains,
							const std::vector<std::deque<int>> & _data_idx_in_shard,
							const Eigen::MatrixXd & _global_adj_matrix);

	// Destructor
	~ShardMerger() = default;

	// Setters
	// void set_shard_partition_list(const std::vector<std::deque<int>>  & _data_idx_in_shard);
	// void set_mcmc_chains(const std::vector<spatialcmc::MCMCChain> & _mcmc_chains);
	// void set_shards(const std::vector<Shard> & _shards);

	// Getters
	unsigned int get_num_iter() const;

	// Merge local clusters at a given iteration
	bayesmix::AlgorithmState merge(size_t iter);

 private:

	// Compute bayes factor for two spatial partitions merging candidates
	double compute_bayes_factor(ShardPartition & lhs, ShardPartition & rhs);

	bool intersects(ShardPartition & lhs, ShardPartition & rhs) {
		// std::cout << "ShardMerger::intersects()" << std::endl;
		// Concatenating vertices
		auto idx = lhs.get_global_cluster_idx();
		int n_row = idx.size();
		auto rhs_idx = rhs.get_global_cluster_idx();
		int n_col = rhs_idx.size();

		// Stack index vectors
		idx.insert(idx.end(), rhs_idx.begin(), rhs_idx.end());

		// Check inersection
		bool condition = global_adj_matrix(idx, idx).topRightCorner(n_row, n_col).sum() > 0;
		// std::cout << "Mat: (" <<  mat.rows() << "," << mat.cols() << ")" << std::endl;
		// std::cout << mat << std::endl;

		// Return
		return condition;
	};

	// Generate spatial partitions form shards at a given iteration
	std::deque<ShardPartition> generate_local_clusters(const size_t & iter);

};

#endif // SPATIALCMC_SHARD_MERGER_H