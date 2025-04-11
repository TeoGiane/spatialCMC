#ifndef SPATIALCMC_SHARD_MERGER_H
#define SPATIALCMC_SHARD_MERGER_H

// STL
#include <deque>

// bayesmix
#include "src/includes.h"

// spatialcmc
#include "shard.h"
#include "shard_partition.h"
#include "spatialcmc_utils.h"

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
	double compute_logBF(ShardPartition & lhs, ShardPartition & rhs);

	// Generate a Mone Carlo sample for the L^2 distance between density functions
	// double sample_l2_distance(ShardPartition & lhs, ShardPartition & rhs, bool prior) {
	// 	bayesmix::AlgorithmState::ClusterState state_lhs, state_rhs;
	// 	state_lhs = prior ? *(lhs.sample_prior()) : *(lhs.sample_full_cond());
	// 	state_rhs = prior ? *(rhs.sample_prior()) : *(rhs.sample_full_cond());
	// 	return compute_l2_distance(state_lhs, state_rhs);
	// }

	// Compute L^2 distance between two densities in closed form, using the associated parameters
	// double compute_l2_distance(bayesmix::AlgorithmState::ClusterState & lhs, bayesmix::AlgorithmState::ClusterState & rhs) {
	// 	// Create buffer
	// 	double out = 0;
	// 	if (lhs.has_uni_ls_state() && rhs.has_uni_ls_state()) {
	// 		// Gaussian case
	// 		double lmean = lhs.uni_ls_state().mean(), lvar = lhs.uni_ls_state().var();
	// 		double rmean = rhs.uni_ls_state().mean(), rvar = rhs.uni_ls_state().var();
	// 		// Compute L2 distance squared
	// 		out = std::exp(-(stan::math::LOG_TWO + stan::math::LOG_SQRT_PI + 0.5*std::log(lvar)));
	// 		out += std::exp(-(stan::math::LOG_TWO + stan::math::LOG_SQRT_PI + 0.5*std::log(rvar)));
	// 		out -= 2 * std::exp(stan::math::normal_lpdf(0, lmean - rmean, std::sqrt(lvar + rvar)));
	// 		// Return
	// 		return out;
	// 	} else if (lhs.has_custom_state() && rhs.has_custom_state()) {
	// 		// Poisson case
	// 		spatialcmc::PoissonState unp_state; 
	// 		lhs.custom_state().UnpackTo(&unp_state); double lrate = unp_state.rate();
	// 		rhs.custom_state().UnpackTo(&unp_state); double rrate = unp_state.rate();
	// 		// Compute L2 distance squared
	// 		out = std::exp(-2*lrate + stan::math::log_modified_bessel_first_kind(0, 2*lrate));
	// 		out += std::exp(-2*rrate + stan::math::log_modified_bessel_first_kind(0, 2*rrate));
	// 		out -= 2*std::exp(-lrate -rrate + stan::math::log_modified_bessel_first_kind(0, 2*std::sqrt(lrate*rrate)));
	// 		// Return
	// 		return out;
	// 	} else {
	// 		// No available cluster state
	// 		throw std::runtime_error("compute_l2_distance() not implemented for this cluster state.");
	// 	}
	// }

	// Generate a Mone Carlo sample for the L^2 distance between density functions
	double sample_qoi(ShardPartition & lhs, ShardPartition & rhs, bool prior) {
		bayesmix::AlgorithmState::ClusterState state_lhs, state_rhs;
		// std::cout << "LHS: prior? " << std::boolalpha << prior << std::endl;
		state_lhs = prior ? *(lhs.sample_prior()) : *(lhs.sample_full_cond());
		// std::cout << "RHS: prior? " << std::boolalpha << prior << std::endl;
		state_rhs = prior ? *(rhs.sample_prior()) : *(rhs.sample_full_cond());
		return compute_qoi(state_lhs, state_rhs);
	}

	// Compute L^2 distance between two densities in closed form, using the associated parameters
	double compute_qoi(bayesmix::AlgorithmState::ClusterState & lhs, bayesmix::AlgorithmState::ClusterState & rhs) {
		// Create buffer
		double out = 0;
		if (lhs.has_uni_ls_state() && rhs.has_uni_ls_state()) {
			// Gaussian case
			double lmean = lhs.uni_ls_state().mean(), lvar = lhs.uni_ls_state().var();
			double rmean = rhs.uni_ls_state().mean(), rvar = rhs.uni_ls_state().var();
			// (Approx) W2 distance
			out = std::abs(lmean-rmean) + std::abs(std::sqrt(lvar)-std::sqrt(rvar));
			// out += (lvar + rvar - 2*std::sqrt(lvar*rvar));
			// Compute sum of absolute distance between each element of parameters' vector
			// out = std::abs(lmean-rmean) + std::abs(std::sqrt(lvar)-std::sqrt(rvar));
			// Return
			return out;
		} else if (lhs.has_custom_state() && rhs.has_custom_state()) {
			// Poisson case
			spatialcmc::PoissonState unp_state; 
			lhs.custom_state().UnpackTo(&unp_state); double lrate = unp_state.rate();
			rhs.custom_state().UnpackTo(&unp_state); double rrate = unp_state.rate();
			// (Approx) W2 distance
			out = std::abs(lrate-rrate); //*(lrate-rrate);
			// out += (lrate + rrate -2*std::sqrt(lrate*rrate));
			// Compute sum of absolute distance between each element of parameters' vector
			// out = std::abs(lrate-rrate) + std::abs(std::sqrt(lrate)-std::sqrt(rrate));
			// Return
			return out;
		} else {
			// No available cluster state
			throw std::runtime_error("compute_qoi() not implemented for this cluster state.");
		}
	}

	bool intersects(ShardPartition & lhs, ShardPartition & rhs) {
		// Concatenating vertices
		auto idx = lhs.get_global_cluster_idx();
		int n_row = idx.size();
		auto rhs_idx = rhs.get_global_cluster_idx();
		int n_col = rhs_idx.size();
		// Stack index vectors
		idx.insert(idx.end(), rhs_idx.begin(), rhs_idx.end());
		// Check inersection
		bool condition = global_adj_matrix(idx, idx).topRightCorner(n_row, n_col).sum() > 0;
		// Return
		return condition;
	};

	// Generate spatial partitions form shards at a given iteration
	std::deque<ShardPartition> generate_local_clusters(const size_t & iter);

};

#endif // SPATIALCMC_SHARD_MERGER_H