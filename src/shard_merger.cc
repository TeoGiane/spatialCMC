#include "shard_merger.h"

ShardMerger::ShardMerger(const std::vector<Shard> & _shards,
												 const std::vector<spatialcmc::MCMCChain> & _mcmc_chains,
												 const std::vector<std::deque<int>> & _data_idx_in_shard,
												 const Eigen::MatrixXd & _global_adj_matrix) : 
	shards(_shards), mcmc_chains(_mcmc_chains), data_idx_in_shard(_data_idx_in_shard), global_adj_matrix(_global_adj_matrix) {
	// Set missing quantities from input
	global_card = global_adj_matrix.rows();
	num_iter = mcmc_chains[0].state_size();
	num_shards = shards.size();
	// for (size_t i = 0; i < num_shards; i++) {
	// 	global_card += shards[i].get_data().rows();
	// }
}

// void ShardMerger::set_shard_partition_list(const std::vector<std::deque<int>>  & _data_idx_in_shard) {
// 	data_idx_in_shard = _data_idx_in_shard;
// }

// void ShardMerger::set_mcmc_chains(const std::vector<spatialcmc::MCMCChain> & _mcmc_chains) {
// 	mcmc_chains = _mcmc_chains;
// }

// void ShardMerger::set_shards(const std::vector<Shard> & _shards) {
// 	shards = _shards; num_shards = shards.size();
// }

unsigned int ShardMerger::get_num_iter() const {
	return num_iter;
}

bayesmix::AlgorithmState ShardMerger::merge(size_t iter) {

	// Generate local clusters 
	auto local_clusters = generate_local_clusters(iter);
	// std::cout << "Generated local clusters" << std::endl;

	// Shuffle local cluster
	std::random_device rd; std::mt19937 g(rd());
	std::shuffle(local_clusters.begin(), local_clusters.end(), g);

	// Check - OK
	// std::cout << "Before Merging: " << std::endl << std::endl;
	// for (auto &&lc : local_clusters) { lc.print(); }
	
	// Merging at a fixed iteration
	for (auto it = local_clusters.begin(); it != local_clusters.end(); ++it) {
		
		// Get lhs as reference (just to make code a little bit more readable)
		auto & lhs = *it;

		// Lambda function that computes the merge condition (with automatic merging)
		auto merge_condition = [this, &lhs/*, &ctx*/] (ShardPartition & rhs) {
			bool merge = (lhs.get_shard_name() != rhs.get_shard_name());
			if(merge) { merge *= intersects(lhs, rhs); }
			if(merge) { merge *= (2*compute_logBF(lhs, rhs) > 10.0); }
			if(merge) { lhs.merge(rhs); }
			return merge;
		};

		// Verify the condition for all the remaining local_clusters and erase elements that have been merged
		auto end = std::remove_if(it, local_clusters.end(), merge_condition);
		local_clusters.erase(end, local_clusters.end());
	}

	// std::cout << "Completed merging at fixed iteration" << std::endl;

	// Check - OK
	// std::cout << "After Merging: " << std::endl << std::endl;
	// for (auto &&lc : local_clusters) { lc.print(); }

	// Build bayesmix::AlgorithmState object with 
	bayesmix::AlgorithmState merged_state = mcmc_chains[0].state(iter);
	merged_state.clear_cluster_states();
	merged_state.clear_cluster_allocs();
	Eigen::VectorXi global_cluster_allocs(global_card);
	unsigned int k = 0;
	for (auto && elem : local_clusters) {
		// Set new cluster states in merged state message
		// bayesmix::AlgorithmState::ClusterState clust_state;
		// clust_state.mutable_general_state()->add_data(elem.sample(1, false)(0));
		// clust_state.set_cardinality(elem.get_card());
		merged_state.add_cluster_states()->CopyFrom(*elem.sample_full_cond());
		// Computing global cluster allocs
		global_cluster_allocs(elem.get_global_cluster_idx()).array() = k;
		k++;
	}

	// Set new cluster allocs in merged state message
	*merged_state.mutable_cluster_allocs() = { global_cluster_allocs.begin(), global_cluster_allocs.end() };

	// Check - OK
	// std::cout << merged_state.DebugString() << std::endl;

	// Return final state
	return merged_state;
}

std::deque<ShardPartition> ShardMerger::generate_local_clusters(const size_t & iter) {
	
	// Shuffle shards parsing order
	// std::vector<unsigned int> parse_order(num_shards);
	// std::iota(parse_order.begin(), parse_order.end(), 0);
	// std::random_device rd; std::mt19937 g(rd());
	// std::shuffle(parse_order.begin(), parse_order.end(), g);

	// Create buffer for output
	std::deque<ShardPartition> local_clusters;
	
	// Global cluster counter is initialized
	size_t clust_counter = 0;
	
	for (size_t s = 0; s < num_shards; s++) {
		
		// get index of shuffled shard
		// unsigned int it = parse_order[s];
		unsigned int it = s;

		// get current state of the MCMC chain
		auto curr_state = mcmc_chains[it].state(iter);

		// Get info from shards
		// GEOSGeometry* geom_in_shard = shards[it].get_geometry();
		auto data_in_shard = shards[it].get_data();
		auto adj_matrix_in_shard = shards[it].get_adjacency_matrix();

		// Get local cluster indexes and partition in each cluster
		std::vector<std::deque<unsigned int>> local_clust_idx(curr_state.cluster_states_size());
		std::vector<std::deque<unsigned int>> global_clust_idx(curr_state.cluster_states_size());
		for (size_t j = 0; j < curr_state.cluster_allocs_size(); j++) {
			unsigned int c = curr_state.cluster_allocs(j);
			local_clust_idx[c].push_back(j);
			global_clust_idx[c].push_back(data_idx_in_shard[it][j]);
		}

		// Building ShardPartition list
		for (size_t k = 0; k < curr_state.cluster_states_size(); k++) {
			// Construct and fill spatial partition
			ShardPartition sp_to_add;
			sp_to_add.set_shard_name(it);
			sp_to_add.set_cluster_name(clust_counter + k);
			sp_to_add.set_data(data_in_shard(local_clust_idx[k], 0));
			sp_to_add.set_global_cluster_idx(global_clust_idx[k]);
			// sp_to_add.set_unique_value(curr_state.cluster_states(k));
			sp_to_add.set_prior_params(curr_state.hierarchy_hypers());
			// sp_to_add.set_geometry(geom_in_clust[k]);
			sp_to_add.set_hierarchy(shards[it].get_hierarchy());
			// Add spatial partition
			local_clusters.push_back(sp_to_add);
		}
		// Update global cluster counter
		clust_counter += curr_state.cluster_states_size();
	}

	// Return
	return local_clusters;
}

double ShardMerger::compute_logBF(ShardPartition & lhs, ShardPartition & rhs) {
	// Number of Monte Carlo samples
	size_t n_sim = 50;
	// Eigen::VectorXd l2_prior(n_sim), l2_post(n_sim); bool prior = true;
  Eigen::VectorXd qoi_prior(n_sim), qoi_post(n_sim); bool prior = true;
  // Compute Monte Carlo samples from prior and posterior
	for (size_t k = 0; k < n_sim; k++) {
		std::cout << "sample: " << k << std::endl;
    // Using L2 distance
	  // l2_prior(k) = sample_l2_distance(lhs, rhs, prior);
    // l2_post(k) = sample_l2_distance(lhs, rhs, !prior);
    // Using L1 distance of the parameter's vector
  	qoi_prior(k) = sample_qoi(lhs, rhs, prior);
    qoi_post(k) = sample_qoi(lhs, rhs, !prior);
	}
	
	// auto lhs_prior = lhs.sample_qoi(n_sim, true);
	// auto rhs_prior = rhs.sample_qoi(n_sim, true);
	// auto lhs_post = lhs.sample_qoi(n_sim, false);
	// auto rhs_post = rhs.sample_qoi(n_sim, false);

	// Compute quantity of interest
	// Eigen::VectorXd qoi_post = (lhs_post - rhs_post).array().abs(); // / (lhs_post.array().abs() + rhs_post.array().abs());
	// Eigen::VectorXd qoi_prior = (lhs_prior - rhs_prior).array().abs(); // / (lhs_prior.array().abs() + rhs_prior.array().abs());

	// Compute epsilon dynamically
	double epsilon = stan::math::quantile(qoi_prior, 0.5);
	// #pragma omp critical
	// std::cout << "l2_prior on thread " << omp_get_thread_num() << ": " << l2_prior.transpose() << std::endl;	
	// double epsilon = stan::math::quantile(l2_prior, 0.5);

	// Compute logBF
	// double l_post_odd = std::log((l2_post.array() < epsilon).count()) - std::log((l2_post.array() >= epsilon).count());
	double l_post_odd = std::log((qoi_post.array() < epsilon).count()) - std::log((qoi_post.array() >= epsilon).count());
	// double l_prior_odd = std::log((l2_prior.array() < epsilon).count()) - std::log((l2_prior.array() >= epsilon).count()); // == 0 by def of epsilon
	double logBF = l_post_odd; // - l_prior_odd;

	// #pragma omp critical
	// std::cout << "ShardMerger::compute_bayes_factor() in thread" << omp_get_thread_num() << ": epsilon = " << epsilon << ", l_post_odd = " << l_post_odd << ", l_prior_odd = " << l_prior_odd << ", logBF = " << logBF << std::endl;
	// double post_odd = (qoi_post.array() < epsilon).count() / static_cast<double>(n);
	// double prior_odd = (qoi_prior.array() < epsilon).count() / static_cast<double>(n);
	// double logBF = (prior_odd == 0) ?  stan::math::NEGATIVE_INFTY : (std::log(post_odd) - std::log(prior_odd));

	// return logBF
	std::cout << "2logBF: " << 2*logBF << std::endl;
	return logBF;
}