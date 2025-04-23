#include "poisson_regression_local_cluster_merger.h"

// Constructor
PoissonRegLocalClusterMerger::PoissonRegLocalClusterMerger(const std::vector<std::vector<bayesmix::AlgorithmState>> & _cluster_chains) : LocalClusterMerger(_cluster_chains) {
    // Resize reg_coeffs_chain
    reg_coeffs_chain.resize(num_iter);
}

// Merge local clusters at a given iteration
bayesmix::AlgorithmState PoissonRegLocalClusterMerger::merge(size_t iter) {  
    // Generate local clusters 
    auto local_clusters = generate_local_clusters(iter);
    // Shuffle local cluster
    auto &rng = bayesmix::Rng::Instance().get();
    std::shuffle(local_clusters.begin(), local_clusters.end(), rng);  
    // Merging at a fixed iteration
    for (auto it = local_clusters.begin(); it != local_clusters.end(); ++it) {
        // Get lhs as reference (just to make code a little bit more readable)
        auto & lhs = *it;
        // Lambda function that computes the merge condition (with automatic merging)
        auto merge_condition = [this, &lhs] (PoissonRegLocalCluster & rhs) {
            bool merge = (lhs.get_shard_label() != rhs.get_shard_label());
            if(merge) { merge *= intersects(lhs, rhs); }
            if(merge) { merge *= (2*compute_logBF(lhs, rhs) > 0.0); }
            if(merge) { lhs.merge_with(rhs); }
            return merge;
        };
        // Verify the condition for all the remaining local_clusters and erase elements that have been merged
        auto end = std::remove_if(it, local_clusters.end(), merge_condition);
        local_clusters.erase(end, local_clusters.end());
    }
    // Build bayesmix::AlgorithmState object with 
    bayesmix::AlgorithmState merged_state = cluster_chains[0][iter];
    merged_state.clear_cluster_states();
    merged_state.clear_cluster_allocs();
    Eigen::VectorXi global_cluster_allocs(global_card);
    unsigned int k = 0;
    for (auto && elem : local_clusters) {
        // Set new cluster states in merged state message
        merged_state.add_cluster_states()->CopyFrom(*elem.sample_full_cond());
        // Computing global cluster allocs
        global_cluster_allocs(elem.get_global_cluster_idx()).array() = k;
        k++;
    }
    // Set new cluster allocs in merged state message
    *merged_state.mutable_cluster_allocs() = { global_cluster_allocs.begin(), global_cluster_allocs.end() };
    // Return final state
    return merged_state;
}

// Generate spatial partitions form shards at a given iteration (DA RIVEDERE)
std::deque<PoissonRegLocalCluster> PoissonRegLocalClusterMerger::generate_local_clusters(const size_t & iter) {
    // Create buffer for output
    std::deque<PoissonRegLocalCluster> local_clusters;
    // Global cluster counter is initialized
    size_t clust_counter = 0;
    for (size_t it = 0; it < num_shards; it++) {
        // Get current state of the MCMC chain
        auto curr_state = cluster_chains[it][iter];
        // Get info from shards
        auto data_in_shard = (*data)(global_numbering[it], Eigen::all);
        auto adj_matrix_in_shard = (*adj_matrix)(global_numbering[it], Eigen::all);
        auto cov_matrix_in_shard = (*cov_matrix)(global_numbering[it], Eigen::all);
        // Get local cluster indexes and partition in each cluster
        std::vector<std::deque<unsigned int>> local_clust_idx(curr_state.cluster_states_size());
        std::vector<std::deque<unsigned int>> global_clust_idx(curr_state.cluster_states_size());
        for (size_t j = 0; j < curr_state.cluster_allocs_size(); j++) {
            unsigned int c = curr_state.cluster_allocs(j);
            local_clust_idx[c].push_back(j);
            global_clust_idx[c].push_back(global_numbering[it][j]);
        }
        // Building ShardPartition list
        for (size_t k = 0; k < curr_state.cluster_states_size(); k++) {
            // Construct and fill spatial partition
            PoissonRegLocalCluster lc(hierarchy->deep_clone());
            lc.set_hierarchy_prior(hier_prior_file);
            lc.set_shard_label(it);
            lc.set_cluster_label(clust_counter + k);
            lc.set_global_cluster_idx(global_clust_idx[k]);
            lc.set_data(data_in_shard(local_clust_idx[k], Eigen::all), cov_matrix_in_shard(local_clust_idx[k], Eigen::all));
            lc.set_reg_coeffs(reg_coeffs_chain[iter]);
            local_clusters.push_back(lc);
        }
        // Update global cluster counter
        clust_counter += curr_state.cluster_states_size();
    }
    // Return
    return local_clusters;
}
