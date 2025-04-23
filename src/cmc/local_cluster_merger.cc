#include "local_cluster_merger.h"

// Constructor
LocalClusterMerger::LocalClusterMerger(const std::vector<std::vector<bayesmix::AlgorithmState>> & _cluster_chains) : cluster_chains(_cluster_chains) {
    // Set members
    num_shards = cluster_chains.size();
    num_iter = cluster_chains[0].size();
}

// Merge local clusters at a given iteration
bayesmix::AlgorithmState LocalClusterMerger::merge(size_t iter) {
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
        auto merge_condition = [this, &lhs] (LocalCluster & rhs) {
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

// Compute bayes factor for two spatial partitions merging candidates
double LocalClusterMerger::compute_logBF(LocalCluster & lhs, LocalCluster & rhs) {
    // Number of Monte Carlo samples
    size_t n_sim = 100;
    Eigen::VectorXd qoi_prior(n_sim), qoi_post(n_sim); bool prior = true;
    // Compute Monte Carlo samples from prior and posterior
    for (size_t k = 0; k < n_sim; k++) {
        qoi_prior(k) = sample_qoi(lhs, rhs, prior);
        qoi_post(k) = sample_qoi(lhs, rhs, !prior);
    }
    // Compute epsilon dynamically
    double epsilon = stan::math::quantile(qoi_prior, 0.5);
    // Compute logBF
    double l_post_odd = std::log((qoi_post.array() < epsilon).count()) - std::log((qoi_post.array() >= epsilon).count());
    double logBF = l_post_odd;
    // Return logBF
    return logBF;
}

double LocalClusterMerger::sample_qoi(LocalCluster & lhs, LocalCluster & rhs, bool prior) {
    bayesmix::AlgorithmState::ClusterState state_lhs, state_rhs;
    state_lhs = prior ? *(lhs.sample_prior()) : *(lhs.sample_full_cond());
    state_rhs = prior ? *(rhs.sample_prior()) : *(rhs.sample_full_cond());
    return compute_qoi(state_lhs, state_rhs);
}

double LocalClusterMerger::compute_qoi(bayesmix::AlgorithmState::ClusterState & lhs, bayesmix::AlgorithmState::ClusterState & rhs) {
    // Create buffer
    double out = 0;
    if (lhs.has_uni_ls_state() && rhs.has_uni_ls_state()) {
        // Gaussian case
        double lmean = lhs.uni_ls_state().mean(), lvar = lhs.uni_ls_state().var();
        double rmean = rhs.uni_ls_state().mean(), rvar = rhs.uni_ls_state().var();
        // (Approx) W2 distance
        out = std::abs(lmean-rmean) + std::abs(std::sqrt(lvar)-std::sqrt(rvar));
        // Return
        return out;
    } else if (lhs.has_custom_state() && rhs.has_custom_state()) {
        // Poisson case
        spatialcmc::PoissonState unp_state; 
        lhs.custom_state().UnpackTo(&unp_state); double lrate = unp_state.rate();
        rhs.custom_state().UnpackTo(&unp_state); double rrate = unp_state.rate();
        // (Approx) W2 distance
        out = std::abs(lrate-rrate);
        // Return
        return out;
    } else {
        // No available cluster state
        throw std::runtime_error("compute_qoi() not implemented for this cluster state.");
    }
}

// Chek if two local clusters intersects each other
bool LocalClusterMerger::intersects(LocalCluster & lhs, LocalCluster & rhs) {
    // Concatenating vertices
    auto idx = lhs.get_global_cluster_idx();
    int n_row = idx.size();
    auto rhs_idx = rhs.get_global_cluster_idx();
    int n_col = rhs_idx.size();
    // Stack index vectors
    idx.insert(idx.end(), rhs_idx.begin(), rhs_idx.end());
    // Check inersection
    bool condition = (*adj_matrix)(idx, idx).topRightCorner(n_row, n_col).sum() > 0;
    // Return
    return condition;
}

// Generate spatial partitions form shards at a given iteration
std::deque<LocalCluster> LocalClusterMerger::generate_local_clusters(const size_t & iter) {
    // Create buffer for output
    std::deque<LocalCluster> local_clusters;
    // Flag for covariates
    bool requires_cov_matrix = hierarchy->is_dependent();
    // Global cluster counter is initialized
    size_t clust_counter = 0;
    for (size_t it = 0; it < num_shards; it++) {
        // Get current state of the MCMC chain
        auto curr_state = cluster_chains[it][iter];
        // Get info from shards
        auto data_in_shard = (*data)(global_numbering[it], Eigen::all);
        auto adj_matrix_in_shard = (*adj_matrix)(global_numbering[it], Eigen::all);
        auto cov_matrix_in_shard = (requires_cov_matrix) ? (*cov_matrix)(global_numbering[it], Eigen::all) : (Eigen::MatrixXd(0,0));
        // Get local cluster indexes and partition in each cluster
        std::vector<std::deque<unsigned int>> local_clust_idx(curr_state.cluster_states_size());
        std::vector<std::deque<unsigned int>> global_clust_idx(curr_state.cluster_states_size());
        for (size_t j = 0; j < curr_state.cluster_allocs_size(); j++) {
            unsigned int c = curr_state.cluster_allocs(j);
            local_clust_idx[c].push_back(j);
            global_clust_idx[c].push_back(global_numbering[it][j]);
        }
        // Build LocalCluster list
        for (size_t k = 0; k < curr_state.cluster_states_size(); k++) {
            // Construct and fill LocalCluster
            LocalCluster lc(hierarchy->deep_clone());
            lc.set_hierarchy_prior(hier_prior_file);
            lc.set_shard_label(it);
            lc.set_cluster_label(clust_counter + k);
            lc.set_global_cluster_idx(global_clust_idx[k]);
            if(requires_cov_matrix){
                lc.set_data(data_in_shard(local_clust_idx[k], Eigen::all), cov_matrix_in_shard(local_clust_idx[k], Eigen::all));
            } else {
                lc.set_data(data_in_shard(local_clust_idx[k], Eigen::all));
            }
            local_clusters.push_back(lc);
        }
        // Update global cluster counter
        clust_counter += curr_state.cluster_states_size();
    }
    // Return
    return local_clusters;
}
