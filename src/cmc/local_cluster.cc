#include "local_cluster.h"

void LocalCluster::set_shard_label(const size_t & _shard_label) {
    shard_label = _shard_label;
}

void LocalCluster::set_cluster_label(const size_t &_cluster_label) {
    cluster_label = _cluster_label;
}

void LocalCluster::set_data(const Eigen::MatrixXd & _data, const Eigen::MatrixXd & _covariates) {
    data = _data;
    covariates = _covariates;
    if(hierarchy->is_dependent() && covariates.size() == 0){
        throw std::runtime_error("Covariance Matrix not yet provided.");
    } else {
        covariates_getter cov_getter(covariates);
        for (int i = 0; i < data.rows(); i++) {
        hierarchy->add_datum(i, data.row(i), false, cov_getter(i));
        }
    }
}

void LocalCluster::set_global_cluster_idx(const std::deque<unsigned int> &_global_clust_idx) {
    global_clust_idx = _global_clust_idx;
}

void LocalCluster::set_hierarchy_prior(const std::string & _hier_prior_file) {
    bayesmix::read_proto_from_file(_hier_prior_file, hierarchy->get_mutable_prior());
    hierarchy->initialize();
}

void LocalCluster::merge_with(const LocalCluster & rhs) {
    // Check if we need covariates
    bool requires_covariates = hierarchy->is_dependent();
    // New cluster name
    shard_label = std::min(shard_label, rhs.get_shard_label());
    cluster_label = std::min(cluster_label, rhs.get_cluster_label());
    // Update data in hierarchy here
    for(size_t i = 0; i < rhs.data.rows(); i++){
        if (requires_covariates) {
            hierarchy->add_datum(data.rows() + i, rhs.data.row(i), false, rhs.covariates.row(i));
        } else {
            hierarchy->add_datum(data.rows() + i, rhs.data.row(i), false);
        }
    }
	// Merge data
	data.conservativeResize(data.rows() + rhs.data.rows(), data.cols());
	data.bottomRightCorner(rhs.data.rows(), data.cols()) = rhs.data;
    // Merge covariates (if present)
    if(requires_covariates){
        covariates.conservativeResize(covariates.rows() + rhs.covariates.rows(), covariates.cols());
        covariates.bottomRightCorner(rhs.covariates.rows(), covariates.cols()) = rhs.covariates;
    }
	// Merge global cluster indexes
	global_clust_idx.insert(global_clust_idx.end(), rhs.global_clust_idx.begin(), rhs.global_clust_idx.end());
}

std::shared_ptr<LocalCluster::ClustState> LocalCluster::sample_prior(void) {
    hierarchy->sample_prior();
    return hierarchy->get_state_proto();
}
  
std::shared_ptr<LocalCluster::ClustState> LocalCluster::sample_full_cond(void){
    hierarchy->sample_full_cond();
    return hierarchy->get_state_proto();
}

// template<typename T>
void LocalCluster::print() {
    std::cout << "Local Cluster: " << cluster_label << " coming from shard " << shard_label << std::endl;
    std::cout << "Data in this partition: " << data.transpose() << std::endl;
    std::cout << "Data size: " << data.size() << std::endl;
    std::cout << "Global data indexes: "; for (auto && elem : global_clust_idx) { std::cout << elem << " "; }; std::cout << std::endl;
    std::cout << "Global data index size: " << global_clust_idx.size() << std::endl;
    std::cout << "Covariates:\n" << covariates << std::endl;
    std::cout << "Unique Value:\n " << sample_full_cond()->DebugString() << std::endl;
    std::cout << std::endl;
}