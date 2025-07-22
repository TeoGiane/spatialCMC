#include "poisson_regression_local_cluster.h"

void PoissonRegLocalCluster::set_data(const Eigen::MatrixXd & _data, const Eigen::MatrixXd & _covariates) {
    LocalCluster::set_data(_data, _covariates);
    std::static_pointer_cast<PoissonRegHierarchy>(hierarchy)->set_covariates(covariates);
}

void PoissonRegLocalCluster::set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) {
    reg_coeffs = _reg_coeffs;
    std::static_pointer_cast<PoissonRegHierarchy>(hierarchy)->set_reg_coeffs(reg_coeffs);
}

void PoissonRegLocalCluster::merge_with(const PoissonRegLocalCluster & rhs) {
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
        std::static_pointer_cast<PoissonRegHierarchy>(hierarchy)->set_covariates(covariates);
    }
	// Merge global cluster indexes
	global_clust_idx.insert(global_clust_idx.end(), rhs.global_clust_idx.begin(), rhs.global_clust_idx.end());
}

void PoissonRegLocalCluster::print() {
    std::cout << "Poisson Regression Local Cluster: " << cluster_label << " coming from shard " << shard_label << std::endl;
    std::cout << "Data in this partition: " << data.transpose() << std::endl;
    std::cout << "Data size: " << data.size() << std::endl;
    std::cout << "Global data indexes: "; for (auto && elem : global_clust_idx) { std::cout << elem << " "; }; std::cout << std::endl;
    std::cout << "Global data index size: " << global_clust_idx.size() << std::endl;
    std::cout << "Covariates:\n" << covariates << std::endl;
    std::cout << "Reg Coeffs:\n" << reg_coeffs.transpose() << std::endl;
    std::cout << "Unique Value:\n " << sample_full_cond()->DebugString() << std::endl;
    std::cout << std::endl;
}