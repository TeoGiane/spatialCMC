#include "local_cluster.h"

template <typename T>
void LocalCluster<T>::set_shard_label(const size_t & _shard_label) {
  shard_label = _shard_label;
}

template <typename T>
void LocalCluster<T>::set_cluster_label(const size_t &_cluster_label) {
  cluster_label = _cluster_label;
}

// template <typename T>
// void LocalCluster<T>::set_data(const Eigen::MatrixXd &_data, const Eigen::MatrixXd & _cov_matrix) {
//   // std::cout << "LocalCluster<T>::set_data()" << std::endl;
//   // std::cout << "reg_coeffs: " << hierarchy->get_reg_coeffs().transpose() << std::endl;
//   data = _data;
//   // Check cov_matrix dimensions
//   if((_cov_matrix.rows() > 0) && (_cov_matrix.rows() == data.rows())) {
//     // std::cout << "Here there is a cov matrix" << std::endl;
//     cov_matrix = _cov_matrix;
//   } else {
//     throw std::runtime_error("Data and covariance matrix dimension mis-match");
//   }
//   // Add data and covariates into hierarchy
//   if (hierarchy->is_dependent()){
// 		for (size_t i = 0; i < data.rows(); i++) {
//       // std::cout << "adding datum" << std::endl;
// 			hierarchy->add_datum(i, data.row(i), false);
// 		}
//   } else {
//     for (size_t i = 0; i < data.rows(); i++) {
//       // std::cout << "adding datum" << std::endl;
// 			hierarchy->add_datum(i, data.row(i), false, cov_matrix.row(i));
// 		}
//   }
// }

template <typename T>
void LocalCluster<T>::set_global_cluster_idx(const std::deque<unsigned int> &_global_clust_idx) {
  global_clust_idx = _global_clust_idx;
}

template<typename T>
void LocalCluster<T>::set_hierarchy_prior(const std::string & _hier_prior_file) {
  bayesmix::read_proto_from_file(_hier_prior_file, hierarchy->get_mutable_prior());
  hierarchy->initialize();
}

template <typename T>
void LocalCluster<T>::merge_with(const LocalCluster<T> & rhs) {
  // Check if we need covariates
  bool requires_covariates = has_dependent_hierarchy();
  // New cluster name
	shard_label = std::min(shard_label, rhs.get_shard_label());
	cluster_label = std::min(cluster_label, rhs.get_cluster_label());
	// Update data in hierarchy here
	for(size_t i = 0; i < rhs.data.rows(); i++){
    if (requires_covariates) {
      hierarchy->add_datum(data.rows() + i, rhs.data.row(i), false, rhs.cov_matrix.row(i));
    } else {
      hierarchy->add_datum(data.rows() + i, rhs.data.row(i), false);
    }
	}
	// Merge data
	data.conservativeResize(data.rows() + rhs.data.rows(), data.cols());
	data.bottomRightCorner(rhs.data.rows(), data.cols()) = rhs.data;
  // Merge covariates (if present)
  if(requires_covariates){
    cov_matrix.conservativeResize(cov_matrix.rows() + rhs.cov_matrix.rows(), cov_matrix.cols());
    cov_matrix.bottomRightCorner(rhs.cov_matrix.rows(), cov_matrix.cols()) = rhs.cov_matrix;
    hierarchy->set_covariates(cov_matrix);
  }
	// Merge global cluster indexes
	global_clust_idx.insert(global_clust_idx.end(), rhs.global_clust_idx.begin(), rhs.global_clust_idx.end());
}

// template <typename T>
// std::shared_ptr<typename LocalCluster<T>::ClustState> LocalCluster<T>::sample_prior() {
//   hierarchy->sample_prior();
//   return hierarchy->get_state_proto();
// }

// template <typename T>
// std::shared_ptr<typename LocalCluster<T>::ClustState> LocalCluster<T>::sample_full_cond() {
//   std::cout << "LocalCluster<T>::sample_full_cond()" << std::endl;
//   std::static_pointer_cast<T>(hierarchy)->sample_full_cond();
//   std::cout << "Here now" << std::endl;
//   return hierarchy->get_state_proto();
// }

template<typename T>
void LocalCluster<T>::print() {
  std::cout << "Local Cluster: " << cluster_label << " coming from shard " << shard_label << std::endl;
	std::cout << "Data in this partition: " << data.transpose() << std::endl;
	std::cout << "Data size: " << data.size() << std::endl;
	std::cout << "Global data indexes: "; for (auto && elem : global_clust_idx) { std::cout << elem << " "; }; std::cout << std::endl;
	std::cout << "Global data index size: " << global_clust_idx.size() << std::endl;
  std::cout << "Cov Matrix:\n" << hierarchy->get_covariates() << std::endl;
  std::cout << "Reg. Coeffs: " << hierarchy->get_reg_coeffs().transpose() << std::endl;
	std::cout << "Unique Value:\n " << sample_full_cond()->DebugString() << std::endl;
	std::cout << std::endl;
}