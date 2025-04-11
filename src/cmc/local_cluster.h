#ifndef SPATIALCMC_CMC_LOCAL_CLUSTER_H
#define SPATIALCMC_CMC_LOCAL_CLUSTER_H

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

template<typename T>
class LocalCluster {

 private:
  // Type aliases
  using HierHypers = bayesmix::AlgorithmState::HierarchyHypers;
  using ClustState = bayesmix::AlgorithmState::ClusterState;
  // Class members
  Eigen::MatrixXd data;
  Eigen::MatrixXd cov_matrix;
  Eigen::VectorXd reg_coeffs;
  // std::shared_ptr<Eigen::VectorXd> reg_coeffs = nullptr;
  size_t shard_label;
  size_t cluster_label;
  std::deque<unsigned int> global_clust_idx;
  std::shared_ptr<T> hierarchy = std::make_shared<T>();

 public:
  // Constructor & Destructor
  LocalCluster() = default;
  ~LocalCluster() = default;
  // Setters
  void set_shard_label(const size_t & _shard_label);
  void set_cluster_label(const size_t & _cluster_label);
  // void set_data(const Eigen::MatrixXd & _data, const Eigen::MatrixXd & _cov_matrix = Eigen::MatrixXd(0,0));
  void set_data(const Eigen::MatrixXd & _data) {
    data = _data;
    if(hierarchy->is_dependent() && cov_matrix.size() == 0){
      throw std::runtime_error("Covariance Matrix not yet provided.");
    } else {
      covariates_getter cov_getter(cov_matrix);
      for (int i = 0; i < data.rows(); i++) {
        hierarchy->add_datum(i, data.row(i), false, cov_getter(i));
      }
    }
  };
  void set_cov_matrix(const Eigen::MatrixXd & _cov_matrix) {
    // std::cout << "LocalCluster<T>::set_cov_matrix()" << std::endl;
    cov_matrix = _cov_matrix;
    // std::cout << "Cov matrix:\n" << cov_matrix << std::endl;
    hierarchy->set_covariates(_cov_matrix);
    // std::cout << "Cov matrix set in hierarchy" << std::endl;
    // std::cout << hierarchy->get_covariates() << std::endl;
  };
  void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { reg_coeffs = _reg_coeffs; hierarchy->set_reg_coeffs(reg_coeffs); };
  void set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx);
  void set_hierarchy_prior(const std::string & _hier_prior_file);
  // void initialize() {hierarchy->initialize();};
    // QUI NON VALE PER TUTTE LE GERARCHIE!!!
    // hierarchy->set_covariates(&cov_matrix);
    // hierarchy->set_reg_coeffs(&reg_coeffs);
    // covariates_getter cov_getter(cov_matrix);
    // for (int i = 0; i < data.rows(); i++) {
    //   hierarchy->add_datum(i, data.row(i), false, cov_getter(i));
    // }
    // std::cout << "Card: " << hierarchy->get_card() << std::endl;
  // }
  // Getters
  size_t get_shard_label() const { return shard_label; };
  Eigen::MatrixXd get_cov_matrix() const { return hierarchy->get_covariates(); };
  size_t get_cluster_label() const { return cluster_label; };
  std::deque<unsigned int> get_global_cluster_idx() const { return global_clust_idx; };
  unsigned int get_card() const { return hierarchy->get_card(); };
  // Merge current spatial partition with the input spatial partition
  void merge_with(const LocalCluster<T> & rhs);
  // Returns a sample from the prior distribution of the hierarchy
  std::shared_ptr<ClustState> sample_prior(){
    // std::cout << "LocalCluster<T>::sample_prior()" << std::endl;
    hierarchy->sample_prior();
    // hierarchy->get_state_proto()->PrintDebugString();
    return hierarchy->get_state_proto();
  };
  // Returns a sample from the full conditional distribution of the hierarchy
  std::shared_ptr<ClustState> sample_full_cond(){
    // std::cout << "LocalCluster<T>::sample_full_cond()" << std::endl;
    hierarchy->sample_full_cond();
    // hierarchy->get_state_proto()->PrintDebugString();
    // std::cout << "Here now" << std::endl;
    return hierarchy->get_state_proto();
  };
  // Print utility for debug
  void print();
  // Check if the underlying hierarchy is dependent or not
  bool has_dependent_hierarchy() const { return hierarchy->is_dependent(); }
};

// Template methods implementations
#include "local_cluster_imp.h"

#endif // SPATIALCMC_CMC_LOCAL_CLUSTER_H