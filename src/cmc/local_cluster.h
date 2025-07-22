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


class LocalCluster {

 protected:
  // Type aliases
  using HierHypers = bayesmix::AlgorithmState::HierarchyHypers;
  using ClustState = bayesmix::AlgorithmState::ClusterState;
  
  // Class members
  Eigen::MatrixXd data;
  Eigen::MatrixXd covariates;
  // Eigen::VectorXd reg_coeffs;
  std::shared_ptr<AbstractHierarchy> hierarchy;

  // Labels
  size_t shard_label;
  size_t cluster_label;
  std::deque<unsigned int> global_clust_idx;

 public:
  
  // Constructor & Destructor
  LocalCluster(std::shared_ptr<AbstractHierarchy> _hierarchy): hierarchy(_hierarchy) {};
  ~LocalCluster() = default;
  
  // Setters
  void set_shard_label(const size_t & _shard_label);
  void set_cluster_label(const size_t & _cluster_label);
  virtual void set_data(const Eigen::MatrixXd & _data, const Eigen::MatrixXd & _covariates = Eigen::MatrixXd(0,0));
  void set_global_cluster_idx(const std::deque<unsigned int> & _global_clust_idx);
  void set_hierarchy_prior(const std::string & _hier_prior_file);
  
  // Getters
  size_t get_shard_label() const { return shard_label; };
  Eigen::MatrixXd get_covariates() const { return covariates; };
  size_t get_cluster_label() const { return cluster_label; };
  std::deque<unsigned int> get_global_cluster_idx() const { return global_clust_idx; };
  unsigned int get_card() const { return hierarchy->get_card(); };
  
  // Merge current spatial partition with the input spatial partition
  virtual void merge_with(const LocalCluster & rhs);
  
  // Returns a sample from the prior distribution of the hierarchy
  std::shared_ptr<ClustState> sample_prior(void);
  
  // Returns a sample from the full conditional distribution of the hierarchy
  std::shared_ptr<ClustState> sample_full_cond(void);
  
  // Print utility for debug
  virtual void print();

};

#endif // SPATIALCMC_CMC_LOCAL_CLUSTER_H