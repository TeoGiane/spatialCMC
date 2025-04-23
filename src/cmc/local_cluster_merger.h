#ifndef SPATIALCMC_CMC_LOCAL_CLUSTER_MERGER_H
#define SPATIALCMC_CMC_LOCAL_CLUSTER_MERGER_H

// STL
#include <deque>

#include <stan/math.hpp>
#include <Eigen/Dense>

// bayesmix
#include "src/includes.h"

// spatialcmc
#include "cmc/local_cluster.h"
#include "cmc/spatialcmc_utils.h"

class LocalClusterMerger {
  protected:
    // Class members
    Eigen::MatrixXd const *data = nullptr;
    Eigen::MatrixXd const *cov_matrix = nullptr;
    Eigen::MatrixXd const *adj_matrix = nullptr;
    std::vector<std::vector<bayesmix::AlgorithmState>> cluster_chains;
    std::vector<std::deque<int>> global_numbering;
    std::shared_ptr<AbstractHierarchy> hierarchy;
    std::string hier_prior_file;
    // Indexes
    unsigned int global_card = 0;
    unsigned int num_shards = 0;
    unsigned int num_iter = 0;

  public:
    // Constructor & Destructor
    LocalClusterMerger(const std::vector<std::vector<bayesmix::AlgorithmState>> & _cluster_chains);
    ~LocalClusterMerger() = default;
    // Setters
    void set_hierarchy(std::shared_ptr<AbstractHierarchy> _hierarchy) { hierarchy = _hierarchy; };
    void set_data(Eigen::MatrixXd const * _data_ptr) { data = _data_ptr; global_card = data->rows(); };
    void set_cov_matrix(Eigen::MatrixXd const * _cov_matrix_ptr) { cov_matrix = _cov_matrix_ptr; };
    void set_adj_matrix(Eigen::MatrixXd const * _adj_matrix_ptr) { adj_matrix = _adj_matrix_ptr; };
    // void set_reg_coeffs(size_t iter, const Eigen::VectorXd & _reg_coeffs) {reg_coeffs_chain[iter] = _reg_coeffs;};
    void set_global_numbering(const std::vector<std::deque<int>> & _global_numbering) { global_numbering = _global_numbering; };
    void set_seed(unsigned int _seed) { auto &rng = bayesmix::Rng::Instance().get(); rng.seed(_seed); };
    void set_hierarchy_prior(const std::string & _hier_prior_file){ hier_prior_file = _hier_prior_file; };
    // Getters
    unsigned int get_num_iter() const { return num_iter; };
    // Merge local clusters at a given iteration
    bayesmix::AlgorithmState merge(size_t iter);

  protected:
    // Compute bayes factor for two spatial partitions merging candidates (+ aux functions)
    double compute_logBF(LocalCluster & lhs, LocalCluster & rhs);
    double sample_qoi(LocalCluster & lhs, LocalCluster & rhs, bool prior);
    double compute_qoi(bayesmix::AlgorithmState::ClusterState & lhs, bayesmix::AlgorithmState::ClusterState & rhs);
    // Chek if two local clusters intersects each other
    bool intersects(LocalCluster & lhs, LocalCluster & rhs);
    // Generate spatial partitions form shards at a given iteration
    std::deque<LocalCluster> generate_local_clusters(const size_t & iter);

};

#endif // SPATIALCMC_CMC_LOCAL_CLUSTER_MERGER_H
