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

template<class T>
class LocalClusterMerger {
 private:

  // Class members
  Eigen::MatrixXd const *data = nullptr;
  Eigen::MatrixXd cov_matrix;
  // Eigen::MatrixXd const *cov_matrix = nullptr;
  Eigen::MatrixXd const *adj_matrix = nullptr;
  std::vector<Eigen::VectorXd> reg_coeffs_chain;
  std::vector<std::vector<bayesmix::AlgorithmState>> cluster_chains;
  std::vector<std::deque<int>> global_numbering;
  std::string hier_prior_file;

  unsigned int global_card = 0;
  unsigned int num_shards = 0;
  unsigned int num_iter = 0;

 public:

	// Constructor & Destructor
  LocalClusterMerger(const std::vector<std::vector<bayesmix::AlgorithmState>> & _cluster_chains);
	~LocalClusterMerger() = default;

	// Setters
  void set_data(Eigen::MatrixXd const * _data_ptr) { data = _data_ptr; global_card = data->rows(); };
  void set_cov_matrix(const Eigen::MatrixXd & _cov_matrix) { cov_matrix = _cov_matrix; };
  void set_adj_matrix(Eigen::MatrixXd const * _adj_matrix_ptr) { adj_matrix = _adj_matrix_ptr; };
  void set_reg_coeffs(size_t iter, const Eigen::VectorXd & _reg_coeffs) {reg_coeffs_chain[iter] = _reg_coeffs;};
  void set_global_numbering(const std::vector<std::deque<int>> & _global_numbering) { global_numbering = _global_numbering; };
  void set_seed(unsigned int _seed) { auto &rng = bayesmix::Rng::Instance().get(); rng.seed(_seed); };
  void set_hierarchy_prior(const std::string & _hier_prior_file){ hier_prior_file = _hier_prior_file; };

	// Getters
	unsigned int get_num_iter() const { return num_iter; };

	// Merge local clusters at a given iteration
	bayesmix::AlgorithmState merge(size_t iter);

 private:

	// Compute bayes factor for two spatial partitions merging candidates (+ aux functions)
	double compute_logBF(LocalCluster<T> & lhs, LocalCluster<T> & rhs);
	double sample_qoi(LocalCluster<T> & lhs, LocalCluster<T> & rhs, bool prior);
	double compute_qoi(bayesmix::AlgorithmState::ClusterState & lhs, bayesmix::AlgorithmState::ClusterState & rhs);

  // Chek if two local clusters intersects each other
	bool intersects(LocalCluster<T> & lhs, LocalCluster<T> & rhs);

	// Generate spatial partitions form shards at a given iteration
	std::deque<LocalCluster<T>> generate_local_clusters(const size_t & iter);

};

// Template methods implementations
#include "local_cluster_merger_imp.h"

#endif // SPATIALCMC_CMC_LOCAL_CLUSTER_MERGER_H


/* OLD STUFF */

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