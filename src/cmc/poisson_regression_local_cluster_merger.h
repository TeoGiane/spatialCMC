#ifndef SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_MERGER_H
#define SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_MERGER_H

// STL
#include <deque>

#include <stan/math.hpp>
#include <Eigen/Dense>

// bayesmix
#include "src/includes.h"

// spatialcmc
#include "cmc/local_cluster_merger.h"
#include "cmc/poisson_regression_local_cluster.h"
#include "cmc/spatialcmc_utils.h"

#include "poisson_regression_algorithm_state.pb.h"


class PoissonRegLocalClusterMerger : public LocalClusterMerger {
  private:
    // Extra Class members
    std::vector<Eigen::VectorXd> reg_coeffs_chain;

  public:
    // Constructor & Destructor
    PoissonRegLocalClusterMerger(const std::vector<std::vector<bayesmix::AlgorithmState>> & _cluster_chains);
	  ~PoissonRegLocalClusterMerger() = default;
    // Setters
    void set_reg_coeffs(size_t iter, const Eigen::VectorXd & _reg_coeffs) {reg_coeffs_chain[iter] = _reg_coeffs;};
	  // Merge local clusters at a given iteration
    bayesmix::AlgorithmState merge(size_t iter);

  private:
	// Generate spatial partitions form shards at a given iteration
	std::deque<PoissonRegLocalCluster> generate_local_clusters(const size_t & iter);
};

#endif // SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_MERGER_H
