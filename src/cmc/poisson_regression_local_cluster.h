#ifndef SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_H
#define SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_H

// STL
#include <iostream>
#include <deque>

// bayesmix
#include <src/includes.h>

// stan-math
#include <stan/math.hpp>

// spatialcmc
#include "hierarchies/poisson_regression_hierarchy.h"
#include "local_cluster.h"
#include "poisson_state.pb.h"
#include "spatialcmc_utils.h"


class PoissonRegLocalCluster : public LocalCluster {
  private:
    // Extra Class members
    Eigen::VectorXd reg_coeffs;

  public:
    // Constructor & Destructor
    PoissonRegLocalCluster(std::shared_ptr<AbstractHierarchy> _hierarchy): LocalCluster(_hierarchy) {};
    ~PoissonRegLocalCluster() = default;
    // Setters
    void set_data(const Eigen::MatrixXd & _data, const Eigen::MatrixXd & _covariates = Eigen::MatrixXd(0,0)) override;
    void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs);  
    // Getters
    Eigen::VectorXd get_reg_coeffs() const { return reg_coeffs; };
    // Merge current spatial partition with the input spatial partition
    void merge_with(const PoissonRegLocalCluster & rhs);
    // Print utility for debug
    void print(void) override;
};

#endif // SPATIALCMC_CMC_POISSON_REGRESSION_LOCAL_CLUSTER_H