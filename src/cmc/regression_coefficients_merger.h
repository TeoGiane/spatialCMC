#ifndef SPATIALCMC_CMC_REGRESSION_COEFFICIENTS_MERGER_H_
#define SPATIALCMC_CMC_REGRESSION_COEFFICIENTS_MERGER_H_

// STL
#include <vector>
#include <stan/math.hpp>
#include <Eigen/Dense>

// bayesmix
#include <src/includes.h>

#pragma omp declare reduction (+: Eigen::MatrixXd: omp_out=omp_out+omp_in) \
        initializer(omp_priv=Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()))

class RegressionCoefficientsMerger {
  private:
    std::vector<std::vector<bayesmix::Vector>> reg_coeffs_chains;
    std::vector<Eigen::MatrixXd> weights;
    Eigen::MatrixXd inv_weights_sum;
    unsigned int num_shards;
    unsigned int num_iter;
    unsigned int cov_size;

  public:
    // Constructor and Destructor
    RegressionCoefficientsMerger(const std::vector<std::vector<bayesmix::Vector>> & _reg_coeffs_chains);
    ~RegressionCoefficientsMerger() = default;
    // Getters and Setters
    unsigned int get_num_iter() const {return num_iter; };
    // Merge at fixed iteration
    bayesmix::Vector merge(size_t iter);
    private:
    // Compute sample variance
    Eigen::MatrixXd compute_sample_variance(const std::vector<bayesmix::Vector> & reg_coeff_chain);
};

#endif // SPATIALCMC_CMC_REGRESSION_COEFFICIENTS_MERGER_H_