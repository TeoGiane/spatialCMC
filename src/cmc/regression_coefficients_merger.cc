#include "regression_coefficients_merger.h"

RegressionCoefficientsMerger::RegressionCoefficientsMerger(const std::vector<std::vector<bayesmix::Vector>> & _reg_coeffs_chains) {
    
    // Compute num_iter, num_shards and cov_size
    num_shards = _reg_coeffs_chains.size();
    num_iter = _reg_coeffs_chains[0].size();
    cov_size = _reg_coeffs_chains[0][0].size();
    
    // Resize containers
    reg_coeff_chains.resize(num_shards);
    weights.resize(num_shards);
    
    // Pre-allocate all matrices in parallel
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < num_shards; i++) {
        reg_coeff_chains[i].resize(num_iter, cov_size);
    }
    
    // Compute weights at construction time in parallel
    weights_sum = Eigen::VectorXd::Zero(cov_size);
    
    // De-serialize regression coefficients and compute weights in a single parallel region
    // This reduces synchronization overhead and improves cache locality
    #pragma omp parallel for reduction(+ : weights_sum)
    for (size_t i = 0; i < num_shards; i++) {
        // De-serialize regression coefficients for shard i
        #pragma omp parallel for schedule(static)
        for (size_t j = 0; j < num_iter; j++) {
            reg_coeff_chains[i].row(j) = bayesmix::to_eigen(_reg_coeffs_chains[i][j]);
        }
        // Compute weights (variance inverse) for shard i as a vector of diagonal elements
        weights[i] = 1.0 / compute_sample_variance(reg_coeff_chains[i]).array();
        weights_sum += weights[i];
    }
    
    // Compute inverse of accumulated weights sum (element-wise)
    // weights = weights_sum.cwiseInverse();
    // inv_weights_sum = weights_sum.cwiseInverse();
}

bayesmix::Vector RegressionCoefficientsMerger::merge(size_t iter) {
    // Prepare buffer
    bayesmix::Vector merged_beta_proto;
    
    // Compute merged regression coefficients' vector
    Eigen::VectorXd merged_beta = Eigen::VectorXd::Zero(cov_size);
    for (size_t i = 0; i < num_shards; i++) {
        Eigen::VectorXd curr_beta = reg_coeff_chains[i].row(iter);
        merged_beta += weights[i].cwiseProduct(curr_beta);
    }
    merged_beta = merged_beta.array() / weights_sum.array();

    bayesmix::to_proto(merged_beta, &merged_beta_proto);
    return merged_beta_proto;
}

Eigen::VectorXd RegressionCoefficientsMerger::compute_sample_variance(const Eigen::MatrixXd & reg_coeff_chain) {
    // Get dimensions from Eigen matrix
    size_t cov_size = reg_coeff_chain.cols();
    // Compute marginal variances for each dimension
    Eigen::VectorXd marginal_variances(cov_size);
    for (size_t j = 0; j < cov_size; j++) {
        marginal_variances(j) = stan::math::variance(reg_coeff_chain.col(j));
    }
    // Return vector of marginal variances (not diagonal matrix)
    return marginal_variances;
}
