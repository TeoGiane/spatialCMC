#include "regression_coefficients_merger.h"

RegressionCoefficientsMerger::RegressionCoefficientsMerger(const std::vector<std::vector<bayesmix::Vector>> & _reg_coeffs_chains)
: reg_coeffs_chains(_reg_coeffs_chains) {
    // Comput num_iterm num_shards and cov_size
    num_shards = reg_coeffs_chains.size();
    num_iter = reg_coeffs_chains[0].size();
    cov_size = reg_coeffs_chains[0][0].size();
    // Compute weights at construction time
    weights.resize(reg_coeffs_chains.size());
    inv_weights_sum.resize(cov_size, cov_size);
    // Compute weights at construction time
    Eigen::MatrixXd weights_sum = Eigen::MatrixXd::Zero(cov_size, cov_size);
    #pragma omp parallel for reduction(+ : weights_sum)
    for (size_t i = 0; i < weights.size(); i++){
    weights[i] = compute_sample_variance(reg_coeffs_chains[i]).inverse();
    weights_sum += weights[i];
    }
    inv_weights_sum = weights_sum.inverse();
}

bayesmix::Vector RegressionCoefficientsMerger::merge(size_t iter) {
    // Prepare buffer
    bayesmix::Vector merged_beta_proto;
    // Compute merged regression coefficients' vector
    Eigen::VectorXd merged_beta(cov_size);
    for (size_t i = 0; i < reg_coeffs_chains.size(); i++) {
        auto curr_beta = bayesmix::to_eigen(reg_coeffs_chains[i][iter]);
        merged_beta += (inv_weights_sum * weights[i]) * curr_beta;
    }
    bayesmix::to_proto(merged_beta, &merged_beta_proto);
    return merged_beta_proto;
}

Eigen::MatrixXd RegressionCoefficientsMerger::compute_sample_variance(const std::vector<bayesmix::Vector> & reg_coeff_chain) {
    Eigen::VectorXd M_curr = bayesmix::to_eigen(reg_coeff_chain[0]);
    Eigen::VectorXd M_old = M_curr;
    Eigen::MatrixXd S_old = Eigen::MatrixXd::Zero(cov_size, cov_size);
    Eigen::MatrixXd S_curr = Eigen::MatrixXd::Zero(cov_size, cov_size);
    for (size_t i = 1; i < num_iter; i++) {
        // Running variance
        auto curr_beta = bayesmix::to_eigen(reg_coeff_chain[i]);
        M_curr = M_old + (curr_beta - M_old) / (i+1);
        S_curr = S_old + (curr_beta - M_old) * (curr_beta - M_curr).transpose();
        // Update
        M_old = M_curr;
        S_old = S_curr;
    }
    // Return
    return S_curr / (num_iter-1);
}
