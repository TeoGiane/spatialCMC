#include "poisson_regression_likelihood.h"

double PoissonRegLikelihood::get_cov_exp_sum() const {
  // std::cout << "PoissonRegLikelihood::get_cov_exp_sum()" << std::endl;
  std::vector<double> adds(card); unsigned int k = 0;
  // std::cout << "reg_coeffs_ptr: " << (*reg_coeffs_ptr).transpose() << std::endl;
  for (auto && i : cluster_data_idx) {
    adds[k++] = log(covariates.row(i)(0)) + (covariates.row(i)(Eigen::seq(1,Eigen::last))) * reg_coeffs;
  }
  return std::exp(stan::math::log_sum_exp(adds));
}

double PoissonRegLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum,
                                          const Eigen::RowVectorXd &covariate) const {
  double offset = covariate(0); Eigen::RowVectorXd cov_vector = covariate(Eigen::seq(1,Eigen::last));
  return stan::math::poisson_log_lpmf(datum(0), log(offset) + log(state.rate) + cov_vector*reg_coeffs);
}

void PoissonRegLikelihood::update_sum_stats(const Eigen::RowVectorXd &datum,
                                            const Eigen::RowVectorXd &covariate, bool add) {
  if (add) {
    data_sum += datum(0);
  } else {
    data_sum -= datum(0);
  }
}

void PoissonRegLikelihood::clear_summary_statistics() {
  data_sum = 0;
}