#include "poisson_regression_likelihood.h"

double PoissonRegLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum,
                                          const Eigen::RowVectorXd &covariate) const {
  double offset = covariate(0); Eigen::RowVectorXd cov_vector = covariate(Eigen::seq(1,Eigen::last));
  return stan::math::poisson_log_lpmf(datum(0), log(offset)+log(state.rate)+cov_vector*(*reg_coeffs_ptr));
}

void PoissonRegLikelihood::update_sum_stats(const Eigen::RowVectorXd &datum,
                                            const Eigen::RowVectorXd &covariate, bool add) {
  // double offset = covariate(0); Eigen::RowVectorXd cov_vector = covariate(Eigen::seq(1,Eigen::last));
  if (add) {
    data_sum += datum(0);
  } else {
    data_sum -= datum(0);
  }
}

void PoissonRegLikelihood::clear_summary_statistics() {
  data_sum = 0;
  // cov_exp_sum = 0;
  // cluster_lpdf = 0;
}