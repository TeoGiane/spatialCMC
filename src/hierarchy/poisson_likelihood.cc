#include "poisson_likelihood.h"

double PoissonLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum) const {
  return stan::math::poisson_log_lpmf(datum(0), log(state.gamma));
}

// double PoissonLikelihood::get_exp_cov_sum() const {
//   std::cout << "PoissonLikelihood::get_exp_cov_sum()" << std::endl;
//   // double out = 0.0;
//   std::cout << "card: " << card << std::endl;
//   std::vector<double> adds(card);
//   unsigned int k = 0;
//   for(int i : cluster_data_idx) {
//     // std::cout << "covariates: " << covariates->row(i) << std::endl;
//     // std::cout << "reg_coeffs: " << (*reg_coeffs).transpose() << std::endl;
//     adds[k++] = covariates->row(i) * (*reg_coeffs);
//     // k++;
//     // adds.push_back(covariates->row(i) * (*reg_coeffs));
//     // out += covariates->row(i) * (*reg_coeffs);
//   }
//   return std::exp(stan::math::log_sum_exp(adds));
// }

int PoissonLikelihood::get_data_sum() const {
  // std::cout << "get_data_sum()" << std::endl;
  // int sum = 0;
  // std::cout << "data:\n" << *dataset_ptr << std::endl;
  // std::vector<int> idx(cluster_data_idx.begin(), cluster_data_idx.end());
  // for(int i : cluster_data_idx) {
  //   sum += (*dataset_ptr)(i,0);
  // }
  // std::cout << "data_in_clust: " << (*dataset_ptr)(idx, 0).transpose() << std::endl;
  // std::cout << "sum: " << sum << std::endl;
  return data_sum;
};

void PoissonLikelihood::update_sum_stats(const Eigen::RowVectorXd &datum, bool add) {
  if (add) {
    data_sum += datum(0);
  } else {
    data_sum -= datum(0);
  }
}

void PoissonLikelihood::clear_summary_statistics() {
  data_sum = 0;
}