#include "poisson_likelihood.h"

double PoissonLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum) const {
  return stan::math::poisson_lpmf(datum(0), state.rate);
}

int PoissonLikelihood::get_data_sum() const {
  return data_sum;
}

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