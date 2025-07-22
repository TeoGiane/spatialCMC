#include "poisson_likelihood.h"

double PoissonLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum, const Eigen::RowVectorXd &covariate) const {
  return stan::math::poisson_lpmf(datum(0), covariate(0) * state.rate);
}

int PoissonLikelihood::get_data_sum() const {
  return data_sum;
}

double PoissonLikelihood::get_offset_sum() const {
  return offset_sum;
}

void PoissonLikelihood::update_sum_stats(const Eigen::RowVectorXd &datum, const Eigen::RowVectorXd &covariate, bool add) {
  if (add) {
    data_sum += datum(0);
    offset_sum += covariate(0);
  } else {
    data_sum -= datum(0);
    offset_sum -= covariate(0);
  }
}

void PoissonLikelihood::clear_summary_statistics() {
  data_sum = 0;
  offset_sum = 0;
}