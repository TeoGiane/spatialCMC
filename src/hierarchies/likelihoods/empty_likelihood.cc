#include "empty_likelihood.h"

double EmptyLikelihood::compute_lpdf(const Eigen::RowVectorXd &datum) const {
  return 0;
}

void EmptyLikelihood::update_sum_stats(const Eigen::RowVectorXd &datum, bool add) {
  return;
}

void EmptyLikelihood::clear_summary_statistics() {
  return;
}
