#ifndef BAYESMIX_HIERARCHIES_LIKELIHOODS_EMPTY_LIKELIHOOD_H_
#define BAYESMIX_HIERARCHIES_LIKELIHOODS_EMPTY_LIKELIHOOD_H_

#include <google/protobuf/stubs/casts.h>

#include <memory>
#include <stan/math/rev.hpp>
#include <vector>

#include "algorithm_state.pb.h"
#include "base_likelihood.h"
#include "states/includes.h"

/**
 * A empty likelihood, using the `State::UniLS` state.
 */

class EmptyLikelihood
    : public BaseLikelihood<EmptyLikelihood, State::Empty> {
 public:
  EmptyLikelihood() = default;
  ~EmptyLikelihood() = default;
  bool is_multivariate() const override { return false; };
  bool is_dependent() const override { return false; };
  void clear_summary_statistics() override;

  // template <typename T>
  // T cluster_lpdf_from_unconstrained(const Eigen::Matrix<T, Eigen::Dynamic, 1> &unconstrained_params) const {
  //   return 0;
  // }

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum, bool add) override;
};

#endif  // BAYESMIX_HIERARCHIES_LIKELIHOODS_UNI_NORM_LIKELIHOOD_H_
