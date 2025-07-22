#ifndef SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_LIKELIHOOD_H_
#define SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_LIKELIHOOD_H_

#include <memory>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>

#include "hierarchies/likelihoods/states/poisson_state.h"
#include "src/hierarchies/likelihoods/base_likelihood.h"

class PoissonLikelihood : public BaseLikelihood<PoissonLikelihood, State::Poisson> {
 public:
  PoissonLikelihood() = default;
  ~PoissonLikelihood() = default;
  bool is_multivariate() const override { return false; };
  bool is_dependent() const override { return true; };
  void clear_summary_statistics() override;

  // Getters and Setters
  int get_data_sum() const;
  double get_offset_sum() const;

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum, const Eigen::RowVectorXd &covariate) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum, const Eigen::RowVectorXd &covariate, bool add) override;
  
  // Class Members
  int data_sum = 0;
  double offset_sum = 0;
};

#endif // SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_LIKELIHOOD_H_