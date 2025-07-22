#ifndef SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_REGRESSION_LIKELIHOOD_H_
#define SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_REGRESSION_LIKELIHOOD_H_

#include <memory>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <google/protobuf/stubs/casts.h>

#include "hierarchies/likelihoods/states/poisson_state.h"
#include "src/hierarchies/likelihoods/base_likelihood.h"

class PoissonRegLikelihood : public BaseLikelihood<PoissonRegLikelihood, State::Poisson> {
 public:
  PoissonRegLikelihood() = default;
  ~PoissonRegLikelihood() = default;
  bool is_multivariate() const override { return false; };
  bool is_dependent() const override { return true; };
  void clear_summary_statistics() override;

  // Setters
  void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { reg_coeffs = _reg_coeffs; };
  void set_covariates(const Eigen::MatrixXd & _covariates) { covariates = _covariates; };
  // Getters
  int get_data_sum(void) const { return data_sum; };
  double get_cov_exp_sum(void) const;
  Eigen::VectorXd get_reg_coeffs() const { return reg_coeffs; };
  Eigen::MatrixXd get_covariates() const { return covariates; };

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum,
                      const Eigen::RowVectorXd &covariate) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum,
                        const Eigen::RowVectorXd &covariate,
                        bool add) override;

  // Common regression coefficient
  Eigen::VectorXd reg_coeffs;
  Eigen::MatrixXd covariates;
  // Summary statistics
  int data_sum = 0;
};

#endif // SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_REGRESSION_LIKELIHOOD_H_