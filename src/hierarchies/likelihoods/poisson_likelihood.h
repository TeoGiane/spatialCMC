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
  bool is_dependent() const override { return false; };
  void clear_summary_statistics() override;

  // Getters and Setters
  // double get_gamma() const { return state.gamma; };
  int get_data_sum() const;
  // double get_exp_cov_sum() const;

  // std::shared_ptr<Eigen::VectorXd> get_reg_coeffs() const { return reg_coeffs; };
  // std::shared_ptr<const Eigen::MatrixXd> get_covariates() const { return covariates; };
  // void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { *reg_coeffs = _reg_coeffs; };
  // void set_covariates(std::shared_ptr<const Eigen::MatrixXd> _covariates) { covariates = _covariates; };

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum, bool add) override; //{ return; };
  
  // Class Members
  double data_sum = 0;
  // std::shared_ptr<const Eigen::MatrixXd> covariates;
  // std::shared_ptr<Eigen::VectorXd> reg_coeffs = std::make_shared<Eigen::VectorXd>();
};

#endif // SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_LIKELIHOOD_H_