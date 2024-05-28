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
  void set_reg_coeffs(Eigen::VectorXd *reg_coeffs) { reg_coeffs_ptr = reg_coeffs; };
  void set_covariates(Eigen::MatrixXd *covariates) { covariates_ptr = covariates; };
  // Getters
  int get_data_sum() const { return data_sum; };
  double get_cov_exp_sum() const {
    std::vector<double> adds(card); unsigned int k = 0;
    for (auto && i : cluster_data_idx) {
      adds[k++] = log(covariates_ptr->row(i)(0))
                + (covariates_ptr->row(i)(Eigen::seq(1,Eigen::last))) * (*reg_coeffs_ptr);
    }
    return std::exp(stan::math::log_sum_exp(adds));
  };
  Eigen::VectorXd get_reg_coeffs() const { return *reg_coeffs_ptr; };
  Eigen::MatrixXd get_covariates() const { return *covariates_ptr; };
  // std::shared_ptr<Eigen::VectorXd> get_offsets() const { return offsets; };
  // std::shared_ptr<Eigen::VectorXd> get_reg_coeffs() const { return reg_coeffs; };
  // std::shared_ptr<const Eigen::MatrixXd> get_covariates() const { return covariates; };
  // void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { *reg_coeffs = _reg_coeffs; };
  // void set_covariates(std::shared_ptr<const Eigen::MatrixXd> _covariates) { covariates = _covariates; };

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum,
                      const Eigen::RowVectorXd &covariate) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum,
                        const Eigen::RowVectorXd &covariate,
                        bool add) override;

  // std::shared_ptr<Eigen::VectorXd> offsets;
  // std::shared_ptr<const Eigen::MatrixXd> covariates;
  // std::shared_ptr<Eigen::VectorXd> reg_coeffs = std::make_shared<Eigen::VectorXd>();

  // Common regression coefficient
  Eigen::VectorXd *reg_coeffs_ptr = nullptr;
  Eigen::MatrixXd *covariates_ptr = nullptr;
  // Summary statistics
  int data_sum = 0;
  // double cov_exp_sum = 0;
  // double cluster_lpdf = 0;
};

#endif // SPATIALCMC_HIERARCHIES_LIKELIHOODS_POISSON_REGRESSION_LIKELIHOOD_H_