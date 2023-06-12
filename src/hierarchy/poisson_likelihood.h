#ifndef POISSON_LIKELIHOOD_H_
#define POISSON_LIKELIHOOD_H_

#include <memory>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <google/protobuf/stubs/casts.h>

#include "src/hierarchies/likelihoods/states/base_state.h"
#include "src/hierarchies/likelihoods/base_likelihood.h"

namespace State {

  class Poisson : public BaseState {
   public:
    double gamma;
    using ProtoState = bayesmix::AlgorithmState::ClusterState;

    ProtoState get_as_proto() const override {
      ProtoState out;
      out.mutable_general_state()->add_data(gamma);
      return out;
    }

    void set_from_proto(const ProtoState &state_, bool update_card) override {
      if (update_card) {card = state_.cardinality();}
      gamma = state_.general_state().data(0);
    }
  };

};  // namespace State

class PoissonLikelihood : public BaseLikelihood<PoissonLikelihood, State::Poisson> {
 public:
  PoissonLikelihood() = default;
  ~PoissonLikelihood() = default;
  bool is_multivariate() const override { return false; };
  bool is_dependent() const override { return false; };
  void clear_summary_statistics() override; // { return; };

  // Getters and Setters
  double get_gamma() const { return state.gamma; };
  int get_data_sum() const;
  // double get_exp_cov_sum() const;

  // std::shared_ptr<Eigen::VectorXd> get_reg_coeffs() const { return reg_coeffs; };
  // std::shared_ptr<const Eigen::MatrixXd> get_covariates() const { return covariates; };
  // void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { *reg_coeffs = _reg_coeffs; };
  // void set_covariates(std::shared_ptr<const Eigen::MatrixXd> _covariates) { covariates = _covariates; };

 protected:
  double compute_lpdf(const Eigen::RowVectorXd &datum) const override;
  void update_sum_stats(const Eigen::RowVectorXd &datum, bool add) override; //{ return; };
  double data_sum = 0;
  // std::shared_ptr<const Eigen::MatrixXd> covariates;
  // std::shared_ptr<Eigen::VectorXd> reg_coeffs = std::make_shared<Eigen::VectorXd>();
};

#endif // POISSON_LIKELIHOOD_H_