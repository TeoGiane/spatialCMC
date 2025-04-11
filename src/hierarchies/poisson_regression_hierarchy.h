#ifndef SPATIALCMC_HIERARCHIES_POISSON_REGRESSION_HIERARCHY_H_
#define SPATIALCMC_HIERARCHIES_POISSON_REGRESSION_HIERARCHY_H_

#include "hierarchy_id.pb.h"

#include "cmc/spatialcmc_utils.h"
#include "hierarchies/likelihoods/poisson_regression_likelihood.h"
#include "hierarchies/priors/gamma_prior_model.h"
#include "hierarchies/updaters/poisson_regression_updater.h"
#include "src/hierarchies/base_hierarchy.h"

class PoissonRegHierarchy : public BaseHierarchy<PoissonRegHierarchy, PoissonRegLikelihood, GammaPriorModel> {
 public:
  PoissonRegHierarchy() = default;
  ~PoissonRegHierarchy() = default;
  
  bayesmix::HierarchyId get_id() const override {
    return bayesmix::HierarchyId::CUSTOM_HIERARCHY;
  };

  void set_default_updater() {
    updater = std::make_shared<PoissonRegUpdater>();
  };

  void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) {
    like->set_reg_coeffs(_reg_coeffs);
  };

  void set_covariates(const Eigen::MatrixXd & _covariates) {
    like->set_covariates(_covariates);
  };

  Eigen::VectorXd get_reg_coeffs() const {
    return like->get_reg_coeffs();
  };

  Eigen::MatrixXd get_covariates(void) const {
    return like->get_covariates();
  }

  double get_data_sum() const {
    return like->get_data_sum();
  };

  double get_cov_exp_sum() const {
    return like->get_cov_exp_sum();
  };

  void initialize_state() override {
    // Get hypers
    auto hypers = prior->get_hypers();
    // Initialize likelihood state
    State::Poisson state;
    state.rate = hypers.shape / hypers.rate;
    like->set_state(state, false);
  };
  
  // Compute marginal lpdf
  virtual double marg_lpdf(ProtoHypersPtr hier_params,
                           const Eigen::RowVectorXd &datum,
                           const Eigen::RowVectorXd &covariate) const override {
    double offset = covariate(0);
    Eigen::RowVectorXd cov_vector = covariate(Eigen::seq(1,Eigen::last));
    Eigen::VectorXd reg_coeffs = like->get_reg_coeffs();
    auto gamma_params = spatialcmc::unpack_to<bayesmix::GammaDistribution>(hier_params->custom_state());
    double rate = gamma_params.rate() / (offset*std::exp(cov_vector*reg_coeffs));
    return stan::math::neg_binomial_lpmf(datum(0), gamma_params.shape(), rate);
  };

};

#endif  // SPATIALCMC_HIERARCHIES_POISSON_GAMMA_HIERARCHY_H_
