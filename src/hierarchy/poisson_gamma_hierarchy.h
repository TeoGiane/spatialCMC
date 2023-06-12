#ifndef POISSON_GAMMA_HIERARCHY_H_
#define POISSON_GAMMA_HIERARCHY_H_

#include "poisson_likelihood.h"
#include "gamma_prior_model.h"
#include "poisson_gamma_updater.h"
#include "hierarchy_id.pb.h"
#include "src/hierarchies/base_hierarchy.h"

class PoissonGammaHierarchy : public BaseHierarchy<PoissonGammaHierarchy, PoissonLikelihood, GammaPriorModel> {
 public:
  // PoissonGammaHierarchy(double shape_, double rate_) {
  //   auto prior = std::make_shared<GammaPriorModel>(shape_, rate_);
  //   set_prior(prior);
  // };
  PoissonGammaHierarchy() = default;
  ~PoissonGammaHierarchy() = default;

  bayesmix::HierarchyId get_id() const override {
    return bayesmix::HierarchyId::UNKNOWN_HIERARCHY;
  }

  void set_default_updater() {
    updater = std::make_shared<PoissonGammaUpdater>();
  }

  void initialize_state() override {
    // Get hypers
    auto hypers = prior->get_hypers();
    // Initialize likelihood state
    State::Poisson state;
    // Questo init è a caso!
    state.gamma = hypers.shape / hypers.rate;
    like->set_state(state, false);
  };

  // void set_covariates(std::shared_ptr<const Eigen::MatrixXd> _covariates) {
  //   like->set_covariates(_covariates);
  // };
  
  // void set_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) { like->set_reg_coeffs(_reg_coeffs); };
  
  double marg_lpdf(ProtoHypersPtr hier_params, const Eigen::RowVectorXd &datum) const override {
    double shape = hier_params->general_state().data(0);
    double rate = hier_params->general_state().data(1);
    // Eigen::VectorXd reg_coeffs = *(like->get_reg_coeffs());
    // double nb_beta = log(rate) - covariate * reg_coeffs;
    // return stan::math::neg_binomial_lpmf(datum(0), shape, std::exp(nb_beta));
    return stan::math::neg_binomial_lpmf(datum(0), shape, rate);
  }
};

#endif  // POISSON_GAMMA_HIERARCHY_H_
