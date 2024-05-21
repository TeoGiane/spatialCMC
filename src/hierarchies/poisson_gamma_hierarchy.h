#ifndef SPATIALCMC_HIERARCHIES_POISSON_GAMMA_HIERARCHY_H_
#define SPATIALCMC_HIERARCHIES_POISSON_GAMMA_HIERARCHY_H_

#include "hierarchy_id.pb.h"

#include "hierarchies/likelihoods/poisson_likelihood.h"
#include "hierarchies/priors/gamma_prior_model.h"
#include "hierarchies/updaters/poisson_gamma_updater.h"
#include "cmc/spatialcmc_utils.h"
#include "src/hierarchies/base_hierarchy.h"

class PoissonGammaHierarchy : public BaseHierarchy<PoissonGammaHierarchy, PoissonLikelihood, GammaPriorModel> {
 public:
  PoissonGammaHierarchy() = default;
  ~PoissonGammaHierarchy() = default;
  bayesmix::HierarchyId get_id() const override {
    return bayesmix::HierarchyId::CUSTOM_HIERARCHY;
  };

  void set_default_updater() {
    updater = std::make_shared<PoissonGammaUpdater>();
  };

  void initialize_state() override {
    // Get hypers
    auto hypers = prior->get_hypers();
    // Initialize likelihood state
    State::Poisson state;
    // Questo init è a caso!
    state.rate = hypers.shape / hypers.rate;
    like->set_state(state, false);
  };
  
  double marg_lpdf(ProtoHypersPtr hier_params, const Eigen::RowVectorXd &datum) const override {
    // Unpack custom state
    bayesmix::GammaDistribution unp_hier_params; hier_params->custom_state().UnpackTo(&unp_hier_params);
    // auto unp_hier_params = spatialcmc::unpack_protobuf_any<bayesmix::GammaDistribution>(hier_params->custom_state());
    return stan::math::neg_binomial_lpmf(datum(0), unp_hier_params.shape(), unp_hier_params.rate());
  };
};

#endif  // SPATIALCMC_HIERARCHIES_POISSON_GAMMA_HIERARCHY_H_
