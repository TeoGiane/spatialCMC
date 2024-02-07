#ifndef BAYESMIX_HIERARCHIES_EMPTY_HIERARCHY_H_
#define BAYESMIX_HIERARCHIES_EMPTY_HIERARCHY_H_

#include "base_hierarchy.h"
#include "hierarchy_id.pb.h"
#include "likelihoods/empty_likelihood.h"
#include "priors/empty_prior_model.h"
#include "updaters/empty_updater.h"

/**
 * Empty hierarchy
 */

class EmptyHierarchy
    : public BaseHierarchy<EmptyHierarchy, EmptyLikelihood, EmptyPriorModel> {
 public:
  EmptyHierarchy() = default;
  ~EmptyHierarchy() = default;

  //! Returns the Protobuf ID associated to this class
  bayesmix::HierarchyId get_id() const override {
    return bayesmix::HierarchyId::Empty;
  }

  //! Sets the default updater algorithm for this hierarchy
  void set_default_updater() { updater = std::make_shared<EmptyUpdater>(); }

  //! Initializes state parameters to appropriate values
  void initialize_state() override {
    // Get hypers
    auto hypers = prior->get_hypers();
    // Initialize likelihood state
    State::Empty state;
    state.fake_field = hypers.fake_field;
    like->set_state(state);
  };

  //! Evaluates the log-marginal distribution of data in a single point
  //! @param hier_params  Pointer to the container of (prior or posterior)
  //! hyperparameter values
  //! @param datum        Point which is to be evaluated
  //! @return             The evaluation of the lpdf
  double marg_lpdf(ProtoHypersPtr hier_params, const Eigen::RowVectorXd &datum) const override {
    return 0;
  }
};

#endif  // BAYESMIX_HIERARCHIES_EMPTY_HIERARCHY_H_
