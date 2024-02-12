#ifndef BAYESMIX_HIERARCHIES_PRIORS_EMPTY_PRIOR_MODEL_H_
#define BAYESMIX_HIERARCHIES_PRIORS_EMPTY_PRIOR_MODEL_H_

#include <memory>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <vector>

#include "empty_hier.pb.h"

#include "hierarchies/likelihoods/states/empty_state.h"
#include "src/hierarchies/priors/base_prior_model.h"
#include "src/utils/rng.h"

namespace Hyperparams {
    struct Empty {};
} // namespace Hyperparams

/**
 * An empty prior model
 */

class EmptyPriorModel
    : public BasePriorModel<EmptyPriorModel, State::Empty, Hyperparams::Empty,
                            spatialcmc::EmptyPrior> {
 public:
  using AbstractPriorModel::ProtoHypers;
  using AbstractPriorModel::ProtoHypersPtr;

  EmptyPriorModel() = default;
  ~EmptyPriorModel() = default;

  double lpdf(const google::protobuf::Message &state_) override;

  State::Empty sample(ProtoHypersPtr hier_hypers = nullptr) override;

  void update_hypers(const std::vector<bayesmix::AlgorithmState::ClusterState>
                         &states) override;

  void set_hypers_from_proto(
      const google::protobuf::Message &hypers_) override;

  std::shared_ptr<bayesmix::AlgorithmState::HierarchyHypers> get_hypers_proto()
      const override;

 protected:
  void initialize_hypers() override;
};

#endif  // BAYESMIX_HIERARCHIES_PRIORS_NIG_PRIOR_MODEL_H_
