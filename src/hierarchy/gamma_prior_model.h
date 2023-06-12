#ifndef GAMMA_PRIOR_MODEL_H_
#define GAMMA_PRIOR_MODEL_H_

#include <memory>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <vector>

#include "algorithm_state.pb.h"
#include "poisson_gamma_prior.pb.h"

#include "src/hierarchies/priors/base_prior_model.h"
#include "poisson_likelihood.h"
#include "src/utils/rng.h"

namespace Hyperparams {
  struct Gamma {
    double shape, rate;
  };
} // namespace Hyperparams

class GammaPriorModel : public BasePriorModel<GammaPriorModel, State::Poisson, Hyperparams::Gamma, spatialcmc::PoissonGammaPrior> {
 public:
  using AbstractPriorModel::ProtoHypers;
  using AbstractPriorModel::ProtoHypersPtr;

  GammaPriorModel() = default;
  // GammaPriorModel(double _shape = -1, double _rate = -1);
  ~GammaPriorModel() = default;
  double lpdf(const google::protobuf::Message &state_) override;
  State::Poisson sample(ProtoHypersPtr hier_hypers = nullptr) override;
  void update_hypers(const std::vector<bayesmix::AlgorithmState::ClusterState> &states) override { return; };
  void set_hypers_from_proto(const google::protobuf::Message &hypers_) override;
  ProtoHypersPtr get_hypers_proto() const override;

 protected:
  // double shape, rate;
  void initialize_hypers() override;
};

#endif // GAMMA_PRIOR_MODEL_H_
