#ifndef SPATIALCMC_HIERARCHIES_LIKELIHOODS_STATES_POISSON_STATE_H
#define SPATIALCMC_HIERARCHIES_LIKELIHOODS_STATES_POISSON_STATE_H

#include <google/protobuf/stubs/casts.h>

#include "poisson_state.pb.h"

#include "spatialcmc_utils.h"
#include "src/hierarchies/likelihoods/states/base_state.h"

namespace State {

  class Poisson : public BaseState {
   public:
    // Members
    double rate;

    // Methods
    using ProtoState = bayesmix::AlgorithmState::ClusterState;

    ProtoState get_as_proto() const override {
      ProtoState out;
      spatialcmc::PoissonState unp_state;
      unp_state.set_rate(rate);
      out.mutable_custom_state()->PackFrom(unp_state);
      return out;
    }

    void set_from_proto(const ProtoState &state_, bool update_card) override {
      if (update_card) {
        card = state_.cardinality();
      }
      auto unp_state = spatialcmc::unpack_protobuf_any<spatialcmc::PoissonState>(state_.custom_state());
      rate = unp_state.rate();
    }
  };

};  // namespace State

#endif // SPATIALCMC_HIERARCHIES_LIKELIHOODS_STATES_POISSON_STATE_H