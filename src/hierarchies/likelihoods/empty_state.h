#ifndef BAYESMIX_HIERARCHIES_LIKELIHOODS_STATES_EMPTY_STATE_H_
#define BAYESMIX_HIERARCHIES_LIKELIHOODS_STATES_EMPTY_STATE_H_

#include <memory>
#include <stan/math/rev.hpp>

#include "algorithm_state.pb.h"
#include "base_state.h"
#include "src/utils/proto_utils.h"

namespace State {

//! A univariate location-scale state with parametrization (mean, var)
//! The unconstrained representation corresponds to (mean, log(var))
class Empty : public BaseState {
 public:
  double fake_field = 0;

  using ProtoState = bayesmix::AlgorithmState::ClusterState;

  void set_from_proto(const ProtoState &state_, bool update_card) override {
    if (update_card) {
      card = state_.cardinality();
    }
    fake_field = state_.general_state().data(0);
  }

  ProtoState get_as_proto() const override {
    ProtoState state;
    state.mutable_general_state()->set_size(1);
    state.mutable_general_state()->mutable_data()->Add(fake_field);
    return state;
  }

};

}  // namespace State

#endif  // BAYESMIX_HIERARCHIES_LIKELIHOODS_STATES_UNI_LS_STATE_H_
