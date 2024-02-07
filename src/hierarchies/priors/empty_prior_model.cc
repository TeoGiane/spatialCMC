#include "empty_prior_model.h"

void EmptyPriorModel::initialize_hypers() {
  hypers->fake_field = prior->fake_field();
}

double EmptyPriorModel::lpdf(const google::protobuf::Message &state_) {
  return 0;
}

State::Empty EmptyPriorModel::sample(ProtoHypersPtr hier_hypers) {
  State::Empty out;
  return out;
}

void EmptyPriorModel::update_hypers(const std::vector<bayesmix::AlgorithmState::ClusterState> &states) {
  return;
}

void EmptyPriorModel::set_hypers_from_proto(const google::protobuf::Message &hypers_) {
  auto &hyperscast = downcast_hypers(hypers_).empty_state();
  hypers->fake_field = hyperscast.fake_field();
  // hypers->fake_field = hyperscast.data(0);
}

std::shared_ptr<bayesmix::AlgorithmState::HierarchyHypers> EmptyPriorModel::get_hypers_proto() const {
  auto out = std::make_shared<bayesmix::AlgorithmState::HierarchyHypers>();
  out->mutable_empty_state()->set_fake_field(hypers->fake_field);
  // out->mutable_general_state()->set_size(1);
  //out->mutable_general_state()->mutable_data()->Add(hypers->fake_field);
  return out;
}
