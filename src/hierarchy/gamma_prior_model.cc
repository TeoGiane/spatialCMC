#include "gamma_prior_model.h"

// GammaPriorModel::GammaPriorModel(double _shape, double _rate) : shape(_shape), rate(_rate) {
//   create_empty_prior();
// };

double GammaPriorModel::lpdf(const google::protobuf::Message &state_) {
  double gamma = downcast_state(state_).general_state().data()[0];
  return stan::math::gamma_lpdf(gamma, hypers->shape, hypers->rate);
}

State::Poisson GammaPriorModel::sample(ProtoHypersPtr hier_hypers) {
  auto &rng = bayesmix::Rng::Instance().get();
  State::Poisson out;

  auto params = (hier_hypers) ? hier_hypers->general_state() : get_hypers_proto()->general_state();
  double shape = params.data()[0];
  double rate = params.data()[1];
  out.gamma = stan::math::gamma_rng(shape, rate, rng);
  // std::cout << "shape: " << shape << std::endl;
  // std::cout << "rate: " << rate << std::endl;
  // std::cout << "gamma: " << out.gamma << std::endl << std::endl;
  return out;
}

void GammaPriorModel::set_hypers_from_proto(
    const google::protobuf::Message &hypers_) {
  auto &hyperscast = downcast_hypers(hypers_).general_state();
  hypers->shape = hyperscast.data()[0];
  hypers->rate = hyperscast.data()[1];
};

GammaPriorModel::ProtoHypersPtr GammaPriorModel::get_hypers_proto() const {
  ProtoHypersPtr out = std::make_shared<ProtoHypers>();
  out->mutable_general_state()->mutable_data()->Add(hypers->shape);
  out->mutable_general_state()->mutable_data()->Add(hypers->rate);
  return out;
};

void GammaPriorModel::initialize_hypers() {

  if (prior->has_fixed_values()) {
    // Set values
    hypers->shape = prior->fixed_values().shape();
    hypers->rate = prior->fixed_values().rate();
    // Checks
    if (hypers->shape <= 0) {
      throw std::runtime_error("shape must be positive");
    }
    if (hypers->rate <= 0) {
      throw std::runtime_error("rate_alpha must be positive");
    }
  }

  else {
    throw std::invalid_argument("Unrecognized hierarchy prior");
  }
}