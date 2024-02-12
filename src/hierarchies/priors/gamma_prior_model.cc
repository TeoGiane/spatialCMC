#include "gamma_prior_model.h"

double GammaPriorModel::lpdf(const google::protobuf::Message &state_) {
  // Downcast state
  auto statecast = downcast_state(state_).custom_state();
  // Unpack custom state
  spatialcmc::PoissonState unp_state; statecast.UnpackTo(&unp_state);
  //auto unp_state = spatialcmc::unpack_protobuf_any<spatialcmc::PoissonState>(statecast);
  return stan::math::gamma_lpdf(unp_state.rate(), hypers->shape, hypers->rate);
}

State::Poisson GammaPriorModel::sample(ProtoHypersPtr hier_hypers) {
  // RNG
  auto &rng = bayesmix::Rng::Instance().get();
  // Output buffer
  State::Poisson out;
  // Get proper HyerarchyHypers
  auto params = (hier_hypers) ? hier_hypers->custom_state() : get_hypers_proto()->custom_state();
  // Unpack custom state
  bayesmix::GammaDistribution unp_params; params.UnpackTo(&unp_params);
  // auto unp_params = spatialcmc::unpack_protobuf_any<bayesmix::GammaDistribution>(params);
  // Sample
  out.rate = stan::math::gamma_rng(unp_params.shape(), unp_params.rate(), rng);
  return out;
}

void GammaPriorModel::set_hypers_from_proto(
    const google::protobuf::Message &hypers_) {
  // Downcast hypers
  auto & hyperscast = downcast_hypers(hypers_).custom_state();
  // Unpack custom state
  bayesmix::GammaDistribution unp_hypers; hyperscast.UnpackTo(&unp_hypers);
  // auto unp_hypers = spatialcmc::unpack_protobuf_any<bayesmix::GammaDistribution>(hyperscast);
  // Set hypers members
  hypers->shape = unp_hypers.shape();
  hypers->rate = unp_hypers.rate();
}

GammaPriorModel::ProtoHypersPtr GammaPriorModel::get_hypers_proto() const {
  ProtoHypersPtr out = std::make_shared<ProtoHypers>();
  bayesmix::GammaDistribution unp_hypers;
  unp_hypers.set_shape(hypers->shape);
  unp_hypers.set_rate(hypers->rate);
  out->mutable_custom_state()->PackFrom(unp_hypers);
  return out;
};

void GammaPriorModel::initialize_hypers() {
  if (prior->has_fixed_values()) {
    hypers->shape = prior->fixed_values().shape();
    hypers->rate = prior->fixed_values().rate();
  }
  else {
    throw std::invalid_argument("Unrecognized hierarchy prior");
  }
  // Checks
  if (hypers->shape <= 0) {
    throw std::runtime_error("shape must be positive");
  }
  if (hypers->rate <= 0) {
    throw std::runtime_error("rate_alpha must be positive");
  }
}