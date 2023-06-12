#include "spp_mixing.h"
#include "src/utils/proto_utils.h"
// #include "utils/beta.h"

void sPPMixing::update_state(const std::vector<std::shared_ptr<AbstractHierarchy>> &unique_values,
                             const std::vector<unsigned int> &allocations) {

  auto priorcast = cast_prior();
  if (priorcast->has_fixed_value()) { return; }
  else { throw std::invalid_argument("Unrecognized mixing prior"); }
}

double sPPMixing::mass_existing_cluster(const unsigned int n, const unsigned int n_clust,
                                        const bool log, const bool propto,
                                        const std::shared_ptr<AbstractHierarchy> hier,
                                        const Eigen::RowVectorXd &covariate) const {
  // Define useful quantities
  std::set<int> tmp(hier->get_data_idx()); // necessario
  std::vector<int> idx(tmp.begin(), tmp.end());
  int added_edges = covariate(idx).sum();

  // Compute mass existing cluster
  double out;
  if (log) {
    out = hier->get_log_card() + state.lambda * added_edges;
    if (!propto)
      throw(std::runtime_error("Normalizing constants are not available. Set propto to true."));
  } else {
    out = 1.0 * hier->get_card() * std::exp(state.lambda * added_edges);
    if (!propto)
      throw(std::runtime_error("Normalizing constants are not available. Set propto to true."));
  }
  return out;
}

double sPPMixing::mass_new_cluster(const unsigned int n, const unsigned int n_clust,
                                   const bool log, const bool propto,
                                   const Eigen::RowVectorXd &covariate) const {
  // Compute mass new cluster
  double out;
  if (log) {
    out = state.logtotmass;
    if (!propto)
      throw(std::runtime_error("Normalizing constants are not available. Set propto to true."));
  } else {
    out = state.totalmass;
    if (!propto)
      throw(std::runtime_error("Normalizing constants are not available. Set propto to true."));
  }
  return out;
}

void sPPMixing::set_state_from_proto(const google::protobuf::Message &state_) {
  auto &statecast = downcast_state(state_);
  state.totalmass = statecast.general_state().data(0);
  state.logtotmass = std::log(state.totalmass);
  state.lambda = statecast.general_state().data(1);
}

std::shared_ptr<bayesmix::MixingState> sPPMixing::get_state_proto() const {
  auto out = std::make_shared<bayesmix::MixingState>();
  out->mutable_general_state()->add_data(state.totalmass);
  out->mutable_general_state()->add_data(state.lambda);
  return out;
}

void sPPMixing::initialize_state() {

  auto priorcast = cast_prior();
  
  if (priorcast->has_fixed_value()) {
    state.totalmass = priorcast->fixed_value().totalmass();
    if (state.totalmass <= 0) {
      throw std::invalid_argument("Total mass parameter must be > 0");
    }
    state.logtotmass = std::log(state.totalmass);
    state.lambda = priorcast->fixed_value().lambda();
  }

  else {
    throw std::invalid_argument("Unrecognized mixing prior");
  }
}
